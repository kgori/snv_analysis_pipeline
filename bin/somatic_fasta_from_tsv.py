#!/usr/bin/env python
from __future__ import print_function
import itertools
import argparse

def parse_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('samplelist')
    parser.add_argument('outfile')
    parser.add_argument('-d', '--debug', action='store_true')
    return parser.parse_args()


def read_sample_list(samplefile):
    hosts = []
    tumours = []
    with open(samplefile) as infile:
        infile.readline()
        for line in infile:
            samplename, sampletype = line.strip().split()
            if sampletype == 'H':
                hosts.append(samplename)
            elif sampletype == 'T':
                tumours.append(samplename)
            else:
                raise ValueError('Unknown sample type: {}'.format(sampletype))
    return hosts, tumours


def process_header(line, hosts, tumours):
    """
    This line is the header line. Return the indices of the tumours and the hosts.
    """
    assert (line[0] == '#' or line[0] == "C")

    host_nr_idx = []
    tumour_nr_idx = []
    host_nv_idx = []
    tumour_nv_idx = []

    labels = line_to_fields(line)

    for i, label in enumerate(labels):

        if label == 'ALT':
            alt_idx = i

        elif label == 'REF':
            ref_idx = i

        elif label.endswith('.nv') or label.endswith('.nr'):
            print(f'Label = {label}')
            samplename, labeltype = label.split('.')

            if labeltype == 'nr':
                if samplename in hosts:
                    host_nr_idx.append(i)
                elif samplename in tumours:
                    tumour_nr_idx.append(i)
                else:
                    raise ValueError('Unknown sample name: {}'.format(samplename))
            elif labeltype == 'nv':
                if samplename in hosts:
                    host_nv_idx.append(i)
                elif samplename in tumours:
                    tumour_nv_idx.append(i)
                else:
                    raise ValueError('Unknown sample name: {}'.format(samplename))

    return {
        'ref_idx' : ref_idx,
        'alt_idx' : alt_idx,
        'host_nr_idx' : host_nr_idx,
        'host_nv_idx' : host_nv_idx,
        'tumour_nr_idx' : tumour_nr_idx,
        'tumour_nv_idx' : tumour_nv_idx,
    }


def to_integer(lst):
    """
    Convert list of strings (representing integer numbers) to integers
    """
    return [int(n) for n in lst]


def line_to_fields(line, sep='\t'):
    return line.rstrip().split(sep)


def get_values_at_index(line, index):
    fields = line_to_fields(line)

    if isinstance(index, int):
        return fields[index]

    return [fields[i] for i in index]


def all_greater_than(n, lst):
    """
    Return True if all values in lst are > n
    """
    return all(value > n for value in lst)


def all_less_than(n, lst):
    """
    Return True if all values in lst are < n
    """
    return all(value < n for value in lst)


def all_equal_to(n, lst):
    """
    Return True if all values in lst == n
    """
    return all(value == n for value in lst)


def get_sample_name(header, index):
    names = line_to_fields(header)
    return names[index].replace('.nv', '')


def grouper(n, iterable):
    iterable = iter(iterable)
    return iter(lambda: list(itertools.islice(iterable, n)), [])


def split_line(s,n=120):
    return '\n'.join(''.join(grp) for grp in grouper(n, s))


_codes = {
    'A': 1,
    'C': 2,
    'G': 4,
    'T': 8,
    'N': 15
}


def char_to_intcode(c):
    return _codes[c]


def is_variant(site):
    """
    Return True if site is variant
    """
    codes = [char_to_intcode(c) for c in site]
    init = codes[0]
    for c in codes[1:]:
        init &= c
    return init == 0


if __name__ == '__main__':
    from collections import defaultdict
    import sys

    args = parse_cli()

    seqs = defaultdict(list)

    hosts, tumours = read_sample_list(args.samplelist)

    with open(args.infile) as infile:
        for linenum, line in enumerate(infile, start=1):
            if line.startswith('#CHROM') or line.startswith("CHROM"):
                header = line
                indices = process_header(header, hosts, tumours)
            else:
                host_nv = to_integer(get_values_at_index(line, indices['host_nv_idx']))
                if all_less_than(3, host_nv):
                    host_nr = to_integer(get_values_at_index(line, indices['host_nr_idx']))
                    tumour_nr = to_integer(get_values_at_index(line, indices['tumour_nr_idx']))
                    if all_greater_than(9, host_nr + tumour_nr):
                        if args.debug:
                            print('Filters pass at line: {}'.format(linenum), file=sys.stderr)
                        # Make a list of tuples of (column index, alt read count value)
                        tumour_nv = list(zip(indices['tumour_nv_idx'], to_integer(get_values_at_index(line, indices['tumour_nv_idx']))))
                        refbase = get_values_at_index(line, indices['ref_idx'])
                        altbase = get_values_at_index(line, indices['alt_idx'])

                        # Check if is variant site
                        tmpseqs = dict()
                        for i, count in tumour_nv:
                            if count >=3:
                                tmpseqs[i] = altbase
                                #seqs[i].append(altbase)
                            else:
                                tmpseqs[i] = refbase
                                #seqs[i].append(refbase)

                        column = list(tmpseqs.values())
                        if is_variant(column):
                            for i, _ in tumour_nv:
                                seqs[i].append(tmpseqs[i])

    print("Writing output to: {}".format(args.outfile), file=sys.stdout)
    with open(args.outfile, 'w') as outfile:
        for i in seqs:
            name = get_sample_name(header, i)
            sequence = ''.join(seqs[i])
            outfile.write('>{}\n{}\n\n'.format(name, split_line(sequence)))
            outfile.flush()

