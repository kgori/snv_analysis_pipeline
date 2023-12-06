#!/usr/bin/env python
from __future__ import print_function

from collections import defaultdict
import gzip
from somatypus_utility_lib import get_variant_from_line
import sys

def is_gzip(filename):
    """
    Check the first two bytes of the file to see if
    it is gzipped
    """
    with open(filename, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'

class FileReader:
    def __init__(self, filename):
        if is_gzip(filename):
            self._file = gzip.open(filename, 'rt')
        else:
            self._file = open(filename, 'r')

    def __enter__(self):
        return self._file

    def __exit__(self, type, value, traceback):
        if not self._file.closed:
            self._file.close()

def parse_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('input_files', nargs='+', help='Input file')
    return parser.parse_args()

def main():
    args = parse_args()
    d = defaultdict(int)
    for input_file in args.input_files:
        print("Processing {}".format(input_file), file=sys.stdout)
        with FileReader(input_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                variant = get_variant_from_line(line)
                nv = line.strip().split('\t')[-1].split(':')[-1]
                try:
                    nv = int(nv)
                except ValueError:
                    continue
                d[variant] += nv
    with open('filter_list.txt', 'w') as outf:
        outf.write('\t'.join(['CHROM', 'POS', 'REF', 'ALT', 'NV_TOTAL\n']))
        for variant in sorted(d):
            nv = d[variant]
            outf.write('{}\t{}\t{}\t{}\t{}\n'.format(variant.chrom, variant.pos, variant.ref, variant.alt, nv))

if __name__ == '__main__':
    main()
