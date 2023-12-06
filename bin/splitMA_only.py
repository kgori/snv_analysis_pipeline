#!/usr/bin/env python

# Splits VCF records describing multiallelic SNVs/indels into individual records

# INPUT
# vcfFile: path to VCF file

 

import sys
import os
import re
import gzip

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


if len(sys.argv) != 2:
    print '\nsplitMA_only.py /path/to/variants.vcf\n'
    sys.exit(0)


script, vcfFile = sys.argv


# Compose paths to output file
if vcfFile.endswith('.gz'):
    outFile = vcfFile[:-7] + '.split.vcf'
else:
    outFile = vcfFile[:-4] + '.split.vcf'
print '\nInput file:  ', vcfFile
print 'Output file: ', outFile


# Split any multi-allelic variant into multiple bi-allelic variants
print '\nSplitting multi-allelic variants...'
with FileReader(vcfFile) as vcf:
    with open(outFile, 'w') as out:
        count = 0
        for line in vcf:
            if line.startswith('#'):
                out.write(line)
            else:
                col = line.strip().split('\t')
                chrom = col[0]
                pos = col[1]
                theId = col[2]
                ref = col[3]
                alts = col[4].split(',')
                qual = col[5]
                filters = col[6]
                info = col[7]
                format = col[8]       
                theRest = list(col[9:])

                # Don't process bi-allellic variants
                if len(alts) == 1:
                    out.write(line)
        
                else:
                    count = count + 1
                    for ind, alt in enumerate(alts):
                        # Create the new info
                        infoElem = info.split(';')
                        fr = infoElem[1][3:].split(',')
                        nf = infoElem[7][3:].split(',')
                        nr = infoElem[8][3:].split(',')
                        pp = infoElem[9][3:].split(',')
                        tr = infoElem[17][3:].split(',')
                        newInfo = ';'.join([ infoElem[0], 'FR='+fr[ind], infoElem[2], infoElem[3], 
                                             infoElem[4], infoElem[5], infoElem[6], 'NF='+nf[ind],
                                             'NR='+nr[ind], 'PP='+pp[ind], infoElem[10], infoElem[11],
                                             infoElem[12], infoElem[13], infoElem[14], infoElem[15], 
                                             infoElem[16], 'TR='+tr[ind], infoElem[18], infoElem[19] ])
                        
                        # Create the new theRest
                        newRest = []
                        for elem in theRest:
                            stats = elem.split(':')
                            nr = stats[4].split(',')
                            nv = stats[5].split(',')
                            # Normal cases
                            if len(nr) > ind and len(nv) > ind:
                                newRest.append(':'.join([ stats[0], stats[1], stats[2], stats[3], nr[ind], nv[ind] ]))
                            # Cases where Platypus gets no data: './.:0,0,0:0:0:0:0'
                            else:
                                newRest.append(':'.join([ stats[0], stats[1], stats[2], stats[3], nr[0], nv[0] ]))
                                
                        newLine = '\t'.join([ chrom, pos, theId, ref, alt, qual, filters, 
                                              newInfo, format, '\t'.join(newRest) ])    
                        out.write(newLine + '\n')

print count, 'multi-allelic variants found'
print 'Done!\n'
