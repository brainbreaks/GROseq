"""
Convert bam file to bigwig file
"""

import sys, os
import argparse
import glob
from misc import call
import multiprocessing

def parse_args():
    parser = argparse.ArgumentParser(description='Convert BAM file to BIGWIG files')
    parser.add_argument("-b", dest = "bamfile", type = str, required = True,
                              help = "bam file" )
    parser.add_argument("-g", dest = "genome", type = str, default = "mm9", required = True,
                              help = "specify the genome: mm9, hg19 or abs path (mm9)" )
    parser.add_argument("-s", "--strandspecific", dest = "strandspecific", 
                              type = str, default = "F", choices={'F', 'T'},
                              help = "Divide into strand plus and minus bigwig files (F)" )
    parser.add_argument("-n", "--normalization", dest = "normalization", 
                              type = str, default = "T", choices={'F', 'T'},
                              help = "Do normalization of bigwig to coverage \
                              of 10 million 100nt reads (T)" )
    parser.add_argument("--threads", dest="threads", default=multiprocessing.cpu_count(), help="Number of threads to use [default is all available - {}]".format(multiprocessing.cpu_count()))
    #parser.add_argument("-q", "MAP_QUAL", dest = "MAP_QUAL", type = int, default = 30,
    #                          help = "Minimum mapping quality for an alignment to be called 'uniquely mapped'. default=30" )
    args = parser.parse_args()
    return args


def parseFILE(args):
    genomefile = args.genome
    normdef = ''
    if args.normalization == 'T':
        normdef = ' -t 1000000000 '

    if args.strandspecific == 'F':
        prefix = args.bamfile.replace('.bam', '')
        call('python2 /bin/bam2wig.py -i %s -o %s -s %s %s' % (args.bamfile, prefix, genomefile, normdef))
        call('rm -rf %s.wig' % prefix)
    else:
        prefix_pos = args.bamfile.replace('.bam', '.pos')
        call('/bin/samtools view -@ %s -F 0x10 -hb -o %s.bam %s' % (args.threads, prefix_pos, args.bamfile))
        call('/bin/samtools sort -@ %s %s.bam -o %s.sort.bam' % (args.threads, prefix_pos, prefix_pos))
        call('mv %s.sort.bam %s.bam' % (prefix_pos, prefix_pos))
        call('/bin/samtools index -@ %s %s.bam' % (args.threads, prefix_pos))
        call('python2 /bin/bam2wig.py -i %s.bam -o %s -s %s %s' % (prefix_pos, prefix_pos, genomefile, normdef))

        # Need to add minus symbol in minus bw file
        prefix_neg = args.bamfile.replace('.bam', '.neg')
        call('/bin/samtools view -@ %s -f 0x10 -hb -o %s.bam %s' % (args.threads, prefix_neg, args.bamfile))
        call('/bin/samtools sort -@ %s %s.bam -o %s.sort.bam' % (args.threads, prefix_neg, prefix_neg))
        call('mv %s.sort.bam %s.bam' % (prefix_neg, prefix_neg))
        call('/bin/samtools index -@ %s %s.bam' % (args.threads, prefix_neg))
        call('python2 /bin/bam2wig.py -i %s.bam -o %s -s %s %s' % (prefix_neg, prefix_neg, genomefile, normdef))
        fout = open('%s.2.wig' % prefix_neg, 'w')
        for line in open('%s.wig' % prefix_neg):
            if line.startswith('variableStep'):
                fout.write(line)
            else:
                fout.write('%s\t-%s\n' % (line.split()[0], line.split()[1]))
        fout.close()
        call('wigToBigWig -clip %s.2.wig %s %s.bw' % (prefix_neg, genomefile, prefix_neg))

        call('rm -rf %s.bam* %s.bam*' % (prefix_pos, prefix_neg))
        call('rm -rf %s.wig %s.wig' % (prefix_pos, prefix_neg))
        call('rm -rf %s.2.wig' % (prefix_neg))

def main():
    args = parse_args()
    parseFILE(args)

if __name__ == '__main__':
    main()
