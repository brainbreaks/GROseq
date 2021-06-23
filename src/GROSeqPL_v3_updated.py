import sys, os, os.path
import argparse
import multiprocessing
from misc import call, AttrDict
import bam2bw
import time
import re
import glob
import pandas as pd
import numpy as np
from collections import OrderedDict

def parse_args():
    parser = argparse.ArgumentParser(description='GRO-Seq Pipeline')
    parser.add_argument("-f", dest = "file", type = str, required = False, help = "Read input file" )
    parser.add_argument("-f1", dest = "file_1", type = str, required = False, help = "Read input file 1 of pair" )
    parser.add_argument("-f2", dest = "file_2", type = str, required = False, help = "Read input file 2 of pair" )
    parser.add_argument("-o", dest = "output", type = str, required = True, help = "Output file name" )
    parser.add_argument("-g", dest = "genome", type = str, required = False, help = "Genome file to align to (do not include .fa)")
    parser.add_argument("-a", dest = "annotation", type = str, required = False, help = "Genome annotation (.bed format)")
    parser.add_argument("-cg", dest = "custom_genome", type = str, required = False, help = "Custom Genome file to align to (includes .fa)" )
    parser.add_argument("--chromInfo", dest = "chromInfo", type = str, required = True, help = "chromInfo file" )
    parser.add_argument("--threads", dest="threads", default=multiprocessing.cpu_count(), help="Number of threads to use [default is all available - {}]".format(multiprocessing.cpu_count()))
    args = parser.parse_args()

    files = {}
    if (args.file_1 or args.file_2) and args.file:
        print >>sys.stderr, "\nCan not use single-ended fasta (-f) file and two pair-end files (-f1, -f2) at the same time\n"
        sys.exit()
    if not ((args.file_1 and args.file_2) or args.file):
        print >>sys.stderr, "\nSingle-ended fasta (-f) file or two pair-end files (-f1, -f2) were not provided\n"
        sys.exit()
    if args.file:
        files["-f"] = args.file
    elif args.file_1 and args.file_2:
        files["-f1"] = args.file_1
        files["-f2"] = args.file_2
    else:
        print >>sys.stderr, "\nBoth files have to be provided in case of pair end sequencing\n"
        sys.exit()

    if args.genome:
        files["-g"] = args.genome + ".1.bt2"
    elif args.custom_genome:
        files["-cg"] = args.custom_genome + ".1.bt2"
    else:
        print >>sys.stderr, "\nPath to genome is not specified (-g or -cg)\n"
        sys.exit()

    if args.annotation:
        files["-a"] = args.annotation

    if args.chromInfo:
        files["--chromInfo"] = args.chromInfo
    else:
        print >>sys.stderr, "\nChromosome sizes are not specified (--chromInfo)\n"
        sys.exit()

    if not args.output:
        print >>sys.stderr, "\nOutput is not specified\n"
        sys.exit()

    f_msg = ""
    for f, p in files.items():
        if not os.path.isfile(p):
            f_msg += "\n{f} = {p}".format(f=f, p=p)
    if f_msg:
        print >>sys.stderr, "\nFollowing files were not found: {f_msg}".format(f_msg=f_msg)
        sys.exit()

    return args

def extract_rpkm(args):
    annotations_cols = OrderedDict([("rpkm_chrom", str), ("rpkm_start", np.float64), ("rpkm_end", np.float64), ("rpkm_symbol", str), ("rpkm_ignore", str), ("rpkm_sense", str)])
    rpkm_cols =  OrderedDict([("rpkm_symbol", str), ("rpkm_ignore", str), ("rpkm_reads", np.float64), ("rpkm_length", np.int64), ("rpkm_rpkm", np.float64)])

    annotation = pd.read_csv(args.annotation, sep='\t', names=annotations_cols.keys(), dtype=annotations_cols, usecols=["rpkm_chrom", "rpkm_start", "rpkm_end", "rpkm_symbol", "rpkm_sense"])

    results = pd.DataFrame()
    for path in [p.format(**vars(args)) for p in ["conv_rpkm/{output}.pos3.cnt", "conv_rpkm/{output}.neg3.cnt"]]:
        path_basename = os.path.basename(path)
        data = pd.read_csv(path, sep="\t", names=rpkm_cols.keys(), dtype=rpkm_cols, usecols=["rpkm_symbol", "rpkm_reads", "rpkm_length", "rpkm_rpkm"])
        data['rpkm_output'] = args.output
        data['rpkm_sense'] = "+" if re.sub(".*\.(.*)3.cnt", r"\1", path_basename)=="pos" else "-"

        results_sense = data.merge(annotation, how='inner', on=["rpkm_symbol", "rpkm_sense"])
        results = results.append(results_sense)

    results.to_csv('conv_rpkm/{output}.rpkm.tsv'.format(**vars(args)), sep="\t", index=False)


def read_concat(args):
    index_1 = args.file_1.find("L001")
    cat_R1 = args.file_1.replace(args.file_1[index_1:],"L00*_R1*")
    index_2 = args.file_2.find("L001")
    cat_R2 = args.file_2.replace(args.file_2[index_2:],"L00*_R2*")
    args.file_1 = args.file_1.replace("L001_","")
    args.file_2 = args.file_2.replace("L001_","")

    call("cat {0}> {1}".format(cat_R1,args.file_1))
    call("cat {0}> {1}".format(cat_R2,args.file_2))
    call("rm -f {0}".format(cat_R1[:-4]))
    return

def fq_to_bam(args):
    if args.custom_genome:
        # path, folder = os.path.split(args.custom_genome)
        args.genome = args.custom_genome

    if args.file_1:
        call("/bin/bowtie2/bowtie2 -p {threads} -x {genome} -1 {file_1} -2 {file_2} --non-deterministic -S alignment/{output}.sam".format(**vars(args)))
    else:
        call("/bin/bowtie2/bowtie2 -p {threads} -x {genome} -U {file} --non-deterministic -S alignment/{output}.sam".format(**vars(args)))
    call("/bin/samtools view -@ {threads} -bhS -F 4 alignment/{output}.sam -o alignment/{output}.bam".format(**vars(args)))
    call("/bin/samtools sort -@ {threads} alignment/{output}.bam -o alignment/{output}.sort.bam".format(**vars(args)))

    call("mv alignment/{output}.sort.bam alignment/{output}.bam".format(**vars(args)))
    call("/bin/samtools index -@ {threads} alignment/{output}.bam alignment/{output}.bam.bai".format(**vars(args)))
    call("/bin/samtools stats -@ {threads} alignment/{output}.bam > alignment/{output}.stats.txt".format(**vars(args)), shell=True)

    return

def wig_convert(args):
    #call("/bin/samtools view -b -F 0x10 alignment/{0}.bam -o alignment/{0}.pos.bam".format(args.output))
    #call("/bin/samtools view -b -f 0x10 alignment/{0}.bam -o alignment/{0}.neg.bam".format(args.output))
    # call("python2 /bin/bam2bw.py -b alignment/{0}.bam -g {1} -s T".format(args.output, args.chromInfo))

    bam2bw_args = AttrDict({
        "threads": args.threads,
        "bamfile":  "alignment/{}.bam".format(args.output),
        "genome": args.chromInfo,
        "strandspecific": "T",
        "normalization": "T"
    })

    bam2bw.parseFILE(bam2bw_args)

    return

def homer(args):
    if args.custom_genome:
        call("/bin/homer/bin/makeTagDirectory conv_rpkm/{0} alignment/{0}.bam -genome {1}.fa".format(args.output,args.custom_genome))
    else:
        call("/bin/homer/bin/makeTagDirectory conv_rpkm/{0} alignment/{0}.bam -genome {1}.fa".format(args.output,args.genome))
    # call("/bin/homer/bin/findPeaks conv_rpkm/{0} -style groseq -o auto -uniqmap /usr/local/homer/data/uniqmap/mm9-uniqmap".format(args.output))
    call("/bin/homer/bin/findPeaks conv_rpkm/{0} -style groseq -o auto".format(args.output))
    return

def convergent(args):
    call("grep -v \"#\" conv_rpkm/{0}/transcripts.txt | awk '{{OFS=\"\t\"}}{{print $2, $3, $4, $1, $6, $5}}' > conv_rpkm/{0}.transcript.bed".format(args.output), shell=True)
    call("/bin/bedtools/bedtools sort -i conv_rpkm/{0}.transcript.bed > see".format(args.output), shell=True)
    call("mv see conv_rpkm/{0}.transcript.bed".format(args.output), shell=True)
    call("grep -v \"+\" conv_rpkm/{0}.transcript.bed > conv_rpkm/{0}.transcripts.minus.bed".format(args.output), shell=True)
    call("grep \"+\" conv_rpkm/{0}.transcript.bed > conv_rpkm/{0}.transcripts.plus.bed".format(args.output), shell=True)
    call("/bin/bedtools/bedtools multiinter -i conv_rpkm/{0}.transcripts.plus.bed conv_rpkm/{0}.transcripts.minus.bed |grep \"1,2\" |awk \'{{OFS=\"\\t\"}}{{ if ($3-$2>100){{print $1, $2, $3}} }}\' > conv_rpkm/{0}.transcripts.convergent.bed".format(args.output), shell=True)
    return

def gmean_cal(args):
    call("awk '{{OFS=\"\\t\"}}{{printf \"%s\\t%s\\t%s\\t%s_%s_%s\\t.\\t.\\n\", $1, $2, $3, $1, $2, $3}}' conv_rpkm/{output}.transcripts.convergent.bed > conv_rpkm/{output}.transcripts.convergent.6.bed".format(**vars(args)), shell=True)
    call("/bin/samtools view -@ {threads} alignment/{output}.pos.bam | gfold count -annf BED -ann conv_rpkm/{output}.transcripts.convergent.6.bed -tag stdin -o conv_rpkm/{output}.pos.cnt".format(**vars(args)), shell=True)
    call("/bin/samtools view -@ {threads} alignment/{output}.neg.bam | gfold count -annf BED -ann conv_rpkm/{output}.transcripts.convergent.6.bed -tag stdin -o conv_rpkm/{output}.neg.cnt".format(**vars(args)), shell=True)
    call("/bin/samtools view -@ {threads} alignment/{output}.pos.bam | gfold count -annf BED -s T -ann {annotation} -tag stdin -o conv_rpkm/{output}.pos3.cnt".format(**vars(args)), shell=True)
    call("/bin/samtools view -@ {threads} alignment/{output}.neg.bam | gfold count -annf BED -s T -ann {annotation} -tag stdin -o conv_rpkm/{output}.neg3.cnt".format(**vars(args)), shell=True)
    call("python2 /bin/COVT.gmean_cal.py conv_rpkm/{output} > conv_rpkm/{output}.gmean.cnt".format(**vars(args)), shell=True)
    call("sort -gr -k 4 conv_rpkm/{output}.gmean.cnt > conv_rpkm/{output}.gmean.cnt.sort".format(**vars(args)), shell=True)
    call("python2 /bin/Add.col.py conv_rpkm/{output}.gmean.cnt.sort > conv_rpkm/{output}.gmean.bed".format(**vars(args)), shell=True)
    return

def GROSeqPL(args):
    call("mkdir -p alignment conv_rpkm")
    #if (os.path.isfile(args.file_1) and "L001" in args.file_1):
    #    read_concat(args)
    fq_to_bam(args)
    wig_convert(args)
    homer(args)
    convergent(args)
    gmean_cal(args)
    extract_rpkm(args)

    call("rm -rf alignment/{output}.sam".format(**vars(args)))


def main():
    args = parse_args()

    start = time.time()
    GROSeqPL(args)
    end = time.time()
    print("Total time: {:.1f} minutes".format((end - start)/60))

if __name__ == "__main__":
    main()
