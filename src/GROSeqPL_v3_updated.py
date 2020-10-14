import sys, os
import argparse
import multiprocessing
from misc import call, AttrDict
import bam2bw
import time

def parse_args():
    parser = argparse.ArgumentParser(description='GRO-Seq Pipeline')
    parser.add_argument("-f", dest = "file", type = str, required = False, help = "Read input file" )
    parser.add_argument("-f1", dest = "file_1", type = str, required = False, help = "Read input file 1 of pair" )
    parser.add_argument("-f2", dest = "file_2", type = str, required = False, help = "Read input file 2 of pair" )
    parser.add_argument("-o", dest = "output", type = str, required = True, help = "Output file name" )
    parser.add_argument("-g", dest = "genome", type = str, required = False, help = "Genome file to align to (do not include .fa)")
    parser.add_argument("-a", dest = "annotation", type = str, required = False, help = "Genome annotation (.gtf/.gtf.tz file)")
    parser.add_argument("-cg", dest = "custom_genome", type = str, required = False, help = "Custom Genome file to align to (includes .fa)" )
    parser.add_argument("--chromInfo", dest = "chromInfo", type = str, required = True, help = "chromInfo file" )
    parser.add_argument("--threads", dest="threads", default=multiprocessing.cpu_count(), help="Number of threads to use [default is all available - {}]".format(multiprocessing.cpu_count()))
    args = parser.parse_args()
    return args

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
    call("/bin/samtools view alignment/{output}.pos.bam | gfold count -annf BED -s T -ann {annotation} -tag stdin -o conv_rpkm/{output}.pos3.cnt".format(**vars(args)), shell=True)
    call("/bin/samtools view alignment/{output}.neg.bam | gfold count -annf BED -s T -ann {annotation} -tag stdin -o conv_rpkm/{output}.neg3.cnt".format(**vars(args)), shell=True)
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

    call("rm -rf alignment/{output}.sam".format(**vars(args)))


def main():
    args = parse_args()

    start = time.time()
    GROSeqPL(args)
    end = time.time()
    print("Total time: {:.1f} minutes".format((end - start)/60))

main()
