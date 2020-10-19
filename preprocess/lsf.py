#!/usr/bin/python -u
import os
import re
import json
import argparse

def find_files(path, pattern):
    for root, dirnames, filenames in os.walk(path):
        for filename in filenames:
            if re.match(pattern, filename):
                yield os.path.join(root, filename)

parser = argparse.ArgumentParser(description='GRO-seq LSF runner')
parser.add_argument('tag', nargs="?", default="no-tag", help="All of the output files are prefixed with this tag. This parameter can be omited altogether or \"no-tag\" can can be specified as a tag name with a same effect (default: no-tag)")
parser.add_argument('--pattern', dest="pattern", help="If set only samples with specified regular expression pattern will be evaluated (default: all fastaq files in subdirectories are evaluated)")
parser.add_argument('--overwrite', dest="overwrite", action='store_true',  help="By default GRO-seq will not overwrite old results and fail. If this flag is specified old results are overwritten")
args = parser.parse_args()

tag = "no-tag" if args.tag=="" else args.tag

fasta_files = list(find_files(".", r'^.*\.fastq$'))
if args.pattern is not None:
    fasta_files = [f for f in fasta_files if re.search(args.pattern, f)]

runtime_limit = 240 * len(fasta_files)

bsub_script = """#!/usr/bin/python -u
#BSUB -J GROseq.{tag}[1-{fasta_files_n}]%2
#BSUB -W {runtime_limit}
#BSUB -R "rusage[mem=8GB]"
#BSUB -R "span[hosts=1]"
#BSUB -n 12,24
#BSUB -e lsf_output.%J.%I.err
#BSUB -o lsf_output.%J.%I.out

import os
import re
import multiprocessing

job_index = int(os.environ.get('LSB_JOBINDEX'))
job_hosts = os.environ.get('LSB_HOSTS').split(" ")


tag = "{tag}"
prefix = tag+"_" if tag != "no-tag" else ""
overwrite = {overwrite}
fasta_files = {fasta_files}   
fasta = fasta_files[job_index-1]
sample = re.sub(r"\.fastq$", "", os.path.basename(fasta))
output = prefix + sample
output_exists = os.path.isfile("conv_rpkm/{{output}}.gmean.bed".format(output=output))
cmd = 'singularity exec --bind {{pwd}}:/mount ./groseq_latest.sif groseq --threads {{threads}} -f {{fasta}} -a mm10.refGene.longest_transcripts.bed -g ./mm10 -o {{output}} --chromInfo mm10.chrom.sizes'.format(pwd=os.getcwd(), fasta=fasta, output=output, threads=len(job_hosts))

print("Running {{job_name}} ({{job_id}}:{{job_index}})...".format(job_name=os.environ.get('LSB_JOBNAME'), job_id=os.environ.get('LSB_JOBID'), job_index=job_index))
print("Host: " + ", ".join(set(job_hosts)))
print("Available CPU: " + str(len(job_hosts)))
print("Sample: " + sample)
print("Output: {{output}} (exists: {{exists}})".format(output=output, exists="yes" if output_exists else "no"))
print('> ' + cmd)

if not overwrite and output_exists:
    print("Output already exists. Use --overwrite flag or specify a different tag")
    exit(1)
else:
    os.system(cmd)
    pass
"""

print(bsub_script.format(fasta_files_n=len(fasta_files), fasta_files=json.dumps(fasta_files),
                         runtime_limit=runtime_limit, tag=tag, overwrite=args.overwrite))
