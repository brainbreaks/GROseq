#!/usr/bin/python -u
import os
import re
import sys
import json
import math
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
parser.add_argument('--max-slots', dest="max_slots", type=int, default=8,  help="Describes how many jobs can run in parallel at most (default: 8)")
args = parser.parse_args()

tag = "no-tag" if args.tag=="" else args.tag

fasta_files = list(find_files(".", r'^.*\.fastq$'))
if args.pattern is not None:
    fasta_files = [f for f in fasta_files if re.search(args.pattern, f)]

#
# Calculate optimal resource requirements
#
max_slots = args.max_slots
fasta_sizes = [os.path.getsize(f)/1073741824. for f in fasta_files]
memory_limit = int(math.ceil(1.5*sum(sorted(fasta_sizes, reverse=True)[:1]) + 4))
runtime_limit = int(math.ceil(sum(fasta_sizes)*15))


prefix = tag+"_" if tag != "no-tag" else ""

any_output_exists = False
for fasta in fasta_files:
    output = prefix + re.sub(r"\.fastq$", "", os.path.basename(fasta))
    output_file = "conv_rpkm/{output}.gmean.bed".format(output=output)
    if os.path.isfile(output_file):
        sys.stderr.write("!!!! OUTPUT '{output}' EXISTS ({output_file}), NOT EXECUTING !!!!\nTry using --overwrite argument or specify a different tag!\n".format(output=output, output_file=output_file))
        any_output_exists = True

if any_output_exists and not args.overwrite:
    exit(2)


bsub_script = """#!/usr/bin/python -u
#BSUB -J groseq-{tag}[1-{fasta_files_n}]%{max_slots}
#BSUB -W {runtime_limit}
#BSUB -R "rusage[mem={memory_limit}GB]"
#BSUB -R "span[hosts=1]"
#BSUB -n 12,24
#BSUB -e lsf_output.groseq-{tag}.%J.%I.err
#BSUB -o lsf_output.groseq-{tag}.%J.%I.out

import os
import re
import multiprocessing
import datetime

# datetime object containing current date and time
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
expected_time1 = os.path.getsize(fasta)/1073741824*7
expected_time2 = os.path.getsize(fasta)/1073741824*9
cmd = 'singularity exec --bind {{pwd}}:/mount ./groseq_latest.sif groseq --threads {{threads}} -f {{fasta}} -a mm10.refGene.longest_transcripts.bed -g ./mm10 -o {{output}} --chromInfo mm10.chrom.sizes'.format(pwd=os.getcwd(), fasta=fasta, output=output, threads=len(job_hosts))

print("Running {{job_name}} ({{job_id}}:{{job_index}})...".format(job_name=os.environ.get('LSB_JOBNAME'), job_id=os.environ.get('LSB_JOBID'), job_index=job_index))
print("Host: " + ", ".join(set(job_hosts)))
print("Available CPU: " + str(len(job_hosts)))
print("Sample: " + sample)
print("Output: {{output}} (exists: {{exists}})".format(output=output, exists="yes" if output_exists else "no"))
print("Started at: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))	
print("Expected time: {{expected_time1}} - {{expected_time2}} minutes".format(expected_time1=expected_time1, expected_time2=expected_time2))
print("Expected finish: {{expected_finish1}} - {{expected_finish2}}".format(expected_finish1=(datetime.datetime.now() + datetime.timedelta(minutes=expected_time1)).strftime("%H:%M:%S"), expected_finish2=(datetime.datetime.now() + datetime.timedelta(minutes=expected_time2)).strftime("%H:%M:%S")))
print('> ' + cmd)

if not overwrite and output_exists:
    print("Output already exists. Use --overwrite flag or specify a different tag")
    exit(1)
else:
    os.system(cmd)
    pass
"""

print(bsub_script.format(fasta_files_n=len(fasta_files), fasta_files=json.dumps(fasta_files),
                         runtime_limit=runtime_limit,
                         memory_limit=memory_limit,
                         tag=tag, overwrite=args.overwrite, max_slots=max_slots))
