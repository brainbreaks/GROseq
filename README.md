GRO-seq pipeline packege (Wei, Pei-Chi group; DKFZ)
====================================================
Global Run-On Sequencing (GRO-Seq) pipeline for analyzing transcription activity of genes from engaged RNA polymerase. 

Download required packages and files
-------------------
To successfully build GRO-seq image first required libraries and files must be downloaded. This can be done by running 
following command
> python download.py

Build docker image
-------------------
To build Docker image you need to execute
> docker build --squash --build-arg http_proxy="http://www.inet.dkfz-heidelberg.de:80" --build-arg http_proxy="http://www.inet.dkfz-heidelberg.de:80" --rm -t groseq .

Import/Export docker image
> docker save groseq > groseq.tar
> docker load < groseq.tar

Run container
-------------------
To run a container first put location of reference genome fasta file (genome.fa) and all of the bowtie2 index files (genome.*.bt2) into separate directory. You will need
to reference them with -g parameter (e.g. ./mm10)
Run container. Current directory will be accessible inside the container
> docker run -v ${PWD}:/mount -u 0:$(id -g ${USER}) --name groseq_test -it groseq -f AS-512172-LR-52456/fastq/AS-512172-LR-52456_R1.fastq -g ./mm10 -o results --chromInfo mm10.chrom.sizes
> docker run -v ${PWD}:/mount -u 0:$(id -g ${USER}) --name <memorable_name> -it groseq -f ./path/to/sample -g ./path/to/genome -o ./path/to/results --chromInfo mm10.chrom.size

Inspect container
-------------------
> docker run -v ${PWD}:/mount -u 0:$(id -g ${USER}) -it --entrypoint bash groseq

 
 
 #>>> mkdir -p alignment conv_rpkm
#docker run -v ${PWD}:/mount -u 0:$(id -u ${USER}) -it --entrypoint bash groseq
#> docker run -v ${PWD}:/mount -u $(id -u ${USER}):$(id -g ${USER}) groseq -f SRR2820488_first100.fastq -g mm9 -o results.txt --chromInfo mm9.chrom.sizes
#>
#>/bin/homer/bin/makeTagDirectory
