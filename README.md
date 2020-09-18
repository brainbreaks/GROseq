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
> docker build -t groseq .

Run container
-------------------
Run container. Current directory will be accessible inside the container
> docker run -v ${PWD}:/mount -u $(id -u ${USER}):$(id -g ${USER}) groseq -f SRR2820488_first100.fastq -g /data/mm9 -o results.txt --chromInfo mm9.chrom.sizes

Inspect container
-------------------
> docker run -v ${PWD}:/mount -u $(id -u ${USER}):$(id -g ${USER}) -it --entrypoint bash groseq
 
