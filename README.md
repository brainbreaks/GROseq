 GRO-seq pipeline package (Wei, Pei-Chi group; DKFZ)
====================================================
Global Run-On Sequencing (GRO-Seq) pipeline for analyzing transcription activity of genes from engaged RNA polymerase. 


Table of Contents
=================

  * Using GRO-seq docker image
    * [Singularity](#singularity) 
      * [Pull GRO-seq image from DockerHUB](#singularity-pull)
      * [Run container](#singularity-run)
    * [Docker](#docker) 
      * [Pull GRO-seq image from DockerHUB](#docker-pull)
      * [Run container](#docker-run)
    * [Build GRO-seq docker image](#build) 
      * [Download required packages and files](#build-download)
      * [Build docker image](#build-build)


<a name="singularity">Singularity</a>
====================================================

<a name="singularity-pull">Pull GRO-seq image from data server</a>
----------------------------------------------------
```console
singularity pull docker://sandrejev/groseq:latest
```

<a name="singularity-run">Run container</a>
----------------------------------------------------
Before running GROseq pipeline you will need to obtain genome(fasta), bowtie index (*.bt2), chromosome sizes (*.chrom.sizes) and annotation. For mm9, mm10 and hg19 these can be downloaded automatically with `download` command. The only exception being
annotation *.bed file.
```console
singularity exec --bind ${PWD}:/mount ./groseq_latest.sif download mm10
```

Run groseq pipeline. Keep in mind that annotation file (`-a` flag is not provided automatically)
```console
singularity exec --bind ${PWD}:/mount ./groseq_latest.sif groseq -f AS-512172-LR-52456/fastq/AS-512172-LR-52456_R1.fastq -a mm10.refGene.longest_transcripts.bed -g ./mm10 -o "singularity1" --chromInfo mm10.chrom.sizes
```




<a name="docker">Docker</a>
====================================================

<a name="docker-pull">Pull GRO-seq image from DockerHUB</a>
----------------------------------------------------
```console
docker pull sandrejev:groseq
```

<a name="docker-run">Run container</a>
----------------------------------------------------
Before running GROseq pipeline you will need to obtain genome(fasta), bowtie index (*.bt2), chromosome sizes (*.chrom.sizes) and annotation. For mm9, mm10 and hg19 these can be downloaded automatically with `download` command. The only exception being
annotation *.bed file.
```console
docker run -v ${PWD}:/mount -u $(id -g ${USER}):$(id -g ${USER}) -it --entrypoint download groseq mm10
```

Run groseq pipeline. Keep in mind that annotation file (`-a` flag is not provided automatically)
```console
docker run -v ${PWD}:/mount -u $(id -g ${USER}):$(id -g ${USER}) -it groseq -f AS-512172-LR-52456/fastq/AS-512172-LR-52456_R1.fastq -a mm10.refGene.longest_transcripts.bed -g ./mm10 -o AS-512172-LR-52456 --chromInfo mm10.chrom.sizes
```

<a name="docker-inspect">Inspect container</a>
----------------------------------------------------
```console
docker run -v ${PWD}:/mount -u $(id -g ${USER}):$(id -g ${USER}) -it --entrypoint bash groseq
```



<a name="build">Building GRO-seq docker image</a>
====================================================

<a name="build-download">Download required packages and files</a>
----------------------------------------------------
To successfully build GRO-seq image first required libraries and files must be downloaded. This can be done by running following command
```console
python download.py
```

<a name="build-build">Build docker image</a>
----------------------------------------------------
To build Docker image you need to execute
```console
docker build --squash --build-arg http_proxy="http://www.inet.dkfz-heidelberg.de:80" --build-arg http_proxy="http://www.inet.dkfz-heidelberg.de:80" --rm -t groseq .
```

Convert cached docker image to singularity (for local testing)
```console
singularity pull docker-daemon:groseq:latest
```

