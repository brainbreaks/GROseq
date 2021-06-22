 GRO-seq pipeline package (Wei, Pei-Chi group; DKFZ)
====================================================
Global Run-On Sequencing (GRO-Seq) pipeline for analyzing transcription activity of genes from engaged RNA polymerase. 


Table of Contents
=================

  * Using GRO-seq docker image
    * [Singularity](#singularity) 
      * [Pull GRO-seq image from DockerHUB](#singularity-pull)
      * [Run container](#singularity-run)
      * [Inspect container](#singularity-inspect)
    * [Docker](#docker) 
      * [Pull GRO-seq image from DockerHUB](#docker-pull)
      * [Run container](#docker-run)
      * [Inspect container](#docker-inspect)
    * [Build GRO-seq docker image](#build) 
      * [Download required packages and files](#build-download)
      * [Build docker image](#build-build)
      * [Push docker image to Docker HUB](#build-build)
      * [Convert cached docker image to singularity (for local testing)](#build-convert)

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
singularity exec -B `pwd` groseq_latest.sif download mm10
```

Run groseq pipeline. Keep in mind that annotation file (`-a` flag) is created automatically for you from geneRef.gtf
```console
singularity exec -B `pwd` groseq_latest.sif groseq -f AS-512172-LR-52456/fastq/AS-512172-LR-52456_R1.fastq -a mm10.refGene.bed -g ./mm10 -o "singularity1" --chromInfo mm10.chrom.sizes
```

You can create annotation file with different clipping at the start or end of the transcript using `longest-transcript` command
```console
singularity exec -B `pwd` groseq_latest.sif longest-transcript mm10.refGene.gtf.gz mm10.refGene.bed --clip-start=50
```

You can extract rpkm using `extract-rpkm` command (done automatically as part of the pipeline)
```console
singularity exec -B `pwd` groseq_latest.sif extract-rpkm -a mm10.refGene.bed -o AS-512172-LR-52456_R1 --clip-start=50
```

For convenience GRO-seq image contains a script that can be used to run the pipeline on LSF cluster 
```console
# Run GRO-seq pipeline on all *.fastq files in the folder
singularity exec -B `pwd` groseq_latest.sif lsf | bsub

# Run GRO-seq pipeline on all *.fastq files in the folder that match the pattern. Additionaly prefix the results with tag PREFIX
singularity exec -B `pwd` groseq_latest.sif lsf PREFIX --pattern "AS-512178" | bsub 
```


<a name="singularity-inspect">Inspect container</a>
----------------------------------------------------
```console
singularity shell groseq_latest.sif
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

Run groseq pipeline. Keep in mind that annotation file (`-a` flag) is created automatically for you from geneRef.gtf
```console
docker run -v ${PWD}:/mount -u $(id -g ${USER}):$(id -g ${USER}) -it groseq -f AS-512172-LR-52456/fastq/AS-512172-LR-52456_R1.fastq -a mm10.refGene.bed -g ./mm10 -o AS-512172-LR-52456 --chromInfo mm10.chrom.sizes
```

You can create annotation file with different clipping at the start or end of the transcript using `longest-transcript` command
```console
docker run -v ${PWD}:/mount -u $(id -g ${USER}):$(id -g ${USER}) -it --entrypoint longest-transcript groseq mm10.refGene.gtf.gz mm10.refGene.bed --clip-start=50
```

You can extract rpkm using `extract-rpkm` command (done automatically as part of the pipeline)
```console
docker run -v ${PWD}:/mount -u $(id -g ${USER}):$(id -g ${USER}) -it --entrypoint extract-rpkm groseq -a mm10.refGene.bed -o AS-512172-LR-52456_R1
```



For convenience GRO-seq image contains a script that can be used to run the pipeline on LSF cluster.
```console
# Run GRO-seq pipeline on all *.fastq files in the folder
docker run -v ${PWD}:/mount -u $(id -g ${USER}):$(id -g ${USER}) -it --entrypoint lsf groseq | bsub

# Run GRO-seq pipeline on all *.fastq files in the folder that match the pattern. Additionaly prefix the results with tag PREFIX
docker run -v ${PWD}:/mount -u $(id -g ${USER}):$(id -g ${USER}) -it --entrypoint lsf groseq PREFIX --pattern "AS-512178" | bsub 
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
python download.py dependencies
```

<a name="build-build">Build docker image</a>
----------------------------------------------------
To build Docker image you need to execute
```console
docker build --squash --build-arg http_proxy="http://www.inet.dkfz-heidelberg.de:80" --build-arg https_proxy="http://www.inet.dkfz-heidelberg.de:80" --rm -t sandrejev/groseq:latest .
```

<a name="build-push">Push docker image to Docker HUB</a>
----------------------------------------------------
```console
docker login
docker push sandrejev/groseq:latest
```

<a name="build-convert">Convert cached docker image to singularity (for local testing)</a>
----------------------------------------------------
```console
singularity pull docker-daemon:sandrejev/groseq:latest
```