FROM python:2.7.18-slim
LABEL version="1"
LABEL software="GRO-seq"
MAINTAINER Sergej Andrejev <sandrejev@gmail.com>

# Install required standard libraries
RUN apt-get update
RUN apt-get install make gcc g++ perl \
    ca-certificates ghostscript \
    zip unzip wget \
    zlib1g-dev libbz2-dev liblzma-dev libopenblas-dev -y && \
    apt-get clean;

# Install python dependencies
#RUN pip2 install --proxy="${http_proxy}" numpy


ENV DESTINATION=/bin
ENV MAKEFLAGS='-j 32'

# Copy GRO-seq files to docker image
#COPY data/hg19/*  /d/
#COPY data/mm9/*  /d/
#COPY data/mm10/*  /d/
COPY data/dependencies/bedtools-2.29.1.tar.gz  $DESTINATION/bedtoos_dir.zip
COPY data/dependencies/RSeQC-2.6.4.tar.gz  $DESTINATION/RSeQC-2.6.4.tar.gz
COPY data/dependencies/bigWigToWig  $DESTINATION/bigWigToWig
COPY data/dependencies/wigToBigWig  $DESTINATION/wigToBigWig
COPY data/dependencies/bowtie2-2.2.9-linux-x86_64.zip $DESTINATION/bowtie2.zip
COPY data/dependencies/bx-python.zip $DESTINATION/bx-python.zip
COPY data/dependencies/gfold-1.1.4-gsl2.2_2.tar.bz2 $DESTINATION/gfold_dir.tar.bz2
COPY data/dependencies/gsl-2.2.1.tar.gz $DESTINATION/gsl-2.2.1.tar.gz
COPY data/dependencies/samtools-1.6.tar.bz2 $DESTINATION/samtools_dir.zip
COPY data/dependencies/numpy-1.16.6.tar.gz $DESTINATION/numpy-1.16.6.tar.gz
COPY data/dependencies/cython-0.29.19.tar.gz $DESTINATION/cython-0.29.19.tar.gz
COPY data/dependencies/nose-1.0.0.tar.gz $DESTINATION/nose-1.0.0.tar.gz
COPY data/dependencies/pysam-0.16.0.tar.gz $DESTINATION/pysam-0.16.0.tar.gz





# Create a web server, It is used to access local copies of packages and libraries
# RUN chmod -R 755 /d/* && echo 'cd /d && python -m SimpleHTTPServer 8000 &> /dev/null &\nsleep 1' > /bin/ws && chmod 755 /bin/ws


# Speed up make buy using multiple cores
#RUN cp /usr/bin/make /usr/bin/make.bak
#RUN echo 'make.bak --jobs=32 $@' > /usr/bin/make

# Install cython
RUN cd $DESTINATION && \
    tar xzf cython-0.29.19.tar.gz && \
    cd cython-0.29.19 && python2 setup.py install && \
    cd $DESTINATION && \
    rm -Rf cython-0.29.19.tar.gz cython-0.29.19

# Install numpy
RUN cd $DESTINATION && \
    tar xzf numpy-1.16.6.tar.gz && \
    cd numpy-1.16.6 && python2 setup.py install && \
    cd $DESTINATION && \
    rm -Rf numpy-1.16.6.tar.gz numpy-1.16.6

# Install nose
RUN cd $DESTINATION && \
    tar xzf nose-1.0.0.tar.gz && \
    cd nose-1.0.0 && python2 setup.py install && \
    cd $DESTINATION && \
    rm -Rf nose-1.0.0.tar.gz nose-1.0.0

# Install pysam
RUN cd $DESTINATION && \
    tar xzf pysam-0.16.0.tar.gz && \
    cd pysam-0.16.0 && python2 setup.py install && \
    cd $DESTINATION && \
    rm -Rf pysam-0.16.0.tar.gz pysam-0.16.0


# Install bowtie2
RUN cd $DESTINATION && \
    unzip $DESTINATION/bowtie2.zip -d $DESTINATION/bowtie2_dir && \
    mkdir $DESTINATION/bowtie2 && \
    cp -r $DESTINATION/bowtie2_dir/*/* $DESTINATION/bowtie2 && \
    rm -Rf $DESTINATION/bowtie2.zip $DESTINATION/bowtie2_dir

# Install bx-python dependencies
#pip uninstall -y -q bx-python && \
RUN unzip $DESTINATION/bx-python.zip -d $DESTINATION/bx-python_dir && \
    cd $DESTINATION/bx-python_dir/* && python2 setup.py install -q  && \
    rm -Rf $DESTINATION/bx-python_dir* -Rf $DESTINATION/bx-python.zip
# git clone https://github.com/bxlab/bx-python.git && cd bx-python && python2 setup.py install
# bam2bw and bam2wig should be used from here

# Install RSeQC
RUN cd $DESTINATION && \
    tar xzf RSeQC-2.6.4.tar.gz && \
    cd RSeQC-2.6.4 && python2 setup.py install && \
    cd $DESTINATION && \
    rm -Rf RSeQC-2.6.4.tar.gz RSeQC-2.6.4


# Install SAMTOOLS
RUN mkdir $DESTINATION/samtools_dir && tar xjf $DESTINATION/samtools_dir.zip -C $DESTINATION/samtools_dir; mv $DESTINATION/samtools_dir/*/* $DESTINATION/samtools_dir; \
    cd $DESTINATION/samtools_dir/ && make; \
    mv $DESTINATION/samtools_dir/samtools $DESTINATION/samtools; \
    rm -rf $DESTINATION/samtools_dir*;

# Install BEDTOOLS
RUN cd $DESTINATION && \
    mkdir $DESTINATION/bedtoos_dir && tar xzf $DESTINATION/bedtoos_dir.zip -C $DESTINATION/bedtoos_dir; mv $DESTINATION/bedtoos_dir/*/* $DESTINATION/bedtoos_dir; \
    cd $DESTINATION/bedtoos_dir/ && make; \
    mv $DESTINATION/bedtoos_dir/bin $DESTINATION/bedtools; \
    rm -rf $DESTINATION/bedtoos_dir* $DESTINATION/example $DESTINATION/doc;

# Install gfold
RUN cd $DESTINATION && \
    mkdir $DESTINATION/gfold_dir && tar xvjf $DESTINATION/gfold_dir.tar.bz2 -C $DESTINATION/gfold_dir && \
    mv $DESTINATION/gfold_dir/bin/gfold $DESTINATION/gfold && chmod 755 $DESTINATION/gfold; \
    rm -rf $DESTINATION/gfold_dir*;

# Install GSL 2.2.1
RUN cd $DESTINATION && \
    tar xzf gsl-2.2.1.tar.gz && \
    cd gsl-2.2.1 && ./configure && make && make install && \
    ldconfig && \
    cd $DESTINATION && \
    rm -Rf gsl-2.2.1.tar.gz gsl-2.2.1

# Install bigWigToWig and wigToBigWig from bigWig
RUN chmod 755 $DESTINATION/bigWigToWig
RUN chmod 755 $DESTINATION/wigToBigWig

# Copy pipeline source files
COPY src/*  $DESTINATION/

# Create shortcut for running GRO-seq
RUN echo '#!/bin/bash \n python2 /bin/GROSeqPL_v3_updated.py "$@"' > $DESTINATION/groseq; chmod 755 $DESTINATION/groseq

# Install HOMER. I update configureHomer.pl because every time a package is downloaded it is overwriten
RUN mkdir $DESTINATION/homer && \
    cd $DESTINATION/homer && \
    http_proxy="${http_proxy}" https_proxy="${https_proxy}" wget http://homer.ucsd.edu/homer/configureHomer.pl && \
    perl ./configureHomer.pl -install homer && \
    rm -Rf $DESTINATION/homer/cpp

RUN apt-get remove -y \
 libcairo2-dev \
 ca-certificates zip unzip git wget libcurl4-openssl-dev \
 tcl-dev tcl8.6-dev tk-dev tk8.6-dev tcl tcl8.6 tk tk8.6 \
 zlib1g-dev libbz2-dev liblzma-dev libbz2-dev \
 libjbig-dev libjpeg-dev libjpeg62-turbo-dev libpng-dev \
 libxml2-dev libxrender-dev libxslt1-dev libxss-dev \
 autotools-dev gcc g++ libgcc-8-dev libstdc++-8-dev libc-dev-bin libc6-dev libyaml-dev \
 libmount-dev dpkg-dev \
 x11proto-core-dev x11proto-dev x11proto-scrnsaver-dev x11proto-xext-dev xtrans-dev libx11-dev \
 libopenblas-dev && \
 apt-get autoremove -y && apt-get clean
# apt-get autoremove -y &&

#COPY /data/dependencies/homer/configureHomer.pl  $DESTINATION/homer/configureHomer.pl
#COPY /data/dependencies/homer/update.txt  $DESTINATION/homer/update.txt
#RUN mkdir $DESTINATION/homer
#COPY /data/homer_hg19_mm9_mm10
#COPY /data/homer_hg19_mm9_mm10 $DESTINATION/homer
#COPY /data/homer_hg19_mm9_mm10 $DESTINATION/homer
#RUN cd $DESTINATION/homer && perl configureHomer.pl -make
#RUN ws && cd $DESTINATION/homer && wget -q $HOMER_URL -O $DESTINATION/homer/configureHomer.pl && perl configureHomer.pl -keepScript -install homer;
#RUN ws && cd $DESTINATION/homer && wget -q $HOMER_URL -O $DESTINATION/homer/configureHomer.pl && perl configureHomer.pl -install mm9;
#RUN ws && cd $DESTINATION/homer && wget -q $HOMER_URL -O $DESTINATION/homer/configureHomer.pl && perl configureHomer.pl -install mm10;
#RUN ws && cd $DESTINATION/homer && wget -q $HOMER_URL -O $DESTINATION/homer/configureHomer.pl && perl configureHomer.pl -install hg19;

## fdisk
## apt
#dpkg
#dash
#coreutils
#bash
#base-passwd
#e2fsprogs

WORKDIR /mount

ENTRYPOINT ["groseq"]
