FROM python:2.7.18
LABEL version="1"
LABEL software="GRO-seq"
MAINTAINER Sergej Andrejev <sandrejev@gmail.com>

# Install required standard libraries
RUN apt-get update
RUN apt-get install make gcc g++ perl lynx \
    ca-certificates ghostscript \
    libncurses5-dev libncursesw5-dev zip unzip git wget \
    zlib1g-dev libbz2-dev liblzma-dev libopenblas-dev -y;

# Install python dependencies
RUN pip2 install RSeQC numpy


ENV DESTINATION=/bin
#ENV HOMER_URL=http://homer.ucsd.edu/homer/configureHomer.pl
#ENV BOWTIE2_URL=https://github.com/BenLangmead/bowtie2/releases/download/v2.2.9/bowtie2-2.2.9-linux-x86_64.zip
#ENV SAMTOOLS_URL=https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2
#ENV BEDTOOLS_URL=https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
#ENV GFOLD_URL=https://anaconda.org/bioconda/gfold/1.1.4/download/linux-64/gfold-1.1.4-gsl2.2_2.tar.bz2
#ENV BXPYTHON_URL=https://github.com/bxlab/bx-python/archive/master.zip
#ENV GSL_URL=http://www.artfiles.org/gnu.org/gsl/gsl-2.2.1.tar.gz
#ENV bigWigToWig_URL=http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToWig
#ENV wigToBigWig_URL=http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
ENV HOMER_URL=http://127.0.0.1:8000/dependencies/homer/configureHomer.pl
ENV BOWTIE2_URL=http://127.0.0.1:8000/dependencies/bowtie2-2.2.9-linux-x86_64.zip
ENV SAMTOOLS_URL=http://127.0.0.1:8000/dependencies/samtools-1.6.tar.bz2
ENV BEDTOOLS_URL=http://127.0.0.1:8000/dependencies/bedtools-2.29.1.tar.gz
ENV GFOLD_URL=http://127.0.0.1:8000/dependencies/gfold-1.1.4-gsl2.2_2.tar.bz2
ENV GSL_URL=http://127.0.0.1:8000/dependencies/gsl-2.2.1.tar.gz
ENV bigWigToWig_URL=http://127.0.0.1:8000/dependencies/bigWigToWig
ENV wigToBigWig_URL=http://127.0.0.1:8000/dependencies/wigToBigWig
ENV BXPYTHON_URL=http://127.0.0.1:8000/dependencies/bx-python.zip

# Copy GRO-seq files to docker image
COPY data/hg19/*  /data/
COPY data/mm9/*  /data/
COPY data/dependencies  /data/dependencies

# Create a web server, It is used to access local copies of packages and libraries
RUN echo 'cd /data && python -m SimpleHTTPServer 8000 &> /dev/null &\nsleep 1' > /bin/ws && chmod 755 /bin/ws


# Speed up make buy using multiple cores
RUN cp /usr/bin/make /usr/bin/make.bak
RUN echo 'make.bak --jobs=32 $@' > /usr/bin/make

# Install bowtie2
RUN cd $DESTINATION && \
    ws && wget -q $BOWTIE2_URL -O $DESTINATION/bowtie2.zip && \
    unzip $DESTINATION/bowtie2.zip -d $DESTINATION/bowtie2_dir && \
    mkdir $DESTINATION/bowtie2 && \
    cp -r $DESTINATION/bowtie2_dir/*/* $DESTINATION/bowtie2 && \
    rm -Rf $DESTINATION/bowtie2.zip $DESTINATION/bowtie2_dir

# Install bx-python dependencies
RUN ws && pip uninstall -y -q bx-python && \
    wget -q $BXPYTHON_URL -O $DESTINATION/bx-python.zip && \
    unzip $DESTINATION/bx-python.zip -d $DESTINATION/bx-python_dir && \
    cd $DESTINATION/bx-python_dir/* && python2 setup.py install -q  && \
    rm -Rf $DESTINATION/bx-python_dir
# git clone https://github.com/bxlab/bx-python.git && cd bx-python && python2 setup.py install
# bam2bw and bam2wig should be used from here

# Install SAMTOOLS
RUN ws && wget -q $SAMTOOLS_URL -O $DESTINATION/samtools_dir.zip;  \
    mkdir $DESTINATION/samtools_dir && tar xjf $DESTINATION/samtools_dir.zip -C $DESTINATION/samtools_dir; mv $DESTINATION/samtools_dir/*/* $DESTINATION/samtools_dir; \
    cd $DESTINATION/samtools_dir/ && make; \
    mv $DESTINATION/samtools_dir/samtools $DESTINATION/samtools; \
    rm -rf $DESTINATION/samtools_dir*;

# Install BEDTOOLS
RUN ws && wget -q $BEDTOOLS_URL -O $DESTINATION/bedtoos_dir.zip;  \
    mkdir $DESTINATION/bedtoos_dir && tar xzf $DESTINATION/bedtoos_dir.zip -C $DESTINATION/bedtoos_dir; mv $DESTINATION/bedtoos_dir/*/* $DESTINATION/bedtoos_dir; \
    cd $DESTINATION/bedtoos_dir/ && make; \
    mv $DESTINATION/bedtoos_dir/bin $DESTINATION/bedtools; \
    rm -rf $DESTINATION/bedtoos_dir;

# Install HOMER. I update configureHomer.pl because every time a package is downloaded it is overwriten
COPY /data/dependencies/homer  $DESTINATION/homer
RUN ws && cd $DESTINATION/homer && wget -q $HOMER_URL -O $DESTINATION/homer/configureHomer.pl && perl configureHomer.pl -keepScript -install homer;
RUN ws && cd $DESTINATION/homer && wget -q $HOMER_URL -O $DESTINATION/homer/configureHomer.pl && perl configureHomer.pl -install mm9;
RUN ws && cd $DESTINATION/homer && wget -q $HOMER_URL -O $DESTINATION/homer/configureHomer.pl && perl configureHomer.pl -install hg19;

# Install bigWigToWig and wigToBigWig from bigWig
RUN ws && wget $bigWigToWig_URL -O /bin/bigWigToWig && chmod 755 $DESTINATION/bigWigToWig
RUN ws && wget $wigToBigWig_URL -O /bin/wigToBigWig && chmod 755 $DESTINATION/wigToBigWig

# Install gfold
RUN ws && wget $GFOLD_URL -O $DESTINATION/gfold_dir.tar.bz2 && \
    mkdir $DESTINATION/gfold_dir && tar xvjf $DESTINATION/gfold_dir.tar.bz2 -C $DESTINATION/gfold_dir && \
    mv $DESTINATION/gfold_dir/bin/gfold $DESTINATION/gfold && chmod 755 $DESTINATION/gfold; \
    rm -rf $DESTINATION/gfold_dir*;

# Install GSL 2.2.1
RUN ws && wget $GSL_URL && \
    tar xzf gsl-2.2.1.tar.gz && \
    cd gsl-2.2.1 && ./configure && make && make install && \
    ldconfig

# Create shortcut for running GRO-seq
RUN echo '#!/bin/bash \n python2 /bin/GROSeqPL_v3_updated.py "$@"' > $DESTINATION/groseq; chmod 755 $DESTINATION/groseq

# Copy pipeline source files
COPY src/*  $DESTINATION/

WORKDIR /mount

ENTRYPOINT ["groseq"]
