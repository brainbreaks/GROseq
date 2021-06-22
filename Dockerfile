FROM python:2.7.18-slim
LABEL version="1"
LABEL software="GRO-seq"
MAINTAINER Sergej Andrejev <sandrejev@gmail.com>

ARG http_proxy
ARG https_proxy
ENV DESTINATION=/bin
ENV MAKEFLAGS='-j 32'

# Install required standard libraries
RUN http_proxy="${http_proxy}" https_proxy="${https_proxy}" apt-get update
RUN http_proxy="${http_proxy}" https_proxy="${https_proxy}" apt-get install \
    locales locales-all \
    make gcc gfortran g++ perl libcurl4-openssl-dev \
    ca-certificates ghostscript \
    zip unzip wget mime-support \
    libopenblas-base libopenblas-dev liblapack-dev liblapack-dev \
    libncurses6 libncursesw6 libncurses-dev \
    zlib1g-dev bzip2 libbz2-dev liblzma-dev -y && \
    apt-get clean && \
    locale-gen en_US.UTF-8 && \
    dpkg-reconfigure locales

ENV LANG=en_US.UTF-8
ENV LC_ADDRESS=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8
ENV LC_COLLATE=en_US.UTF-8
ENV LC_CTYPE=en_US.UTF-8
ENV LC_IDENTIFICATION=en_US.UTF-8
ENV LC_MEASUREMENT=en_US.UTF-8
ENV LC_MESSAGES=en_US.UTF-8
ENV LC_MONETARY=en_US.UTF-8
ENV LC_NAME=en_US.UTF-8
ENV LC_NUMERIC=en_US.UTF-8
ENV LC_PAPER=en_US.UTF-8
ENV LC_TELEPHONE=en_US.UTF-8
ENV LC_TIME=en_US.UTF-8



# Copy dependencies into docker image
COPY dependencies/bedtools-2.29.1.tar.gz $DESTINATION/bedtools-2.29.1.tar.gz
COPY dependencies/bowtie2-2.2.9-linux-x86_64.zip $DESTINATION/bowtie2-2.2.9-linux-x86_64.zip
COPY dependencies/bx-python-0.8.9.tar.gz $DESTINATION/bx-python-0.8.9.tar.gz
COPY dependencies/cython-0.29.19.tar.gz $DESTINATION/cython-0.29.19.tar.gz
COPY dependencies/gfold-1.1.4-gsl2.2_2.tar.bz2 $DESTINATION/gfold-1.1.4-gsl2.2_2.tar.bz2
COPY dependencies/gsl-2.2.1.tar.gz $DESTINATION/gsl-2.2.1.tar.gz
COPY dependencies/nose-1.0.0.tar.gz $DESTINATION/nose-1.0.0.tar.gz
COPY dependencies/numpy-1.16.6.tar.gz $DESTINATION/numpy-1.16.6.tar.gz
COPY dependencies/pysam-0.16.0.tar.gz $DESTINATION/pysam-0.16.0.tar.gz
COPY dependencies/RSeQC-2.6.4.tar.gz $DESTINATION/RSeQC-2.6.4.tar.gz
COPY dependencies/samtools-1.6.tar.bz2 $DESTINATION/samtools-1.6.tar.bz2
COPY dependencies/six-1.15.0.tar.gz $DESTINATION/six-1.15.0.tar.gz
COPY dependencies/requests-2.24.0.tar.gz $DESTINATION/requests-2.24.0.tar.gz
COPY dependencies/tqdm-4.49.0.tar.gz $DESTINATION/tqdm-4.49.0.tar.gz
COPY dependencies/pydevtools-0.8.0.tar.gz $DESTINATION/pydevtools-0.8.0.tar.gz
COPY dependencies/pandas-0.24.2.tar.gz $DESTINATION/pandas-0.24.2.tar.gz
COPY dependencies/dateutil-2.8.1.tar.gz $DESTINATION/dateutil-2.8.1.tar.gz
COPY dependencies/wigToBigWig $DESTINATION/wigToBigWig
COPY dependencies/bigWigToWig $DESTINATION/bigWigToWig

# Install cython
RUN cd $DESTINATION && \
    tar xzf cython-0.29.19.tar.gz && \
    cd cython-0.29.19 && python2 setup.py install && \
    cd $DESTINATION && \
    rm -Rf cython-0.29.19.tar.gz cython-0.29.19

# Install six
RUN cd $DESTINATION && \
    tar xzf six-1.15.0.tar.gz && \
    cd six-1.15.0 && python2 setup.py install && \
    cd $DESTINATION && \
    rm -Rf six-1.15.0.tar.gz six-1.15.0

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

# Install requests
RUN cd $DESTINATION && \
    tar xzf requests-2.24.0.tar.gz && \
    cd requests-2.24.0 && python2 setup.py install && \
    cd $DESTINATION && \
    rm -Rf requests-2.24.0.tar.gz requests-2.24.0


# Install tqdm
RUN cd $DESTINATION && \
    tar xzf tqdm-4.49.0.tar.gz && \
    cd tqdm-4.49.0 && python2 setup.py install && \
    cd $DESTINATION && \
    rm -Rf tqdm-4.49.0.tar.gz tqdm-4.49.0

# Install bowtie2
RUN cd $DESTINATION && \
    unzip $DESTINATION/bowtie2-2.2.9-linux-x86_64.zip -d $DESTINATION/bowtie2_dir && \
    mkdir $DESTINATION/bowtie2 && \
    cp -r $DESTINATION/bowtie2_dir/*/* $DESTINATION/bowtie2 && \
    rm -Rf $DESTINATION/bowtie2-2.2.9-linux-x86_64.zip $DESTINATION/bowtie2_dir


# Install bx-python dependencies
# bam2bw and bam2wig could be used from here
RUN cd $DESTINATION && \
    tar xzf bx-python-0.8.9.tar.gz && \
    cd bx-python-0.8.9 && python2 setup.py install && \
    cd $DESTINATION && \
    rm -Rf bx-python-0.8.9.tar.gz bx-python-0.8.9


# Install RSeQC
RUN cd $DESTINATION && \
    tar xzf RSeQC-2.6.4.tar.gz && \
    cd $DESTINATION/RSeQC-2.6.4 && python2 setup.py install && \
    cd $DESTINATION && \
    rm -Rf RSeQC-2.6.4.tar.gz RSeQC-2.6.4


# Install SAMTOOLS
RUN cd $DESTINATION && \
    tar xjf samtools-1.6.tar.bz2 && \
    cd $DESTINATION/samtools-1.6 && make && \
    cp $DESTINATION/samtools-1.6/samtools $DESTINATION/samtools && \
    rm -rf $DESTINATION/samtools-1.6.tar.bz2 $DESTINATION/samtools-1.6


# Install BEDTOOLS
RUN cd $DESTINATION && \
    tar xzf bedtools-2.29.1.tar.gz && \
    cd $DESTINATION/bedtools2 && make && \
    mv $DESTINATION/bedtools2/bin $DESTINATION/bedtools && \
    rm -Rf $DESTINATION/bedtools-2.29.1.tar.gz $DESTINATION/bedtools2

# Install gfold
RUN cd $DESTINATION && \
    tar xjf $DESTINATION/gfold-1.1.4-gsl2.2_2.tar.bz2 && \
    mv $DESTINATION/bin/gfold $DESTINATION/gfold && chmod 755 $DESTINATION/gfold && \
    rm -rf $DESTINATION/gfold-1.1.4-gsl2.2_2.tar.bz2

# Install GSL 2.2.1
RUN cd $DESTINATION && \
    tar xzf gsl-2.2.1.tar.gz && \
    cd gsl-2.2.1 && ./configure && make && make install && \
    ldconfig && \
    cd $DESTINATION && \
    rm -Rf gsl-2.2.1.tar.gz gsl-2.2.1


# Install dateutil 2.8.1
RUN cd $DESTINATION && \
    tar xzf dateutil-2.8.1.tar.gz && \
    cd python-dateutil-2.8.1 && python2 setup.py install  && \
    cd $DESTINATION && \
    rm -Rf dateutil-2.8.1.tar.gz python-dateutil-2.8.1

# Install pandas 0.24.2
RUN cd $DESTINATION && \
    tar xzf pandas-0.24.2.tar.gz && \
    cd pandas-0.24.2 && python2 setup.py install  && \
    cd $DESTINATION && \
    rm -Rf pandas-0.24.2.tar.gz pandas-0.24.2


# Install pybedtools 0.8.0
RUN cd $DESTINATION && \
    tar xzf pydevtools-0.8.0.tar.gz && \
    cd pybedtools-0.8.0 && python2 setup.py cythonize && python2 setup.py install && \
    ldconfig && \
    cd $DESTINATION && \
    rm -Rf pydevtools-0.8.0.tar.gz pydevtools-0.8.0



# Install bigWigToWig and wigToBigWig from bigWig
RUN chmod 755 $DESTINATION/bigWigToWig $DESTINATION/wigToBigWig

# Create shortcut for running GRO-seq
RUN echo '#!/bin/bash \n python2 /bin/GROSeqPL_v3_updated.py "$@"' > $DESTINATION/groseq; chmod 755 $DESTINATION/groseq

# Install HOMER. I update configureHomer.pl because every time a package is downloaded it is overwriten
RUN mkdir $DESTINATION/homer && \
    cd $DESTINATION/homer && \
    http_proxy="${http_proxy}" https_proxy="${https_proxy}" wget http://homer.ucsd.edu/homer/configureHomer.pl && \
    perl ./configureHomer.pl -install homer && \
    rm -Rf $DESTINATION/homer/cpp

# Remove unused packages (TODO: can probably remove more)
RUN apt-get remove -y \
 libcairo2-dev \
 ca-certificates zip unzip git wget \
 tcl-dev tcl8.6-dev tk-dev tk8.6-dev tcl tcl8.6 tk tk8.6 \
 zlib1g-dev libbz2-dev liblzma-dev libbz2-dev \
 libjbig-dev libjpeg-dev libjpeg62-turbo-dev libpng-dev \
 libxml2-dev libxrender-dev libxslt1-dev libxss-dev \
 autotools-dev gcc g++ libgcc-8-dev libstdc++-8-dev libc-dev-bin libc6-dev libyaml-dev \
 libmount-dev dpkg-dev \
 x11proto-core-dev x11proto-dev x11proto-scrnsaver-dev x11proto-xext-dev xtrans-dev libx11-dev && \
 apt-get autoremove -y && apt-get clean

# Copy pipeline source files. All files are mentioneded separately so that docker build is rerun when source code is changing
COPY src/Add.col.py  $DESTINATION/Add.col.py
COPY src/bam2bw.py  $DESTINATION/bam2bw.py
COPY src/bam2wig.py  $DESTINATION/bam2wig.py
COPY src/COVT.gmean_cal.py  $DESTINATION/COVT.gmean_cal.py
COPY src/GROSeqPL_v3_updated.py  $DESTINATION/GROSeqPL_v3_updated.py
COPY src/misc.py  $DESTINATION/misc.py

# Copy preprocessing script into image
COPY preprocess/download.py  $DESTINATION/download.py
RUN echo '#!'"/bin/bash \n python2 $DESTINATION/download.py "'"$@"' > $DESTINATION/download; chmod 755 $DESTINATION/download

# Copy extract_rpkm script into image
COPY src/extract_rpkm.py  $DESTINATION/extract_rpkm.py
RUN echo '#!'"/bin/bash \n python2 $DESTINATION/extract_rpkm.py "'"$@"' > $DESTINATION/extract-rpkm; chmod 755 $DESTINATION/extract-rpkm

# Copy LSF script into image
COPY preprocess/lsf.py  $DESTINATION/lsf.py
RUN echo '#!'"/bin/bash \n python2 $DESTINATION/lsf.py "'"$@"'> $DESTINATION/lsf; chmod 755 $DESTINATION/lsf


# Copy LSF script into image
COPY preprocess/gff_longest_transcript.py  $DESTINATION/gff_longest_transcript.py
RUN echo '#!'"/bin/bash \n python2 /$DESTINATION/gff_longest_transcript.py "'"$@"' > $DESTINATION/longest-transcript; chmod 755 $DESTINATION/longest-transcript

ENV PATH="/bin/bedtools:${PATH}"

WORKDIR /mount

ENTRYPOINT ["groseq"]
