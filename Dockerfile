FROM python:2.7.18
LABEL version="1"
LABEL software="GRO-seq"
MAINTAINER Sergej Andrejev <sandrejev@gmail.com>

ENV DESTINATION=/bin
ENV BOWTIE2_URL=https://github.com/BenLangmead/bowtie2/releases/download/v2.2.9/bowtie2-2.2.9-linux-x86_64.zip
ENV SAMTOOLS_URL=https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2
ENV BEDTOOLS_URL=https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz

# Install required standard libraries
RUN apt-get update
RUN apt-get install make gcc g++ perl \
    ca-certificates ghostscript \
    libncurses5-dev libncursesw5-dev zip unzip git wget \
    zlib1g-dev libbz2-dev liblzma-dev -y;

# Install bowtie2
RUN wget -q $BOWTIE2_URL -O $DESTINATION/bowtie2.zip && \
    unzip $DESTINATION/bowtie2.zip -d $DESTINATION/bowtie2 && \
    rm $DESTINATION/bowtie2.zip

# Install SAMTOOLS
RUN wget -q $SAMTOOLS_URL -O $DESTINATION/samtools_dir.zip;  \
    mkdir $DESTINATION/samtools_dir && tar xjf $DESTINATION/samtools_dir.zip -C $DESTINATION/samtools_dir; mv $DESTINATION/samtools_dir/*/* $DESTINATION/samtools_dir; \
    cd $DESTINATION/samtools_dir/ && make; \
    mv $DESTINATION/samtools_dir/samtools $DESTINATION/samtools; \
    rm -rf $DESTINATION/samtools_dir*;

# Install BEDTOOLS
RUN wget -q $BEDTOOLS_URL -O $DESTINATION/bedtoos_dir.zip;  \
    mkdir $DESTINATION/bedtoos_dir && tar xzf $DESTINATION/bedtoos_dir.zip -C $DESTINATION/bedtoos_dir; mv $DESTINATION/bedtoos_dir/*/* $DESTINATION/bedtoos_dir; \
    cd $DESTINATION/bedtoos_dir/ && make; \
    mv $DESTINATION/bedtoos_dir/bin $DESTINATION/bedtools; \
    rm -rf $DESTINATION/bedtoos_dir;

# Install HOMER
RUN mkdir $DESTINATION/homer && cd $DESTINATION/homer && wget -q http://homer.ucsd.edu/homer/configureHomer.pl && \
    perl configureHomer.pl -install homer mm9;

# Copy GRO-seq files to docker image
COPY src/*  $DESTINATION/groseq/
COPY data/*  ./



# Install BEDTOOLS
# Remove unneeded packages
#    apk del build-base zlib-dev ca-certificates wget

WORKDIR /data

#CMD ["samtools"]
#ENTRYPOINT ["samtools"]