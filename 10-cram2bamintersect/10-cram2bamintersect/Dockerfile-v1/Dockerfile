FROM phusion/baseimage:0.11 as builder

LABEL maintainer="Achal Neupane <achal@wustl.edu>"
LABEL org.label-schema.schema-version="1.0"
# LABEL org.label-schema.build-date=$BUILD_DATE
LABEL org.label-schema.name="achalneupane/cram2bam4intersect"
LABEL org.label-schema.description="Image for samtools and bedtools"

### Define versions of programs to install
ENV SAMTOOLS_VERSION 1.9
ENV BEDTOOLS_VERSION 2.27.1

#set the working directory to be used when the docker gets run
WORKDIR /tmp

# Do a few updates of the base system
RUN apt-get update -y \
  && apt-get install --no-install-recommends -y \
       wget \
       make \
       g++ \
       ncurses-dev \
       libz-dev \
       libbz2-dev \
       liblzma-dev \
	   libnss-sss \
	   python \
       time \
       tzdata \
  && wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
  && tar xjvf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
  && cd /tmp/samtools-${SAMTOOLS_VERSION}/ \
  && make \
  && cp -av /tmp/samtools-${SAMTOOLS_VERSION}/samtools /usr/bin/ \
  && cd /tmp \
  && wget -q https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz \
  && tar xzvf bedtools-${BEDTOOLS_VERSION}.tar.gz \
  && cd /tmp/bedtools${BEDTOOLS_VERSION%%.*}/ \
  && make \ 
  && cp -av /tmp/bedtools${BEDTOOLS_VERSION%%.*}/bin/bedtools /usr/bin/ \
  && apt-get clean all \
  && rm -rfv /var/lib/apt/lists/* /tmp/* /var/tmp/*   


#FROM phusion/baseimage:0.11
ENV SAMTOOLS_VE1SION 1.9
ENV BEDTOOLS_VERSION 2.27.1
ENV TZ "America/Chicago"
COPY ./Capture_*.GRCh38.bed /usr/local/lib/
ENV COVERED_BED /usr/local/lib/Capture_Covered.GRCh38.bed
ENV PADDED_BED /usr/local/lib/Capture_Padded.GRCh38.bed

COPY ./entrypoint.sh /usr/local/bin/

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]





