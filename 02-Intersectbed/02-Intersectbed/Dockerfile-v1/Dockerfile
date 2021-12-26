FROM phusion/baseimage:0.11 as builder
LABEL maintainer="Achal Neupane <achal@wustl.edu>"

LABEL org.label-schema.schema-version="1.0"
# LABEL org.label-schema.build-date=$BUILD_DATE
LABEL org.label-schema.name="achalneupane/intersectbed"
LABEL org.label-schema.description="Image for bedtools"


ENV BEDTOOLS_VERSION 2.27.1

WORKDIR /tmp

RUN apt-get update -y \
  && apt-get install --no-install-recommends -y \
       wget \
       make \
       g++ \
       libz-dev \
       python \
  && wget -q https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz \
  && tar xzvf bedtools-${BEDTOOLS_VERSION}.tar.gz \
  && cd /tmp/bedtools${BEDTOOLS_VERSION%%.*}/ \
  && make

FROM phusion/baseimage:0.11


RUN apt-get update -y \
  && apt-get install --no-install-recommends -y \
       libnss-sss \
       time \
       tzdata \
  && apt-get clean all \
  && rm -rfv /var/lib/apt/lists/* /tmp/* /var/tmp/* 

ENV TZ "America/Chicago"
COPY ./Capture_*.GRCh38.bed /usr/local/lib/
ENV COVERED_BED /usr/local/lib/Capture_Covered.GRCh38.bed
ENV PADDED_BED /usr/local/lib/Capture_Padded.GRCh38.bed
ENV BEDTOOLS_VERSION 2.27.1

COPY --from=builder /tmp/bedtools*/bin/bedtools /usr/bin/
COPY ./entrypoint.sh /usr/local/bin/

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]