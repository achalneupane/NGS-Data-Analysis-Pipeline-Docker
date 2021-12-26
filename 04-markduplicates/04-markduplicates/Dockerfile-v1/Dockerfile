FROM openjdk:8-jre
LABEL maintainer="Achal Neupane <achal@wustl.edu>"

LABEL org.label-schema.schema-version="1.0"
# LABEL org.label-schema.build-date=$BUILD_DATE
LABEL org.label-schema.name="achalneupane/markduplicates"
LABEL org.label-schema.description="Image for picard"

ENV PICARD_VERSION 2.18.21

WORKDIR /tmp

RUN apt-get update -y \
  && apt-get install --no-install-recommends -y \
       libnss-sss \
       time \
  && wget -q -O /usr/bin/picard.jar https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar \
  && ln -sf /usr/share/zoneinfo/America/Chicago /etc/localtime \
  && echo "America/Chicago" > /etc/timezone \
  && dpkg-reconfigure --frontend noninteractive tzdata \
  && apt-get clean all \
  && rm -rfv /var/lib/apt/lists/* /tmp/* /var/tmp/* 

ENV PICARD /usr/bin/picard.jar

COPY ./entrypoint.sh /usr/local/bin/


ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]