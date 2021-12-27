FROM openjdk:8-jre

# wget https://github.com/statgen/verifyBamID/releases/download/v1.1.3/verifyBamID

LABEL maintainer="Achal Neupane <achal@wustl.edu>"

LABEL org.label-schema.name="achalneupane/verifybamid"
LABEL org.label-schema.description="Image for verifybamid"


RUN apt-get update
COPY verifyBamID /usr/bin/

ENV PATH=/usr/bin/:$PATH

## To run
# docker run -ti achalneupane/verifybamid verifyBamID --help
