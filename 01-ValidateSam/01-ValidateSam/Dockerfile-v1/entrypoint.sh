#!/bin/bash

VERSION="0.0.1"
# Entrypoint script for docker container achalneupane/validatesam:${VERSION}

# Parameters must be defined so the script knows which file(s) to process

# Required parameters (must provide or container will quit)
#   BASE       location of reference files (needed for picard NM validation)
#   BAMFILE    .sam or .bam input file to validate

# Optional parameters
#   SHELLDROP   Drop to shell instead of running anything (used with docker)
#   MEM         Memory Limit in GB (e.g. 32)
#   TIMING      Do /usr/bin/time -v timing of steps

DATE="/bin/date +%s"
display_date () {
    /bin/date -d @${1} +%Y%m%d_%H%M%S
}

date_diff () {
    earlier=${1}
    later=${2}
    diff=$((${later}-${earlier}))
    if [ ${diff} -gt 86400 ]; then
        date -u -d @$((${diff}-86400)) +"%jd%-Hh%-Mm%-Ss"
    else
        date -u -d @${diff} +"%-Hh%-Mm%-Ss"
    fi
}

quit () {
  echo "[$(display_date $(${DATE}))] Run failed at ${1} step, exit code: 1"
  # /bin/mail -r "docker@${HOSTNAME}" -s "${SYSTEM}: ${BAMFILE} FAIL at ${1} (${PBS_JOBID}${SLURM_JOBID}${LSB_JOBID})" "${EMAIL}" < /dev/null > /dev/null
  if [ "${SHELLDROP}" -eq 1 ]; then echo "Dropping to shell"; exec "/bin/bash"; else exit 1; fi
}

# trap "{ echo \"[$(display_date $(${DATE}))] Terminated by SIGTERM \"; quit \"${CUR_STEP}\"; }" SIGTERM
trap "echo \"[$(display_date $(${DATE}))] Terminated by SIGTERM \" && sleep 10s" SIGTERM

# Option for usage in docker
if [ "${SHELLDROP:=0}" -eq 1 ]; then echo "Dropping to shell"; exec "/bin/bash"; fi

if [ -z "${MEM}" ]; then echo "WARNING, memory limit (in GB) not provided in variable \${MEM}; defaulting to 4G"; MEM=4; fi

if [ -z "${BAMFILE}" ]; then
    echo "ERROR, must provide .sam or .bam file in variable \${BAMFILE}"
    quit "Input File"
fi
if [ ! -r "${BAMFILE}" ]; then
    echo "ERROR, cannot read ${BAMFILE}"
    echo "Please make sure you provide the full path"
    quit "Input File"
fi

# Set locations based on ${BASE}, which is set via env variable
# Maybe consider dropping this requirement and skipping NM validation
if [ -z "${BASE}" ]; then echo "Error, BASE not specified"; quit "Job Config"; fi
echo "Running on system ${SYSTEM:=UNDEFINED}, BASE set to ${BASE}"

if [ ! -z "${TIMING}" ]; then TIMING=(/usr/bin/time -v); fi

# Just in case of symlinks, which cause problems for some programs
REF="$(readlink -f "${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta")"

JAVAOPTS="-Xms2g -Xmx${MEM}g -XX:+UseSerialGC -Dpicard.useLegacyParser=false"

# ValidateSamFile
start=$(${DATE}); echo "[$(display_date ${start})] ValidateSamFile starting"
CUR_STEP="ValidateSam"
"${TIMING[@]}" /usr/bin/java ${JAVAOPTS} -jar "${PICARD}" ValidateSamFile -R "${REF}" -I "${BAMFILE}" "$@"
exitcode=$?
end=$(${DATE}); echo "[$(display_date ${end})] ValidateSamFile finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"

exit ${exitcode}