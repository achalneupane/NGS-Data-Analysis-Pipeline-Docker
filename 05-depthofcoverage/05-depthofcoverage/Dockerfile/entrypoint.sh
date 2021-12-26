#!/bin/bash

VERSION="0.1.0"

# Entrypoint script for Docker container achalneupane/depthofcoverage:${VERSION}

#### Required parameters (must provide or container will quit)
#   BASE        location of parent directories for input_files, running_jobs
#   FULLSM      FULLSM (SM^DNA^PR 3-part ID) for the individual whose .bam need to be marked / merged
# 	BAMFILE		*.aln.srt.isec-${RUN_TYPE}.markd.bam input file to check base quality
# 	OUTFILE		*.${RUN_TYPE}_recal.table
# 	REF			Provide link to reference location
# 	RUN_TYPE		covered or padded exome
# 	GATK_THREADS	4 
#	COVERED_BED
#	PADDED_BED


### Optional parameters
# SHELLDROP 	Drop to shell isntead of running anything (used with docker)
# TIMING 		Do /usr/bin/time -v timing of steps

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
if ! [ -r "${BAMFILE}" ]; then
    echo "ERROR, cannot read ${BAMFILE}"
    echo "Please make sure you provide the full path"
    quit "Input File"
fi
if [ -z "${RUN_TYPE}" ]; then
    echo "ERROR, must provide run type (genome, paddedexome, exome) in variable \${RUN_TYPE}"
    # TODO: Should validate RUN_TYPE against genome, exome, paddedexome
    quit "Job Config"
fi

if [ ! -z "${TIMING}" ]; then TIMING=(/usr/bin/time -v); fi

# Set locations based on ${BASE}, which is set via env variable
# Maybe consider dropping this requirement and skipping NM validation
if [ -z "${BASE}" ]; then echo "Error, BASE not specified"; quit "Job Config"; fi
echo "Running on system ${SYSTEM:=UNDEFINED}, BASE set to ${BASE}"

if [ "${RUN_TYPE}" = "genome" ]; then
    GATK_PARAMS=""
elif [ "${RUN_TYPE}" = "exome" ]; then
    GATK_PARAMS="-L ${COVERED_BED}"
elif [ "${RUN_TYPE}" = "paddedexome" ]; then
    GATK_PARAMS="-L ${PADDED_BED}"
else
    echo "Error, invalid RUN_TYPE specified"
    quit "Job Config"
fi

JAVAOPTS="-Xms2g -Xmx${MEM}g -XX:+UseSerialGC -Dpicard.useLegacyParser=false"

# DepthOfCoverage
CUR_STEP="DepthOfCoverage"
start=$(${DATE}); echo "[$(display_date ${start})] ${CUR_STEP} starting"
if [ "${RUN_TYPE}" = "genome" ]; then
  "${TIMING[@]}" java ${JAVAOPTS} -jar "${GATK360}" \
    -T "${CUR_STEP}" \
    -R "${REF}" \
    -nt ${GATK_THREADS} \
    -ct 10 -ct 15 -ct 20 -ct 30 -ct 40 -ct 50 -ct 60 -ct 70 -ct 80 -ct 90 -ct 100 \
    --omitDepthOutputAtEachBase \
    --omitIntervalStatistics \
    --omitLocusTable \
    -I "${BAMFILE}" \
    -o "${BAMFILE%.bam}.genome-coverage"
else
  "${TIMING[@]}" java ${JAVAOPTS} -jar "${GATK360}" \
    -T "${CUR_STEP}" \
    -R "${REF}" \
    -nt ${GATK_THREADS} \
    -ct 10 -ct 15 -ct 20 -ct 30 -ct 40 -ct 50 -ct 60 -ct 70 -ct 80 -ct 90 -ct 100 \
    --omitDepthOutputAtEachBase \
    --omitIntervalStatistics \
    --omitLocusTable \
    -L "${COVERED_BED}" \
    -I "${BAMFILE}" \
    -o "${BAMFILE%.bam}.exome-coverage"
fi
exitcode=$?
end=$(${DATE}); echo "[$(display_date ${end})] ${CUR_STEP} finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"

exit ${exitcode}