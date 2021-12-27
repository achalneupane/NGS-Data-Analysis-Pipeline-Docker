#!/bin/bash

VERSION="0.1.0"
# Entrypoint script for docker container achalneupane/markduplicates:${VERSION}

# Parameters must be defined so the script knows which file(s) to process

# Required parameters (must provide or container will quit)
#   BASE        location of parent directories for input_files, running_jobs
#   FULLSM      FULLSM (SM^DNA^PR 3-part ID) for the individual whose .bam need to be marked / merged
#   RUN_TYPE    required for now, may be dropped in the future
#   $@          full paths to .bam files to mark / merge (command-line parameter), at least 1 required

# Optional parameters
#   WORKDIR     parent directory to create JOB_TMP directory to hold all files. Often uses ${PBS_JOBID} or ${SLURM_JOBID} or ${LSB_JOBID}. Required on fenix
#   OUT_DIR     directory to hold output file from this script; defaults to ${BASE}/variant_calling/running_jobs/${FULLSM}
#   MEM         Memory Limit in GB (e.g. 32), defaults to 4
#   SHELLDROP   Drop to shell instead of running anything (used with docker)
#   TIMING      Do /usr/bin/time -v timing of steps

# Possibly added in the future to make the image more flexible (any full path inputs, any full path output)
#   OUTFILE     name of output marked / merged .bam; defaults to ${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.bam

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
  echo "Leftover files may be in ${JOB_TMP}"
  if [ -d "${JOB_TMP}" ]; then rmdir -v "${JOB_TMP}"; fi  # this might fail if files are leftover, but is handy if the job fails early
  # /bin/mail -r "docker@${HOSTNAME}" -s "${SYSTEM}: ${BAMFILE} FAIL at ${1} (${PBS_JOBID}${SLURM_JOBID}${LSB_JOBID})" "${EMAIL}" < /dev/null > /dev/null
  if [ "${SHELLDROP}" -eq 1 ]; then echo "Dropping to shell"; exec "/bin/bash"; else exit 1; fi
}

# trap "{ echo \"[$(display_date $(${DATE}))] Terminated by SIGTERM \"; quit \"${CUR_STEP}\"; }" SIGTERM
trap "echo \"[$(display_date $(${DATE}))] Terminated by SIGTERM \" && sleep 10s" SIGTERM

# Option for usage in docker
if [ "${SHELLDROP:=0}" -eq 1 ]; then echo "Dropping to shell"; exec "/bin/bash"; fi

if [ -z "${FULLSM}" ]; then echo "ERROR, must provide 3-part ID in variable \${FULLSM}"; quit "Job Config"; fi
if [ -z "${MEM}" ]; then echo "WARNING, memory limit (in GB) not provided in variable \${MEM}; defaulting to 4G"; MEM=4; fi

# Set locations based on ${BASE}, which is set via env variable
if [ -z "${BASE}" ]; then echo "Error, BASE not specified"; quit "Job Config"; fi
echo "Running on system ${SYSTEM:=UNDEFINED}, BASE set to ${BASE}"
if [ -z "${WORKDIR}" ]; then 
  if [ "${SYSTEM}" = "MGI" ]; then
    WORKDIR="${BASE}/tmp/${LSB_JOBID}.tmpdir"
    mkdir -p "${WORKDIR}" || { echo "Error, cannot create ${WORKDIR}"; quit "Setup WORKDIR"; }
  else 
    echo "Error, WORKDIR not specified, refusing to guess an appropriate location on this unknown filesystem"
    quit "Job Config"
  fi
fi

if [ $# -lt 1 ]; then echo "ERROR, must provide full paths to .bam files as command-line parameters"; quit "Job Config"; fi
INPUT_BAM_LIST=()
MARKD_ARGS=()
for BAMFILE in "${@}"; do
  if [ ! -r "${BAMFILE}" ]; then echo -e "ERROR, cannot read ${BAMFILE}\nPlease make sure you provide the full path"; quit "Input File"; fi
  INPUT_BAM_LIST+=("${BAMFILE}")
  MARKD_ARGS+=("-I" "${BAMFILE}")
done

if [ ! -z "${TIMING}" ]; then TIMING=(/usr/bin/time -v); fi

if [ -z "${OUT_DIR}" ]; then OUT_DIR="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}"; fi
mkdir -p "${OUT_DIR}" || { echo "Error, cannot create ${OUT_DIR}"; quit "Setup OUT_DIR"; }

# Create a directory to contain all of the files generated while running the script
JOB_TMP="$(mktemp -d -p "${WORKDIR}/" "${LSB_JOBID:-0}.XXXXXXXXXXXXXXXX")" \
  || { echo "Error, cannot create ${JOB_TMP}"; quit "Setup JOB_TMP"; }

JAVAOPTS="-Xms2g -Xmx${MEM}g -XX:+UseSerialGC -Dpicard.useLegacyParser=false"

# MarkDuplicates
start=$(${DATE}); echo "[$(display_date ${start})] MarkDuplicates starting"
CUR_STEP="MarkDuplicates"
# Would be better to create -O name from consensus commonalities in INPUT_BAM_LIST
"${TIMING[@]}" java ${JAVAOPTS} -jar "${PICARD}" \
  MarkDuplicates \
  "${MARKD_ARGS[@]}" \
  -O "${OUT_DIR}/${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.bam" \
  -METRICS_FILE "${OUT_DIR}/${FULLSM}_metrics.txt" \
  -QUIET true \
  -MAX_RECORDS_IN_RAM $((MEM*250000)) \
  -TMP_DIR "${JOB_TMP}" \
  -ASSUME_SORTED TRUE \
  -CREATE_INDEX TRUE
exitcode=$?
end=$(${DATE}); echo "[$(display_date ${end})] MarkDuplicates finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"

rmdir -v "${JOB_TMP}"  # this might fail if files are leftover, but is handy if the job fails early

echo "Contents of ${OUT_DIR}:"
ls -ltr "${OUT_DIR}" | sed '1d'
if [ ${exitcode} -eq 0 ]; then  # assume CLEANUP=1
    rm -fv "${INPUT_BAM_LIST[@]}"
    # rm -fv "${BASE}/input_files/${FULLSM}/"*.rgfile  # Should this be input_files to replicate same structure? bam_to_fastq would have to put them there
    # rm -f "${BASE}/input_files/${FULLSM}/".*.startrun "${BASE}/input_files/${FULLSM}/".*.runfailed
    # rmdir -v "${BASE}/input_files/${FULLSM}"
else
    echo "Something failed"
    # /bin/mail -s "${SYSTEM}: ${FULLSM} FAIL at MarkDuplicates (${PBS_JOBID}${SLURM_JOBID})" "${EMAIL}" < /dev/null > /dev/null
fi

exit ${exitcode}