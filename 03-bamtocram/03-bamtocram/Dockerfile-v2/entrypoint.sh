#!/bin/bash

VERSION="0.2.0"
# Entrypoint script for docker container achalneupane/bamtocram:${VERSION}
# modified fromprevius version so that this one requests user to specify directory where reference genome is placed

# Parameters must be defined so the script knows which file(s) to process

# Required parameters (must provide or container will quit)
#   BASE       location of reference files (needed for picard NM validation)
#   BAMFILE    .sam or .bam input file to convert to .cram
#   REF			location of reference genome file

# Optional parameters
#   CRAMFILE   .cram output file; defaults to ${BAMFILE%.bam}.cram
#   OUT_DIR     directory to place completed .cram file; defaults to ${BASE}/variant_calling/completed_jobs
#   WORKDIR     parent directory to create JOB_TMP directory to hold all files. Often uses ${PBS_JOBID} or ${SLURM_JOBID} or ${LSB_JOBID}. Required on fenix
#   SHELLDROP   Drop to shell instead of running anything (used with docker)
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
  echo "Leftover files may be in ${JOB_TMP}"
  # /bin/mail -r "docker@${HOSTNAME}" -s "${SYSTEM}: ${BAMFILE} FAIL at ${1} (${PBS_JOBID}${SLURM_JOBID}${LSB_JOBID})" "${EMAIL}" < /dev/null > /dev/null
  if [ "${SHELLDROP}" -eq 1 ]; then echo "Dropping to shell"; exec "/bin/bash"; else exit 1; fi
}

# trap "{ echo \"[$(display_date $(${DATE}))] Terminated by SIGTERM \"; quit \"${CUR_STEP}\"; }" SIGTERM
trap "echo \"[$(display_date $(${DATE}))] Terminated by SIGTERM \" && sleep 10s" SIGTERM

# Option for usage in docker
if [ "${SHELLDROP:=0}" -eq 1 ]; then echo "Dropping to shell"; exec "/bin/bash"; fi

if [ -z "${BAMFILE}" ]; then
    echo "ERROR, must provide .sam or .bam file in variable \${BAMFILE}"
    quit "Input File"
fi
if ! [ -r "${BAMFILE}" ]; then
    echo "ERROR, cannot read ${BAMFILE}"
    echo "Please make sure you provide the full path"
    quit "Input File"
fi

# Set locations based on ${BASE}, which is set via env variable
# Maybe consider dropping this requirement and skipping NM validation
if [ -z "${BASE}" ]; then echo "Error, BASE not specified"; quit "Job Config"; fi
echo "Running on system ${SYSTEM:=UNDEFINED}, BASE set to ${BASE}"

if [ -z "${WORKDIR}" ]; then 
  if [ "${SYSTEM}" = "MGI" ]; then
    WORKDIR="${BASE}/tmp/"
  else
    echo "Error, WORKDIR not specified, refusing to guess an appropriate location on this unknown filesystem"
    quit "Job Config"
  fi
fi
mkdir -p "${WORKDIR}" || { echo "Error, cannot create ${WORKDIR}"; quit "Setup WORKDIR"; }

# Create a directory to contain all of the files generated while running the script
JOB_TMP="$(mktemp -d -p "${WORKDIR}/" "${LSB_JOBID:-0}.XXXXXXXXXXXXXXXX")" \
  || { echo "Error, cannot create ${JOB_TMP}"; quit "Setup JOB_TMP"; }

if ! [ -z "${TIMING}" ]; then TIMING=(/usr/bin/time -v); fi
if [ -z "${OUT_DIR}" ]; then OUT_DIR="${BASE}/WXS_Aquilla/03-FINAL/${PR}/${FULLSM}/"; fi
mkdir -p "${OUT_DIR}" || { echo "Error, cannot create ${OUT_DIR}"; quit "Setup OUT_DIR"; }


# Just in case of symlinks, which cause problems for some programs
# REF="$(readlink -f "${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta")"

BAMFILENAME="${BAMFILE##*/}"
if [ -z "${CRAMFILE}" ]; then
  CRAMFILENAME="${BAMFILENAME%.bam}.cram"
else
  if ! [[ "${CRAMFILE}" =~ \.cram$ ]]; then
    echo "ERROR, output CRAMFILE must end in .cram"
    quit "Job Config";
  fi
  CRAMFILENAME="${CRAMFILE##*/}"
fi


# Convert bam to cram
start=$(${DATE}); echo "[$(display_date ${start})] bamtocram starting"
CUR_STEP="bamtocram"
# Output is staged in ${JOB_TMP} and moved to ${OUT_DIR} only once complete (to complement with other automation workflows)
mv -vi "${BAMFILE}" "${JOB_TMP}" \
  && "${TIMING[@]}" samtools view -T "${REF}" -C -o "${JOB_TMP}/${CRAMFILENAME}" "${JOB_TMP}/${BAMFILENAME}"
exitcode="$?"
end=$(${DATE}); echo "[$(display_date ${end})] bamtocram finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"

if [ ${exitcode} -eq 0 ]; then
  mv -vi "${JOB_TMP}/${CRAMFILENAME}" "${OUT_DIR}/" \
    && rm -fv "${JOB_TMP}/${BAMFILENAME}" \
    && rmdir -v "${JOB_TMP}"
else
  echo "ERROR, convertion to .cram failed"
  quit "samtools view"
fi

exit ${exitcode}