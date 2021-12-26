#!/bin/bash

VERSION="0.1.0"
# Entrypoint script for docker container achalneupane/cram2bam:${VERSION}

# Parameters must be defined so the script knows which file(s) to process

# Required parameters (must provide or container will quit)
#   BASE       location of reference files (needed for picard NM validation)
#   RAWDIR     location of raw files
#   CRAMFILE   .cram input file to convert to .bam
#   RUN_TYPE   genome, paddedexome, or exome -- the type of intersect to perform
# # This cram to bam operation requires a specific GR38 reference, we need touse the same one MGI used -- REF="GRCh38_MGI_20190813/all_sequences.fa"
#	REF="$(readlink -f "${BASE}/Genome_Ref/GRCh38/GRCh38_MGI_20190813/all_sequences.fa")"


# Optional parameters
#   OUTDIR     directory to place completed .cram file; defaults to ${BASE}/variant_calling/completed_jobs
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

if [ -z "${CRAMFILE}" ]; then
    echo "ERROR, must provide .cram file in variable \${CRAMFILE}"
    quit "Input File"
fi
if ! [ -r "${CRAMFILE}" ]; then
    echo "ERROR, cannot read ${CRAMFILE}"
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
if [ -z "${OUTDIR}" ]; then OUTDIR="${CRAMFILE%/*}/"; echo "WARNING, no OUTDIR was given so it has been set as same dir as CRAMFILES"; fi

CRAMNAME="${CRAMFILE##*/}"
if [ -z "${BAMFILE}" ]; then
  BAMFILENAME="${CRAMNAME%.final.cram}.bam"
else
  if ! [[ "${BAMFILE}" =~ \.bam$ ]]; then
    echo "ERROR, output BAMFILE must end in .bam"
    quit "Job Config";
  fi
  BAMFILENAME="${BAMFILE##*/}"
fi

if [ "${RUN_TYPE}" = "genome" ]; then
    # cp to preserve the ${FULLSM_RGID}.aln.srt.bam, which might be reused for a WGS analysis later
        CUR_STEP="cramtobam"
		start=$(${DATE}); echo "[$(display_date ${start})] ${CUR_STEP} on a ${RUN_TYPE} file"; echo "sample being processed is ${CRAMFILE}"
        "${TIMING[@]}" samtools view -b -T "${REF}" -o "${OUTDIR}/${BAMFILENAME%.*}.${RUN_TYPE}.bam" "${CRAMFILE}"
        exitcode="$?"
        end=$(${DATE}); echo "[$(display_date ${end})] cramtobam for ${RUN_TYPE} finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"
else
        # Convert cram to bam
        CUR_STEP="cramtobam"
		start=$(${DATE}); echo "[$(display_date ${start})] ${CUR_STEP} starting on a ${RUN_TYPE} file"; echo "sample being processed is ${CRAMFILE}"
        "${TIMING[@]}" samtools view -b -T "${REF}" -o "${OUTDIR}/${BAMFILENAME%.*}.${RUN_TYPE}.bam" "${CRAMFILE}"
        exitcode="$?"
        end=$(${DATE}); echo "[$(display_date ${end})] ${CUR_STEP} for ${RUN_TYPE} finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"
fi

exit ${exitcode}
