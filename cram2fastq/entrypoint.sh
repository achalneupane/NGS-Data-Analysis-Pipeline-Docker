#!/bin/bash

VERSION="1.1"
# Entrypoint script for docker container achalneupane/cram2bam:${VERSION}

# Parameters must be defined so the script knows which file(s) to process

# Required parameters (must provide or container will quit)
#   BASE       location of reference files (needed for picard NM validation)
#   CRAMFILE    .sam or .cram input file to revert
#   FULLSM     FULLSM (SM^DNA^PR 3-part ID) for the individual whose .cram need to be reverted (needed for naming)

# Optional parameters
#   OUT_DIR       Directory to place reverted .fastq files; defaults to ${CRAMFILE%/*}
#   SHELLDROP     Drop to shell instead of running anything (used with docker)
#   MEM           Memory Limit in GB (e.g. 32)
#   TIMING        Do /usr/bin/time -v timing of steps
#   CREATE_RGFILE Create .rgfile for each .fastq after reverting
#   REMOVE_INPUT  Delete the input .cram upon successful completion

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
  # /bin/mail -r "docker@${HOSTNAME}" -s "${SYSTEM}: ${CRAMFILE} FAIL at ${1} (${PBS_JOBID}${SLURM_JOBID}${LSB_JOBID})" "${EMAIL}" < /dev/null > /dev/null
  if [ "${SHELLDROP}" -eq 1 ]; then echo "Dropping to shell"; exec "/bin/bash"; else exit 1; fi
}

# trap "{ echo \"[$(display_date $(${DATE}))] Terminated by SIGTERM \"; quit \"${CUR_STEP}\"; }" SIGTERM
trap "echo \"[$(display_date $(${DATE}))] Terminated by SIGTERM \" && sleep 10s" SIGTERM

# Option for usage in docker
if [ "${SHELLDROP:=0}" -eq 1 ]; then echo "Dropping to shell"; exec "/bin/bash"; fi

if [ -z "${FULLSM}" ]; then echo "ERROR, must provide 3-part ID in variable \${FULLSM}"; quit "Job Config"; fi
if [ -z "${MEM}" ]; then echo "WARNING, memory limit (in GB) not provided in variable \${MEM}; defaulting to 4G"; MEM=4; fi

if [ -z "${CRAMFILE}" ]; then
    echo "ERROR, must provide cram file in variable \${CRAMFILE}"
    quit "Input File"
fi
if [ ! -r "${CRAMFILE}" ]; then
    echo "ERROR, cannot read ${CRAMFILE}"
    echo "Please make sure you provide the full path"
    quit "Input File"
fi


if [ ! -r "${REF}" ]; then
    echo "ERROR, cannot read ${REF}"
    echo "Please make sure you provide the full path"
    quit "Input File"
fi


# Remove files as you go. Set to 0 for testing. Only set if CLEANUP has not already been declared at the command line.
if [ -z "${REMOVE_INPUT}" ]; then REMOVE_INPUT=0; fi

if [ -z "${OUT_DIR}" ]; then OUT_DIR="${CRAMFILE%/*}"; fi
mkdir -p "${OUT_DIR}" || { echo "Error, cannot create ${OUT_DIR}"; quit "Setup OUT_DIR"; }

# Set locations based on ${BASE}, which is set via env variable
# Maybe consider dropping this requirement and skipping NM validation
if [ -z "${BASE}" ]; then echo "Error, BASE not specified"; quit "Job Config"; fi
echo "Running on system ${SYSTEM:=UNDEFINED}, BASE set to ${BASE}"

if [ ! -z "${TIMING}" ]; then TIMING=(/usr/bin/time -v); fi

JAVAOPTS="-Xms2g -Xmx${MEM}g -XX:+UseSerialGC -Dpicard.useLegacyParser=false"

# RevertSam
CUR_STEP="RevertSam"
start=$(${DATE}); echo "[$(display_date ${start})] ${CUR_STEP} starting"
"${TIMING[@]}" /usr/local/openjdk-8/bin/java ${JAVAOPTS} -jar "${PICARD}" \
  "${CUR_STEP}" \
  -I "${CRAMFILE}" \
  -R ${REF} \
  -O /dev/stdout \
  -SORT_ORDER queryname \
  -COMPRESSION_LEVEL 0 \
  -VALIDATION_STRINGENCY SILENT \
  | /usr/local/openjdk-8/bin/java ${JAVAOPTS} -jar "${PICARD}" \
      SamToFastq \
      -I /dev/stdin \
      -R ${REF} \
      -OUTPUT_PER_RG true \
      -RG_TAG ID \
      -OUTPUT_DIR "${OUT_DIR}" \
      -VALIDATION_STRINGENCY SILENT
exitcode=$?
end=$(${DATE}); echo "[$(display_date ${end})] ${CUR_STEP} finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"

# Transfer of files
CUR_STEP="Transfer of files"
start=$(${DATE}); echo "[$(display_date ${start})] ${CUR_STEP} starting"
CRAMLINES=$(samtools idxstats "${CRAMFILE}" | awk '{s+=$3+$4} END {print s*4}')
FQLINES=$(cat ${OUT_DIR}/*.fastq | wc -l)
if [ ! -z ${DEBUG} ]; then
  echo "${CRAMLINES} lines in .cram"
  echo "${FQLINES} lines in all .fastq files"
  if [ $(echo "scale=2;${FQLINES}/${CRAMLINES} > 0.90" | bc) -eq 0 ]; then
    echo "Warning, .fastq files contain less than 90% of the number of reads of .cram file"
  fi
fi
exitcode=$?
end=$(${DATE}); echo "[$(display_date ${end})] ${CUR_STEP} finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"


## Create RGFILE
CUR_STEP="Create RGFILE"
start=$(${DATE}); echo "[$(display_date ${start})] ${CUR_STEP} starting"
if [ ! -z "${CREATE_RGFILE}" ]; then
  SM="$(echo "${FULLSM}" | cut -d^ -f1)"
  DNA="$(echo "${FULLSM}" | cut -d^ -f2)"
  PR="$(echo "${FULLSM}" | cut -d^ -f3)"
  IFS=$'\n' RGS=($(samtools view -H "${CRAMFILE}" | grep "^@RG"))
  echo "Creating ${#RGS[@]} .rgfiles for newly created .fastq files" 
  echo "Moving ${#RGS[@]} .fastq files to ${OUT_DIR}/"
  for RG in ${RGS[@]}; do
    RGID="$(echo ${RG} | grep -oP "(?<=ID:)[^[:space:]]*")"
    RGID_NEW="$(echo ${RGID} | cut -d: -f2- | sed 's/:/^/g')"
    mv -vf "${OUT_DIR}/${RGID//:/_}_1.fastq" "${OUT_DIR}/${FULLSM}.${RGID_NEW}.r1.fastq"
    if [ -f "${OUT_DIR}/${RGID//:/_}_2.fastq" ]; then mv -vf "${OUT_DIR}/${RGID//:/_}_2.fastq" "${OUT_DIR}/${FULLSM}.${RGID_NEW}.r2.fastq"; fi
    RGPU="$(echo ${RG} | grep -oP "(?<=PU:)[^[:space:]]*")"
    RGLB="${SM}.${PR}"
    echo "@RG\tID:${RGID}\tPL:illumina\tPU:${RGPU}\tLB:${RGLB}\tSM:${SM}\tDS:${SM}^${DNA}^${PR}" > "${OUT_DIR}/${FULLSM}.${RGID_NEW}.rgfile"
  done
fi
exitcode=$?
end=$(${DATE}); echo "[$(display_date ${end})] ${CUR_STEP} finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"

## Cleaning up files
CUR_STEP="Cleaning up files"
start=$(${DATE}); echo "[$(display_date ${start})] ${CUR_STEP} starting"
if [ ${exitcode} -eq 0 ] && [ ${REMOVE_INPUT} -eq 1 ]; then
  rm -fv "${CRAMFILE}" "${CRAMFILE%.cram}.crai" "${CRAMFILE}.crai" 2>/dev/null
  echo -e "REMOVE_INPUT was set to ${REMOVE_INPUT} then CRAMFILE is removed"
  else
  echo -e "REMOVE_INPUT was set to ${REMOVE_INPUT} then CRAMFILE is NOT removed"
fi
exitcode=$?
end=$(${DATE}); echo "[$(display_date ${end})] ${CUR_STEP} finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"

exit ${exitcode}