#!/bin/bash

VERSION="1.0"
# Entrypoint script for docker container achalneupane/bam2fq:${VERSION} 

# Parameters must be defined so the script knows which file(s) to process


# Required parameters (must provide or container will quit)
#   BASE       location of reference files (needed for picard NM validation)
#   BAMFILE    .sam or .bam input file to revert
#   FULLSM     FULLSM (SM^DNA^PR 3-part ID) for the individual whose .bam need to be reverted (needed for naming)

# Optional parameters
#   OUT_DIR       Directory to place reverted .fastq files; defaults to ${BAMFILE%/*}
#   SHELLDROP     Drop to shell instead of running anything (used with docker)
#   MEM           Memory Limit in GB (e.g. 32)
#   TIMING        Do /usr/bin/time -v timing of steps
#   CREATE_RGFILE Create .rgfile for each .fastq after reverting
#   REMOVE_INPUT  Delete the input .bam upon successful completion

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

if [ -z "${FULLSM}" ]; then echo "ERROR, must provide 3-part ID in variable \${FULLSM}"; quit "Job Config"; fi
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

# Remove files as you go. Set to 0 for testing. Only set if CLEANUP has not already been declared at the command line.
if [ -z "${REMOVE_INPUT}" ]; then REMOVE_INPUT=0; fi

# if [ -z "${TMP_DIR}" ]; then TMP_DIR="${BAMFILE%/*}"; fi
# mkdir -p "${TMP_DIR}" || { echo "Error, cannot create ${TMP_DIR}"; quit "Setup TMP_DIR"; }

# TMP_DIR="/tmp/${LSB_JOBID}.tmpdir"
# TMP_DIR=${TMP_DIR}
# mkdir -p "/tmp/${LSB_JOBID}.tmpdir"
echo "Using TMP_DIR: ${TMP_DIR}"


# Set locations based on ${BASE}, which is set via env variable
# Maybe consider dropping this requirement and skipping NM validation
if [ -z "${BASE}" ]; then echo "Error, BASE not specified"; quit "Job Config"; fi
echo "Running on system ${SYSTEM:=UNDEFINED}, BASE set to ${BASE}"

if [ ! -z "${TIMING}" ]; then TIMING=(/usr/bin/time -v); fi




if [ -z "${DNA}" ]; then
    echo "ERROR, must provide DNA in variable \${DNA}"
    exit 1
fi

### Create a test bam-little bam
## samtools view -@ ${THREADS} -h ${BAMFILE} | head -n 1000000 | samtools view -bS - > "${OUT_DIR}/little.bam"

IFS=$'\n' RGS=($(samtools view -@ ${THREADS} -h ${BAMFILE} | head -n 10000000 | grep ^HS2000 | cut -d$'\t' -f1| cut -d: -f1,2 | sort -V | uniq | grep ^HS2000))

echo "Readgroups are ${RGS[@]}"
unset IFS


args=(tee)
for RG in ${RGS[@]}; do
  args+=(\>\(grep -A3 --no-group-separator \"^@${RG/^/:}:\" \| gzip \> ${OUT_DIR}/${SM}^${DNA}^${PR}.${RG/:/.}.fq.gz\))
done
args+=(\>/dev/null)

JAVA="/usr/local/openjdk-8/bin/java"
# JAVAOPTS="-Xms2g -Xmx${MEM}g -XX:+UseSerialGC -Dpicard.useLegacyParser=false"
JAVAOPTS="-Xms4g -Xmx${MEM}g -XX:ParallelGCThreads=${THREADS} -Djava.io.tmpdir=${TMP_DIR}"


CUR_STEP="RevertSam"
start=$(${DATE}); echo "[$(display_date ${start})] ${CUR_STEP} starting"
"${TIMING[@]}" ${JAVA} ${JAVAOPTS} -jar "${PICARD}" \
  "${CUR_STEP}" \
  I="${BAMFILE}" \
  O=/dev/stdout \
  SORT_ORDER=queryname \
  COMPRESSION_LEVEL=0 \
  VALIDATION_STRINGENCY=SILENT \
  TMP_DIR=${TMP_DIR} \
  | ${JAVA} ${JAVAOPTS} -jar "${PICARD}" \
  SamToFastq \
  I=/dev/stdin \
  FASTQ=/dev/stdout \
  INTERLEAVE=TRUE \
  VALIDATION_STRINGENCY=SILENT \
  TMP_DIR=${TMP_DIR} | eval ${args[@]}
arr=(${PIPESTATUS[@]}); exitcode=0; for i in ${arr[@]}; do ((exitcode+=i)); done
exitcode=$?
end=$(${DATE}); echo "[$(display_date ${end})] ${CUR_STEP} finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"

# TMP_DIR="${DEST}/ADNI_WGS/${FULLSM}"
# /tmp/${LSB_JOBID}.tmpdir

# Check reads counts
CUR_STEP="Doing Sanity check.."
start=$(${DATE}); echo "[$(display_date ${start})] ${CUR_STEP} starting"
BAMLINES=$(samtools idxstats "${BAMFILE}" -@ ${THREADS} | awk '{s+=$3+$4} END {print s*4}')
FQLINES=$(zcat ${OUT_DIR}/${FULLSM}*.fq.gz | wc -l)
if [[ ! -z ${DEBUG} ]]; then
echo "${BAMLINES} lines in .bam"
echo "${FQLINES} lines in all .fastq.gz files"
if [[ $(echo "scale=2;${FQLINES}/${BAMLINES} > 0.90" | bc) -eq 0 ]]; then
echo "Warning, .fastq.gz files contain less than 90% of the number of reads of .bam file"
fi
fi
exitcode=$?
end=$(${DATE}); echo "[$(display_date ${end})] ${CUR_STEP} finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"


# ## Create RGFILE
CUR_STEP="Create RGFILE"
start=$(${DATE}); echo "[$(display_date ${start})] ${CUR_STEP} starting"
if [ ! -z "${CREATE_RGFILE}" ]; then
SM="$(echo "${FULLSM}" | cut -d^ -f1)"
DNA="$(echo "${FULLSM}" | cut -d^ -f2)"
PR="$(echo "${FULLSM}" | cut -d^ -f3)"
echo "Creating ${#RGS[@]} .rgfiles for newly created .fastq.gz files" 
for RG in ${RGS[@]}; do
RGID="${RG}"
RGID_NEW="$(echo ${RG/:/.})"
RGPU="${RG}"
RGLB=${DNA}
echo "@RG\tID:${RGID}\tDS:${SM}^${DNA}^${PR}\tLB:${RGLB}\tPL:illumina\tPU:${RGPU}\tSM:${SM}" > "${OUT_DIR}/${FULLSM}.${RGID_NEW}.rgfile"
done
fi
exitcode=$?
end=$(${DATE}); echo "[$(display_date ${end})] ${CUR_STEP} finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"

## Cleaning up files
CUR_STEP="Cleaning up files"
start=$(${DATE}); echo "[$(display_date ${start})] ${CUR_STEP} starting"
if [ ${exitcode} -eq 0 ] && [ ${REMOVE_INPUT} -eq 1 ]; then
  rm -fv "${BAMFILE}" "${BAMFILE%.bam}.bai" "${BAMFILE}.bai" 2>/dev/null
  echo -e "REMOVE_INPUT was set to ${REMOVE_INPUT} then BAMFILE is removed"
  else
  echo -e "REMOVE_INPUT was set to ${REMOVE_INPUT} then BAMFILE is NOT removed"
fi
exitcode=$?
end=$(${DATE}); echo "[$(display_date ${end})] ${CUR_STEP} finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"

exit ${exitcode}
