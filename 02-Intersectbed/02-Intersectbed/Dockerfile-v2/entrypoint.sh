#!/bin/bash

VERSION="2.0"
# Entrypoint script for docker container achalneupane/intersectbedref:${VERSION} 
# modified from version 1.0 so that thsi ones asks the user to provide refenrece and intermediate files

# Parameters must be defined so the script knows which file(s) to process

# Required parameters (must provide or container will quit)
#   BAMFILE    .sam or .bam input file to intersect
#   RUN_TYPE   genome, paddedexome, or exome -- the type of intersect to perform
#   BEDFILE		provide bed file to intersect according to previosuly defined RUN_TYPE


# Optional parameters
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
  # /bin/mail -r "docker@${HOSTNAME}" -s "${SYSTEM}: ${BAMFILE} FAIL at ${1} (${PBS_JOBID}${SLURM_JOBID}${LSB_JOBID})" "${EMAIL}" < /dev/null > /dev/null
  if [ "${SHELLDROP}" -eq 1 ]; then echo "Dropping to shell"; exec "/bin/bash"; else exit 1; fi
}

# trap "{ echo \"[$(display_date $(${DATE}))] Terminated by SIGTERM \"; quit \"${CUR_STEP}\"; }" SIGTERM
trap "echo \"[$(display_date $(${DATE}))] Terminated by SIGTERM \" && sleep 10s" SIGTERM

# Option for usage in docker
if [ "${SHELLDROP:=0}" -eq 1 ]; then echo "Dropping to shell"; exec "/bin/bash"; fi

if [ -z "${BAMFILE}" ]; then
    echo "ERROR, must provide .bam file in variable \${BAMFILE}"
    quit "Input File"
fi
if [ ! -r "${BAMFILE}" ]; then
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

if [ "${RUN_TYPE}" = "genome" ]; then
    # cp to preserve the ${FULLSM_RGID}.aln.srt.bam, which might be reused for a WGS analysis later
    cp -a "${BAMFILE}" "${BAMFILE%.bam}.isec-${RUN_TYPE}.bam"
else
    # INPUT: _sorted.bam; OUTPUT: _${RUN_TYPE}_sorted.bam
    start=$(${DATE}); echo "[$(display_date ${start})] intersectBed starting"
    CUR_STEP="intersectBed"
    if [ "${RUN_TYPE}" = "exome" ]; then
      BEDFILE="${COVERED_BED}"
    elif [ "${RUN_TYPE}" = "paddedexome" ]; then
      BEDFILE="${PADDED_BED}"
    fi
    "${TIMING[@]}" bedtools intersect -u -a "${BAMFILE}" -b "${BEDFILE}" > "${BAMFILE%.bam}.isec-${RUN_TYPE}.bam"
    exitcode=$?
    end=$(${DATE}); echo "[$(display_date ${end})] intersectBed finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"
fi

exit ${exitcode}