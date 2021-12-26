#!/bin/bash

VERSION="0.1.0"
# Entrypoint script for docker container achalneupane/carm2bam4intersect:${VERSION}

# Parameters must be defined so the script knows which file(s) to process

# Required parameters (must provide or container will quit)
#   BASE       location of reference files (needed for picard NM validation)
#	RAWDIR		location of raw files   
#	CRAMFILE    .cram input file to convert to .bam
#   RUN_TYPE   genome, paddedexome, or exome -- the type of intersect to perform
		

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
JOB_TMP="$(mktemp -d -p "${WORKDIR}/" "${LSB_JOBID:-0}.XXXXXXXXXXXXXXXX")"; echo "The JOB_TMP is ${JOB_TMP}" \
  || { echo "Error, cannot create ${JOB_TMP}"; quit "Setup JOB_TMP"; }

if ! [ -z "${TIMING}" ]; then TIMING=(/usr/bin/time -v); fi
if [ -z "${OUTDIR}" ]; then OUTDIR="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/"; fi
mkdir -p "${OUTDIR}" || { echo "Error, cannot create ${OUTDIR}"; quit "Setup OUTDIR"; }


# This cram to bam operation requires a specific GR38 reference, we need touse the same one MGI used -- REF="GRCh38_MGI_20190813/all_sequences.fa"
REF="$(readlink -f "${BASE}/Genome_Ref/GRCh38/GRCh38_MGI_20190813/all_sequences.fa")"

CRAMFILENAME="${CRAMFILE##*/}"
echo "INFO: CRAMFILENAME provided is ${CRAMFILENAME}"
if [ -z "${BAMFILE}" ]; then
  BAMFILENAME="${CRAMFILENAME%.*.*}"
  echo "$BAMFILENAME"
else
  if ! [[ "${BAMFILE}" =~ \.bam$ ]]; then
    echo "ERROR, output BAMFILE must end in .bam"
    quit "Job Config";
  fi
  BAMFILENAME="${BAMFILE##*/}"
  echo "INFO: BAMFILENAME provided is ${BAMFILENAME}"
fi


if [ "${RUN_TYPE}" = "genome" ]; then
    # cp to preserve the ${FULLSM_RGID}.aln.srt.bam, which might be reused for a WGS analysis later
    # cp -a "${BAMFILE}" "${BAMFILE%.bam}.isec-${RUN_TYPE}.bam"
	start=$(${DATE}); echo "[$(display_date ${start})] cramtobam"
	CUR_STEP="cramtobam"
	
	# Output is staged in ${JOB_TMP} and moved to ${OUTDIR} only once complete (to complement with other automation workflows)
	mv -vi "${CRAMFILE}" "${JOB_TMP}" \
	&& "${TIMING[@]}" samtools view -b -T ${REF} -o "${JOB_TMP}/${BAMFILENAME}" "${JOB_TMP}/${CRAMFILENAME}"
	exitcode="$?"
	end=$(${DATE}); echo "[$(display_date ${end})] cramtobam for ${RUN_TYPE} finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"
	
	if [ ${exitcode} -eq 0 ]; then
	mv -vi "${JOB_TMP}/${BAMFILENAME}" "${OUTDIR}/" \
		&& mv -fv "${JOB_TMP}/${CRAMFILENAME}" "${RAWDIR}/" \
		&& rmdir -v "${JOB_TMP}" \
		&& echo "After processing, CRAMFILE is located at ${RAWDIR} and BAMFILE is located at ${OUTDIR} and ready for Validate Sam"
	else
		echo "ERROR, convertion to .bam failed"
	quit "samtools view"
	fi
else
	# Convert cram to bam & intersect
	start=$(${DATE}); echo "[$(display_date ${start})] cramtobam and intersect starting"
	CUR_STEP="cramtobam4intersect"
    if [ "${RUN_TYPE}" = "exome" ]; then
      BEDFILE="${COVERED_BED}"
    elif [ "${RUN_TYPE}" = "paddedexome" ]; then
      BEDFILE="${PADDED_BED}"
    fi
	
	# Output is staged in ${JOB_TMP} and moved to ${OUTDIR} only once complete (to complement with other automation workflows)
	mv -vi "${CRAMFILE}" "${JOB_TMP}" \
	&& "${TIMING[@]}" samtools view -b -T "${REF}" "${JOB_TMP}/${CRAMFILENAME}" | bedtools intersect -abam stdin -b "${PADDED_BED}" > "${JOB_TMP}/${BAMFILENAME}"
	exitcode="$?"
	end=$(${DATE}); echo "[$(display_date ${end})] cramtobam4intersect for ${RUN_TYPE} finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"

	if [ ${exitcode} -eq 0 ]; then
	mv -vi "${JOB_TMP}/${BAMFILENAME}" "${OUTDIR}/" \
		&& mv -fv "${JOB_TMP}/${CRAMFILENAME}" "${RAWDIR}/" \
		&& echo "After processing, CRAMFILE ${CRAMFILENAME} is located at ${RAWDIR} and BAMFILE ${BAMFILENAME} is located at ${OUTDIR} and ready for Validate Sam" \
		#&& rmdir -v "${JOB_TMP}"
		
	else
		mv -vi "${JOB_TMP}/${CRAMFILENAME}" "${RAWDIR}/" \
		&& echo "ERROR, convertion to .bam and / or intersect failed; CRAMFILE will be placed back to their original directory"
	quit "samtools view"
	fi
fi


exit ${exitcode}



