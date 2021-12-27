#!/bin/bash

VERSION="0.0.1"
# Entrypoint script for docker container achalneupane/samorbam2fastq:${VERSION}

## This version has the modification that the path to the reference must be given

# Parameters must be defined so the script knows which file(s) to process
# Required parameters (must provide or container will quit)
#   BASE       location of reference files (needed for picard NM validation)
#   BAMFILE    .sam or .bam input file to validate
#	FASTQNAME	Name to be given to the FASTQ if different from the BAMFILE
# 	RAWDIR		location to place the FASTQ if differnet from teh BAMfile location
# 	REF			needed for cram to bam conversion

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

# Set locations based on ${BASE}, which is set via env variable
# Maybe consider dropping this requirement and skipping NM validation
if [ -z "${BASE}" ]; then echo "Error, BASE not specified"; quit "Job Config"; fi
echo "Running on system ${SYSTEM:=UNDEFINED}, BASE set to ${BASE}"
if [ ! -z "${TIMING}" ]; then TIMING=(/usr/bin/time -v); fi
if [ -z "${MEM}" ]; then echo "WARNING, memory limit (in GB) not provided in variable \${MEM}; defaulting to 4G"; MEM=4; fi

# Check presence of all required variables
if [ -z "${RAWDIR}" ]; then
    echo "WARNING, no RAWDIR was given - location to place the FASTQ if differnet from the BAMfile"
    echo "Newly generated FQ file will be stored along with the bamfile on ${BAMFILE%/*}"
    RAWDIR="${BAMFILE%/*}"
	else
	echo "the newly enerated FQ files will be stored in ${RAWDIR}"
fi

if [ -z "${INFILE}" ]; then echo "ERROR, must provide .cram or .bam file in variable \${INFILE}"; quit "Input File"; fi
if [ ! -r "${INFILE}" ]; then echo "ERROR, cannot read ${INFILE}"; echo "Please make sure you provide the full path"; quit "Input File"; fi	

if [ -z "${FQNAME}" ]; then
    echo "WARNING, no FQNAME given"
    INFILENAME="${INFILE##*/}"
	FQNAME="${INFILENAME%.*}"
	FQ1="${FQNAME}_r1.fastq" ; FQ2="${FQNAME}_r2.fastq"
	echo "The newly generated FQfiles will be named after the INFILE, as ${FQ1} and ${FQ2}"
	else
	FQ1="${FQNAME}_r1.fastq" ; FQ2="${FQNAME}_r2.fastq"
	echo "the newly enerated FQ files will be named as ${FQ1} and ${FQ2}"
fi

JAVAOPTS="-Xms2g -Xmx${MEM}g -XX:+UseSerialGC -Dpicard.useLegacyParser=false"
# BAM or CRAM to Fastq
if [ "${INFILE} ~ INFILE.bam"]; then
	echo "The INFILE provided is a BAMFILE ${INFILE}"
	BAMFILENAME="${INFILE##*/}"
	
	CUR_STEP="BamtoFastq"
	start=$(${DATE}); echo "[$(display_date ${start})] ${CUR_STEP} starting"
	"${TIMING[@]}" /usr/bin/java ${JAVAOPTS} -jar "${PICARD}" SamToFastq VALIDATION_STRINGENCY=SILENT INPUT=${BAMFILE} FASTQ=${RAWDIR}/${FQ1} SECOND_END_FASTQ=${RAWDIR}/${FQ2} "$@"
	exitcode=$?
	end=$(${DATE}); echo "[$(display_date ${end})] ${CUR_STEP} finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"
	exit ${exitcode}
		
else
	# Check presence of REF file
	if [ -z ${REF} ]; then 
		[ ${INFILE} ~ ${INFILE.cram} ]; then
		echo "The INFILE provided is a CRAMFILE ${CRAMFILE}"
		CRAMFILENAME="${INFILE##*/}"
		CUR_STEP="CramtoBam"
		start=$(${DATE}); echo "[$(display_date ${start})] ${CUR_STEP} starting"
		REF=""
		# cram to bam
		# bam to fastq
		exitcode=$?
		end=$(${DATE}); echo "[$(display_date ${end})] ${CUR_STEP} finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"
		exit ${exitcode}
	else
		echo "No reference provided

fi 



