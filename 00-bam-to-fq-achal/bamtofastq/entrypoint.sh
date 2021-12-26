#!/bin/bash

VERSION="2.0"
# Entrypoint script for docker container achalneupane/bamtofastq:${VERSION} 


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

if [ -z "${INBAM}" ]; then
    echo "ERROR, must provide .sam or .bam file in variable \${INBAM}"
    quit "Input File"
fi
if [ ! -r "${INBAM}" ]; then
    echo "ERROR, cannot read ${INBAM}"
    echo "Please make sure you provide the full path"
    quit "Input File"
fi

# Remove files as you go. Set to 0 for testing. Only set if CLEANUP has not already been declared at the command line.
if [ -z "${REMOVE_INPUT}" ]; then REMOVE_INPUT=0; fi


# Set locations based on ${BASE}, which is set via env variable
# Maybe consider dropping this requirement and skipping NM validation
if [ -z "${BASE}" ]; then echo "Error, BASE not specified"; quit "Job Config"; fi
echo "Running on system ${SYSTEM:=UNDEFINED}, BASE set to ${BASE}"

if [ ! -z "${TIMING}" ]; then TIMING=(/usr/bin/time -v); fi



# RevertSam
CUR_STEP="samtools"
echo "RUNNING Current step for SM::${SM}, DNA::${DNA}, PR::${PR}, FULLSM::${FULLSM}"


start=$(${DATE}); echo "[$(display_date ${start})] ${CUR_STEP} starting"
"${TIMING[@]}" ${CUR_STEP} collate -uO ${INBAM} | ${CUR_STEP} fastq - -@ ${THREADS} -N -1 ${FQ_OUT1} -2 ${FQ_OUT2} -0 /dev/null -s /dev/null -n 
exitcode=$?
end=$(${DATE}); echo "[$(display_date ${end})] ${CUR_STEP} finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"

# BAM lines
BAMLINES=$(samtools idxstats "${INBAM}" | awk '{s+=$3+$4} END {print s*4}')
# FQ lines
FQ1LINES=$(cat ${FQ_OUT1} | wc -l)
FQ2LINES=$(cat ${FQ_OUT2} | wc -l)


# if the DEBUG string length is not zero; then
set -x
if [ ! -z ${DEBUG} ]; then
echo "${BAMLINES} lines in .bam"
echo "${FQ1LINES} lines in all ${FQ_OUT1} files"
echo "${FQ2LINES} lines in all ${FQ_OUT2} files"

if [ $(echo "scale=2;${FQ1LINES}/${BAMLINES} > 0.40" | bc) -eq 0 ]; then
    echo "Warning, ${FQ_OUT1} file contains ${FQ1LINES} lines - less than 40% of the number of reads of .bam file"
  fi
if [ $(echo "scale=2;${FQ2LINES}/${BAMLINES} > 0.4" | bc) -eq 0 ]; then
    echo "Warning, ${FQ_OUT1} file contains ${FQ1LINES} lines - less than 40% of the number of reads of .bam file"
  fi  
fi


## Now create RG files
# in place of 'tab' can also use \\t
echo "DOING paired FASTQ lane splitting !!!"

if [ ${FQ1LINES} -eq ${FQ2LINES} ] && [ ${FQ2LINES} -ne 0 ]; then
# FQ1
echo "Splitting lanes from FASTQ: ${FQ_OUT1}"
awk -v sm="$SM" -v dna="$DNA" -v pr="$PR" -v fullsm="$FULLSM" -v out_dir="${OUT_DIR}" -v tab="\\\\t" '
BEGIN{FS = ":"} 
/^@HS2000/{
  arr=$1
  sub(/^@+/,"",arr)
  lane=$2
  close(outputFile)  
  outputFile=out_dir"/"fullsm"."arr"."lane".r1.fastq"
  outputFile2=out_dir"/"fullsm"."arr"."lane".rgfile"
  print "@RG" tab "ID:"arr"."lane tab "DS:"fullsm tab "LB:"dna tab "PL:illumina" tab "PU:"arr"."lane tab "SM:"sm > (outputFile2)
  close(outputFile2)
}
{
  print >> (outputFile)
}' ${FQ_OUT1}

# FQ2
echo "Splitting lanes from FASTQ: ${FQ_OUT2}"
awk -v sm="$SM" -v dna="$DNA" -v pr="$PR" -v fullsm="$FULLSM" -v out_dir="${OUT_DIR}" -v tab="\\\\t" '
BEGIN{FS = ":"} 
/^@HS2000/{
  arr=$1
  sub(/^@+/,"",arr)
  lane=$2
  close(outputFile)  
  outputFile=out_dir"/"fullsm"."arr"."lane".r2.fastq"
}
{
  print >> (outputFile)
}' ${FQ_OUT2}

fi

if [ ${exitcode} -eq 0 ] && [ ${REMOVE_INPUT} -eq 1 ]; then
 rm -fv "${INBAM}" 2>/dev/null
 rm -fv "${FQ_OUT1}" 2>/dev/null
 rm -fv "${FQ_OUT2}" 2>/dev/null
fi

# Need to add CLEANUP variable as well

exit ${exitcode}