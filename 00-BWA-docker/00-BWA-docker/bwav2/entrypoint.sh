#!/bin/bash

VERSION="2.1"
## Entrypoint script for docker container achalneupane/bwaref:${VERSION}.
## Different from bwa:1.0, this entrypoiunt requires the user to specify the REF file to use, instead of pointing to a directory in the system
## v2.1 update - correction so the script performs cleanup of entryfiles; I change CACHING from 0 to 1 at line #9

GC_THREADS=2  # benchmarking found this to be the fastest for most jobs, vs. 4 or 8
CACHING=0
if [ -z "${THREADS}" ]; then THREADS=8; fi  # Would be nice to figure out how many cores we have
if [ -z "${BWA_PIPE_SORT}" ]; then BWA_PIPE_SORT=1; fi  # Assuming pipe directly from bwa to samtools sort

# Parameters must be defined so the script knows which file(s) to process

# Required parameters (must provide or container will quit)
#   BASE        location of reference files for alignment
# 	REF

# Required parameters (choose BAMFILE or FQ1+FQ2+RGFILE or FQ+RGFILE)
#   BAMFILE     starting .bam file to operate on (not required if FQ1 and FQ2 are given)
#   FQ1, FQ2    starting .fastq files to operate on (not required if BAMFILE is given)
#   FQ          starting from a single interleaved .fastq, .fq, .fastq.gz, or .fq.gz file
#   RGFILE      required if providing FQ1+FQ2 or FQ
#   CLEANUP     set to 0 (don't cleanup) or 1 (do cleanup) to remove temporary files as the script runs
#   REMOVE_INPUT  Delete the input .bam upon successful completion

# Optional parameters
#   WORKDIR     parent directory to create JOB_TMP directory to hold all files (which often uses ${PBS_JOBID} or ${SLURM_JOBID} or ${LSB_JOBID}). Required on fenix
#   OUT_DIR     directory to hold output file from this script; defaults to ${INDIR}
#   FULLSM      sample ID, used for determining output filenames. Will be auto-detected from BAMFILE or FQ1 if not supplied
#   FULLSM_RGID sample ID_RGID, used for determining output filenames. Will be set to FULLSM if not supplied
#   CLEANUP     set to 0 (don't cleanup) or 1 (do cleanup) to remove temporary files as the script runs
#   THREADS     by default set to 4 on ewok, 8 on CHPC
#   SHELLDROP   Drop to shell instead of running anything (used with docker)
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
  echo "Leftover files may be in ${JOB_TMP}"
  echo "Cleaning up ${SORT_TMP}"
  rm -rfv "${SORT_TMP}"
  if [ -d "${JOB_TMP}" ]; then rmdir -v "${JOB_TMP}"; fi  # this might fail if files are leftover, but is handy if the job fails early
  touch "${INDIR}/.${FULLSM_RGID}.runfailed"
  if [ ${CACHING} -eq 1 ]; then
    echo "Cleaning up cached files in /tmp/${PBS_JOBID}"
#    cd "/tmp/${PBS_JOBID}"
#    rm -fv "${REF}"* "${FQ}" "${FQ1}" "${FQ2}" 2>/dev/null
  fi
  # /bin/mail -r "docker@${HOSTNAME}" -s "${SYSTEM}: ${FULLSM_RGID} FAIL at ${1} (${PBS_JOBID}${SLURM_JOBID}${LSB_JOBID})" "${EMAIL}" < /dev/null > /dev/null
  if [ "${SHELLDROP}" -eq 1 ]; then echo "Dropping to shell"; exec "/bin/bash"; else exit 1; fi
}

# trap "{ echo \"[$(display_date $(${DATE}))] Terminated by SIGTERM \"; quit \"${CUR_STEP}\"; }" SIGTERM
trap "echo \"[$(display_date $(${DATE}))] Terminated by SIGTERM \" && sleep 10s" SIGTERM

while getopts ":1:2:f:b:s:n:w:r:t:cgpei" opt; do
  case $opt in
    1)
      FQ1="$OPTARG"
      ;;
    2)
      FQ2="$OPTARG"
      ;;
    f)
      FQ="$OPTARG"
      ;;
    b)
      BAMFILE="$OPTARG"
      ;;
    s)
      FULLSM="$OPTARG"
      ;;
    n)
      FULLSM_RGID="$OPTARG"
      ;;
    w)
      WORKDIR="$OPTARG"
      ;;
    r)
      RGFILE="$OPTARG"
      ;;
    t)
      THREADS=$OPTARG
      ;;
    c)
      CLEANUP=$OPTARG
      ;;
    g)
      RUN_TYPE="genome"
      ;;
    p)
      RUN_TYPE="paddedexome"
      ;;
    e)
      RUN_TYPE="exome"
      ;;
    i)
      FASTQ_TYPE="interleaved"
      ;;
    \?)
      echo "Invalid option: -$OPTARG"
      exit 2
      ;;
    :)
      echo "Option -$OPTARG requires an argument"
      exit 2
  esac
done
shift $((OPTIND-1))

# Option for usage in docker
if [ "${SHELLDROP:=0}" -eq 1 ]; then echo "Dropping to shell"; exec "/bin/bash"; fi

# Remove files as you go. Set to 0 for testing. Only set if CLEANUP has not already been declared at the command line.
if [ -z "${CLEANUP}" ]; then CLEANUP=1; fi

# Set locations based on ${BASE}, which is set via env variable
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

if [[ -n "${BAMFILE}" && (-n "${FQ}" || -n "${FQ1}" || -n "${FQ2}") ]]; then
  echo "Error, BAMFILE and FASTQ files specified as input, this is usually an error with ENV variables"
  quit "Input File List"
fi

# Do all file does-it-exist error-checking up front
if [ -n "${BAMFILE}" ]; then  # If starting with a .bam file
  MODE="bam"
  if [ ! -r "${BAMFILE}" ]; then echo "Error, file \"${BAMFILE}\" does not exist or cannot be read"; quit "Reading Input File"; fi
  INDIR="${BAMFILE%/*}"  # Extract directory from full path
  INFILE="${BAMFILE##*/}"  # Extract filename from full path
  if [ "${INDIR}" = "${BAMFILE}" ]; then  # If the BAMFILE is in the current directory, there will be no slashes in BAMFILE
    INDIR="."
  fi
elif [ -n "${FQ1}" ] && [ -n "${FQ2}" ]; then  # If starting with 2 .fq files
  MODE="2xfq"
  if ! [[ -r "${FQ1}" && -r "${FQ2}" ]]; then echo "Error, ${FQ1} and/or ${FQ2} do not exist or cannot be read"; quit "Reading Input File"; fi
  INDIR="${FQ1%/*}"  # Extract directory from full path
  INFILE="${FQ1##*/}"  # Extract filename from full FQ1 path
  if [ "${INDIR}" = "${FQ1}" ]; then  # If the file is in the current directory, there will be no slashes in FQ1
    INDIR="."
  fi
  if [ -z "${FULLSM}" ]; then  # Extract ID from filename if not already defined
    FULLSM="${INFILE%_[rR][12].*}"
  fi
  if [ -z "${RGFILE}" ]; then echo "[ERROR] Must provide RGFILE when running in FQ1+FQ2 mode"; quit "Reading Input File"; fi
  # DEBUG -- this may cause a problem if we ever receive large .fq.bz2 files, but for WES it is okay to uncompress
  if [[ ${FQ1} =~ bz2$ ]]; then echo "Decompressing .fq.bz2 file"; if bunzip2 "${FQ1}"; then FQ1="${FQ1%.bz2}"; else quit "Decompressing Input File"; fi; fi
  if [[ ${FQ2} =~ bz2$ ]]; then echo "Decompressing .fq.bz2 file"; if bunzip2 "${FQ2}"; then FQ2="${FQ2%.bz2}"; else quit "Decompressing Input File"; fi; fi
elif [ -n "${FQ}" ]; then  # If starting with 1 interleaved .fq file
  MODE="fq"
  if [ ! -r "${FQ}" ]; then echo "Error, input file ${FQ} does not exist or cannot be read"; quit "Reading Input File"; fi
  INDIR="${FQ%/*}"  # Extract filename from full FQ path
  INFILE="${FQ##*/}"  # Extract filename from full FQ path
  if [ "${INDIR}" = "${FQ}" ]; then  # If the FQ is in the current directory, there will be no slashes in FQ
    INDIR="."
  fi
  if [ -z "${FULLSM}" ]; then  # Extract ID from filename if not already defined
    FULLSM="${INFILE%_[rR][12].*}"
  fi
  if [ -z "${RGFILE}" ]; then echo "[ERROR] Must provide RGFILE when running in FQ mode"; quit "Reading Input File"; fi
  if [[ ${FQ} =~ bz2$ ]]; then echo "Decompressing .fq.bz2 file"; if bunzip2 "${FQ}"; then FQ="${FQ%.bz2}"; else quit "Decompressing Input File"; fi; fi
else
  echo "[ERROR] Must provide 1 .bam or 2 .fastq files or 1 interleaved .fastq file" 1>&2
  echo "[ERROR] Cannot determine input file type(s), executing arguments as command ($@)" 1>&2
  exec "$@"
  exit 1
fi

# Remove files as you go. Set to 0 for testing. Only set if CLEANUP has not already been declared at the command line.
if [ -z "${REMOVE_INPUT}" ]; then REMOVE_INPUT=0; fi

if [ -z "${OUT_DIR}" ]; then OUT_DIR="${INDIR}"; fi
mkdir -p "${OUT_DIR}" || { echo "Error, cannot create ${OUT_DIR}"; quit "Setup OUT_DIR"; }

# Create a directory to contain all of the files generated while running the script
JOB_TMP="$(mktemp -d -p "${WORKDIR}/" "${LSB_JOBID:-0}.XXXXXXXXXXXXXXXX")" \
  || { echo "Error, cannot create ${JOB_TMP}"; quit "Setup JOB_TMP"; }
# Create a directory to contain all of the temporary .bam files
SORT_TMP="$(mktemp -d -p "${JOB_TMP}")" \
  || { echo "Error, cannot create temporary directory for sorting ${SORT_TMP}"; quit "Setup SORT_TMP"; }

# Setup java Garbage Collection, depending on the number of threads
if [ "${THREADS:=1}" -eq 1 ]; then
  JAVAOPTS="-Xms4g -Xmx32g -XX:+UseSerialGC"
else
  JAVAOPTS="-Xms4g -Xmx32g -XX:+UseParallelGC -XX:ParallelGCThreads=${GC_THREADS}"
fi
MEM=$((32/${THREADS}))

if [ -z "${FULLSM_RGID}" ]; then
  FULLSM_RGID="${FULLSM}"  # Workaround for now. In the future figure out how many RGs a sample has
fi

run_start=$(${DATE})
echo "[$(display_date $(${DATE}))] Starting run on ${SYSTEM}, host ${HOSTNAME}; Run Parameters:"
touch "${INDIR}/.${FULLSM_RGID}.startrun"
for i in VERSION SYSTEM SLURM_JOB_NODELIST SLURM_NTASKS REF FULLSM FULLSM_RGID INDIR BAMFILE FQ FQ1 FQ2 RGFILE WORKDIR JOB_TMP SORT_TMP OUT_DIR BWA_PIPE_SORT THREADS CLEANUP RUN_TYPE MODE FASTQ_TYPE CACHING; do
  if [ ! -z ${!i} ]; then  # bash indirection: not the value of $i, but the value of the parameter name in $i
    printf "%26s ${!i}\n" "${i}"
  fi
done
echo
echo "[$(display_date $(${DATE}))] Contents of working directory ${JOB_TMP}:"
ls -ltr "${JOB_TMP}" | sed '1d' | grep -vE "${FULLSM_RGID}.out|${FULLSM_RGID}.err"
echo

# Align the input file (.bam or .fq or .r1.fq + .r2.fq)
if [ "${MODE}" = "bam" ] || [ "${MODE}" = "2xfq" ] || [ "${MODE}" = "fq" ]; then
  # Cache files if enabled, regardless of start (bam, 2xfastq, fastq, GATKready)
  ### Would need to update this to use ${TMPDIR} or something similar in docker [20190107] ###
  if [ ${CACHING} -eq 1 ]; then
    echo "Caching REF to /tmp/${PBS_JOBID}"
    cp -av "${REF}" "/tmp/${PBS_JOBID}/" \
      && cp -av "${REF}.fai" "/tmp/${PBS_JOBID}/" \
      && cp -av "${REF%.fasta}.dict" "/tmp/${PBS_JOBID}/" \
      && cp -av "${REF}.amb" "/tmp/${PBS_JOBID}/" \
      && cp -av "${REF}.ann" "/tmp/${PBS_JOBID}/" \
      && cp -av "${REF}.bwt" "/tmp/${PBS_JOBID}/" \
      && cp -av "${REF}.pac" "/tmp/${PBS_JOBID}/" \
      && cp -av "${REF}.sa" "/tmp/${PBS_JOBID}/" \
      && REF="/tmp/${PBS_JOBID}/${REF##*/}"
    if [ "${MODE}" = "2xfq" ]; then
      echo "Caching FQ1, and FQ2 to /tmp/${PBS_JOBID}"
      OLD_FQ1="${FQ1}"
      OLD_FQ2="${FQ2}"
      cp -av "${FQ1}" "/tmp/${PBS_JOBID}/" && FQ1="/tmp/${PBS_JOBID}/${FQ1##*/}"
      cp -av "${FQ2}" "/tmp/${PBS_JOBID}/" && FQ2="/tmp/${PBS_JOBID}/${FQ2##*/}"
    elif [ "${MODE}" = "fq" ]; then
      echo "Caching FQ to /tmp/${PBS_JOBID}"
      OLD_FQ="${FQ}"
      cp -av "${FQ}" "/tmp/${PBS_JOBID}/" && FQ="/tmp/${PBS_JOBID}/${FQ##*/}"
    fi
  echo
  fi
  
  # Input: ( .bam || .fq.gz || ( .fq.gz && .fq.gz ) ); OUTPUT: .aln.bam or .srt.bam
  bwa_done=0; bwa_failcount=0
  while [ ${bwa_done} -eq 0 ]; do
    CUR_STEP="bwa+samtools"
    start=$(${DATE}); echo "[$(display_date ${start})] ${CUR_STEP} starting"
    # Put some checks in to validate ${RG}? or just let bwa handle it? Minimally need @RG RGPU RGID RGPL RGLB SM
    if [ ! -z ${RGFILE} ] && [ -r ${RGFILE} ]; then
       echo "[$(display_date $(${DATE}))] Loading ReadGroup info from ${RGFILE}"
       RG=$(<${RGFILE})
       if [ ${BWA_PIPE_SORT} -eq 1 ]; then
         SAMTOOLS_CMD=(samtools sort -@ "${THREADS}" -m "${MEM}G" -o "${OUT_DIR}/${FULLSM_RGID}.aln.srt.bam" -T "${SORT_TMP}/")
       else
         SAMTOOLS_CMD=(samtools view -b -1 -o "${JOB_TMP}/${FULLSM_RGID}.aln.bam")
         SAMTOOLS_CMD+=(\&\& samtools sort -@ "${THREADS}" -m "${MEM}G" -o "${OUT_DIR}/${FULLSM_RGID}.aln.srt.bam" -T "${SORT_TMP}/" "${JOB_TMP}/${FULLSM_RGID}.aln.bam")
         SAMTOOLS_CMD+=(\&\& rm -fv "${JOB_TMP}/${FULLSM_RGID}.aln.bam")
       fi
       if [ ! -z "${TIMING}" ]; then TIMING=(/usr/bin/time -v); SAMTOOLS_CMD=("${TIMING[@]}" "${SAMTOOLS_CMD[@]}"); fi
       if [ "${MODE}" = "bam" ]; then
         "${TIMING[@]}" bwa mem -t ${THREADS} -R "${RG}" -M -p "${REF}" \
           <(/usr/bin/java \
               -Dsamjdk.buffer_size=131072 \
               -Dsamjdk.compression_level=1 \
               -XX:GCTimeLimit=50 \
               -XX:GCHeapFreeLimit=10 \
               -Xmx128m \
               -jar "${PICARD}" SamToFastq INTERLEAVE=TRUE I="${BAMFILE}" FASTQ=/dev/stdout) \
           | eval "${SAMTOOLS_CMD[@]}"
       elif [ "${MODE}" = "2xfq" ]; then
         "${TIMING[@]}" bwa mem -t ${THREADS} -R "${RG}" -M "${REF}" "${FQ1}" "${FQ2}" \
           | eval "${SAMTOOLS_CMD[@]}"
       elif [ "${MODE}" = "fq" ]; then
         # Support for single-ended reads or interleaved reads in a single .fq file
         if [ "${FASTQ_TYPE}" = "interleaved" ]; then INTERLEAVE_OPTION="-p"; else INTERLEAVE_OPTION=""; fi
         "${TIMING[@]}" bwa mem -t ${THREADS} -R "${RG}" -M "${INTERLEAVE_OPTION}" "${REF}" "${FQ}" \
           | eval "${SAMTOOLS_CMD[@]}"
       else  # This should never be triggered with the current script
         echo "ERROR, how did you get here?"
         quit "Program Flow Error"
       fi
    else
      echo "[$(display_date $(${DATE}))] Warning, no RGFILE not specified or cannot be read. Not adding RG info to .bam"
        "${TIMING[@]}" bwa mem -t ${THREADS} -M "${REF}" "${FQ1}" "${FQ2}" \
          | /usr/bin/java -Djava.io.tmpdir=${JOB_TMP}/${FULLSM_RGID} ${JAVAOPTS} \
              -jar "${PICARD}" SortSam \
              I=/dev/stdin \
              O="${FULLSM_RGID}.srt.bam" \
              SO=coordinate \
              MAX_RECORDS_IN_RAM=2000000
    fi
    arr=(${PIPESTATUS[@]})
    exitcode=0; for i in ${arr[@]}; do ((exitcode+=i)); done
    end=$(${DATE}); echo "[$(display_date ${end})] ${CUR_STEP} finished, exit codes: ${arr[@]}, step time $(date_diff ${start} ${end})"
    ls -ltr "${JOB_TMP}" | sed '1d' | grep -vE "${FULLSM_RGID}.out|${FULLSM_RGID}.err"
    if [ ${exitcode} -eq 0 ]; then
      bwa_done=1;
      if [ ${CLEANUP} -eq 1 ]; then
        # Remove temporary directory created by SortSam
        if [ -d "${JOB_TMP}/${FULLSM_RGID}/${USER}" ]; then rm -rfv "${JOB_TMP}/${FULLSM_RGID}/${USER}"; fi
        if [ -d "${SORT_TMP}" ]; then rm -rfv "${SORT_TMP}"; fi
        if [ -d "${JOB_TMP}" ]; then rmdir -v "${JOB_TMP}"; fi
        if [ "${MODE}" = "bam" ] || [ "${SYSTEM}" = "CHPC2" ]; then rm -fv "${BAMFILE}"; fi
        if { [ "${MODE}" = "2xfq" ] || [ "${MODE}" = "fq" ]; } && { [ "${SYSTEM}" = "CHPC2" ] || [ "${SYSTEM}" = "FSL" ]; }; then
          if [ ${CACHING} -eq 1 ]; then
		   # If CACHING was enabled, the FQ FQ1 FQ2 variables no longer point to the original input files
            rm -fv "${OLD_FQ}" "${OLD_FQ1}" "${OLD_FQ2}" 2>/dev/null
          else
            if [ ${REMOVE_INPUT} -eq 1 ]; then
              echo "Removing input files ${FQ} ${FQ1} ${FQ2}"
              rm -fv "${FQ}" "${FQ1}" "${FQ2}" 2>/dev/null
			 fi
          fi
        fi
        # Would normally CLEANUP AND DELETE ${FQ} and ${RGFILE}
        # Just delete ${BAMFILE} or ${FQ1}+${FQ2} or ${FQ} now; delete all ${RGFILE} and rmdir after MarkDuplicates (assume input file location)
      fi
    else
      # So much for the while loop -- fix this later to retry bwa
      quit "${CUR_STEP}"
    fi
    echo
  done
fi
end=$(${DATE}); echo "[$(display_date ${end})] Run finished, overall time $(date_diff ${run_start} ${end})"

exit 0