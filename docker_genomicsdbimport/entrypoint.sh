#!/bin/bash

VERSION="0.1"

# Entrypoint script for docker container achalneupane/genomicsdbimport:${VERSION}

# Parameters must be defined so the script knows which file(s) to process
# Required parameters (must provide or container will quit)
#   BASE        location of parent directories for input_files, running_jobs
#   BAMFILE     .bam input file on which to perform Haplotype Variant Calling
#   RUN_TYPE    required for now, may be dropped in the future
# Just in case of symlinks, which cause problems for some programs
# REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
# CCDS="${BASE}/GATK_pipeline/Reference/GRCh83/GRCh38.CCDS.exons.sorted.bed" ;\
# DBSNP1="${BASE}/Genome_Ref/GRCh38/20190522_bundle/dbsnp_138.hg38.vcf.gz" ;\
# PADDED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Padded.GRCh38.bed" ;\
# COVERED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Covered.GRCh38.bed" ;\

# Optional parameters
#   MEM         Memory Limit in GB (e.g. 32), defaults to 4
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

# Mem
if [ -z "${MEM}" ]; then echo "WARNING, memory limit (in GB) not provided in variable \${MEM}; defaulting to 4G"; MEM=4; fi
# gatk threads
if [ -z "${GATK_THREADS}" ]; then echo "WARNING, Number of GATK threads not provided in variable \${GATK_THREADS}; defaulting to 5 threads"; GATK_THREADS=5; fi
# batch size
if [ -z "${batch_size}" ]; then echo "WARNING, Sample batch size not provided in variable \${batch_size}; defaulting to 50 samples"; batch_size=50; fi
# Max interval
if [ -z "${MaxIntervals}" ]; then echo "WARNING, Max intervals for parallelization not provided in variable \${MaxIntervals}; defaulting to 1 variant"; MaxIntervals=1; fi
# temporary directory path
if [ -z "${tmpPATH}" ]; then echo "WARNING, temporary directory not provided in variable \${tmpPATH}; defaulting to tmp directory variant"; tmpPATH="/tmp"; fi
# Merge rule
if [ -z "${mergeTF}" ]; then echo "WARNING, merge intervals boolean not provided in variable \${mergeTF}; defaulting to TRUE"; mergeTF="true"; fi
#consolidate rule
if [ -z "${consolidateTF}" ]; then echo "WARNING, consolidate-batch boolean not provided in variable \${consolidateTF}; defaulting to FALSE"; consolidateTF="false"; fi



# GATK intervals and BED; either -L or --interval
if [ -z "${GATKinverval}" ]; then echo "WARNING, GATKinverval not defined in variable \${GATKinverval}; defaulting to NULL"; GATKinverval=""; fi
# Coveered intervals, either chromosome number or bed regions (could be multiple intervals)   
if [ -z "${covered_bed}" ]; then echo "WARNING, covered_bed not defined in variable \${covered_bed}; defaulting to NULL"; covered_bed=""; fi

wantedIntervals="${GATKinverval} ${covered_bed}"


if [ -z "${mylist}" ]; then
    echo "ERROR, must provide sample list file in variable \${mylist}"
    quit "Input File"
fi
if ! [ -r "${mylist}" ]; then
    echo "ERROR, cannot read ${mylist}"
    echo "Please make sure you provide the full path"
    quit "Input File"
fi



# if [ -z "${RUN_TYPE}" ]; then
#     echo "ERROR, must provide run type (genome, paddedexome, exome) in variable \${RUN_TYPE}"
#     # TODO: Should validate RUN_TYPE against genome, exome, paddedexome
#     quit "Job Config"
# fi

if [ ! -z "${TIMING}" ]; then TIMING=(/usr/bin/time -v); fi

# Set locations based on ${BASE}, which is set via env variable
# Maybe consider dropping this requirement and skipping NM validation
if [ -z "${BASE}" ]; then echo "Error, BASE not specified"; quit "Job Config"; fi
echo "Running on system ${SYSTEM:=UNDEFINED}, BASE set to ${BASE}"

# if [ "${RUN_TYPE}" = "genome" ]; then
#     GATK_PARAMS=""
# elif [ "${RUN_TYPE}" = "exome" ]; then
#     GATK_PARAMS="-L ${COVERED_BED}"
# elif [ "${RUN_TYPE}" = "paddedexome" ]; then
#     GATK_PARAMS="-L ${PADDED_BED}"
# else
#     echo "Error, invalid RUN_TYPE specified"
#     quit "Job Config"
# fi



################# 
#GenomicsDBImport
#################
CUR_STEP="GenomicsDBImport"
start=$(${DATE}); echo "[$(display_date ${start})] ${CUR_STEP} starting"
"${TIMING[@]}" ${GATKv4} --java-options "-Xms4G -Xmx32G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true"  \
 GenomicsDBImport \
    --genomicsdb-workspace-path ${WORKDIR} \
    --batch-size ${batch_size} \
    --sample-name-map ${mylist} \
    ${wantedIntervals} \
    --reader-threads ${GATK_THREADS} \
    --tmp-dir=${tmpPATH} \
    --consolidate ${consolidateTF} \
    --merge-input-intervals ${mergeTF} \
    --max-num-intervals-to-import-in-parallel ${MaxIntervals}
exitcode=$?
end=$(${DATE}); echo "[$(display_date ${end})] ${CUR_STEP} finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"
################# 




# if [ "${SYSTEM}" = "MGI" ]; then
#   mv -vf "${BAMFILE}" "${BASE}/tmp/"
# fi

exit ${exitcode}