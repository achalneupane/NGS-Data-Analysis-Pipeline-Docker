#!/bin/bash

VERSION="0.1"

# GenotypeGVCFs
# Entrypoint script for docker container achalneupane/genotypegvcfs:${VERSION}

# Parameters must be defined so the script knows which file(s) to process
# Required parameters (must provide or container will quit)
#   BASE        location of parent directories for input_files, running_jobs

# # To run:
# EMAIL="achal@wustl.edu"
#   export BASE="/gscmnt/gc2645/wgs"; \
#   export TILEDB_DISABLE_FILE_LOCKING=1; \
#   export MEM=32; \
#   export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
#   export DB="-V gendb:///${BASE}/WXS_Aquilla/gvcfTest/output"; \
#   export tmpPATH="--tmp-dir=${BASE}/WXS_Aquilla/gvcfTest/temp"; \
#   export OUTvcf="${BASE}/WXS_Aquilla/gvcfTest/SNP/test.vcf"; \
#   bsub \
#   -J "Allchr" \
#     -u "${EMAIL}" \
#     -n1 -W 2880 \
#     -M 46000000 \
#     -R "rusage[mem=49152]" \
#     -o "${BASE}/WXS_Aquilla/gvcfTest/SNP" \
#     -q research-hpc \
#     -a "docker(achalneupane/genotypegvcfs)" \
#   entrypoint.sh;

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

# temporary directory path
if [ -z "${tmpPATH}" ]; then echo "WARNING, temporary directory not provided in variable \${tmpPATH}; defaulting to tmp directory variant"; tmpPATH="/tmp"; fi

# REF
if [ -z "${REF}" ]; then
    echo "ERROR, must provide REFERENCE.fa file in variable \${REF}"
    quit "Input File"
fi
if ! [ -r "${REF}" ]; then
    echo "ERROR, cannot read ${REF}"
    echo "Please make sure you provide the full path"
    quit "Input File"
fi


# Database
if [ -z "${DB}" ]; then
    echo "ERROR, must provide genomicsDB database DIR in variable \${DB}"
    quit "Input Dir"
fi


# OUTvcf
if [ -z "${OUTvcf}" ]; then
    echo "ERROR, must provide VCF output file name in variable \${OUTvcf}"
    quit "Output File"
fi



if [ ! -z "${TIMING}" ]; then TIMING=(/usr/bin/time -v); fi

# Set locations based on ${BASE}, which is set via env variable
# Maybe consider dropping this requirement and skipping NM validation
if [ -z "${BASE}" ]; then echo "Error, BASE not specified"; quit "Job Config"; fi
echo "Running on system ${SYSTEM:=UNDEFINED}, BASE set to ${BASE}"

# gendb://my_database

################# 
#GenotypeGVCFs
#################
CUR_STEP="GenotypeGVCF"
start=$(${DATE}); echo "[$(display_date ${start})] ${CUR_STEP} starting"
"${TIMING[@]}" ${GATKv4} --java-options "-Xms4G -Xmx32G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true"  \
GenotypeGVCFs \
   -R ${REF} \
   ${DB} \
   ${tmpPATH} \
   -O ${OUTvcf}
exitcode=$?
end=$(${DATE}); echo "[$(display_date ${end})] ${CUR_STEP} finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"
################# 




exit ${exitcode}
