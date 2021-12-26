#!/bin/bash

VERSION="0.2"
GATK_THREADS=1
# Entrypoint script for docker container achalneupane/variantevalref:${VERSION}

# instead of defining DBSNP144 or DBSNP138, we will name them DBSNP1 and then we specify which DNSP database we use
# this runs on GATV3.6

# Required parameters (must provide or container will quit)
#   BASE        location of parent directories for input_files, running_jobs
#   GVCF        .g.vcf.gz input file on which to perform Variant Evaluation
#   RUN_TYPE    required for now, may be dropped in the future
# 	CCDS 
# Just in case of symlinks, which cause problems for some programs
# REF="${BASE}/Genome_Ref/GRCh38/human_g1k_v38_decoy.fasta"
# DBSNP1="${BASE}/Genome_Ref/GRCh38/20190522_bundle/dbsnp_138.hg38.vcf.gz"
 

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

if [ -z "${MEM}" ]; then echo "WARNING, memory limit (in GB) not provided in variable \${MEM}; defaulting to 4G"; MEM=4; fi

if [ -z "${GVCF}" ]; then echo "ERROR, must provide .g.vcf.gz file in variable \${GVCF}"; quit "Input File"; fi

if ! [ -r "${GVCF}" ]; then echo "ERROR, cannot read ${GVCF}"; echo "Please make sure you provide the full path"; quit "Input File"; fi

if [ -z "${RUN_TYPE}" ]; then echo "ERROR, must provide run type (genome, paddedexome, exome) in variable \${RUN_TYPE}"; 
# TODO: Should validate RUN_TYPE against genome, exome, paddedexome
    quit "Job Config"
fi

if [ ! -z "${TIMING}" ]; then TIMING=(/usr/bin/time -v); fi

# Set locations based on ${BASE}, which is set via env variable
# Maybe consider dropping this requirement and skipping NM validation
if [ -z "${BASE}" ]; then echo "Error, BASE not specified"; quit "Job Config"; fi
echo "Running on system ${SYSTEM:=UNDEFINED}, BASE set to ${BASE}"

SAMPLEID="$(echo "${GVCF##*/}" | sed "s/.aln.srt.isec-${RUN_TYPE}.markd.recal.raw.snps.indels.g.vcf.gz//g")"
if [ "${RUN_TYPE}" = "paddedexome" ]; then RUN_TYPE="exome"; fi

if [ "${RUN_TYPE}" = "genome" ]; then
    GATK_PARAMS=""
elif [ "${RUN_TYPE}" = "exome" ]; then
    GATK_PARAMS="-L ${CCDS}"
else
    echo "Error, invalid RUN_TYPE specified"
    quit "Job Config"
fi

JAVAOPTS="-Xms2g -Xmx${MEM}g -XX:+UseSerialGC -Dpicard.useLegacyParser=false"

# VariantEval
CUR_STEP="VariantEval"
start=$(${DATE}); echo "[$(display_date ${start})] ${CUR_STEP} starting"
# --eval: does not support carets in the evaluation name as of GATK-3.6.0
# Change carets to hyphens so the .gatkreport has the ${SAMPLEID} instead of just "eval"
SAMPLEID_VE="${SAMPLEID//^/-}"
echo ${SAMPLEID_VE}
"${TIMING[@]}" java ${JAVAOPTS} -jar "${GATK360}" \
  -T "${CUR_STEP}" \
  -R "${REF}" \
  ${GATK_PARAMS} \
  -nt "${GATK_THREADS}" \
  --dbsnp "${DBSNP1}" \
  --eval:"${SAMPLEID_VE}" "${GVCF}" \
  -o "${GVCF%.g.vcf.gz}.${RUN_TYPE}-varianteval.gatkreport"
exitcode=$?
end=$(${DATE}); echo "[$(display_date ${end})] ${CUR_STEP} finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"

if [ "${SYSTEM}" = "MGI" ]; then
  mkdir -p "${BASE}/variant_calling/completed_jobs/" \
    && mv -vi "${GVCF%/*}" "${BASE}/variant_calling/completed_jobs/"
fi

exit ${exitcode}