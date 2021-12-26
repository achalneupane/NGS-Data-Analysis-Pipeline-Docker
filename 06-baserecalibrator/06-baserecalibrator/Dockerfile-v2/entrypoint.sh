#!/bin/bash

VERSION="0.1.0"

# Entrypoint script for Docker container achalneupane/baserecalibrator:${VERSION}

#### PARAMETERS that MUST be defined:
# Required parameters (must provide or container will quit)
#   BASE        location of parent directories for input_files, running_jobs
#   BAMFILE     .bam input file on which to perform Base Recalibration == *.GATKready.bam input file to check base quality; alternatively "${WORKDIR}/${SAMPLEID}_GATKready.bam"
#   RUN_TYPE    required for now, may be dropped in the futur
# 	OUTFILE		*.${RUN_TYPE}_recal.table1 
#	REF="${BASE}/Genome_ref/GRCh38Genome_Ref/Homo_sapiens_assembly38.fasta"
# 	MILLS_GOLD="${BASE}/Genome_ref/GRCh38/20190522_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
# 	DBSNP1		DBSNP138="${BASE}/Genome_ref/GRCh38/20190522_bundle/dbsnp_138.hg38.vcf.gz"
# 	ONEKGP1="${BASE}/Genome_ref/GRCh38Genome_Ref/20190522_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
# 	CCDS="${BASE}/GATK_pipeline/Reference/GRCh37.CCDS.exons.sorted.bed" ## NEED TO UPDATE  - Used only at Variant Eval
#   $@          full paths to .bam files to mark / merge (command-line parameter), at least 1 required

### Optional parameters
# SHELLDROP 	Drop to shell isntead of running anything (used with docker)
# TIMING 		Do /usr/bin/time -v timing of steps
# MEM         Memory Limit in GB (e.g. 32), defaults to 4

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

if [ -z "${MEM}" ]; then echo "WARNING, memory limit (in GB) not provided in variable \${MEM}; defaulting to 4G"; MEM=4; fi

if [ -z "${BAMFILE}" ]; then
    echo "ERROR, must provide .sam or .bam file in variable \${BAMFILE}"
    quit "Input File"
fi
if ! [ -r "${BAMFILE}" ]; then
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

# Set locations based on ${BASE}, which is set via env variable
# Maybe consider dropping this requirement and skipping NM validation
if [ -z "${BASE}" ]; then echo "Error, BASE not specified"; quit "Job Config"; fi
echo "Running on system ${SYSTEM:=UNDEFINED}, BASE set to ${BASE}"


if [ "${RUN_TYPE}" = "genome" ]; then
    GATK_PARAMS=""
elif [ "${RUN_TYPE}" = "exome" ]; then
    GATK_PARAMS="-L ${COVERED_BED}"
elif [ "${RUN_TYPE}" = "paddedexome" ]; then
    GATK_PARAMS="-L ${PADDED_BED}"
else
    echo "Error, invalid RUN_TYPE specified"
    quit "Job Config"
fi

GATK_THREADS=1
## FOLLOWING GATK best practices https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165
# BaseRecalibrator #1 https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_BaseRecalibrator.php
# usage:  gatk BaseRecalibrator -I my_reads.bam -R reference.fasta --known-sites sites_of_variation.vcf --known-sites another/optional/setOfSitesToMask.vcf -O recal_data.table
#   INPUT: _GATKready.bam; OUTPUT: _${RUN_TYPE}_recal.table
CUR_STEP="BaseRecalibrator I"
start=$(${DATE}); echo "[$(display_date ${start})] ${CUR_STEP} starting"
${GATKv4} --java-options "-Xms4G -Xmx32G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	BaseRecalibrator \
	-I "${BAMFILE}" \
	-R "${REF}" \
	--known-sites "${MILLS_GOLD}" \
	--known-sites "${DBSNP1}" \
	--known-sites "${ONEKGP1}" \
	-O "${BAMFILE%.bam}.recal.table1"
exitcode=$?
end=$(${DATE}); echo "[$(display_date ${end})] ${CUR_STEP} finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"
#exit ${exitcode}

## ApplyRecalibrator BQSR - https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_ApplyBQSR.php
# usage: gatk ApplyBQSR -R reference.fasta -I input.bam --bqsr-recal-file recalibration.table -O output.bam
CUR_STEP="Apply BQSR"
start=$(${DATE}); echo "[$(display_date ${start})] ${CUR_STEP} starting"
${GATKv4} --java-options "-Xms4G -Xmx32G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	ApplyBQSR \
	-R "${REF}" \
	-I "${BAMFILE}" \
	-bqsr-recal-file "${BAMFILE%.bam}.recal.table1" ${GATK_PARAMS} \
	-O "${BAMFILE%.bam}.recal.bam"
exitcode=$?
end=$(${DATE}); echo "[$(display_date ${end})] ${CUR_STEP} finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"
#exit ${exitcode}

##Analyze covariates - https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_AnalyzeCovariates.php
## 1 - plot a single recalibration table
CUR_STEP="AnalyzeCovariates I - one plot"
start=$(${DATE}); echo "[$(display_date ${start})] ${CUR_STEP} starting"
${GATKv4} --java-options "-Xms4G -Xmx32G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	AnalyzeCovariates \
	-bqsr "${BAMFILE%.bam}.recal.table1" \
	-plots "${BAMFILE%.bam}_AnalyzeCovariates.pdf"
exitcode=$?
end=$(${DATE}); echo "[$(display_date ${end})] ${CUR_STEP} finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"
#exit ${exitcode}

## 2 - plot "before" (first pass) and "after" (second pass) recalibration tables 
## 2.1 - generate the post recal table (table2) by using the recalibrated bam file as input file
CUR_STEP="BaseRecalibrator II"
start=$(${DATE}); echo "[$(display_date ${start})] ${CUR_STEP} starting"
${GATKv4} --java-options "-Xms4G -Xmx32G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true"  \
	BaseRecalibrator \
	-I "${BAMFILE}" \
	-R "${REF}" \
	--known-sites "${MILLS_GOLD}" \
	--known-sites "${DBSNP1}" \
	--known-sites "${ONEKGP1}" \
	-O "${BAMFILE%.bam}.recal.table2"
exitcode=$?
end=$(${DATE}); echo "[$(display_date ${end})] ${CUR_STEP} finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"
#exit ${exitcode}

## 2.2 plot before and after
CUR_STEP="AnalyzeCovariates II - before after plots"
start=$(${DATE}); echo "[$(display_date ${start})] ${CUR_STEP} starting"
${GATKv4} --java-options "-Xms4G -Xmx32G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true"  \
	AnalyzeCovariates \
	-before "${BAMFILE%.bam}.recal.table1"\
	-after "${BAMFILE%.bam}.recal.table2" \
	-plots "${BAMFILE%.bam}_before-after-plots.pdf"
exitcode=$?
end=$(${DATE}); echo "[$(display_date ${end})] ${CUR_STEP} finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"
#exit ${exitcode}






