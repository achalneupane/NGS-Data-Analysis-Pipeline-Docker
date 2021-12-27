#!/bin/bash
#make sure vcf file is unzipped for this software
#BAMFILE="MAP_86661^8038741470^202004_USUHS_EOAD-WGS_gDNA_EOLUS.aln.srt.isec-paddedexome.markd.recal.bam"
#VCFFILE="MAP_86661^8038741470^202004_USUHS_EOAD-WGS_gDNA_EOLUS.aln.srt.isec-paddedexome.markd.recal.raw.snps.indels.g.vcf"
#INDEX_BAMFILE="MAP_86661^8038741470^202004_USUHS_EOAD-WGS_gDNA_EOLUS.aln.srt.isec-paddedexome.markd.recal.bai"
#OUT_DIRECTORY="/home/blenny/verifybamid_output_MAP_86661.log"


VERSION="2.0"

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


if [ -z "${BAMFILE}" ]; then echo "ERROR, must provide .bam file in variable \${BAMFILE}"; exit
fi
if ! [ -r "${BAMFILE}" ]; then echo "ERROR, cannot read ${BAMFILE}. Please make sure you provide the full path."; exit
fi

if [ -z "${VCFFILE}" ]; then echo "ERROR, must provide .vcf file in variable \${VCFFILE}"; exit
fi
if ! [ -r "${VCFFILE}" ]; then echo "ERROR, cannot read ${VCFFILE}. Please make sure you provide the full path."; exit
fi

if [ -z "${INDEX_BAMFILE}" ]; then echo "ERROR, must provide .bai file in variable \${INDEX_BAMFILE}"; exit
fi
if ! [ -r "${INDEX_BAMFILE}" ]; then echo "ERROR, cannot read ${INDEX_BAMFILE}. Please make sure you provide the full path."; exit
fi

if [ ! -z "${TIMING}" ]; then TIMING=(/usr/bin/time -v); fi

CUR_STEP="verifyBamID"
start=$(${DATE}); echo "[$(display_date ${start})] ${CUR_STEP} starting"
"${TIMING[@]}" "${CUR_STEP}" --vcf "${VCFFILE}" --bam "${BAMFILE}" --bai "${INDEX_BAMFILE}" --out "${WORKDIR}" --verbose
exitcode=$?
end=$(${DATE}); echo "[$(display_date ${end})] ${CUR_STEP} finished, exit code: ${exitcode}, step time $(date_diff ${start} ${end})"
exit ${exitcode}