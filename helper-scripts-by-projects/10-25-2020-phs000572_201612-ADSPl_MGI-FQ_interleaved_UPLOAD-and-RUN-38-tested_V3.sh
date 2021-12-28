# Updated on 10/08/2020 by Achal. This script was used to process ADSP samples from $PR: phs000572_201612
# First use phs000572_201612.sh script to download ADSP SRAs, then use Create_metafile.sh to create metafile for this project. 


unset all ENV variables
for c in $(env | cut -d '=' -f 1); do unset $c; done

cd /30/dbGaP/6109/sra/phs000572_201612/fqgz
# ls -l /80/AD/AD_Seq_Data/02.-BWA_GATK/02-gVCFs/GRCh37/singlerun/PPMI_3432^multi^PPMI_WES
# locate PPMI_3432 | grep ".g.vcf.gz$"
# cd /80/AD/AD_Seq_Data/02.-BWA_GATK/02-gVCFs/GRCh37/singlerun


###################On MGI
# remove lines matching patterns
# sed '/SRR1391349\|SRR1411264/d' phs000572_201612_sm_dna-pr_rgbase_FQext_fullsm_1-250.csv > "${PR}_sm_dna-pr_rgbase_FQext_fullsm.csv"



## Log onto MGI
PR="phs000572_201612"
USER="achal"
BASE="/gscmnt/gc2645/wgs"
DEST="${BASE}/WXS_Aquilla/01-RAW"
mkdir -p ${PR}
EMAIL="achal@wustl.edu"
USER="achal"


#################
# missing samples
cd /gscmnt/gc2645/wgs/WXS_Aquilla/02-TRANSIT/phs000572_201612
ls ./*/*vcf.gz| cut -d/ -f2 > /gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/${PR}_sm_dna-pr_rgbase_FQext_fullsm_vcf_done_samples.csv
cd /gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/phs000572_201612
ls | grep phs000572_201612$ | cut -d/ -f2 > /gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/${PR}_sm_dna-pr_rgbase_FQext_fullsm_all_samples.csv

VCF="/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/${PR}_sm_dna-pr_rgbase_FQext_fullsm_vcf_done_samples.csv"
ALL_Samples="/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/${PR}_sm_dna-pr_rgbase_FQext_fullsm_all_samples.csv"
comm -23 <(sort ${ALL_Samples}) <(sort ${VCF}) >  /gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/${PR}_sm_dna-pr_rgbase_FQext_missing_6_samples.csv

MISSING="/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/${PR}_sm_dna-pr_rgbase_FQext_missing_6_samples.csv"


FILE="${DEST}/${PR}_sm_dna-pr_rgbase_FQext_fullsm_6_samples_missing.csv" 

rm  ${FILE} 
####################################################
for dir in $(cat ${MISSING})/; do
    filename=${dir}
    SM="$(echo "$filename" | cut -f 1 -d '^')"
    DNA="$(echo "$filename" | cut -d "^" -f 2 )"
    # R1="r1.fq.bz2"
    # R2="r2.fq.bz2"
    FULLSM=${dir%/}
    # RGBASE=( $(find $dir -name "*.rgfile" -printf "%f\n" | cut -d "." -f 1,2 | sort -V) )
    RGBASE=( $(find $dir -name "*.rgfile" -printf "%f\n" | cut -d "." -f 1,2 | sort -t"^" -n -k1) )
    for RG in ${RGBASE[@]}; do
    R=( $(find $dir -name "$RG.fq.gz" -printf "%f\n" | cut -d "." -f 3,4,5) )    
    echo "${SM},${DNA},${PR},${RG},${R},${FULLSM}"  >> ${FILE}
done      
done

################

# mv ${PR}_sm_dna-pr_rgbase_FQext_fullsm_21_40_samples.csv "${PR}_sm_dna-pr_rgbase_FQext_fullsm.csv"
# mv ${PR}_sm_dna-pr_rgbase_FQext_fullsm_61_80_samples.csv "${PR}_sm_dna-pr_rgbase_FQext_fullsm.csv"
# mv ${PR}_sm_dna-pr_rgbase_FQext_fullsm_81_100_samples.csv "${PR}_sm_dna-pr_rgbase_FQext_fullsm.csv"
# mv ${PR}_sm_dna-pr_rgbase_FQext_fullsm_101_116_samples.csv "${PR}_sm_dna-pr_rgbase_FQext_fullsm.csv"
# mv ${PR}_sm_dna-pr_rgbase_FQext_fullsm_117_120_samples.csv "${PR}_sm_dna-pr_rgbase_FQext_fullsm.csv"
# mv ${PR}_sm_dna-pr_rgbase_FQext_fullsm_121_140_samples.csv "${PR}_sm_dna-pr_rgbase_FQext_fullsm.csv"
# mv ${PR}_sm_dna-pr_rgbase_FQext_fullsm_141_160_samples.csv "${PR}_sm_dna-pr_rgbase_FQext_fullsm.csv"
# mv ${PR}_sm_dna-pr_rgbase_FQext_fullsm_161_176_samples.csv "${PR}_sm_dna-pr_rgbase_FQext_fullsm.csv"
# mv ${PR}_sm_dna-pr_rgbase_FQext_fullsm_177_180_samples.csv "${PR}_sm_dna-pr_rgbase_FQext_fullsm.csv"
# mv ${PR}_sm_dna-pr_rgbase_FQext_fullsm_final_11_samples.csv "${PR}_sm_dna-pr_rgbase_FQext_fullsm.csv"
mv $FILE "${PR}_sm_dna-pr_rgbase_FQext_fullsm.csv"


# ## for  missing  samples that have FASTQ on MGI, get this file from create_metafile.sh
# cp /gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/${PR}_sm_dna-pr_rgbase_FQext_fullsm_missing_samples.csv "/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/${PR}_sm_dna-pr_rgbase_FQext_fullsm.csv"


# ## get only lines matching patterns
# cat ${PR}_sm_dna-pr_rgbase_FQext_fullsm_1-250.csv | grep "SRR1391349\|SRR1411264"
##############
## 2019-08-02
##############

### In Fenix

metafile="${PR}_sm_dna-pr_rgbase_FQext_fullsm.csv"
# vim ${metafile}   ## copy here the metadata for the DIAN project from the excel file
sed -i 's/\t/,/g' ${metafile}

################################

LOOKUP="${BASE}/WXS_Aquilla/01-RAW/${PR}_sm_dna-pr_rgbase_FQext_fullsm.csv"

#Generate the workfile
cut -d, -f1-3 $LOOKUP | tr "," "^" | sort -u > ${DEST}/${PR}-worklist_fullsm.csv
# WORKLISTtmp="${DEST}/${PR}-worklist_fullsm.csv"
WORKLIST="${DEST}/${PR}-worklist_fullsm.csv"

#################################################################################


bjobs | grep RUN | cut -d\  -f1 | xargs bkill 
bjobs | grep PEN | cut -d\  -f1 | xargs bkill 
bjobs | grep SSUSP | cut -d\  -f1 | xargs bkill 

# 9, 11, 12, 16


START=1;
END=6;
SHELLDROP=1
# for FULLSM in $(sed -n "${START},${END}p" "${WORKLIST}"); do \
Snumber=1;
export BASE="/gscmnt/gc2645/wgs"; \
export WORKDIR="${BASE}/tmp"; \
export THREADS=4; \
export BWA_PIPE_SORT=1; \
export PR="phs000572_201612"
export BASE="/gscmnt/gc2645/wgs"; \
export TIMING=1; \
LOOKUP_COL_SM=1; \
LOOKUP_COL_DNA=2; \
export REMOVE_INPUT=1; \
LOOKUP_COL_PR=3; \
LOOKUP_COL_RGBASE=4; \
LOOKUP_COL_FQEXT=5; \
# LOOKUP_COL_FQ2EXT=6; \
unset FBASE; \
export RUN_TYPE="paddedexome"; \
export FASTQ_TYPE="interleaved"; \
export MODE="fq"; \
for FULLSM in $(sed -n "${START},${END}p" "${WORKLIST}"); do \
echo "*****************************************************"; \
echo "*****************************************************"; \
echo "Doing sample number**********: " $Snumber; \
echo "*****************************************************"; \
echo "*****************************************************"; \
((Snumber=${Snumber}+1)); \
SM="$(echo "${FULLSM}" | cut -d^ -f1)"; \
IFS=$'\n' export DNA=($(awk -F, "\$${LOOKUP_COL_SM} == \"${SM}\"" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_DNA} | sort -u)); \
if [ ${#DNA[@]} -gt 1 ]; then echo "Warning, \${DNA} not unique for ${SM} (n=${#DNA[@]}: ${DNA[@]})"; fi; \
DNA="$(echo "${FULLSM}" | cut -d^ -f2)"; \
PR="$(echo "${FULLSM}" | cut -d^ -f3)"; \
export OUT_DIR="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}"; \
echo -e "00 - Starting jobs per sample FULLSM ${FULLSM}"; \
for RGBASE in $(grep "${SM},${DNA},${PR}" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_RGBASE}); do \
SM="$(echo "${FULLSM}" | cut -d^ -f1)"; \
DNA="$(echo "${FULLSM}" | cut -d^ -f2)"; \
PR="$(echo "${FULLSM}" | cut -d^ -f3)"; \
## change this to where fastq files are located
FBASE="/home/achal/SRA/phs000572_201612/fqgz/${FULLSM}"; \
FQEXT=($(awk -F, "\$${LOOKUP_COL_RGBASE} == \"${RGBASE}\"" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_FQEXT})); \
# FQ2EXT=($(awk -F, "\$${LOOKUP_COL_RGBASE} == \"${RGBASE}\"" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_FQ2EXT})); \
RGFILE="${FBASE}/${RGBASE}.rgfile";\
FQ="${FBASE}/${RGBASE}.${FQEXT}";\
# FQ2="${FBASE}/${RGBASE}.${FQ2EXT}";\
DEST="${BASE}/WXS_Aquilla/01-RAW";\
echo -e "00a - Uploading FASTQ and rgfiles for sample ${FULLSM} and RGBASE ${RGBASE}\nFQ:${FQ}\nRGFILE:${RGFILE}\nDNA:${DNA}";\
mkdir ${DEST}/${PR}/${FULLSM}; \
# ### rsync with copy referent files as we are copying symlinks for this PPMI data from fenix
# rsync -avh -L ${USER}@fenix.psych.wucon.wustl.edu:${FQ} ${DEST}/${PR}/${FULLSM}/ ;\
# rsync -avh ${USER}@fenix.psych.wucon.wustl.edu:${RGFILE} ${DEST}/${PR}/${FULLSM}/ ;\
# # done
# # done
# # done
export RGFILE="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.rgfile"; \
export FULLSM_RGID="${RGBASE}"; \
unset BAMFILE; \
unset FQ1; \
unset FQ2; \
IFS=$'\n' FQEXT=($(awk -F, "\$${LOOKUP_COL_RGBASE} == \"${RGBASE}\"" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_FQEXT})); \
# IFS=$'\n' FQ2EXT=($(awk -F, "\$${LOOKUP_COL_RGBASE} == \"${RGBASE}\"" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_FQ2EXT})); \
echo -e "01 - Starting bwa per sample ${FULLSM} and RGBASE ${RGBASE}";\
export RGFILE="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.rgfile"; \
export FQ="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.${FQEXT}"; \
# export FQ2="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.${FQ2EXT}"; \
export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
echo "######### RGBASE is:  $RGBASE"; \
echo "######### RGFILE is: $(ls -lht $RGFILE)"; \
echo "######### FASTQ is : $(ls -lht $FQ)"; \
# done; \
# done
export MEM=76; \
bsub \
     -J "${RGBASE}_s01alnsrt" \
     -o "${BASE}/WXS_Aquilla/04-LOGS/${RGBASE}_s01alnsrt.%J" \
     -u "${EMAIL}" \
     -n ${THREADS} -W 2880 \
     -M 86000000 \
     -R "rusage[mem=86152]" \
     -q research-hpc \
     -a 'docker(vifehe/bwaref)' \
     entrypoint.sh; \
echo -e "02 - Starting validatesam per sample ${FULLSM} and RGBASE ${RGBASE}\nBAMFILE:${BAMFILE}"; \
export MEM=42; \
export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.bam"; \
echo ${BAMFILE}; \
export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
bsub \
     -w "done(\"${RGBASE}_s01alnsrt\")" \
     -J "${RGBASE}_s02vldate" \
     -o "${BASE}/WXS_Aquilla/04-LOGS/${RGBASE}_s02vldate.%J" \
     -u "${EMAIL}" \
     -n1 -W 1360 \
     -R "rusage[mem=48000]" \
     -q research-hpc \
     -a "docker(vifehe/validatesamref)" \
     entrypoint.sh -IGNORE INVALID_VERSION_NUMBER -IGNORE INVALID_TAG_NM; \
echo -e "03 - Starting intersectbed per sample ${FULLSM} and RGBASE ${RGBASE}\nBAMFILE:${BAMFILE}";\
export RUN_TYPE="paddedexome"; \
export BEDFILE="${BASE}/Genome_Ref/GRCh38/Capture_Padded.GRCh38.bed"; \
export COVERED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Covered.GRCh38.bed" ;\
export PADDED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Padded.GRCh38.bed" ;\
export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.bam"; \
bsub \
     -w "done(\"${RGBASE}_s02vldate\")" \
     -J "${RGBASE}_s03intsct" \
     -o "${BASE}/WXS_Aquilla/04-LOGS/${RGBASE}_s03intsct.%J"\
     -u "${EMAIL}" \
     -n1 -W 1360 \
     -q research-hpc \
     -a "docker(vifehe/intersectbedref)" \
     entrypoint.sh; \
echo -e "04 - Starting validatesam per sample ${FULLSM} and RGBASE ${RGBASE}\nBAMFILE:${BAMFILE}"; \
export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.isec-${RUN_TYPE}.bam"; \
bsub \
     -w "done(\"${RGBASE}_s03intsct\")" \
     -J "${RGBASE}_s04vldate" \
     -o "${BASE}/WXS_Aquilla/04-LOGS/${RGBASE}_s04vldate.%J"\
     -u "${EMAIL}" \
     -n1 -W 1360 \
     -R "rusage[mem=48000]" \
     -q research-hpc \
     -a "docker(vifehe/validatesamref)" \
     entrypoint.sh -IGNORE INVALID_VERSION_NUMBER -IGNORE INVALID_TAG_NM -IGNORE MATE_NOT_FOUND; \
echo -e "05 - Starting bamtocram per sample ${FULLSM} and RGBASE ${RGBASE}\nBAMFILE:${BAMFILE}\nOUT_DIR:${OUT_DIR}"; \
# export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.bam"; \
export OUT_DIR="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/"; \
sleep ${DELAY2:-10}s; \
done; \
echo -e "DONE ANALYZING RGBASES per sample ${FULLSM}"; \
IFS=$'\n' RGBASES=($(grep "${FULLSM}" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_RGBASE})); \
INPUT_LIST=(); \
WAIT_LIST=(); \
CLEANUP_LIST=();\
for RGBASE in ${RGBASES[@]}; do \
   INPUT_LIST+=("${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.isec-${RUN_TYPE}.bam"); \
   WAIT_LIST+=("&&" "done(\"${RGBASE}_s04vldate\")"); \
    CLEANUP_LIST+=("${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.${FQEXT}"); \
  # CLEANUP_LIST+=("${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.${FQ2EXT}");\
done; \
echo -e "06 - Starting markduplicates for FULLSM ${FULLSM} with INPUT_LIST:${INPUT_LIST[@]}"; \
export FULLSM="${FULLSM}"; \
export MEM=38; \
export OUT_DIR="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/"; \
bsub \
   -w ${WAIT_LIST[@]:1} \
   -J "${FULLSM}_s06mrkdup" \
   -o "${BASE}/WXS_Aquilla/04-LOGS/${FULLSM}_s06mrkdup.%J" \
   -u "${EMAIL}" \
   -n1 -W 1440 \
   -M 46000000 \
   -R "rusage[mem=49152]" \
   -q research-hpc \
   -a "docker(vifehe/markduplicates)" \
   ${INPUT_LIST[@]}; \
  CLEANUP_LIST+=("${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.isec-${RUN_TYPE}.markd.bam"); \
echo -e "07 - Starting validatesam for FULLSM ${FULLSM}\nBAMFILE=${BAMFILE}"; \
export MEM=32; \
export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.bam"; \
export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
bsub \
   -w "done(\"${FULLSM}_s06mrkdup\")" \
   -J "${FULLSM}_s07vldate" \
   -o "${BASE}/WXS_Aquilla/04-LOGS/${FULLSM}_s07vldate.%J" \
   -u "${EMAIL}" \
   -n1 -W 1440 \
   -R "rusage[mem=46000]" \
   -q research-hpc \
   -a "docker(vifehe/validatesamref)" \
   entrypoint.sh -IGNORE INVALID_TAG_NM -IGNORE MATE_NOT_FOUND; \
echo -e "08 - Starting depthofcoverage for FULLSM ${FULLSM}\nBAMFILE=${BAMFILE}"; \
unset MEM; \
export RUN_TYPE="paddedexome"; \
export BEDFILE="${BASE}/Genome_Ref/GRCh38/Capture_Padded.GRCh38.bed"; \
export COVERED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Covered.GRCh38.bed"; \
export PADDED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Padded.GRCh38.bed"; \
export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
export GATK_THREADS="4"; \
bsub \
   -w "done(\"${FULLSM}_s07vldate\")" \
   -J "${FULLSM}_s08dcvrge" \
   -o "${BASE}/WXS_Aquilla/04-LOGS/${FULLSM}_s08dcvrge.%J" \
   -u "${EMAIL}" \
   -n1 -W 1440 \
   -q research-hpc \
   -a "docker(vifehe/depthofcoverage)" \
   entrypoint.sh; \
echo -e "09 - Starting baserecal for FULLSM ${FULLSM}\nBAMFILE=${BAMFILE}"; \
export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.bam"; \
export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
export DBSNP1="${BASE}/Genome_Ref/GRCh38/20190522_bundle/dbsnp_138.hg38.vcf.gz" ;\
export ONEKGP1="${BASE}/Genome_Ref/GRCh38/20190522_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz" ;\
export MILLS_GOLD="${BASE}/Genome_Ref/GRCh38/20190522_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" ;\
bsub \
   -w "done(\"${FULLSM}_s07vldate\")" \
   -J "${FULLSM}_s09bsercl" \
   -o "${BASE}/WXS_Aquilla/04-LOGS/${FULLSM}_s09bsercl.%J" \
   -u "${EMAIL}" \
   -n1 -W 1440 \
   -q research-hpc \
   -a "docker(vifehe/baserecalv2)" \
   entrypoint.sh; \
CLEANUP_LIST+=("${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.isec-${RUN_TYPE}.markd.recal.bam"); \
echo -e "10 - Starting validatesam for FULLSM ${FULLSM}\nBAMFILE ${BAMFILE}"; \
export MEM=32;\
export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.recal.bam"; \
export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
bsub \
   -w "done(\"${FULLSM}_s09bsercl\")" \
   -J "${FULLSM}_s10vldate" \
   -o "${BASE}/WXS_Aquilla/04-LOGS/${FULLSM}_s10vldate.%J" \
   -u "${EMAIL}" \
   -n1 -W 1360 \
   -R "rusage[mem=46000]" \
   -q research-hpc \
   -a "docker(vifehe/validatesamref)" \
   entrypoint.sh -IGNORE INVALID_VERSION_NUMBER -IGNORE INVALID_TAG_NM -IGNORE MATE_NOT_FOUND; \
echo -e "11 - Starting haplotype caller for FULLSM ${FULLSM}\nBAMFILE=${BAMFILE}"; \
export MEM=32; \
export RUN_TYPE="paddedexome"; \
export DBSNP1="${BASE}/Genome_Ref/GRCh38/20190522_bundle/dbsnp_138.hg38.vcf.gz"; \
export OUTDIR="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/";\
export PADDED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Padded.GRCh38.bed"; \
export COVERED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Covered.GRCh38.bed"; \
export BEDFILE="${BASE}/Genome_Ref/GRCh38/Capture_Padded.GRCh38.bed"; \
export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
bsub \
   -w "done(\"${FULLSM}_s10vldate\")" \
   -J "${FULLSM}_s11hcallr" \
   -o "${BASE}/WXS_Aquilla/04-LOGS/${FULLSM}_s11hcallr.%J" \
   -u "${EMAIL}" \
   -n1 -W 2880 \
   -M 46000000 \
   -R "rusage[mem=49152]" \
   -q research-hpc \
   -a "docker(vifehe/haplocallerv2)" \
   entrypoint.sh; \
echo -e "12 - Starting varianteval caller for FULLSM ${FULLSM}\GVCF=${GVCF}"; \
export MEM=32; \
export GVCF="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.recal.raw.snps.indels.g.vcf.gz"; \
export DBSNP1="${BASE}/Genome_Ref/GRCh38/20190522_bundle/dbsnp_138.hg38.vcf.gz"; \
export CCDS="/gscmnt/gc2645/wgs/Genome_Ref/GRCh38/GRCh38.CCDS.exons.sorted.bed"; \
export CLEANUP="${CLEANUP_LIST[@]}"; \
bsub \
   -w "done(\"${FULLSM}_s11hcallr\")" \
   -J "${FULLSM}_s12vnteval" \
   -o "${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}_s12vnteval.%J" \
   -u "${EMAIL}" \
   -n1 -W 2880 \
   -q research-hpc \
   -a "docker(vifehe/variantevalref)" \
   entrypoint.sh; \
sleep ${DELAY3:-10}s; \
echo "Finished submitting jobs for ${FULLSM}"; \
done; \
echo "Finished all samples within ${PR}"


WORKLIST="${BASE}/WXS_Aquilla/01-RAW/${PR}-worklist_fullsm.csv"
for x in $(cat ${WORKLIST}); do ls -lh ../02-TRANSIT/${PR}/${x}/*.vcf.gz; done

# kill jobs that faile at the intermediate steps
# bjobs -l| grep phs000572_201612 | grep Job | cut -d'<' -f2| cut -d'>' -f1 |  xargs bkill 

### CLEANING UP FASTQ files for ${PR}
# 1 - generate a list of gvcf files finished, just keep the FULLSMID info
ls ../02-TRANSIT/${PR}/*/*.g.vcf.gz | cut -d "/" -f 4 | sort | uniq > ${PR}-fqtodelete.list
ls ../02-TRANSIT/${PR}/*/*.g.vcf.gz | cut -d "/" -f 4 | sort | uniq > ${PR}-fqtodelete.list2

cat ${PR}-fqtodelete.list2 ${PR}-fqtodelete.list | sort | uniq -c | awk '{if ($1==1) print $2}' > ${PR}-fqtodelete.LIST

## delete all fastq & cram within that FULLSMID directory for the samples that gvcf is completed
for x in $(cat ${PR}-fqtodelete.list); do rm ${PR}/${x}/*.fq.gz; done
for x in $(cat ${PR}-fqtodelete.list); do rm ../02-TRANSIT/${PR}/${x}/*.bam; done




ls ../02-TRANSIT/${PR}/*/*.bam | wc -l
ls ./${PR}/*/*.fq.gz | wc -l

ls ../02-TRANSIT/${PR}/*/*.aln.srt.bam | wc -l
265
ls ../02-TRANSIT/${PR}/*/*.aln.srt.isec-paddedexome.markd.bam | wc -l
32
ls ../02-TRANSIT/${PR}/*/*.aln.srt.isec-*.markd.bam | wc -l
32
ls ../02-TRANSIT/${PR}/*/*.aln.srt.isec-*.markd.recal.table1 | wc -l
32
ls ../02-TRANSIT/${PR}/*/*.aln.srt.isec-*.markd.exome-coverage.sample_statistics | wc -l
31
ls ../02-TRANSIT/${PR}/*/*.aln.srt.isec-*.markd.exome-coverage.sample_summary | wc -l
31
ls ../02-TRANSIT/${PR}/*/*.aln.srt.isec-*.markd.recal.raw.snps.indels.g.vcf.gz | wc -l
82

### CLEAR UP SOME FILES:
fernandezv@virtual-workstation1:/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW$ rm ./${PR}/*/*.fastq


