# Updated on 9/20/2020 by Achal. This script was used to process ADSP samples from $PR: phs000572_201802
# First use phs000572_201802.sh script to download ADSP SRAs, then use Create_metafile.sh to create metafile for this project. 


unset all ENV variables
for c in $(env | cut -d '=' -f 1); do unset $c; done

cd /30/dbGaP/6109/sra/phs000572_201802/fqgz2
# ls -l /80/AD/AD_Seq_Data/02.-BWA_GATK/02-gVCFs/GRCh37/singlerun/PPMI_3432^multi^PPMI_WES
# locate PPMI_3432 | grep ".g.vcf.gz$"
# cd /80/AD/AD_Seq_Data/02.-BWA_GATK/02-gVCFs/GRCh37/singlerun


###################On MGI
HOMEDIR="/40/AD/AD_Seq_Data/01.-RawData"
PR="phs000572_201802"
# remove lines matching patterns
# sed '/SRR1391349\|SRR1411264/d' phs000572_201802_sm_dna-pr_rgbase_FQext_fullsm_1-250.csv > "${PR}_sm_dna-pr_rgbase_FQext_fullsm.csv"

# cat phs000572_201802_sm_dna-pr_rgbase_FQext_fullsm_1-250.csv > "${PR}_sm_dna-pr_rgbase_FQext_fullsm.csv"
# cat phs000572_201802_sm_dna-pr_rgbase_FQext_fullsm_251-500.csv > "${PR}_sm_dna-pr_rgbase_FQext_fullsm.csv"
# cat phs000572_201802_sm_dna-pr_rgbase_FQext_fullsm_501-750.csv > "${PR}_sm_dna-pr_rgbase_FQext_fullsm.csv"
# cat phs000572_201802_sm_dna-pr_rgbase_FQext_fullsm_751-1000.csv > "${PR}_sm_dna-pr_rgbase_FQext_fullsm.csv"
# cat phs000572_201802_sm_dna-pr_rgbase_FQext_fullsm_1001-1250.csv > "${PR}_sm_dna-pr_rgbase_FQext_fullsm.csv"
# cat phs000572_201802_sm_dna-pr_rgbase_FQext_fullsm_1251-1500.csv > "${PR}_sm_dna-pr_rgbase_FQext_fullsm.csv"
# cat phs000572_201802_sm_dna-pr_rgbase_FQext_fullsm_1501-1540_plus_missing.csv > "${PR}_sm_dna-pr_rgbase_FQext_fullsm.csv"
cat phs000572_201802_sm_dna-pr_rgbase_FQext_fullsm_last_4_missing_samples.csv > "${PR}_sm_dna-pr_rgbase_FQext_fullsm.csv"
# ## get only lines matching patterns
# cat phs000572_201802_sm_dna-pr_rgbase_FQext_fullsm_1-250.csv | grep "SRR1391349\|SRR1411264"
##############
## 2019-08-02
##############

### In Fenix

metafile="${PR}_sm_dna-pr_rgbase_FQext_fullsm.csv"
# vim ${metafile}   ## copy here the metadata for the DIAN project from the excel file
# sed -i 's/\t/,/g' ${metafile}


## Log onto MGI
USER="achal"
BASE="/gscmnt/gc2645/wgs"
PR="phs000572_201802"
DEST="${BASE}/WXS_Aquilla/01-RAW"
mkdir -p ${PR}

# Copy over the metafile
# scp -r ${USER}@fenix.psych.wucon.wustl.edu:${HOMEDIR}/${metafile} ${DEST}
dir="/home/achal/01.-RawData"
# scp -r ${USER}@fenix.psych.wucon.wustl.edu:${dir}/${metafile} ${DEST}
EMAIL="achal@wustl.edu"
USER="achal"
## UPLOAD FQ RGFILE, run samples and delete
BASE="/gscmnt/gc2645/wgs"
PR="phs000572_201802"

DEST="${BASE}/WXS_Aquilla/01-RAW"

# ##############################

LOOKUP="${BASE}/WXS_Aquilla/01-RAW/${PR}_sm_dna-pr_rgbase_FQext_fullsm.csv"

#Generate the workfile
cut -d, -f1-3 $LOOKUP | tr "," "^" | sort -u > ${DEST}/${PR}-worklist_fullsm.csv
# WORKLISTtmp="${DEST}/${PR}-worklist_fullsm.csv"
WORKLIST="${DEST}/${PR}-worklist_fullsm.csv"

#################################################################################


# bjobs | grep RUN | cut -d\  -f1 | xargs bkill 
# bjobs | grep PEN | cut -d\  -f1 | xargs bkill 
# bjobs | grep SSUSP | cut -d\  -f1 | xargs bkill 

#1-250;  120 position missing A-ADC-AD001631^SRR1411947^phs000572_201802
# 251-500; A-ADC-AD003216^SRR1428657^phs000572_201802; A-ADC-AD004072^SRR1358891^phs000572_201802; A-ADC-AD005075^SRR1412783^phs000572_201802; 
# A-ADC-AD005199^SRR1439166^phs000572_201802; A-ADC-AD006720^SRR1578753^phs000572_201802; A-ADC-AD006969^SRR1347639^phs000572_201802
#501-750; A-ADC-AD007682^SRR1451532^phs000572_201802; A-LOAD-LD001077^SRR1423602^phs000572_201802; 
# 751-1000; A-VAN-VU000741^SRR1210956^phs000572_201802

# # 13 missing
# A-VAN-VU000741^SRR1210956^phs000572_201802
# C-ASPS-30098^SRR1533483^phs000572_201802
# C-FHS-51057^SRR1391127^phs000572_201802
# C-RS-50510^SRR1127287^phs000572_201802
# C-RS-50759^SRR1450256^phs000572_201802
# A-ADC-AD007682^SRR1451532^phs000572_201802
# A-ADC-AD005199^SRR1439166^phs000572_201802
# A-ADC-AD006720^SRR1578753^phs000572_201802
# A-ADC-AD006969^SRR1347639^phs000572_201802
# A-ADC-AD003216^SRR1428657^phs000572_201802
# A-ADC-AD004072^SRR1358891^phs000572_201802
# A-ADC-AD005075^SRR1412783^phs000572_201802
# A-ADC-AD001631^SRR1411947^phs000572_201802


EMAIL="achal@wustl.edu";


LOOKUP="${BASE}/WXS_Aquilla/01-RAW/${PR}_sm_dna-pr_rgbase_FQext_fullsm.csv";
WORKLIST="${DEST}/${PR}-worklist_fullsm.csv";


START=1;
END=4;
SHELLDROP=1
# for FULLSM in $(sed -n "${START},${END}p" "${WORKLIST}"); do \
Snumber=1;
export BASE="/gscmnt/gc2645/wgs"; \
export WORKDIR="${BASE}/tmp"; \
export THREADS=4; \
export BWA_PIPE_SORT=1; \
export TIMING=1; \
LOOKUP_COL_SM=1; \
LOOKUP_COL_DNA=2; \
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
FBASE="/30/dbGaP/6109/sra/phs000572_201802/fqgz/${FULLSM}"; \
FQEXT=($(awk -F, "\$${LOOKUP_COL_RGBASE} == \"${RGBASE}\"" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_FQEXT})); \
# FQ2EXT=($(awk -F, "\$${LOOKUP_COL_RGBASE} == \"${RGBASE}\"" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_FQ2EXT})); \
RGFILE="${FBASE}/${RGBASE}.rgfile";\
FQ="${FBASE}/${RGBASE}.${FQEXT}";\
# FQ2="${FBASE}/${RGBASE}.${FQ2EXT}";\
DEST="${BASE}/WXS_Aquilla/01-RAW";\
echo -e "00a - Uploading FASTQ and rgfiles for sample ${FULLSM} and RGBASE ${RGBASE}\nFQ:${FQ}\nRGFILE:${RGFILE}\nDNA:${DNA}";\
# echo $FILE_RG
# echo $FILE_FQ1
# echo $FILE_FQ2
# done
# done
# done
mkdir ${DEST}/${PR}/${FULLSM}; \
# rsync with copy referent files as we are copying symlinks for this PPMI data from fenix
rsync -avh -L ${USER}@fenix.psych.wucon.wustl.edu:${FQ} ${DEST}/${PR}/${FULLSM}/ ;\
rsync -avh ${USER}@fenix.psych.wucon.wustl.edu:${RGFILE} ${DEST}/${PR}/${FULLSM}/ ;\
# done
# done
# done
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
export MEM=50; \
bsub \
     -J "${RGBASE}_s01alnsrt" \
     -o "${BASE}/WXS_Aquilla/04-LOGS/${RGBASE}_s01alnsrt.%J" \
     -u "${EMAIL}" \
     -n ${THREADS} -W 2160 \
     -M 66000000 \
     -R "rusage[mem=69152]" \
     -q research-hpc \
     -a 'docker(vifehe/bwaref)' \
     entrypoint.sh; \
echo -e "02 - Starting validatesam per sample ${FULLSM} and RGBASE ${RGBASE}\nBAMFILE:${BAMFILE}";
export MEM=16; \
export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.bam"; \
echo ${BAMFILE}; \
export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
bsub \
     -w "done(\"${RGBASE}_s01alnsrt\")" \
     -J "${RGBASE}_s02vldate" \
     -o "${BASE}/WXS_Aquilla/04-LOGS/${RGBASE}_s02vldate.%J"\
     -u "${EMAIL}" \
     -n1 -W 1360 \
     -R "rusage[mem=18192]" \
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
     -R "rusage[mem=18192]" \
     -q research-hpc \
     -a "docker(vifehe/validatesamref)" \
     entrypoint.sh -IGNORE INVALID_VERSION_NUMBER -IGNORE INVALID_TAG_NM -IGNORE MATE_NOT_FOUND; \
# echo -e "05 - Starting bamtocram per sample ${FULLSM} and RGBASE ${RGBASE}\nBAMFILE:${BAMFILE}\nOUT_DIR:${OUT_DIR}";\
# export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.bam"; \
export OUT_DIR="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/"; \
sleep ${DELAY2:-10}s; \
done; \
echo -e "DONE ANALYZING RGBASES per sample ${FULLSM}";\
IFS=$'\n' RGBASES=($(grep "${FULLSM}" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_RGBASE})); \
INPUT_LIST=(); \
WAIT_LIST=(); \
CLEANUP_LIST=();\
for RGBASE in ${RGBASES[@]}; do \
   INPUT_LIST+=("${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.isec-${RUN_TYPE}.bam"); \
   WAIT_LIST+=("&&" "done(\"${RGBASE}_s04vldate\")"); \
    CLEANUP_LIST+=("${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.${FQEXT}");\
  # CLEANUP_LIST+=("${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.${FQ2EXT}");\
done; \
echo -e "06 - Starting markduplicates for FULLSM ${FULLSM} with INPUT_LIST:${INPUT_LIST[@]}"; \
export FULLSM="${FULLSM}"; \
export MEM=32; \
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
  CLEANUP_LIST+=("${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.isec-${RUN_TYPE}.markd.bam");\
echo -e "07 - Starting validatesam for FULLSM ${FULLSM}\nBAMFILE=${BAMFILE}"; \
export MEM=16; \
export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.bam"; \
export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
bsub \
   -w "done(\"${FULLSM}_s06mrkdup\")" \
   -J "${FULLSM}_s07vldate" \
   -o "${BASE}/WXS_Aquilla/04-LOGS/${FULLSM}_s07vldate.%J" \
   -u "${EMAIL}" \
   -n1 -W 1440 \
   -R "rusage[mem=18192]" \
   -q research-hpc \
   -a "docker(vifehe/validatesamref)" \
   entrypoint.sh -IGNORE INVALID_TAG_NM -IGNORE MATE_NOT_FOUND; \
echo -e "08 - Starting depthofcoverage for FULLSM ${FULLSM}\nBAMFILE=${BAMFILE}"; \
unset MEM; \
export RUN_TYPE="paddedexome"; \
export BEDFILE="${BASE}/Genome_Ref/GRCh38/Capture_Padded.GRCh38.bed"; \
export COVERED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Covered.GRCh38.bed" ;\
export PADDED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Padded.GRCh38.bed" ;\
export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
export GATK_THREADS="4";\
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
CLEANUP_LIST+=("${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.isec-${RUN_TYPE}.markd.recal.bam");\
echo -e "10 - Starting validatesam for FULLSM ${FULLSM}\nBAMFILE ${BAMFILE}";\
export MEM=16;\
export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.recal.bam"; \
export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
bsub \
   -w "done(\"${FULLSM}_s09bsercl\")" \
   -J "${FULLSM}_s10vldate" \
   -o "${BASE}/WXS_Aquilla/04-LOGS/${FULLSM}_s10vldate.%J" \
   -u "${EMAIL}" \
   -n1 -W 1360 \
   -R "rusage[mem=18192]" \
   -q research-hpc \
   -a "docker(vifehe/validatesamref)" \
   entrypoint.sh -IGNORE INVALID_VERSION_NUMBER -IGNORE INVALID_TAG_NM -IGNORE MATE_NOT_FOUND; \
echo -e "11 - Starting haplotype caller for FULLSM ${FULLSM}\nBAMFILE=${BAMFILE}"; \
export MEM=32; \
export RUN_TYPE="paddedexome"; \
export DBSNP1="${BASE}/Genome_Ref/GRCh38/20190522_bundle/dbsnp_138.hg38.vcf.gz" ;\
export OUTDIR="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/";\
export PADDED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Padded.GRCh38.bed" ;\
export COVERED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Covered.GRCh38.bed" ;\
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
export MEM=2; \
export GVCF="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.recal.raw.snps.indels.g.vcf.gz"; \
export DBSNP1="${BASE}/Genome_Ref/GRCh38/20190522_bundle/dbsnp_138.hg38.vcf.gz" ;\
export CCDS="/gscmnt/gc2645/wgs/Genome_Ref/GRCh38/GRCh38.CCDS.exons.sorted.bed";\
export CLEANUP="${CLEANUP_LIST[@]}"
bsub \
   -w "done(\"${FULLSM}_s11hcallr\")" \
   -J "${FULLSM}_s12vnteval" \
   -o "${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}_s12vnteval.%J" \
   -u "${EMAIL}" \
   -n1 -W 1360\
   -q research-hpc \
   -a "docker(vifehe/variantevalref)" \
   entrypoint.sh; \
sleep ${DELAY3:-10}s; \
echo "Finished submitting jobs for ${FULLSM}"; \
done; \
echo "Finished all samples within ${PR}"


WORKLIST="${BASE}/WXS_Aquilla/01-RAW/${PR}-worklist_fullsm.csv"
for x in $(cat ${WORKLIST}); do ls -lh ../02-TRANSIT/phs000572_201802/${x}/*.vcf.gz; done

# missing; 1, 2, 4, 8
A-ADC-AD007682^SRR1451532^phs000572_201802
A-ADC-AD005075^SRR1412783^phs000572_201802
A-ADC-AD003216^SRR1428657^phs000572_201802
A-ADC-AD001631^SRR1411947^phs000572_201802
### CLEANING UP FASTQ files for phs000572_201802
# 1 - generate a list of gvcf files finished, just keep the FULLSMID info
ls ../02-TRANSIT/phs000572_201802/*/*.g.vcf.gz | cut -d "/" -f 4 | sort | uniq > phs000572_201802-fqtodelete.list
ls ../02-TRANSIT/phs000572_201802/*/*.g.vcf.gz | cut -d "/" -f 4 | sort | uniq > phs000572_201802-fqtodelete.list2

cat phs000572_201802-fqtodelete.list2 phs000572_201802-fqtodelete.list | sort | uniq -c | awk '{if ($1==1) print $2}' > phs000572_201802-fqtodelete.LIST

## delete all fastq & cram within that FULLSMID directory for the samples that gvcf is completed
for x in $(cat phs000572_201802-fqtodelete.list); do rm phs000572_201802/${x}/*.fq.gz; done
for x in $(cat phs000572_201802-fqtodelete.list); do rm ../02-TRANSIT/phs000572_201802/${x}/*.bam; done




ls ../02-TRANSIT/phs000572_201802/*/*.bam | wc -l
ls ./phs000572_201802/*/*.fq.gz | wc -l

ls ../02-TRANSIT/phs000572_201802/*/*.aln.srt.bam | wc -l
265
ls ../02-TRANSIT/phs000572_201802/*/*.aln.srt.isec-paddedexome.markd.bam | wc -l
32
ls ../02-TRANSIT/phs000572_201802/*/*.aln.srt.isec-*.markd.bam | wc -l
32
ls ../02-TRANSIT/phs000572_201802/*/*.aln.srt.isec-*.markd.recal.table1 | wc -l
32
ls ../02-TRANSIT/phs000572_201802/*/*.aln.srt.isec-*.markd.exome-coverage.sample_statistics | wc -l
31
ls ../02-TRANSIT/phs000572_201802/*/*.aln.srt.isec-*.markd.exome-coverage.sample_summary | wc -l
31
ls ../02-TRANSIT/phs000572_201802/*/*.aln.srt.isec-*.markd.recal.raw.snps.indels.g.vcf.gz | wc -l
82

### CLEAR UP SOME FILES:
fernandezv@virtual-workstation1:/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW$ rm ./phs000572_201802/*/*.fastq





###################################################################################################
###################################################################################################
###################################################################################################
################ 4 samples keep failing at Mark duplicate step, so completing them here ###########
###################################################################################################
###################################################################################################
###################################################################################################

bjobs | grep RUN | cut -d\  -f1 | xargs bkill 
bjobs | grep PEN | cut -d\  -f1 | xargs bkill 
# bjobs | grep SSUSP | cut -d\  -f1 | xargs bkill 

# missing; 1, 2, 4, 8
A-ADC-AD007682^SRR1451532^phs000572_201802
A-ADC-AD005075^SRR1412783^phs000572_201802
A-ADC-AD003216^SRR1428657^phs000572_201802
A-ADC-AD001631^SRR1411947^phs000572_201802


START=1;
END=4;
SHELLDROP=0
# for FULLSM in $(sed -n "${START},${END}p" "${WORKLIST}"); do \
Snumber=1;
export BASE="/gscmnt/gc2645/wgs"; \
export WORKDIR="${BASE}/tmp"; \
export THREADS=4; \
export BWA_PIPE_SORT=1; \
export TIMING=1; \
LOOKUP_COL_SM=1; \
LOOKUP_COL_DNA=2; \
LOOKUP_COL_PR=3; \
LOOKUP_COL_RGBASE=4; \
LOOKUP_COL_FQEXT=5; \
unset FBASE; \
export RUN_TYPE="paddedexome"; \
export FASTQ_TYPE="interleaved"; \
export MODE="fq"; \
for FULLSM in $(sed -n "${START},${END}p" "${WORKLIST}"); do \
  echo "*****************************************************"; \
  echo "*****************************************************"; \
  echo "Doing sample number**********: " ${Snumber}; \
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
FBASE="/30/dbGaP/6109/sra/phs000572_201802/fqgz/${FULLSM}"; \
FQEXT=($(awk -F, "\$${LOOKUP_COL_RGBASE} == \"${RGBASE}\"" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_FQEXT})); \
RGFILE="${FBASE}/${RGBASE}.rgfile";\
FQ="${FBASE}/${RGBASE}.${FQEXT}";\
DEST="${BASE}/WXS_Aquilla/01-RAW";\
echo -e "00a - Uploading FASTQ and rgfiles for sample ${FULLSM} and RGBASE ${RGBASE}\nFQ:${FQ}\nRGFILE:${RGFILE}\nDNA:${DNA}";\
done; \
echo -e "DONE ANALYZING RGBASES per sample ${FULLSM}";\
IFS=$'\n' RGBASES=($(grep "${FULLSM}" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_RGBASE})); \
INPUT_LIST=(); \
WAIT_LIST=(); \
CLEANUP_LIST=();\
for RGBASE in ${RGBASES[@]}; do \
INPUT_LIST+=("${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.isec-${RUN_TYPE}.bam"); \
CLEANUP_LIST+=("${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.${FQEXT}");\
sleep ${DELAY2:-1}s; \
done; \
sleep ${DELAY2:-10}s; \
echo -e "06 - Starting markduplicates for FULLSM ${FULLSM} with INPUT_LIST:${INPUT_LIST[@]}"; \
export FULLSM="${FULLSM}"; \
export MEM=32; \
export OUT_DIR="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/"; \
bsub \
  -J "${FULLSM}_s06mrkdup" \
  -o "${BASE}/WXS_Aquilla/04-LOGS/${FULLSM}_s06mrkdup.%J" \
  -u "${EMAIL}" \
  -n1 -W 1440 \
  -M 46000000 \
  -R "rusage[mem=49152]" \
  -q research-hpc \
  -a "docker(vifehe/markduplicates)" \
${INPUT_LIST[@]}; \
CLEANUP_LIST+=("${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.isec-${RUN_TYPE}.markd.bam");\
echo -e "07 - Starting validatesam for FULLSM ${FULLSM}\nBAMFILE=${BAMFILE}"; \
export MEM=6; \
export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.bam"; \
export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
bsub \
   -w "done(\"${FULLSM}_s06mrkdup\")" \
   -J "${FULLSM}_s07vldate" \
   -o "${BASE}/WXS_Aquilla/04-LOGS/${FULLSM}_s07vldate.%J" \
   -u "${EMAIL}" \
   -n1 -W 1440 \
   -R "rusage[mem=8192]" \
   -q research-hpc \
   -a "docker(vifehe/validatesamref)" \
   entrypoint.sh -IGNORE INVALID_TAG_NM -IGNORE MATE_NOT_FOUND; \
echo -e "08 - Starting depthofcoverage for FULLSM ${FULLSM}\nBAMFILE=${BAMFILE}"; \
unset MEM; \
export RUN_TYPE="paddedexome"; \
export BEDFILE="${BASE}/Genome_Ref/GRCh38/Capture_Padded.GRCh38.bed"; \
export COVERED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Covered.GRCh38.bed" ;\
export PADDED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Padded.GRCh38.bed" ;\
export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
export GATK_THREADS="4";\
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
CLEANUP_LIST+=("${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.isec-${RUN_TYPE}.markd.recal.bam");\
echo -e "10 - Starting validatesam for FULLSM ${FULLSM}\nBAMFILE ${BAMFILE}";\
export MEM=6;\
export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.recal.bam"; \
export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
bsub \
   -w "done(\"${FULLSM}_s09bsercl\")" \
   -J "${FULLSM}_s10vldate" \
   -o "${BASE}/WXS_Aquilla/04-LOGS/${FULLSM}_s10vldate.%J" \
   -u "${EMAIL}" \
   -n1 -W 360 \
   -R "rusage[mem=8192]" \
   -q research-hpc \
   -a "docker(vifehe/validatesamref)" \
   entrypoint.sh -IGNORE INVALID_VERSION_NUMBER -IGNORE INVALID_TAG_NM -IGNORE MATE_NOT_FOUND; \
echo -e "11 - Starting haplotype caller for FULLSM ${FULLSM}\nBAMFILE=${BAMFILE}"; \
export MEM=32; \
export RUN_TYPE="paddedexome"; \
export DBSNP1="${BASE}/Genome_Ref/GRCh38/20190522_bundle/dbsnp_138.hg38.vcf.gz" ;\
export OUTDIR="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/";\
export PADDED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Padded.GRCh38.bed" ;\
export COVERED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Covered.GRCh38.bed" ;\
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
export MEM=2; \
export GVCF="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.recal.raw.snps.indels.g.vcf.gz"; \
export DBSNP1="${BASE}/Genome_Ref/GRCh38/20190522_bundle/dbsnp_138.hg38.vcf.gz" ;\
export CCDS="/gscmnt/gc2645/wgs/Genome_Ref/GRCh38/GRCh38.CCDS.exons.sorted.bed";\
export CLEANUP="${CLEANUP_LIST[@]}"
bsub \
   -w "done(\"${FULLSM}_s11hcallr\")" \
   -J "${FULLSM}_s12vnteval" \
   -o "${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}_s12vnteval.%J" \
   -u "${EMAIL}" \
   -n1 -W 360\
   -q research-hpc \
   -a "docker(vifehe/variantevalref)" \
   entrypoint.sh; \
sleep ${DELAY3:-10}s; \
echo "Finished submitting jobs for ${FULLSM}"; \
done; \
echo "Finished all samples within ${PR}"