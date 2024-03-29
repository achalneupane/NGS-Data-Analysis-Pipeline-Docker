### 2019-08 - README file for Genentech-WGS sequences processing over at MGI

##############
## 2019-08-05
##############

### In Fenix

# generate metadata file for the PROJECT in question
HOMEDIR="/40/AD/AD_Seq_Data/01.-RawData"
PR="Genentech_WGS"

##List directories within the PR directory
## Macrogen like DIAN have several bams that need to be merged
metafile="${PR}-sm_dna_pr_rgbase_fq1ext_fq2ext.csv"
vim ${metafile}   ## copy here the metadata for the DIAN project from the excel file
sed -i 's/\t/,/g' ${metafile}

## Log onto MGI
USER="achal"

`RANDOM=$$
NUMBER=$(((${RANDOM}%5)+1))
alias sshmgi="ssh ${USER}@virtual-workstation${NUMBER}.gsc.wustl.edu"`

sshmgi


## move to directory where all processing samples are:
cd /gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW
#   This step can be automated for dealing with 100s of samples
#   Ideally harmonized symlink structure is already setup for each project, which can be the src for scp -rp
BASE="/gscmnt/gc2645/wgs"
PR="Genentech_WGS"
DEST="${BASE}/WXS_Aquilla/01-RAW"
cd ${DEST}
mkdir -p ${PR}

# Copy over the metafile
# metafile="${PR}-sm_dna_pr_rgbase_fq1ext_fq2ext.csv"
metafile="${PR}-MGIID-SMID-DNA-PR-metafile_achal.csv"
HOMEDIR="/40/AD/AD_Seq_Data/01.-RawData"
scp -r ${USER}@fenix.psych.wucon.wustl.edu:${HOMEDIR}/${metafile} ${DEST}

#Generate the workfile
cut -d, -f1-3 ${DEST}/${metafile} | tr "," "^" | sort -u > ${DEST}/${PR}-worklist_fullsm.csv
workfile="${DEST}/${PR}-worklist_fullsm.csv"

# generate a directory for each one of the samples within the project 
for FULLSM in $(cat ${workfile}); do mkdir -p ${PR}/${FULLSM}; done

## transfer RGFILES from original directory to directory at MGI
while read -r line;  do 
COUNTER=$((COUNTER+1))
FULLSM="$(echo ${line} | cut -d, -f1-3 | tr "," "^")"
SM="$(echo "${FULLSM}" | cut -d^ -f1)"
DNA="$(echo "${FULLSM}" | cut -d^ -f2)"
PR="$(echo "${FULLSM}" | cut -d^ -f3)"
RGBASE="$(echo ${line} | cut -d, -f4)"
FQ1EXT="$(echo ${line} | cut -d, -f5)"
FQ2EXT="$(echo ${line} | cut -d, -f6)"
FBASE="/80/AD/AD_Seq_Data/01.-RawData/Genentech/${PR}/${FULLSM}"
FILE_RG="${FBASE}/${RGBASE}.rgfile"
FILE_FQ1="${FBASE}/${RGBASE}.${FQ1EXT}"
FILE_FQ2="${FBASE}/${RGBASE}.${FQ2EXT}"
echo -e "$FULLSM\n$SM\n$DNA\n$PR\n$RGBASE\n$FQ1EXT\n$FQ2EXT\n$FBASE\n$FILE_RG\n$FILE_FQ1\n$FILE_FQ2"
#done < ${metafile}
#scp -p ${USER}@fenix.psych.wucon.wustl.edu:${FILE_RG} ${DEST}/${PR}/${FULLSM}/
# rsync -avh ${USER}@fenix.psych.wucon.wustl.edu:${FILE_FQ1} ${DEST}/${PR}/${FULLSM}/
# rsync -avh ${USER}@fenix.psych.wucon.wustl.edu:${FILE_FQ2} ${DEST}/${PR}/${FULLSM}/
done < ${metafile}

# already done gvcf
  ./*/*vcf.gz | cut -d/ -f2 > /gscmnt/gc2645/wgs/WXS_Aquilla/02-TRANSIT/Genentech_WGS/already_done.txt
/gscmnt/gc2645/wgs/WXS_Aquilla/02-TRANSIT/Genentech_WGS/already_done.txt


grep -wv -F -f /gscmnt/gc2645/wgs/WXS_Aquilla/02-TRANSIT/Genentech_WGS/already_done.txt  "${BASE}/WXS_Aquilla/01-RAW/${PR}-worklist_fullsm.csv" > "${BASE}/WXS_Aquilla/01-RAW/${PR}-worklist_fullsm_missing.csv"


# difference between WORKLIS and already done ..., we create a new worklist

##############
## 2019-10-22
##############

## 1 - check that the number of uploaded files matches what we have on the metafile
'fernandezv@virtual-workstation1:/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW$ ls Genentech_WGS/*/*.rgfile | wc -l
332
fernandezv@virtual-workstation1:/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW$ ls Genentech_WGS/*/*.gz | wc -l
664
fernandezv@virtual-workstation1:/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW$ ls Genentech_WGS/*/*r1.fq.gz | wc -l
332'
## 2 - set the correct STAR Tand END for the workfile
fernandezv@virtual-workstation1:/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW$ wc -l ${WORKLIST}
47 /gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/Genentech_WGS-worklist_fullsm.csv


bjobs | grep RUN | cut -d\  -f1 | xargs bkill 
bjobs | grep PEN | cut -d\  -f1 | xargs bkill 
bjobs | grep SSUSP | cut -d\  -f1 | xargs bkill 

## We'll make 4 batches of 83 - batch 1/4

#### Start processign fq files
START=1
END=9
DELAY=10
EMAIL="fernandezv@wustl.edu"
SHELLDROP=0
Snumber=1
PR="Genentech_WGS"
export BASE="/gscmnt/gc2645/wgs"; \
export WORKDIR="${BASE}/tmp"; \
export THREADS=4; \
export BWA_PIPE_SORT=1; \
export TIMING=1; \
WORKLIST="${BASE}/WXS_Aquilla/01-RAW/${PR}-worklist_fullsm_missing.csv"
# LOOKUP="${BASE}/WXS_Aquilla/01-RAW/${PR}-sm_dna_pr_rgbase_fq1ext_fq2ext.csv"
LOOKUP="${BASE}/WXS_Aquilla/01-RAW/${PR}-MGIID-SMID-DNA-PR-metafile_achal.csv"
export BASE="/gscmnt/gc2645/wgs"; \
export WORKDIR="${BASE}/tmp"; \
export THREADS=16; \
export BWA_PIPE_SORT=1; \
export TIMING=1; \
LOOKUP_COL_SM=1; \
LOOKUP_COL_DNA=2; \
LOOKUP_COL_PR=3; \
LOOKUP_COL_RGBASE=4; \
LOOKUP_COL_FQ1EXT=5; \
LOOKUP_COL_FQ2EXT=6; \
export RUN_TYPE="paddedexome"; \
for FULLSM in $(sed -n "${START},${END}p" "${WORKLIST}"); do \
  echo "*****************************************************"; \
  echo "*****************************************************"; \
  echo "Doing sample number**********: " $Snumber; \
  echo "*****************************************************"; \
  echo "*****************************************************"; \
  ((Snumber=${Snumber}+1)); \
FBASE="/40/AD/AD_Seq_Data/01.-RawData/Genentech/${PR}/${FULLSM}"; \
SM="$(echo "${FULLSM}" | cut -d^ -f1)"; \
IFS=$'\n' export DNA=($(awk -F, "\$${LOOKUP_COL_SM} == \"${SM}\"" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_DNA} | sort -u)); \
if [ ${#DNA[@]} -gt 1 ]; then echo "Warning, \${DNA} not unique for ${SM} (n=${#DNA[@]}: ${DNA[@]})"; fi; \
DNA="$(echo "${FULLSM}" | cut -d^ -f2)"; \
PR="$(echo "${FULLSM}" | cut -d^ -f3)"; \
export OUT_DIR="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}"; \
echo -e "00 - Starting jobs per sample FULLSM ${FULLSM}"; \
for RGBASE in $(grep "${SM},${DNA},${PR}" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_RGBASE}); do \
SM="$(echo "${FULLSM}" | cut -d^ -f1)";\
DNA="$(echo "${FULLSM}" | cut -d^ -f2)";\
PR="$(echo "${FULLSM}" | cut -d^ -f3)";\
FQ1EXT=($(awk -F, "\$${LOOKUP_COL_RGBASE} == \"${RGBASE}\"" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_FQ1EXT})); \
FQ2EXT=($(awk -F, "\$${LOOKUP_COL_RGBASE} == \"${RGBASE}\"" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_FQ2EXT})); \
RGFILE="${FBASE}/${RGBASE}.rgfile"; \
FQ1="${FBASE}/${RGBASE}.${FQ1EXT}"; \
FQ2="${FBASE}/${RGBASE}.${FQ2EXT}"; \
DEST="${BASE}/WXS_Aquilla/01-RAW"; \
echo -e "00a - Uploading FASTQ and rgfiles for sample ${FULLSM} and RGBASE ${RGBASE}\nFQ1:${FQ1}\nFQ2:${FQ2}\nRGFILE:${RGFILE}\nDNA:${DNA}"; \
mkdir ${DEST}/${PR}/${FULLSM}; \
# rsync -avh -L ${USER}@fenix.psych.wucon.wustl.edu:${FQ1} ${DEST}/${PR}/${FULLSM}/; \
# rsync -avh -L ${USER}@fenix.psych.wucon.wustl.edu:${FQ2} ${DEST}/${PR}/${FULLSM}/; \
# rsync -avh ${USER}@fenix.psych.wucon.wustl.edu:${RGFILE} ${DEST}/${PR}/${FULLSM}/; \
export RGFILE="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.rgfile"; \
 if [ "$(echo -e $(cat ${RGFILE}) | grep -oP "(?<=SM:)[^[:space:]]*")" != "${SM}" ]; then \
  echo -e "Error, SM ${SM} mismatch between worklist and ${RGFILE}"; continue; \
 fi; \
 if [ "$(echo -e $(cat ${RGFILE}) | grep -oP "(?<=DS:)[^[:space:]]*")" != "${FULLSM}" ]; then \
  echo -e "Error, FULLSM ${FULLSM} mismatch between worklist and ${RGFILE}"; continue; \
 fi; \
 if [ "${FULLSM}.$(echo -e $(cat ${RGFILE}) | grep -oP "(?<=PU:)[^[:space:]]*" | tr ":" "^")" != "${RGBASE}" ]; then \
  echo -e "Error, RGBASE ${RGBASE} mismatch between LOOKUP and ${RGFILE}"; continue; \
 fi; \
export FULLSM_RGID="${RGBASE}"; \
unset BAMFILE; \
IFS=$'\n' FQ1EXT=($(awk -F, "\$${LOOKUP_COL_RGBASE} == \"${RGBASE}\"" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_FQ1EXT})); \
IFS=$'\n' FQ2EXT=($(awk -F, "\$${LOOKUP_COL_RGBASE} == \"${RGBASE}\"" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_FQ2EXT})); \
echo -e "01 - Starting bwa per sample ${FULLSM} and RGBASE ${RGBASE}";\
export RGFILE="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.rgfile"; \
export FQ1="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.${FQ1EXT}"; \
export FQ2="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.${FQ2EXT}"; \
export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
echo "FQ1 FILE present here: ######## $(ls -lht $FQ1)"; \
echo "FQ2 FILE present here: ######## $(ls -lht $FQ2)"; \
echo "FILE_RG FILE present here: ######## $(ls -lht $RGFILE)"; \
export MEM=50; \
# done; \
# done
   bsub \
     -J "${RGBASE}_s01alnsrt" \
     -o "${BASE}/WXS_Aquilla/04-LOGS/${RGBASE}_s01alnsrt.%J" \
     -u "${EMAIL}" \
     -n ${THREADS} -W 2880 \
     -M 69000000 \
     -R "rusage[mem=69152]" \
     -q research-hpc \
     -a 'docker(vifehe/bwaref)' \
     entrypoint.sh; \
export MEM=16; \
export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.bam"; \
echo -e "02 - Starting validatesam per sample ${FULLSM} and RGBASE ${RGBASE}\nBAMFILE:${BAMFILE}"; \
echo ${BAMFILE}; \
export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
bsub \
     -w "done(\"${RGBASE}_s01alnsrt\")" \
     -J "${RGBASE}_s02vldate" \
     -o "${BASE}/WXS_Aquilla/04-LOGS/${RGBASE}_s02vldate.%J"\
     -u "${EMAIL}" \
     -n1 -W 1360 \
     -R "rusage[mem=24000]" \
     -q research-hpc \
     -a "docker(vifehe/validatesamref)" \
     entrypoint.sh -IGNORE INVALID_VERSION_NUMBER -IGNORE INVALID_TAG_NM; \
echo -e "03 - Starting intersectbed per sample ${FULLSM} and RGBASE ${RGBASE}\nBAMFILE:${BAMFILE}"; \
export RUN_TYPE="paddedexome"; \
export BEDFILE="${BASE}/Genome_Ref/GRCh38/Capture_Padded.GRCh38.bed"; \
export COVERED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Covered.GRCh38.bed" ; \
export PADDED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Padded.GRCh38.bed" ; \
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
echo -e "04 - Starting validatesam per sample ${FULLSM} and RGBASE ${RGBASE}\nBAMFILE:${BAMFILE}";\
export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.isec-${RUN_TYPE}.bam"; \
bsub \
     -w "done(\"${RGBASE}_s03intsct\")" \
     -J "${RGBASE}_s04vldate" \
     -o "${BASE}/WXS_Aquilla/04-LOGS/${RGBASE}_s04vldate.%J"\
     -u "${EMAIL}" \
     -n1 -W 1360 \
     -R "rusage[mem=24000]" \
     -q research-hpc \
     -a "docker(vifehe/validatesamref)" \
     entrypoint.sh -IGNORE INVALID_VERSION_NUMBER -IGNORE INVALID_TAG_NM -IGNORE MATE_NOT_FOUND; \
echo -e "05 - Starting bamtocram per sample ${FULLSM} and RGBASE ${RGBASE}\nBAMFILE:${BAMFILE}\nOUT_DIR:${OUT_DIR}"; \
# export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.bam"; \
export OUT_DIR="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/"; \
done; \
echo -e "DONE ANALYZING RGBASES per sample ${FULLSM}";\
IFS=$'\n' RGBASES=($(grep "${FULLSM}" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_RGBASE})); \
INPUT_LIST=(); \
WAIT_LIST=(); \
CLEANUP_LIST=();\
for RGBASE in ${RGBASES[@]}; do \
   INPUT_LIST+=("${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.isec-${RUN_TYPE}.bam"); \
   WAIT_LIST+=("&&" "done(\"${RGBASE}_s04vldate\")"); \
    CLEANUP_LIST+=("${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}${FQ1EXT}"); \
  CLEANUP_LIST+=("${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}${FQ2EXT}"); \
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
  CLEANUP_LIST+=("${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.isec-${RUN_TYPE}.markd.bam"); \
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
   -R "rusage[mem=16000]" \
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
export DBSNP1="${BASE}/Genome_Ref/GRCh38/20190522_bundle/dbsnp_138.hg38.vcf.gz"; \
export ONEKGP1="${BASE}/Genome_Ref/GRCh38/20190522_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz"; \
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
export MEM=16; \
export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.recal.bam"; \
export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
bsub \
   -w "done(\"${FULLSM}_s09bsercl\")" \
   -J "${FULLSM}_s10vldate" \
   -o "${BASE}/WXS_Aquilla/04-LOGS/${FULLSM}_s10vldate.%J" \
   -u "${EMAIL}" \
   -n1 -W 360 \
   -R "rusage[mem=16000]" \
   -q research-hpc \
   -a "docker(vifehe/validatesamref)" \
   entrypoint.sh -IGNORE INVALID_VERSION_NUMBER -IGNORE INVALID_TAG_NM -IGNORE MATE_NOT_FOUND; \
echo -e "11 - Starting haplotype caller for FULLSM ${FULLSM}\nBAMFILE=${BAMFILE}"; \
export MEM=32; \
export RUN_TYPE="paddedexome"; \
export DBSNP1="${BASE}/Genome_Ref/GRCh38/20190522_bundle/dbsnp_138.hg38.vcf.gz"; \
export OUTDIR="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/"; \
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
export MEM=16; \
export GVCF="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.recal.raw.snps.indels.g.vcf.gz"; \
export DBSNP1="${BASE}/Genome_Ref/GRCh38/20190522_bundle/dbsnp_138.hg38.vcf.gz"; \
export CCDS="/gscmnt/gc2645/wgs/Genome_Ref/GRCh38/GRCh38.CCDS.exons.sorted.bed"; \
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
sleep ${DELAY2:-10}s; \
echo "Finished submitting jobs for ${FULLSM}"; \
done; \
echo "Finished all samples within ${PR}"; 




### CLEANING UP FASTQ files for DIAN_REDCLOUD
# 1 - generate a list of gvcf files finished, just keep the FULLSMID info
PR="Genentech_WGS"
ls ../02-TRANSIT/${PR}/*/*.g.vcf.gz | cut -d "/" -f 4 | sort | uniq > ${PR}-fqtodelete.list

#ls ../02-TRANSIT/${PR}/*/*.g.vcf.gz | cut -d "/" -f 4 | sort | uniq > ${PR}-fqtodelete.list2
#cat 201812_MGI_DIANWGS_REDCLOUD-fqtodelete.list2 201812_MGI_DIANWGS_REDCLOUD-fqtodelete.list | sort | uniq -c | awk '{if ($1==1) print $2}' > 201812_MGI_DIANWGS_REDCLOUD-fqtodelete.LIST

## delete all fastq & cram within that FULLSMID directory for the samples that gvcf is completed
for x in $(cat ${PR}-fqtodelete.list); do rm ${PR}/${x}/*.fq.gz;  done
for x in $(cat ${PR}-fqtodelete.LIST); do rm ${PR}/${x}/*.final.cram; done


for FULLSM in $(cat ${WORKLIST}); do \
#for FULLSM in $(sed -n "${START},${END}p" "${WORKLIST}"); do \
if [ -f ${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}.*.g.vcf.gz ]; then \
echo -e "GVCF file for sample ${FULLSM} is complete"; \
rm ${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/*.fq.gz; \
else echo -e "GVCF file for sample ${FULLSM} is NOT complete"; fi; done;

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@ PENDING
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
##############
## 2019-09-11
##############

PR="Genentech_WGS"
### Check how many ouitputs there are of each kind, and whether that is the expected:
ls ../02-TRANSIT/${PR}/*/*.isec-paddedexome.bam  | wc -l

ls ../02-TRANSIT/${PR}/*/*.cram | wc -l

ls ../02-TRANSIT/${PR}/*/*.bam | wc -l

ls ../02-TRANSIT/${PR}/*/*aln.srt.isec-paddedexome.bam | wc -l

ls ../03-FINAL/${PR}/*/*.cram | wc -l

ls 02-TRANSIT/${PR}/*/*.vcf.gz | wc -l

###Also, originally there are  134 RGfiles, suggesting that 4 smaples did not process

##Do a loop of some sort that for each GRFILE chekcs whether the expected bam, padded-exome and cram files re three

PR="Genentech_WES"
echo -e "${PR}_FULLSM_RG_sample\tPADDEDM-BAM\tCRAM\tBAM" > OUT1a
while read -r line ; do \
FULLSM=`echo ${line} | cut -d "," -f 7`;\
RG=`echo ${line} | cut -d "," -f 4`;\
echo -e "${FULLSM}/${RG}" >> SAMPLES ;   \
if [ -f "/gscmnt/gc2645/wgs/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RG}.aln.srt.isec-paddedexome.bam" ]; then \
echo -e "YES" >> PADDEDBAM; else \
echo -e "NO" >> PADDEDBAM; fi; \
if [ -f "/gscmnt/gc2645/wgs/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RG}.aln.srt.cram" ]; then \
echo -e "YES" >> CRAM; else \
echo -e "NO" >> CRAM; fi; \
if [ -f "/gscmnt/gc2645/wgs/WXS_Aquilla/03-FINAL/${PR}/${FULLSM}/${RG}.aln.srt.bam" ]; then \
echo -e "YES" >> BAM ; else \
echo -e "NO">> BAM; fi;
done < ${PR}-sm_dna_pr_rgbase_fq1ext_fq2ext.csv
paste SAMPLES PADDEDBAM CRAM BAM > OUT1b
cat OUT1a OUT1b > ${PR}-sm_dna_pr_rgbase_fq1ext_fq2ext.OUT1	


###################### MERGE gVCFs togeher
  export CCDS="${BASE}/GATK_pipeline/Reference/GRCh83/GRCh38.CCDS.exons.sorted.bed" ;\