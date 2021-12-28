### Oct-5-2020 - README file for PPMI_WES sequences processing over at MGI 

##############
## Oct-5-2020
##############

### In Fenix
# These samples are not propely labelled- have multiple fastq and 246 samples have multiple runs. Some samples are also not concurrent in terms of IBD with GWAS samples. Multi-run samples also have cryptic IBD so they need to be checked before we merge them.    
# for f in *R1*; do
# ls | grep ${f/R1/R2}
# done

cd /40/Cruchaga_Data/PD/3-GeneticData/4-WES/WholeExome_RawData_March2017
# ls -l /80/AD/AD_Seq_Data/02.-BWA_GATK/02-gVCFs/GRCh37/singlerun/PPMI_3432^multi^PPMI_WES
# locate PPMI_3432 | grep ".g.vcf.gz$"
# cd /80/AD/AD_Seq_Data/02.-BWA_GATK/02-gVCFs/GRCh37/singlerun
# PPMI IBD check
# cd /30/PD/3-GeneticData/4-WES/WholeExome_RawData_March2017/ready/ppmi_ibd
# cat ppmi_ibd-PPMI_full_imputed_gwas_201605_GRCh37_IBD.genomels -lht
# cat /30/PD/3-GeneticData/4-WES/WholeExome_RawData_March2017/ready/ppmi_ibd/merge_run_files.list

### 2019-08 - README file for DIAN_WGS sequences processing over at MGI

##############
## 2019-08-02
##############

### In Fenix
HOMEDIR="/40/AD/AD_Seq_Data/01.-RawData"
PR="PPMI_WES"
metafile="${PR}_sm_dna-pr_rgbase_fq1ext_fq2ext_fullsm.csv"
vim ${metafile}   ## copy here the metadata for the DIAN project from the excel file
sed -i 's/\t/,/g' ${metafile}

## Log onto MGI
USER="achal"
## move to directory where all processing samples are:
cd /gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW
#   This step can be automated for dealing with 100s of samples
#   Ideally harmonized symlink structure is already setup for each project, which can be the src for scp -rp
BASE="/gscmnt/gc2645/wgs"
PR="PPMI_WES"
DEST="${BASE}/WXS_Aquilla/01-RAW"
mkdir -p ${PR}

USER="achal"
# Copy over the metafile
# scp -r ${USER}@fenix.psych.wucon.wustl.edu:${HOMEDIR}/${metafile} ${DEST}
dir="/home/achal/01.-RawData"
scp -r ${USER}@fenix.psych.wucon.wustl.edu:${dir}/${metafile} ${DEST}


EMAIL="achal@wustl.edu"
USER="achal"
## UPLOAD FQ1 FQ2 RGFILE, run samples and delete
BASE="/gscmnt/gc2645/wgs"
PR="PPMI_WES"

DEST="${BASE}/WXS_Aquilla/01-RAW"


LOOKUP="${BASE}/WXS_Aquilla/01-RAW/${PR}_sm_dna-pr_rgbase_fq1ext_fq2ext_fullsm.csv"
#Generate the workfile
cut -d, -f1-3 ${LOOKUP} | tr "," "^" | sort -u > ${DEST}/${PR}-worklist_fullsm.csv
WORKLIST="${DEST}/${PR}-worklist_fullsm.csv"


START=1
END=918

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
FBASE="/40/Cruchaga_Data/PD/3-GeneticData/4-WES/WholeExome_RawData_March2017/ready/${FULLSM}" \
FQ1EXT=($(awk -F, "\$${LOOKUP_COL_RGBASE} == \"${RGBASE}\"" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_FQ1EXT})); \
FQ2EXT=($(awk -F, "\$${LOOKUP_COL_RGBASE} == \"${RGBASE}\"" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_FQ2EXT})); \
RGFILE="${FBASE}/${RGBASE}.rgfile";\
FQ1="${FBASE}/${RGBASE}.${FQ1EXT}";\
FQ2="${FBASE}/${RGBASE}.${FQ2EXT}";\
DEST="${BASE}/WXS_Aquilla/01-RAW";\
echo -e "00a - Uploading FASTQ and rgfiles for sample ${FULLSM} and RGBASE ${RGBASE}\nFQ1:${FQ1}\nFQ2:${FQ2}\nRGFILE:${RGFILE}\nDNA:${DNA}";\
mkdir ${DEST}/${PR}/${FULLSM};\
# scp -p ${USER}@fenix.psych.wucon.wustl.edu:${FILE_RG} ${DEST}/${PR}/${FULLSM}/ ;\
# rsync with copy referent files as we are copying symlinks for this PPMI data from fenix
rsync -avh -L ${USER}@fenix.psych.wucon.wustl.edu:${FQ1} ${DEST}/${PR}/${FULLSM}/ ;\
rsync -avh -L ${USER}@fenix.psych.wucon.wustl.edu:${FQ2} ${DEST}/${PR}/${FULLSM}/ ;\
rsync -avh ${USER}@fenix.psych.wucon.wustl.edu:${RGFILE} ${DEST}/${PR}/${FULLSM}/ ;\
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
export MEM=32; \
   bsub \
     -J "${RGBASE}_s01alnsrt" \
     -o "${BASE}/WXS_Aquilla/04-LOGS/${RGBASE}_s01alnsrt.%J" \
     -u "${EMAIL}" \
     -n ${THREADS} -W 2160 \
     -M 46000000 \
     -R "rusage[mem=49152]" \
     -q research-hpc \
     -a 'docker(vifehe/bwaref)' \
     entrypoint.sh; \
echo -e "02 - Starting validatesam per sample ${FULLSM} and RGBASE ${RGBASE}\nBAMFILE:${BAMFILE}";
export MEM=6; \
export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.bam"; \
echo ${BAMFILE}; \
export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
bsub \
     -w "done(\"${RGBASE}_s01alnsrt\")" \
     -J "${RGBASE}_s02vldate" \
     -o "${BASE}/WXS_Aquilla/04-LOGS/${RGBASE}_s02vldate.%J"\
     -u "${EMAIL}" \
     -n1 -W 360 \
     -R "rusage[mem=8192]" \
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
     -n1 -W 360 \
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
     -n1 -W 360 \
     -R "rusage[mem=8192]" \
     -q research-hpc \
     -a "docker(vifehe/validatesamref)" \
     entrypoint.sh -IGNORE INVALID_VERSION_NUMBER -IGNORE INVALID_TAG_NM -IGNORE MATE_NOT_FOUND; \
echo -e "05 - Starting bamtocram per sample ${FULLSM} and RGBASE ${RGBASE}\nBAMFILE:${BAMFILE}\nOUT_DIR:${OUT_DIR}";\
export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.bam"; \
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
    CLEANUP_LIST+=("${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.${FQ1EXT}");\
  CLEANUP_LIST+=("${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.${FQ2EXT}");\
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
sleep ${DELAY2:-10}s; \
echo "Finished submitting jobs for ${FULLSM}"; \
done; \
echo "Finished all samples within ${PR}"; 


WORKLIST="${BASE}/WXS_Aquilla/01-RAW/${PR}-worklist_fullsm.csv"
for x in $(cat ${WORKLIST}); do ls -lh ../02-TRANSIT/PPMI_WES/${x}/*.vcf.gz; done



####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
###################################### Merge multi run samples: WORK from BAM LIST for merging  #####################################################
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
# # Low IBD (multi-run samples; did not merge them)
# IID1  IID2  RT  EZ  Z0  Z1  Z2  PI_HAT  PHE DST PPC RATIO IBS0  IBS1  IBS2  HOMHOM  HETHET  Not to be merged
# PPMI_3205^SI500-SI504 PPMI_3205^SI610-SI613 UN  NA  0 0.3406  0.6594  0.8297  -1  0.926585  1 NA  0 132 767 0 119 
# PPMI_3362^SI1045-SI1063 PPMI_3362^SI1657-SI1687 UN  NA  0.1274  0.233 0.6396  0.7561  -1  0.9069  1 18.1429 11  148 754 7 127 
# PPMI_3362^SI1657-SI1687 PPMI_3362^SI604-SI616 UN  NA  0 0.3601  0.6399  0.8199  -1  0.922367  1 NA  0 143 778 0 141 PPMI_3362^SI1657-SI1687
# PPMI_3364^SI1674-SI1579 PPMI_3364^SI916-SI887 UN  NA  0 0.2765  0.7235  0.8618  -1  0.940401  1 NA  0 113 835 0 146 PPMI_3364^SI1674-SI1579
# PPMI_3431^SI1652-SI1701 PPMI_3431^SI602-SI620 UN  NA  0.0118  0.3665  0.6217  0.805 -1  0.917038  1 125 1 147 750 1 125 
# PPMI_3523^SI1418-SI1317 PPMI_3523^SI1831-SI1805 UN  NA  0 0.5094  0.4906  0.7453  -1  0.890199  1 NA  0 177 629 0 94  PPMI_3523^SI1418-SI1317
# PPMI_3621^SI1340-SI1392 PPMI_3621^SI599-SI615 UN  NA  0.0123  0.4236  0.5641  0.7759  -1  0.90454 1 112 1 162 696 1 112 

# #samples with two runs; Removed the conflicts and were not required to merge after removing one run
# PPMI_3803^SI699-SI681 GWAS_PPMI_3904  UN  NA  0.0856  0.6504  0.2639  0.5891  -1  0.830972  0.9882  7.3333  4 159 331 3 22  Conflict
# PPMI_3807^SI715-SI710 GWAS_PPMI_3702  UN  NA  0.0685  0.8098  0.1217  0.5266  -1  0.802376  0.9848  7 3 177 283 3 21  Conflict
# PPMI_3812^SI697-SI698 GWAS_PPMI_3781  UN  NA  0.031 0.8767  0.0923  0.5306  -1  0.800587  0.956 10  1 134 206 1 10  Conflict
# PPMI_3857^SI677-SI704 GWAS_PPMI_3775  UN  NA  0 0.775 0.225 0.6125  -1  0.832944  0.9982  NA  0 143 285 0 17  Conflict

# # Sample Conflicts (Single runs; left them as is); Need to discuss if we should really eliminate them before doing downstream analysis
# PPMI_3773 PPMI_3784 UN  NA  0.2311  0.3778  0.3911  0.58  -1  0.840791  1 11.25 21  264 676 12  135 Conflict
# PPMI_3784 GWAS_PPMI_3773  UN  NA  0.2751  0.3163  0.4086  0.5667  -1  0.83923 1 8.9333  25  259 677 15  134 Conflict



cd /gscmnt/gc2645/wgs/WXS_Aquilla/02-TRANSIT/PPMI_WES
ls | grep PPMI_ > /Merged_PPMI_WES/sample_runs.list
cat sample_runs.list | cut -d^ -f1 | uniq -c | sort > sample_runs_list_counts.txt
sed -e 's/[\t ]//g;/^$/d'  sample_runs_list_counts.txt | sed -e 's/PPMI/,PPMI/g' > clean_sample_runs_list_counts.txt
grep -v ^1 clean_sample_runs_list_counts.txt > clean_multi_run_list.txt
# single run
grep ^1 clean_sample_runs_list_counts.txt > clean_single_run_list.txt
cat clean_multi_run_list.txt |  cut -d, -f2 > multi_run_list_to_be_merged.txt
# extract bam list to be merged
ls /gscmnt/gc2645/wgs/WXS_Aquilla/02-TRANSIT/PPMI_WES/*/*.isec-paddedexome.markd.bam > PPMI_WES_ALL_unmerged_paddedexome.markd.BAM_FILES.txt


grep -Ff multi_run_list_to_be_merged.txt PPMI_WES_ALL_unmerged_paddedexome.markd.BAM_FILES.txt > PPMI_WES_filtered_to_be_merged_paddedexome.markd_BAM_FILES.txt
# grep -v -Ff  PPMI_WES_filtered_to_be_merged_paddedexome.markd_BAM_FILES.txt
# remove IBD unmatched list
grep -v -Ff PPMI_Do_NOT_merge.txt  PPMI_WES_filtered_to_be_merged_paddedexome.markd_BAM_FILES.txt > PPMI_WES_filtered_to_be_merged_paddedexome.markd.BAM_FILES_Final.txt
# create metafile
cd /gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW
rm PPMI_WES_filtered_to_be_merged_paddedexome.markd.BAM_FILES_Final_FULLSM.csv
for lines in $(cat /gscmnt/gc2645/wgs/WXS_Aquilla/02-TRANSIT/PPMI_WES/Merged_PPMI_WES/PPMI_WES_filtered_to_be_merged_paddedexome.markd.BAM_FILES_Final.txt); do
  samples="$(basename ${lines}| cut -d^ -f1)";
  # echo "Samples: ${samples}"
  DNA="$(basename ${lines}| cut -d^ -f2)";
  BAMNAME="$(basename ${lines})";
  # echo "DNA: ${DNA}"
  PR="PPMI_WES";
  # MULTI is going to be DNA here for sample labels
  MULTI="multi"
  FULLSM="${samples}^multi^${PR}";
  RGBASE="${samples}^$DNA^${PR}";
  echo "${FULLSM},${samples},${MULTI},${PR},${DNA},${lines},${BAMNAME},${RGBASE}" >> PPMI_WES_filtered_to_be_merged_paddedexome.markd.BAM_FILES_Final_FULLSM.csv
done

# remove lines matching patterns (removing IBD conflicts from BAM-metafile to merge the right sample runs only).
sed '/PPMI_3803\|PPMI_3807\|PPMI_3812\|PPMI_3857/d' PPMI_WES_filtered_to_be_merged_paddedexome.markd.BAM_FILES_Final_FULLSM.csv > temp_file
mv temp_file PPMI_WES_filtered_to_be_merged_paddedexome.markd.BAM_FILES_Final_FULLSM.csv


EMAIL="achal@wustl.edu"
USER="achal"
## UPLOAD FQ1 FQ2 RGFILE, run samples and delete
BASE="/gscmnt/gc2645/wgs"
PR="PPMI_WES"
DEST="${BASE}/WXS_Aquilla/01-RAW"

LOOKUP="${DEST}/PPMI_WES_filtered_to_be_merged_paddedexome.markd.BAM_FILES_Final_FULLSM.csv"
cut -d, -f1 ${LOOKUP} | tr "," "^" | sort -u > ${DEST}/${PR}-worklist_fullsm.csv
WORKLIST="${DEST}/${PR}-worklist_fullsm.csv"


#####save all files just in case we need them; Will delete SAVECOPY later
# cp -R PPMI_300*/ SAVECOPY/


START=1
END=242
Snumber=1
# for FULLSM in $(sed -n "${START},${END}p" "${WORKLIST}"); do \
export BASE="/gscmnt/gc2645/wgs"; \
export WORKDIR="${BASE}/tmp"; \
export THREADS=16; \
export BWA_PIPE_SORT=1; \
export DEST="${BASE}/WXS_Aquilla/01-RAW"; \
export TIMING=1; \
LOOKUP_COL_SM=2; \
LOOKUP_COL_DNA=3; \
LOOKUP_COL_PR=4; \
LOOKUP_COL_BAR=5; \
LOOKUP_COL_BAM=7; \
LOOKUP_COL_RGBASE=8; \
export RUN_TYPE="paddedexome"; \
for FULLSM in $(sed -n "${START},${END}p" "${WORKLIST}"); do \
  echo "*****************************************************"; \
  echo "*****************************************************"; \
  echo "Doing sample number**********: " $Snumber; \
  echo "*****************************************************"; \
  echo "*****************************************************"; \
  ((Snumber=${Snumber}+1)); \
SM="$(echo "${FULLSM}" | cut -d^ -f1)"; \
IFS=$'\n' export DNA=($(awk -F, "\$${LOOKUP_COL_SM} == \"${SM}\"" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_DNA} | sort -u)); \
# if [ ${#DNA[@]} -gt 1 ]; then echo "Warning, \${DNA} not unique for ${SM} (n=${#DNA[@]}: ${DNA[@]})"; fi; \
DNA="$(echo "${FULLSM}" | cut -d^ -f2)"; \
PR="PPMI_WES"; \
export OUT_DIR="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}"; \
echo -e "00 - Starting jobs per sample FULLSM ${FULLSM}"; \
for RGBASE in $(grep "${SM},${DNA},${PR}" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_RGBASE}); do \
DEST="${BASE}/WXS_Aquilla/01-RAW";\
echo ${RGBASE}; \
IFS=$'\n' RGBASES=($(grep "${FULLSM}" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_RGBASE})); \
INPUT_LIST=(); \
WAIT_LIST=(); \
CLEANUP_LIST=();\
for RGBASE in ${RGBASES[@]}; do \
INPUT_LIST+=("${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${RGBASE/_MERGED/}/${RGBASE/_MERGED/}.aln.srt.isec-${RUN_TYPE}.markd.bam"); \
sleep ${DELAY1:-1}s; \
done; done; \
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
  CLEANUP_LIST+=("${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${RGBASE/_MERGED/}/${RGBASE/_MERGED/}.aln.srt.isec-${RUN_TYPE}.markd.bam");\
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
sleep ${DELAY2:-10}s; \
echo "Finished submitting jobs for ${FULLSM}"; \
done; \
echo "Finished all samples within ${PR}"; 




cd /gscmnt/gc2645/wgs/WXS_Aquilla/02-TRANSIT/PPMI_WES
# Now remove unwanted folders
for file in $(cat /gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/PPMI_WES-worklist_fullsm.csv); do 
getFile="$(echo $file | cut -d^ -f1)"
DIR="$(ls | grep ${getFile}^SI)"
echo $DIR
rm -rf ${DIR}
done


# remove all PPMI Bams and fastqs
rm ./*PPMI_WES/*.bam

cd /gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/PPMI_WES
rm ./*PPMI_WES/*.fq.gz 
ls ./*PPMI_WES/*.fq.gz | wc -l

ls ../02-TRANSIT/PPMI_WES/*/*.bam | wc -l
0
ls ./PPMI_WES/*/*.fq.gz | wc -l
0
ls ../02-TRANSIT/PPMI_WES/*/*.aln.srt.bam | wc -l
0
ls ../02-TRANSIT/PPMI_WES/*/*.aln.srt.isec-paddedexome.markd.bam | wc -l
0
ls ../02-TRANSIT/PPMI_WES/*/*.aln.srt.isec-*.markd.bam | wc -l
0
ls ../02-TRANSIT/PPMI_WES/*/*.aln.srt.isec-*.markd.recal.table1 | wc -l
591
ls ../02-TRANSIT/PPMI_WES/*/*.aln.srt.isec-*.markd.exome-coverage.sample_statistics | wc -l
591
ls ../02-TRANSIT/PPMI_WES/*/*.aln.srt.isec-*.markd.exome-coverage.sample_summary | wc -l
591
ls ../02-TRANSIT/PPMI_WES/*/*.aln.srt.isec-*.markd.recal.raw.snps.indels.g.vcf.gz | wc -l
591


# # Check the difference in folders on Fenix (John's) and MGI (Achal's)
# # On Fenix
# cd /home/achal/01.-RawData/PPMI_WES
# sdiff -l all_samples_after_merging_john.txt all_samples_after_merging_Achal.txt

# PPMI_3751^SI1566-SI1536^PPMI_WES                              | PPMI_3751^multi^PPMI_WES
# PPMI_3752^SI1554-SI1563^PPMI_WES                              | PPMI_3752^multi^PPMI_WES
# PPMI_3754^SI1573-SI1543^PPMI_WES                              | PPMI_3754^multi^PPMI_WES
# PPMI_3804^SI1568-SI1533^PPMI_WES                              | PPMI_3804^multi^PPMI_WES
# PPMI_3806^SI1542-SI1540^PPMI_WES                              | PPMI_3806^multi^PPMI_WES
# PPMI_3808^SI1570-SI1560^PPMI_WES                              | PPMI_3808^multi^PPMI_WES
# PPMI_3811^SI1537-SI1562^PPMI_WES                              | PPMI_3811^multi^PPMI_WES