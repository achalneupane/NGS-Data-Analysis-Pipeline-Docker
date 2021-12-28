# grep -wv -F -f //gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/Otogenetics_WES/Otogenetics_GVCFs_done.txt  "${BASE}/WXS_Aquilla/01-RAW/Otogenetics_WES-worklist_fullsm.csv" > "/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/${PR}-worklist_missing_achal.csv"


14_6053_206^8008287680^Otogenetics_WES
25_50_1420^8006302353^Otogenetics_WES
25_48_1356^8006302372^Otogenetics_WES
27_99_84605^8006302217^Otogenetics_WES
27_99_85407^8006302338^Otogenetics_WES
27_97_85580^8008288474^Otogenetics_WES
28_5_49^8008287082^Otogenetics_WES
4_653_3^8008285098^Otogenetics_WES
8_64079_3^8008286339^Otogenetics_WES
27_97_84596^8008288217^Otogenetics_WES
# # remove any white space
# sed -i 's/[\t ]//g;/^$/d' 

3,4,6,16,18,19

# First transfer fastq files to active directory
cp * /storage1/fs1/cruchagac/Active/achal/Otogenetics_WGS/ /storage1/fs1/cruchagac/Active/achal/Otogenetics_WGS

4H_158_24^8008291747^Otogenetics_WES
8_64079_3^8008286339^Otogenetics_WES
4_653_3^8008285098^Otogenetics_WES
28_5_49^8008287082^Otogenetics_WES
25_48_1356^8006302372^Otogenetics_WES
14_6053_206^8008287680^Otogenetics_WES
14_6053_204^8008287644^Otogenetics_WES



####
# These samples also needs to be processed again
25_57_1671^8008286673^Otogenetics_WES
25_9_172^8008287262^Otogenetics_WES
26_BDU_BDU98201^8008286227^Otogenetics_WES
25_60_1760^8008287272^Otogenetics_WES
25_64_1831^8008286156^Otogenetics_WES
26_BCR_BCR18850^8011737168^Otogenetics_WES
26_ARH_ARH05003^8008286085^Otogenetics_WES
25_53_1534z^8008286166^Otogenetics_WES
25_66_1907^09-024_E_E07^Otogenetics_WES
25_66_1901^8011737144^Otogenetics_WES
26_ARH_ARH05004^8008286129^Otogenetics_WES
25_9_171^8008287251^Otogenetics_WES
25_57_1678^8008286629^Otogenetics_WES
26_BCR_BCR18853^8011737102^Otogenetics_WES
25_64_1830^8008286177^Otogenetics_WES
25_53_1525^8008286165^Otogenetics_WES


# Adding  some more samples to the redo list
25_60_1761^8008286232^Otogenetics_WES
25_6_105^8008286210^Otogenetics_WES
25_6_107^8008286189^Otogenetics_WES
25_64_1832^8008286167^Otogenetics_WES
25_66_1902^8011737313^Otogenetics_WES
25_66_1903^8011737370^Otogenetics_WES
25_9_171C1^8008285636^Otogenetics_WES
26_ARH_ARH05002^8008286063^Otogenetics_WES
26_BCR_BCR18802^09-024_E_B07^Otogenetics_WES
26_BCR_BCR18844^8011736908^Otogenetics_WES
26_BCR_BCR18846^8011737074^Otogenetics_WES
26_BCR_BCR18864^8011737097^Otogenetics_WES
26_BDU_BDU98201^8008286227^Otogenetics_WES
26_CIH_CIH73203^8008293399^Otogenetics_WES
26_MTJ_MTJ88801^8008285325^Otogenetics_WES


# On Fenix
cd /storage1/fs1/cruchagac/Archive/Otogenetics_WES/fqgz
for line in $(cat /home/achal/01.-RawData/Otogenetics_WES/redo.txt); do
cd $line
cp *.fq.gz /storage1/fs1/cruchagac/Active/achal/Otogenetics_WGS/
cp *.rgfile /storage1/fs1/cruchagac/Active/achal/Otogenetics_WGS/
cd ..
done



# 26_BDU_BDU98201^8008286227^Otogenetics_WES
# [mem_sam_pe] paired reads have different names: "HISEQ2000B:573:C6027ACXX:1:2314:18", "HISEQ2000B:573:C6027ACXX:1:2314:18428:92556"

# Command exited with non-zero status 1
#         Command being timed: "bwa mem -t 4 -R @RG\tID:C6027ACXX:1\tPL:illumina\tPU:HISEQ2000B:573:C6027ACXX:1\tLB:8008286227\tSM:26_BDU_BDU98201\tDS:26_BDU_BDU98201^8008286227^Otogenetics_WES -M /tmp//Homo_sapiens_assembly38.fasta /tmp//26_BDU_BDU98201^8008286227^Otogenetics_WES.HISEQ2000B^573^C6027ACXX^1.r1.fq.gz /tmp//26_BDU_BDU98201^8008286227^Otogenetics_WES.HISEQ2000B^573^C6027ACXX^1.r2.fq.gz"




START=1
END=15
USER="achal"
export BASE="/gscmnt/gc2645/wgs"
export PR="Otogenetics_WES"
WORKLIST="${BASE}/WXS_Aquilla/01-RAW/${PR}-worklist_missing_achal.csv"; \
LOOKUP="${BASE}/WXS_Aquilla/01-RAW/${PR}-sm_dna_pr_rgbase_fq1ext_fq2ext.csv"; \
# sed -i 's/[\t ]//g;/^$/d' $LOOKUP; \
Snumber=1
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
RGBASE="$(echo ${RGBASE/ /})"; \
SM="$(echo "${FULLSM}" | cut -d^ -f1)"; \
DNA="$(echo "${FULLSM}" | cut -d^ -f2)"; \
PR="$(echo "${FULLSM}" | cut -d^ -f3)"; \
FBASE="/storage1/fs1/cruchagac/Active/achal/Otogenetics_WGS"; \
FQ1EXT=($(awk -F, "\$${LOOKUP_COL_RGBASE} == \"${RGBASE}\"" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_FQ1EXT})); \
FQ2EXT=($(awk -F, "\$${LOOKUP_COL_RGBASE} == \"${RGBASE}\"" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_FQ2EXT})); \  
# FQ1EXT="r1.fq.gz"; \
# FQ1EXT="r2.fq.gz"; \
unset FQ1; \
unset FQ2; \
unset RGFILE; \
RGFILE="${FBASE}/${RGBASE}.rgfile";\
FQ1="${FBASE}/${RGBASE}.${FQ1EXT}";\
FQ2="${FBASE}/${RGBASE}.${FQ2EXT}";\
DEST="${BASE}/WXS_Aquilla/01-RAW";\
echo -e "00a - Uploading FASTQ and rgfiles for sample ${FULLSM} and RGBASE ${RGBASE}\nFQ1:${FQ1}\nFQ2:${FQ2}\nRGFILE:${RGFILE}\nDNA:${DNA}";\
mkdir ${DEST}/${PR}/${FULLSM};\
rsync -avh -L ${USER}@compute1-client-1.ris.wustl.edu:${RGFILE} ${DEST}/${PR}/${FULLSM}/; \
rsync -avh -L ${USER}@compute1-client-1.ris.wustl.edu:${FQ1} ${DEST}/${PR}/${FULLSM}/; \
rsync -avh -L ${USER}@compute1-client-1.ris.wustl.edu:${FQ2} ${DEST}/${PR}/${FULLSM}/; \
# export RGFILE="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.rgfile"; \
export FULLSM_RGID="${RGBASE}"; \
unset BAMFILE; \
IFS=$'\n' FQ1EXT=($(awk -F, "\$${LOOKUP_COL_RGBASE} == \"${RGBASE}\"" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_FQ1EXT})); \
IFS=$'\n' FQ2EXT=($(awk -F, "\$${LOOKUP_COL_RGBASE} == \"${RGBASE}\"" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_FQ2EXT})); \
# FQ1EXT="r1.fq.gz"; \
# FQ1EXT="r2.fq.gz"; \
echo -e "01 - Starting bwa per sample ${FULLSM} and RGBASE ${RGBASE}";\
export RGFILE="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.rgfile"; \
export FQ1="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.${FQ1EXT}"; \
export FQ2="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.${FQ2EXT}"; \
 if [ "$(echo -e $(cat ${RGFILE}) | grep -oP "(?<=SM:)[^[:space:]]*")" != "${SM}" ]; then \
  echo -e "Error, SM ${SM} mismatch between worklist and ${RGFILE}"; continue; \
 fi; \
 if [ "$(echo -e $(cat ${RGFILE}) | grep -oP "(?<=DS:)[^[:space:]]*")" != "${FULLSM}" ]; then \
  echo -e "Error, FULLSM ${FULLSM} mismatch between worklist and ${RGFILE}"; continue; \
 fi; \
 if [ "${FULLSM}.$(echo -e $(cat ${RGFILE}) | grep -oP "(?<=PU:)[^[:space:]]*" | tr ":" "^")" != "${RGBASE}" ]; then \
  echo -e "Error, RGBASE ${RGBASE} mismatch between LOOKUP and ${RGFILE}"; continue; \
 fi; \
export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
export MEM=50; \
# done; \
# done
   bsub \
     -J "${RGBASE}_s01alnsrt" \
     -o "${BASE}/WXS_Aquilla/04-LOGS/${RGBASE}_s01alnsrt.%J" \
     -u "${EMAIL}" \
     -n ${THREADS} -W 2860 \
     -M 66000000 \
     -R "rusage[mem=69152]" \
     -q research-hpc \
     -a 'docker(vifehe/bwaref)' \
     entrypoint.sh; \
export MEM=32; \
export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.bam"; \
echo -e "02 - Starting validatesam per sample ${FULLSM} and RGBASE ${RGBASE}\nBAMFILE:${BAMFILE}";
echo ${BAMFILE}; \
export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
bsub \
     -w "done(\"${RGBASE}_s01alnsrt\")" \
     -J "${RGBASE}_s02vldate" \
     -o "${BASE}/WXS_Aquilla/04-LOGS/${RGBASE}_s02vldate.%J"\
     -u "${EMAIL}" \
     -n1 -W 360 \
     -R "rusage[mem=49000]" \
     -q research-hpc \
     -a "docker(vifehe/validatesamref)" \
     entrypoint.sh -IGNORE INVALID_VERSION_NUMBER -IGNORE INVALID_TAG_NM; \
echo -e "03 - Starting intersectbed per sample ${FULLSM} and RGBASE ${RGBASE}\nBAMFILE:${BAMFILE}";\
export RUN_TYPE="paddedexome"; \
export MEM=16; \
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
echo -e "04 - Starting validatesam per sample ${FULLSM} and RGBASE ${RGBASE}\nBAMFILE:${BAMFILE}";\
export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.isec-${RUN_TYPE}.bam"; \
bsub \
     -w "done(\"${RGBASE}_s03intsct\")" \
     -J "${RGBASE}_s04vldate" \
     -o "${BASE}/WXS_Aquilla/04-LOGS/${RGBASE}_s04vldate.%J"\
     -u "${EMAIL}" \
     -n1 -W 1360 \
     -R "rusage[mem=20000]" \
     -q research-hpc \
     -a "docker(vifehe/validatesamref)" \
     entrypoint.sh -IGNORE INVALID_VERSION_NUMBER -IGNORE INVALID_TAG_NM -IGNORE MATE_NOT_FOUND; \
echo -e "05 - Starting bamtocram per sample ${FULLSM} and RGBASE ${RGBASE}\nBAMFILE:${BAMFILE}\nOUT_DIR:${OUT_DIR}";\
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
    CLEANUP_LIST+=("${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}${FQ1EXT}");\
  CLEANUP_LIST+=("${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}${FQ2EXT}");\
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
   -R "rusage[mem=16000]" \
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
export MEM=16; \
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
sleep ${DELAY2:-10}s; \
echo "Finished submitting jobs for ${FULLSM}"; \
done; \
echo "Finished all samples within ${PR}"; 




# START=1;
# END=15
# for line in $(sed -n "${START},${END}p" "${WORKLIST}"); do \
# ls -lht $line
# done




























