### 11-16-2020 - README file for ADNI_WGS on compute1 

##############
## 11-16-2020; Achal Neupane
############## 
# unset all ENV variables
# for c in $(env | cut -d '=' -f 1); do unset $c; done


# First we will unwrap BAM to FASTQ (inter-leaved) on Compute1. If you prefer, you can run everything on MGI, but you'll have to change the BASE and variables. I was working on both MGI and Compute1 so that I could divide the workload between clusters and finish processing these samples relatively fast. 

# # Samples left to process
# ADNI_4643^LP6005123-DNA_E04^ADNI_WGS
# ADNI_4645^LP6005123-DNA_C02^ADNI_WGS

## Transfer BAM files from Fenix to Compute1

START=796
END=797
USER="achal"
BASE="/storage1/fs1/cruchagac/Active/achal"
DEST="${BASE}/WXS_Aquilla/01-RAW"
EMAIL="achal@wustl.edu"
WORKLIST="${DEST}/${PR}-worklist_fullsm_809.csv"
cd /storage1/fs1/cruchagac/Archive/ADNI_WGS

Snumber=1;
START=796; 
END=797;

s1_id="01f0ac4c-9570-11ea-b3c4-0ae144191ee3"
s1_path="storage1/fs1/cruchagac"
for FULLSM in $(sed -n "${START},${END}p" "${WORKLIST}"); do 
echo "Doing Number: ${Snumber}"  
((Snumber=${Snumber}+1))
DNA="$(echo "${FULLSM}" | cut -d^ -f2).bam"
globus transfer "${s1_id}:${s1_path}/Archive/ADNI_WGS/${DNA}" "${s1_id}:${s1_path}/Active/achal/WXS_Aquilla/01-RAW/all_bam_files/${DNA}"
done

# Check the transfer status on Globus site: https://app.globus.org/activity

#################################################################
########## Run the pipeline from bam back to fastq ##################
#################################################################
# bjobs | grep RUN | cut -d\  -f1 | xargs bkill 
# bjobs | grep PEN | cut -d\  -f1 | xargs bkill 



WORKLIST="ADNI_WGS-worklist_fullsm_809.csv"
NUM="$(cat $WORKLIST| wc -l)"


START=796
END=797
USER="achal"
BASE="/storage1/fs1/cruchagac/Active/achal"
DEST="${BASE}/WXS_Aquilla/01-RAW"
EMAIL="achal@wustl.edu"
$Snumber=1
export PR="ADNI_WGS"
WORKLIST="${DEST}/${PR}-worklist_fullsm_809.csv"
for FULLSM in $(sed -n "${START},${END}p" "${WORKLIST}"); do \
echo "*****************************************************"; \
echo "*****************************************************"; \
echo "Doing sample number***********************************: " $Snumber; \
echo "*****************************************************"; \
echo "*****************************************************"; \
export SCRATCH1=/scratch1/fs1/cruchagac; \
export STORAGE1=/storage1/fs1/cruchagac/Active; \
export BAMDIR=/storage1/fs1/cruchagac/Active/achal/WXS_Aquilla/01-RAW; \
export LSF_DOCKER_VOLUMES="$HOME:$HOME $STORAGE1:$STORAGE1 $SCRATCH1:$SCRATCH1"; \
((Snumber=${Snumber}+1)); \
export TMP_DIR="${DEST}/ADNI_WGS/${FULLSM}"
export DEBUG=1; \
export REMOVE_INPUT=1; \
export BASE="/storage1/fs1/cruchagac/Active/achal"; \
export THREADS=32; \
export CREATE_RGFILE=1; \
export BWA_PIPE_SORT=1; \
export TIMING=1; \
export DEST="${BASE}/WXS_Aquilla/01-RAW"; \
mkdir -p ${DEST}/${PR}/${FULLSM}; \
export SM="$(echo "${FULLSM}" | cut -d^ -f1)"; \
export FULLSM=${FULLSM}; \
unset DNA; \
export DNA="$(echo "${FULLSM}" | cut -d^ -f2)"; \
export PR="$(echo "${FULLSM}" | cut -d^ -f3)"; \
export OUT_DIR="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}"; \
export BAMBASE=${DNA}.bam; \
export MEM=32; \
unset BAMFILE; \
export BAMFILE="/storage1/fs1/cruchagac/Active/achal/WXS_Aquilla/01-RAW/all_bam_files/${BAMBASE}"; \
echo -e "00a - In this round we will analyze ${FULLSM}\n${BAMFILE}\n${SM}\n${DNA}\n${PR}"; \
bsub \
  -J "${FULLSM}_s00bam2fqv2" \
  -o "${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${FULLSM}.${BAMBASE}_s00bam2fq.%J" \
  -u "${EMAIL}" \
  -n ${THREADS} -W 5760 \
  -G compute-cruchagac \
  -R "rusage[mem=49152]" \
  -g /achal \
  -q general \
  -a 'docker(achalneupane/bam2fqv2:latest)' \
  entrypoint.sh; \
done



# Create Metafile
FILE="${DEST}/ADNI_WGS_SM_DNA_PR_RG_FQ1ext_FQ2ext_fullsm.csv"
####################################################
for FULLSM in $(sed -n "${START},${END}p" "${WORKLIST}"); do \
echo $FULLSM
ls -lht ${FULLSM}/*fq.gz
SM="$(echo "$FULLSM" | cut -d^ -f1)"
DNA="$(echo "$FULLSM" | cut -d^ -f2 )"
FULLSM=${FULLSM}
for base in `ls ${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM} | grep .rgfile| sort -V`; do
RGBASE=${base/.rgfile/}
RG="$(echo ${RGBASE})"
R="fq.gz"
echo "${SM},${DNA},${PR},${RG},${R},${FULLSM}"  >>  ${FILE}
done
done

# transfer Metafile to MGI
scp  ${FILE} fernandezv@virtual-workstation1.gsc.wustl.edu:/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/

# Rest of the work, we will do on MGI. If you prefer, you can run everything on MGI, but you'll have to change the BASE and variables. I was working on both MGI and Compute1 so I could divide the workload between clusters and finish processing these samples relatively fast. 

## Log onto MGI
USER="achal"
BASE="/gscmnt/gc2645/wgs"
PR="ADNI_WGS"
DEST="${BASE}/WXS_Aquilla/01-RAW"
# ##############################
LOOKUP="${BASE}/WXS_Aquilla/01-RAW/${PR}_SM_DNA_PR_RG_FQ1ext_FQ2ext_fullsm.csv"
# sed -i 's/\t/,/g' ${LOOKUP}

WORKLIST="${DEST}/ADNI_WGS-worklist_fullsm_809.csv"

#################################################################################
# bjobs | grep RUN | cut -d\  -f1 | xargs bkill 
# bjobs | grep PEN | cut -d\  -f1 | xargs bkill 
# bjobs | grep SSUSP | cut -d\  -f1 | xargs bkill 

USER="achal"
START=796;
END=796;
# END="$(cat $WORKLIST| wc -l)";
SHELLDROP=1
# for FULLSM in $(sed -n "${START},${END}p" "${WORKLIST}"); do \
Snumber=1;
export BASE="/gscmnt/gc2645/wgs"; \
export WORKDIR="${BASE}/tmp"; \
export THREADS=16; \
export BWA_PIPE_SORT=1; \
export REMOVE_INPUT=1; \
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
FBASE="${dir}/${FULLSM}"; \
FQEXT=($(awk -F, "\$${LOOKUP_COL_RGBASE} == \"${RGBASE}\"" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_FQEXT})); \
RGFILE="${FBASE}/${RGBASE}.rgfile";\
FQ="${FBASE}/${RGBASE}.${FQEXT}";\
DEST="${BASE}/WXS_Aquilla/01-RAW";\
echo -e "00a - Uploading FASTQ and rgfiles for sample ${FULLSM} and RGBASE ${RGBASE}\nFQ:${FQ}\nRGFILE:${RGFILE}\nDNA:${DNA}";\
mkdir -p ${DEST}/${PR}/${FULLSM}; \
rsync -avh ${USER}@compute1-client-1.ris.wustl.edu:${RGFILE} ${DEST}/${PR}/${FULLSM}/; \
rsync -avh -L ${USER}@compute1-client-1.ris.wustl.edu:${FQ} ${DEST}/${PR}/${FULLSM}/; \
## delete transferred files on compute 1
# ssh ${USER}@compute1-client-1.ris.wustl.edu rm "${FQ}"; \
export RGFILE="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.rgfile"; \
export FULLSM_RGID="${RGBASE}"; \
unset BAMFILE; \
unset FQ1; \
unset FQ2; \
IFS=$'\n' FQEXT=($(awk -F, "\$${LOOKUP_COL_RGBASE} == \"${RGBASE}\"" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_FQEXT})); \
echo -e "01 - Starting bwa per sample ${FULLSM} and RGBASE ${RGBASE}";\
export RGFILE="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.rgfile"; \
export FQ="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.${FQEXT}"; \
export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
export MEM=45; \
bsub \
     -J "${RGBASE}_s01alnsrt" \
     -o "${BASE}/WXS_Aquilla/04-LOGS/${RGBASE}_s01alnsrt.%J" \
     -u "${EMAIL}" \
     -n ${THREADS} -W 2880 \
     -M 56000000 \
     -R "rusage[mem=59152]" \
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
     -R "rusage[mem=8192]" \
     -q research-hpc \
     -a "docker(vifehe/validatesamref)" \
     entrypoint.sh -IGNORE INVALID_VERSION_NUMBER -IGNORE INVALID_TAG_NM -IGNORE MATE_NOT_FOUND; \
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
done; \
echo -e "06 - Starting markduplicates for FULLSM ${FULLSM} with INPUT_LIST:${INPUT_LIST[@]}"; \
export FULLSM="${FULLSM}"; \
export MEM=45; \
export OUT_DIR="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/"; \
bsub \
   -w ${WAIT_LIST[@]:1} \
   -J "${FULLSM}_s06mrkdup" \
   -o "${BASE}/WXS_Aquilla/04-LOGS/${FULLSM}_s06mrkdup.%J" \
   -u "${EMAIL}" \
   -n1 -W 1440 \
   -M 66000000 \
   -R "rusage[mem=69152]" \
   -q research-hpc \
   -a "docker(vifehe/markduplicates)" \
   ${INPUT_LIST[@]}; \
  CLEANUP_LIST+=("${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.isec-${RUN_TYPE}.markd.bam");\
echo -e "07 - Starting validatesam for FULLSM ${FULLSM}\nBAMFILE=${BAMFILE}"; \
export MEM=20; \
export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.bam"; \
export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
bsub \
   -w "done(\"${FULLSM}_s06mrkdup\")" \
   -J "${FULLSM}_s07vldate" \
   -o "${BASE}/WXS_Aquilla/04-LOGS/${FULLSM}_s07vldate.%J" \
   -u "${EMAIL}" \
   -n1 -W 1440 \
   -R "rusage[mem=24192]" \
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
export MEM=20;\
export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.recal.bam"; \
export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
bsub \
   -w "done(\"${FULLSM}_s09bsercl\")" \
   -J "${FULLSM}_s10vldate" \
   -o "${BASE}/WXS_Aquilla/04-LOGS/${FULLSM}_s10vldate.%J" \
   -u "${EMAIL}" \
   -n1 -W 360 \
   -R "rusage[mem=24192]" \
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


# for x in $(cat ${WORKLIST}); do ls -lh ../02-TRANSIT/${PR}/${x}/*.vcf.gz; done

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