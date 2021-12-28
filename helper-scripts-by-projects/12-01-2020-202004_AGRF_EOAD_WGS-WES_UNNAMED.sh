### 12-01-2020 - README file for "202004_AGRF_EOAD_WGS-WES_UNNAMED" sequences processing over at MGI



# generate metadata file for the PROJECT in question
HOMEDIR="/40/AD/AD_Seq_Data/01.-RawData"
PR="202004_AGRF_EOAD_WGS-WES_UNNAMED"

##List directories within the PR directory
metafile="${PR}-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv"
# vim ${metafile}   ## copy here the metadata 
sed -i 's/\t/,/g' ${metafile}

BASE="/gscmnt/gc2645/wgs/"



#####################################
### 3 - Start pipeline from fastq to variant calling
#####################################

for LOOKUP in $(ls ${metafile}); do \
awk -F"," '{print $1","$2","$4","$3","$5}' ${LOOKUP} > ${PR}-MGIID-SMID-DNA-PR-BAMFILE.${LOOKUP##*.}; done


## Create WORKLIST
for LOOKUP in $(ls ${PR}-MGIID-SMID-DNA-PR-BAMFILE.csv*); do
cut -d, -f1,2,3 "${LOOKUP}" | tr "," "^" | sort -u > ${PR}-worklist_fullsm.${LOOKUP##*.};
done



#################################################################
##### Run pipeline from bam back to fastq and then call variants
#################################################################

bjobs | grep RUN | cut -d\  -f1 | xargs bkill 
bjobs | grep PEN | cut -d\  -f1 | xargs bkill 


# If needed incorporate this at line 159, after declaring ${FULLSM}
#if [ -d ${BASE}/WXS_Aquilla/03-FINAL/${PR}/${FULLSM} ] || [ -d ${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM} ]; then echo "Skipping ${FULLSM}, already run(ning)"; continue; else echo "RUN ${FULLSM}"; fi; \

PR="202004_AGRF_EOAD_WGS-WES_UNNAMED"
BASE="/gscmnt/gc2645/wgs/"

# -w "done(\"${FULLSM}_s00rvtbam\")" \

# Revert BAM
Snumber=1
START=1; \
END=2; \
DELAY=10
EMAIL="achal@wustl.edu"
SHELLDROP=0
export BASE="/gscmnt/gc2645/wgs"; \
export WORKDIR="${BASE}/tmp"; \
export THREADS=16; \
export BWA_PIPE_SORT=1; \
export TIMING=1;\
PR="202004_AGRF_EOAD_WGS-WES_UNNAMED"; \
WORKLIST="${BASE}/WXS_Aquilla/01-RAW/${PR}-worklist_fullsm.csv"; \
LOOKUP="${BASE}/WXS_Aquilla/01-RAW/${PR}-MGIID-SMID-DNA-PR-BAMFILE.csv"; \
LOOKUP_COL_SM=1; \
LOOKUP_COL_DNA=2; \
LOOKUP_COL_PR=3; \
LOOKUP_COL_BAMFILE=4; \
LOOKUP_COL_BAMPATH=5; \
export RUN_TYPE="paddedexome"; \
for FULLSM in $(sed -n "${START},${END}p" "${WORKLIST}"); do \
echo "*****************************************************"; \
echo "*****************************************************"; \
echo "Doing sample number***********************************: " $Snumber; \
echo "*****************************************************"; \
echo "*****************************************************"; \
((Snumber=Snumber+1)); \
	export FULLSM="${FULLSM}"; \
	SM="$(echo "${FULLSM}" | cut -d^ -f1)"; \
	IFS=$'\n' export DNA=($(awk -F, "\$${LOOKUP_COL_SM} == \"${SM}\"" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_DNA} | sort -u)); \
	if [ ${#DNA[@]} -gt 1 ]; then echo "Warning, \${DNA} not unique for ${SM} (n=${#DNA[@]}: ${DNA[@]})"; fi; \
	DNA="$(echo "${FULLSM}" | cut -d^ -f2)"; \
	PR="$(echo "${FULLSM}" | cut -d^ -f3)"; \
  mkdir -p "${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/"; \
	export OUT_DIR="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}"; \
	INPUT_LIST=(); \
	WAIT_LIST=(); \
	for BAMBASE in $(grep "${SM},${DNA},${PR}" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_BAMFILE}); do \
	echo -e "00a - In this round we will analyze ${FULLSM}\n${BAMFILE}\n${SM}\n${DNA}\n${PR}"; \
	FENIX_BAM_PATH="$(grep "${SM},${DNA},${PR}" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_BAMPATH})"; \
  FENIX_BAMFILE="${FENIX_BAM_PATH}/${BAMBASE}"; \
	rsync -avh -L ${USER}@fenix.psych.wucon.wustl.edu:${FENIX_BAMFILE} ${OUT_DIR}/ ;\
	export BAMFILE="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${BAMBASE}"; \
	export MEM=36; \
	export CREATE_RGFILE=1; \
	export DEBUG=1; \
	export REMOVE_INPUT=0; \
# done; \
# done
	bsub \
	-J "${FULLSM}_s00rvtbam" \
	-o "${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${FULLSM}.${BAMBASE}_s00rvtbam.%J" \
	-u "${EMAIL}" \
	-n1 -W 2440 \
	-R "rusage[mem=48192]" \
	-q research-hpc \
	-a 'docker(vifehe/revertbam:1.0)' \
	/bin/bash; \
	done;\
  echo -e "00b - Starting revertbam for ${FULLSM}\n${BAMFILE}\n${SM}\n${DNA}\n${PR}"; \
  SAMTOOLS="/gscmnt/gc2645/wgs/variant_calling/samtools"; \
  IFS=$'\n' RGS=($(${SAMTOOLS} view -H "${BAMFILE}" | grep "^@RG")); \
  for RG in ${RGS[@]}; do RGID="$(echo ${RG} | grep -oP "(?<=ID:)[^[:space:]]*")"; \
  RGID_NEW="$(echo ${RGID} | cut -d: -f2- | sed 's/:/^/g')"; \
  RGBASE="${FULLSM}.${RGID_NEW}"; \
  echo -e "01 - Starting bwa per sample ${FULLSM} and RGBASE ${RGBASE}\nFQ1:${FQ1}\nFQ2:${FQ2}\nRGFILE:${RGFILE}\nDNA:${DNA}\nPR:${PR}\nIFS:${IFS}";\
  export RGFILE="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.rgfile"; \
  export FULLSM_RGID="${RGBASE}"; \
  unset BAMFILE; \
  export OUT_DIR="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}"; \
  FQ1EXT="r1.fastq"; \
  FQ2EXT="r2.fastq"; \
echo -e "01 - Starting bwa per sample ${FULLSM} and RGBASE ${RGBASE}"; \
export RGFILE="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.rgfile"; \
export FQ1="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.${FQ1EXT}"; \
export FQ2="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.${FQ2EXT}"; \
export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
echo -e "00a - Doing FASTQ and rgfiles for sample ${FULLSM} and RGBASE ${RGBASE}\nFQ1:${FQ1}\nFQ2:${FQ2}\nRGFILE:${RGFILE}\nDNA:${DNA}"; \
# done; \ 
# done 
export MEM=56; \
   bsub \
     -w "done(\"${FULLSM}_s00rvtbam\")" \
     -J "${RGBASE}_s01alnsrt" \
     -o "${BASE}/WXS_Aquilla/04-LOGS/${RGBASE}_s01alnsrt.%J" \
     -u "${EMAIL}" \
     -n ${THREADS} -W 4160 \
     -M 76000000 \
     -R "rusage[mem=79152]" \
     -q research-hpc \
     -a 'docker(vifehe/bwaref)' \
     entrypoint.sh; \
echo -e "02 - Starting validatesam per sample ${FULLSM} and RGBASE ${RGBASE}\nBAMFILE:${BAMFILE}"; \
export MEM=16; \
export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.bam"; \
echo ${BAMFILE}; \
export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
bsub \
     -w "done(\"${RGBASE}_s01alnsrt\")" \
     -J "${RGBASE}_s02vldate" \
     -o "${BASE}/WXS_Aquilla/04-LOGS/${RGBASE}_s02vldate.%J" \
     -u "${EMAIL}" \
     -n1 -W 1360 \
     -R "rusage[mem=18192]" \
     -q research-hpc \
     -a "docker(vifehe/validatesamref)" \
     entrypoint.sh -IGNORE INVALID_VERSION_NUMBER -IGNORE INVALID_TAG_NM; \
echo -e "03 - Starting intersectbed per sample ${FULLSM} and RGBASE ${RGBASE}\nBAMFILE:${BAMFILE}"; \
export RUN_TYPE="paddedexome"; \
export BEDFILE="${BASE}/Genome_Ref/GRCh38/Capture_Padded.GRCh38.bed"; \
export COVERED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Covered.GRCh38.bed"; \
export PADDED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Padded.GRCh38.bed";\
export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.bam"; \
bsub \
     -w "done(\"${RGBASE}_s02vldate\")" \
     -J "${RGBASE}_s03intsct" \
     -o "${BASE}/WXS_Aquilla/04-LOGS/${RGBASE}_s03intsct.%J" \
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
     -o "${BASE}/WXS_Aquilla/04-LOGS/${RGBASE}_s04vldate.%J" \
     -u "${EMAIL}" \
     -n1 -W 1360 \
     -R "rusage[mem=8192]" \
     -q research-hpc \
     -a "docker(vifehe/validatesamref)" \
     entrypoint.sh -IGNORE INVALID_VERSION_NUMBER -IGNORE INVALID_TAG_NM -IGNORE MATE_NOT_FOUND; \
export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.bam"; \
export OUT_DIR="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/"; \
INPUT_LIST+=("${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.isec-${RUN_TYPE}.bam"); \
WAIT_LIST+=("&&" "done(\"${RGBASE}_s04vldate\")"); \
CLEANUP_LIST+=("${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.${FQ1EXT}");\
CLEANUP_LIST+=("${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.${FQ2EXT}");\
sleep ${DELAY1:-1}s; \
done; \
echo "Finished submitting per-RG jobs for ${BAMBASE}"; \
echo -e "06 - Starting markduplicates for FULLSM ${FULLSM} with INPUT_LIST:${INPUT_LIST[@]}"; \
export FULLSM="${FULLSM}"; \
export MEM=52; \
export OUT_DIR="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/"; \
bsub \
   -w ${WAIT_LIST[@]:1} \
   -J "${FULLSM}_s06mrkdup" \
   -o "${BASE}/WXS_Aquilla/04-LOGS/${FULLSM}_s06mrkdup.%J" \
   -u "${EMAIL}" \
   -n1 -W 2440 \
   -M 66000000 \
   -R "rusage[mem=69152]" \
   -q research-hpc \
   -a "docker(vifehe/markduplicates)" \
   ${INPUT_LIST[@]}; \
  CLEANUP_LIST+=("${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.isec-${RUN_TYPE}.markd.bam"); \
echo -e "07 - Starting validatesam for FULLSM ${FULLSM}\nBAMFILE=${BAMFILE}"; \
export MEM=26; \
export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.bam"; \
export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
bsub \
   -w "done(\"${FULLSM}_s06mrkdup\")" \
   -J "${FULLSM}_s07vldate" \
   -o "${BASE}/WXS_Aquilla/04-LOGS/${FULLSM}_s07vldate.%J" \
   -u "${EMAIL}" \
   -n1 -W 2440 \
   -R "rusage[mem=28192]" \
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
   -n1 -W 2440 \
   -q research-hpc \
   -a "docker(vifehe/depthofcoverage)" \
   entrypoint.sh; \
echo -e "09 - Starting baserecal for FULLSM ${FULLSM}\nBAMFILE=${BAMFILE}"; \
export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.bam"; \
export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
export DBSNP1="${BASE}/Genome_Ref/GRCh38/20190522_bundle/dbsnp_138.hg38.vcf.gz"; \
export ONEKGP1="${BASE}/Genome_Ref/GRCh38/20190522_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz"; \
export MILLS_GOLD="${BASE}/Genome_Ref/GRCh38/20190522_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"; \
bsub \
   -w "done(\"${FULLSM}_s07vldate\")" \
   -J "${FULLSM}_s09bsercl" \
   -o "${BASE}/WXS_Aquilla/04-LOGS/${FULLSM}_s09bsercl.%J" \
   -u "${EMAIL}" \
   -n1 -W 2440 \
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
   -n1 -W 1360 \
   -R "rusage[mem=18192]" \
   -q research-hpc \
   -a "docker(vifehe/validatesamref)" \
   entrypoint.sh -IGNORE INVALID_VERSION_NUMBER -IGNORE INVALID_TAG_NM -IGNORE MATE_NOT_FOUND; \
echo -e "11 - Starting haplotype caller for FULLSM ${FULLSM}\nBAMFILE=${BAMFILE}"; \
export MEM=32; \
export RUN_TYPE="paddedexome"; \
export DBSNP1="${BASE}/Genome_Ref/GRCh38/20190522_bundle/dbsnp_138.hg38.vcf.gz" ;\
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
   -n1 -W 3880 \
   -M 46000000 \
   -R "rusage[mem=49152]" \
   -q research-hpc \
   -a "docker(vifehe/haplocallerv2)" \
   entrypoint.sh; \
echo -e "12 - Starting varianteval caller for FULLSM ${FULLSM}\GVCF=${GVCF}"; \
export MEM=2; \
export GVCF="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.recal.raw.snps.indels.g.vcf.gz"; \
export DBSNP1="${BASE}/Genome_Ref/GRCh38/20190522_bundle/dbsnp_138.hg38.vcf.gz"; \
export CCDS="/gscmnt/gc2645/wgs/Genome_Ref/GRCh38/GRCh38.CCDS.exons.sorted.bed"; \
export CLEANUP="${CLEANUP_LIST[@]}"; \
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

### CHECK UP if ALL gvcfs have been created:
WORKLIST="${BASE}/WXS_Aquilla/01-RAW/${PR}-worklist_fullsm.csv"
for x in $(cat ${WORKLIST}); do ls -lh ../02-TRANSIT/${PR}/${x}/*.vcf.gz; done

### CLEANING UP FASTQ files for DIAN_REDCLOUD
# 1 - generate a list of gvcf files finished, just keep the FULLSMID info
ls ../02-TRANSIT/${PR}/*/*.g.vcf.gz | cut -d "/" -f 4 | sort | uniq > ${PR}-fqtodelete.list
ls ../02-TRANSIT/${PR}/*/*.g.vcf.gz | cut -d "/" -f 4 | sort | uniq > ${PR}-fqtodelete.list2

cat ${PR}-fqtodelete.list2 ${PR}-fqtodelete.list | sort | uniq -c | awk '{if ($1==1) print $2}' > ${PR}-fqtodelete.LIST

## delete all fastq & cram within that FULLSMID directory for the samples that gvcf is completed
for x in $(cat ${PR}-fqtodelete.list); do rm ${PR}/${x}/*.fastq; done
for x in $(cat ${PR}-fqtodelete.list); do rm ../02-TRANSIT/${PR}/${x}/*.bam; done


### Check poresentce of files
### Test if we got Bam files for all samples
ls ./${PR}/*/*.cram | wc -l
34
ls ./${PR}/*/*.genome.bam | wc -l
1 --> 01-RAW/${PR}/955-0004-0287652^0010007593^${PR}/abf5010305ac466eb0d6c4d9c6e1e045.genome.bam
ls ./${PR}/*/*.r1.fastq | wc -l
538
ls ./${PR}/*/*.r2.fastq | wc -l
538
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
