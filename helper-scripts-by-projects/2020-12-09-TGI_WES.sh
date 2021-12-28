################################
# 2020-12-09 - TGI_WES 
################################

SOURCEDIR="/40/AD/AD_Seq_Data/01.-RawData/TGI_WES/bams"

## This project has 174 FULLSMs and 174 bams

## IN Fenix, generate the metafile with the following fields: SM, DNA, PR, BAMFILE as follows: 
WORKDIR="/40/AD/AD_Seq_Data/01.-RawData"
ls TGI_WES/bam/*/*.bam | cut -d "/" -f 4 | cut -d "." -f1 > TGI_WES-worklist-fullsm.csv
ls TGI_WES/bam/*/*.bam | cut -d "/" -f 4 | awk -F "." '{print $1"."$2}' > TGI_WES-bamfiles.list

ls TGI_WES/bam/*/*.bam | cut -d "/" -f 4 | cut -d "." -f1 | awk -F "^" '{print $1","$2","$3","}' > TGI_WES-meatdata-head
paste TGI_WES-meatdata-head TGI_WES-bamfiles.list| sed 's/\t//g' > TGI_WES-sm_dna_pr_bamfile.csv


##At MGI:
mkdir TGI_WES
SOURCE="/40/AD/AD_Seq_Data/01.-RawData"
PR="TGI_WES"
metafile="${PR}-sm_dna_pr_bamfile_achal.csv"
DEST="/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW"
scp -p ${USER}@fenix.psych.wucon.wustl.edu:${SOURCE}/${metafile} ${DEST}/;



## on MGI - Generate the workfile
cut -d, -f1-3 ${DEST}/${metafile} | tr "," "^" | sort -u > ${DEST}/${PR}-worklist_fullsm.csv
# WORKLIST="${DEST}/${PR}-worklist_fullsm.csv"


# Find the missing samples
ls ./*/*vcf.gz| cut -d/ -f2 > /gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/TGI_WES-worklist_VCF_done.csv

VCF="/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/TGI_WES-worklist_VCF_done.csv"
ALL_samples="/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/TGI_WES-worklist_fullsm.csv"
comm -23 $ALL_samples $VCF

# This sample found missing
echo "10R_R95_13^Original_38A_F10^TGI_WES" > "${BASE}/WXS_Aquilla/01-RAW/${PR}-worklist_fullsm-missing.csv" 



## Run these samples using pipeline as defined in  MGI-BAM-UPLOAD-FASTQ-RUN-v1.sh
## Sets to define:
START=1; \
END=1; \
DELAY=10
EMAIL="fernandezv@wustl.edu"
SHELLDROP=0

export BASE="/gscmnt/gc2645/wgs"; \
export WORKDIR="${BASE}/tmp"; \
export THREADS=16; \
export BWA_PIPE_SORT=1; \
export TIMING=1;\
USER="achal";\
PR="TGI_WES"
WORKLIST="${BASE}/WXS_Aquilla/01-RAW/${PR}-worklist_fullsm-missing.csv" 
LOOKUP="${BASE}/WXS_Aquilla/01-RAW/${PR}-sm_dna_pr_bamfile_achal.csv"; \
LOOKUP_COL_SM=1; \
LOOKUP_COL_DNA=2; \
LOOKUP_COL_PR=3; \
LOOKUP_COL_BAMFILE=4; \
export RUN_TYPE="paddedexome";


for FULLSM in $(sed -n "${START},${END}p" "${WORKLIST}"); do \
	export FULLSM="${FULLSM}"; \
	SOURCE="/40/AD/AD_Seq_Data/01.-RawData/TGI_WES/bams/";\
	SM="$(echo "${FULLSM}" | cut -d^ -f1)"; \
	IFS=$'\n' export DNA=($(awk -F, "\$${LOOKUP_COL_SM} == \"${SM}\"" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_DNA} | sort -u)); \
	if [ ${#DNA[@]} -gt 1 ]; then echo "Warning, \${DNA} not unique for ${SM} (n=${#DNA[@]}: ${DNA[@]})"; fi; \
	DNA="$(echo "${FULLSM}" | cut -d^ -f2)"; \
	PR="$(echo "${FULLSM}" | cut -d^ -f3)"; \
	export OUT_DIR="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}"; \
	INPUT_LIST=(); \
	WAIT_LIST=(); \
	CLEANUP_LIST=();\
	DEST="/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/${PR}/${FULLSM}";\
	echo -e "Create directory ${DEST}";\
	mkdir ${DEST};
	for BAMBASE in $(grep "${SM},${DNA},${PR}" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_BAMFILE}); do \
	BAMFILE="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${BAMBASE}"; \
	echo -e "00a - Uploading BAMFILE ${BAMFILE} for sample ${FULLSM}";\
	rsync -avh ${USER}@fenix.psych.wucon.wustl.edu:${SOURCE}/${BAMBASE%.*}.* ${DEST}/;
	echo -e "00b - Starting revertbam for ${FULLSM}\n${BAMFILE}\n${SM}\n${DNA}\n${PR}"; \
	BAMBASE="${BAMBASE/.bam/}"; \
	export BAMFILE="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${BAMBASE}.bam"; \
	export MEM=16; \
	export CREATE_RGFILE=1; \
	export DEBUG=1; \
	export REMOVE_INPUT=0; \
	export OUT_DIR="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}" ;\
	echo -e "OUTPUT is in "${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${BAMBASE}_s00rvtbam.%J" \n OUTDIR is ${OUT_DIR}";\
	# done; \
 #    done
	bsub \
	-J "${FULLSM}_s00rvtbam" \
	-o "${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${BAMBASE}_s00rvtbam.%J" \
	-u "${EMAIL}" \
	-n1 -W 1440 \
	-R "rusage[mem=18192]" \
	-q research-hpc \
	-a 'docker(vifehe/revertbam:1.1)' \
	/bin/bash; \
	echo -e "00c - reading readgroups RG of BAMfile"
	SAMTOOLS="/gscmnt/gc2645/wgs/variant_calling/samtools"; \
	BAMFILE="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${BAMBASE}"; \
	IFS=$'\n' RGS=($(${SAMTOOLS} view -H "${BAMFILE}.bam" | grep "^@RG")); \
	for RG in ${RGS[@]}; do RGID="$(echo ${RG} | grep -oP "(?<=ID:)[^[:space:]]*")"; \
	RGID_NEW="$(echo ${RGID} | cut -d: -f2- | sed 's/:/^/g')"; \
	RGBASE="${FULLSM}.${RGID_NEW}"; \
	export RGFILE="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.rgfile"; \
	export FULLSM_RGID="${RGBASE}"; \
	unset BAMFILE; \
	export OUT_DIR="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}"; \
	FQ1EXT="r1.fastq"; \
	FQ2EXT="r2.fastq"; \
	echo -e "01 - Starting bwa per sample ${FULLSM} and RGBASE ${RGBASE}\nFQ1:${FQ1}\nFQ2:${FQ2}\nRGFILE:${RGFILE}\nDNA:{$DNA}\nPR:${PR}";\
	export FQ1="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.${FQ1EXT}"; \
	export FQ2="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.${FQ2EXT}"; \
	export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
	export MEM=42; \
	export CLEANUP=1; 
	export REMOVE_INPUT=0; \
	bsub \
	-w "done(\"${FULLSM}_s00rvtbam\")" \
	-J "${RGBASE}_s01alnsrt" \
	-o "${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}_s01alnsrt.%J" \
	-u "${EMAIL}" \
	-n ${THREADS} -W 2880 \
	-M 66000000 \
	-R "rusage[mem=69152]" \
	-q research-hpc \
	-a 'docker(vifehe/bwaref)' \
	/bin/bash; \
	echo -e "02 - Perform ValidateSam-I on newly aligned bamfile ${BAMFILE}"
	export MEM=16; \
	export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.bam"; \
	echo ${BAMFILE}; \
	export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
	bsub \
	-w "done(\"${RGBASE}_s01alnsrt\")" \
	-J "${RGBASE}_s02vldate" \
	-o "${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}_s02vldate.%J" \
	-u "${EMAIL}" \
	-n1 -W 1360 \
	-R "rusage[mem=18192]" \
	-q research-hpc \
	-a "docker(vifehe/validatesamref)" \
	entrypoint.sh -IGNORE INVALID_VERSION_NUMBER -IGNORE INVALID_TAG_NM; \
	echo -e "03 - perform intersectbed"
	export RUN_TYPE="paddedexome"; \
	export BEDFILE="${BASE}/Genome_Ref/GRCh38/Capture_Padded.GRCh38.bed"; \
	bsub \
	-w "done(\"${RGBASE}_s02vldate\")" \
	-J "${RGBASE}_s03intsct" \
	-o "${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}_s03intsct.%J" \
	-u "${EMAIL}" \
	-n1 -W 1360 \
	-q research-hpc \
	-a "docker(vifehe/intersectbedref)" \
	entrypoint.sh; \
	echo -e "04 - Perform ValidateSam-II on ${RGBASE}.aln.srt.isec-${RUN_TYPE}.bam"
	export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.isec-${RUN_TYPE}.bam"; \
	export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
	bsub \
	-w "done(\"${RGBASE}_s03intsct\")" \
	-J "${RGBASE}_s04vldate" \
	-o "${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}_s04vldate.%J" \
	-u "${EMAIL}" \
	-n1 -W 360 \
	-R "rusage[mem=18192]" \
	-q research-hpc \
	-a "docker(vifehe/validatesamref)" \
	entrypoint.sh -IGNORE INVALID_VERSION_NUMBER -IGNORE INVALID_TAG_NM -IGNORE MATE_NOT_FOUND; \
	unset MEM; \
	export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.bam"; \
	export OUT_DIR="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/"; \
	INPUT_LIST+=("${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.isec-${RUN_TYPE}.bam"); \
	WAIT_LIST+=("&&" "done(\"${RGBASE}_s04vldate\")"); \
	CLEANUP_LIST+=("${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${BAMBASE}");\
	CLEANUP_LIST+=("${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.${FQ1EXT}");\
	CLEANUP_LIST+=("${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.${FQ2EXT}");\
	sleep ${DELAY1:-1}s; \
	done; done; \
	echo "Finished submitting per-RG jobs for ${BAMBASE}"; \
	echo -e "05 - MarkDuplicates for ${INPUT_LIST[@]}"
	export MEM=32; \
	export OUT_DIR="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}"; \
	bsub \
	-w ${WAIT_LIST[@]:1} \
    -J "${FULLSM}_s05mrkdup" \
    -o "${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}_s05mrkdup.%J" \
    -u "${EMAIL}" \
    -n1 -W 1440 \
    -M 46000000 \
    -R "rusage[mem=49152]" \
    -q research-hpc \
    -a "docker(vifehe/markduplicates)" \
    ${INPUT_LIST[@]}; \
	CLEANUP_LIST+=("${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.isec-${RUN_TYPE}.markd.bam");\
	echo -e "06 - ValidateSam-III on ${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.bam"
	export MEM=6; \
	export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.bam"; \
	export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
	bsub \
    -w "done(\"${FULLSM}_s05mrkdup\")" \
    -J "${FULLSM}_s06vldate" \
    -o "${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}_s06vldate.%J" \
    -u "${EMAIL}" \
    -n1 -W 1440 \
    -R "rusage[mem=8192]" \
    -q research-hpc \
    -a "docker(vifehe/validatesamref)" \
    entrypoint.sh -IGNORE INVALID_TAG_NM -IGNORE MATE_NOT_FOUND; \
	echo -e "07 - Depth of Coverage"
	unset MEM; \
  	export RUN_TYPE="paddedexome" ;\
	export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
	export PADDED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Padded.GRCh38.bed" ;\
	export COVERED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Covered.GRCh38.bed" ;\
	export GATK_THREADS="4"
  bsub \
    -w "done(\"${FULLSM}_s06vldate\")" \
    -J "${FULLSM}_s07dcvrge" \
    -o "${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}_s07dcvrge.%J" \
    -u "${EMAIL}" \
    -n1 -W 1440 \
    -q research-hpc \
    -a "docker(vifehe/depthofcoverage)" \
    entrypoint.sh; \
	echo -e "08 - Base Recalibrator"
	export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.bam"; \
	export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
	export DBSNP1="${BASE}/Genome_Ref/GRCh38/20190522_bundle/dbsnp_138.hg38.vcf.gz" ;\
	export ONEKGP1="${BASE}/Genome_Ref/GRCh38/20190522_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz" ;\
	export MILLS_GOLD="${BASE}/Genome_Ref/GRCh38/20190522_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" ;\
  bsub \
    -w "done(\"${FULLSM}_s06vldate\")" \
    -J "${FULLSM}_s08bsercl" \
    -o "${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}_s08bsercl.%J" \
    -u "${EMAIL}" \
    -n1 -W 1440 \
    -q research-hpc \
    -a "docker(vifehe/baserecalv2)" \
    entrypoint.sh; \
	CLEANUP_LIST+=("${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.isec-${RUN_TYPE}.markd.recal.bam");\
	echo -e "09 - ValidateSam IV on ${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.recal.bam"
	export MEM=16; \
	export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.recal.bam"; \
	export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
	bsub \
    -w "done(\"${FULLSM}_s08bsercl\")" \
    -J "${FULLSM}_s09vldate" \
    -o "${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}_s09vldate.%J" \
    -u "${EMAIL}" \
    -n1 -W 1360 \
    -R "rusage[mem=18192]" \
    -q research-hpc \
    -a "docker(vifehe/validatesamref)" \
    entrypoint.sh -IGNORE INVALID_VERSION_NUMBER -IGNORE INVALID_TAG_NM -IGNORE MATE_NOT_FOUND; \
	echo -e "10 - Haplotype Caller"
	export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
	export RUN_TYPE="paddedexome" ;\
	export MEM=32; \
	export DBSNP1="${BASE}/Genome_Ref/GRCh38/20190522_bundle/dbsnp_138.hg38.vcf.gz" ;\
	export OUTDIR="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/"
	bsub \
    -w "done(\"${FULLSM}_s09vldate\")" \
    -J "${FULLSM}_s10hcallr" \
    -o "${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}_s10hcallr.%J" \
    -u "${EMAIL}" \
    -n1 -W 2880 \
    -M 46000000 \
    -R "rusage[mem=49152]" \
    -q research-hpc \
    -a "docker(vifehe/haplocallerv2)" \
    entrypoint.sh; \
	echo -e "11 - VariantEval"
	export MEM=8; \
	export GVCF="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.recal.raw.snps.indels.g.vcf.gz"; \
	export CCDS="/gscmnt/gc2645/wgs/Genome_Ref/GRCh38/GRCh38.CCDS.exons.sorted.bed"
	export CLEANUP="${CLEANUP_LIST[@]}"
	bsub \
    -w "done(\"${FULLSM}_s07dcvrge\") && done(\"${FULLSM}_s10hcallr\")" \
    -J "${FULLSM}_s12vnteval" \
    -o "${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}_s12vnteval.%J" \
    -u "${EMAIL}" \
    -n1 -W 1360\
    -q research-hpc \
    -a "docker(vifehe/variantevalref:1.2)" \
    entrypoint.sh; \
	echo "Finished submitting jobs for ${FULLSM}"; \
	sleep ${DELAY2:-10}s; \
	echo "Finished submitting per-bam jobs for ${FULLSM}"; \
	done; \
echo "Finished all samples within ${PR} from ${START} to ${END}"

### 2020-01-03 - Launch first sample to test the pipeline runs ok.
##  2020-01-04 - jobs didn't go through so I cancelled all of them and rethrown again
###JOBS:
53258   fernand RUN   research-h virtual-wor blade17-1-1 *s00rvtbam Jan  4 14:32
53259   fernand PEND  research-h virtual-wor             *s01alnsrt Jan  4 14:32
53260   fernand PEND  research-h virtual-wor             *s02vldate Jan  4 14:32
53261   fernand PEND  research-h virtual-wor             *s03intsct Jan  4 14:32
53262   fernand PEND  research-h virtual-wor             *s04vldate Jan  4 14:32
53263   fernand PEND  research-h virtual-wor             *4b-bamcrm Jan  4 14:32
53264   fernand PEND  research-h virtual-wor             *s05mrkdup Jan  4 14:32
53265   fernand PEND  research-h virtual-wor             *s06vldate Jan  4 14:32
53266   fernand PEND  research-h virtual-wor             *s07dcvrge Jan  4 14:32
53267   fernand PEND  research-h virtual-wor             *s08bsercl Jan  4 14:32
53268   fernand PEND  research-h virtual-wor             *s09vldate Jan  4 14:32
53269   fernand PEND  research-h virtual-wor             *s10hcallr Jan  4 14:32
53270   fernand PEND  research-h virtual-wor             *12vnteval Jan  4 14:32

### 2020-01-06
## checked previous sample and it ran OK
## Launch first set of smaples: start=2; END=100;

### 2020-01-12
## check how many samples have been processed
ls ../02-TRANSIT/TGI_WES/*/*.vcf.gz | wc -l
98
## launch second set

### 2020-01-16 
## check how many samples have run
ls ../02-TRANSIT/TGI_WES/*/*.vcf.gz | wc -l
163

## we have a few that were pending for several days so I killed them ..let's see which samples are missing to process

BASE="/gscmnt/gc2645/wgs"
PR="TGI_WES"
WORKLIST="${PR}-worklist_fullsm.csv"
for FULLSM in $(cat ${WORKLIST}); do \
if [ -f ${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}.*.g.vcf.gz ]; then echo "YES"; else echo -e "GVCF file for sample ${FULLSM} is NOT complete"; fi; done

GVCF file for sample 11_15_15.201^8006302878^TGI_WES is NOT complete
GVCF file for sample 11_8_8.226^8008285859^TGI_WES    11_8_8.226^8008285859^TGI_WES.H22CV^4 is NOT complete
GVCF file for sample 65467^8006305441^TGI_WES is NOT complete
GVCF file for sample 65533^8005873588^TGI_WES is NOT complete
GVCF file for sample ADNI_0814^8006304109^TGI_WES is NOT complete
GVCF file for sample ADNI_1035^8006303895^TGI_WES is NOT complete
GVCF file for sample ADNI_1263^8006303805^TGI_WES is NOT complete
GVCF file for sample ADNI_1295^8006303449^TGI_WES is NOT complete
GVCF file for sample ADNI_1334^8006303427^TGI_WES is NOT complete
GVCF file for sample ADNI_1373^8006303758^TGI_WES is NOT complete
GVCF file for sample ADNI_1409^8006303784^TGI_WES is NOT complete

### Let's launch these samples again. Place them on a new file:
11_15_15.201^8006302878^TGI_WES
11_8_8.226^8008285859^TGI_WES
65467^8006305441^TGI_WES
65533^8005873588^TGI_WES
ADNI_0814^8006304109^TGI_WES
ADNI_1035^8006303895^TGI_WES
ADNI_1263^8006303805^TGI_WES
ADNI_1295^8006303449^TGI_WES
ADNI_1334^8006303427^TGI_WES
ADNI_1373^8006303758^TGI_WES
ADNI_1409^8006303784^TGI_WES

START=1; \
END=11; \
DELAY=10
EMAIL="fernandezv@wustl.edu"
SHELLDROP=0

export BASE="/gscmnt/gc2645/wgs"; \
export WORKDIR="${BASE}/tmp"; \
export THREADS=16; \
export BWA_PIPE_SORT=1; \
export TIMING=1;\
USER="fernandezv";\
PR="TGI_WES"
WORKLIST="${BASE}/WXS_Aquilla/01-RAW/${PR}-worklist_fullsm.csv-missing"; \
LOOKUP="${BASE}/WXS_Aquilla/01-RAW/${PR}-sm_dna_pr_bamfile.csv"; \
LOOKUP_COL_SM=1; \
LOOKUP_COL_DNA=2; \
LOOKUP_COL_PR=3; \
LOOKUP_COL_BAMFILE=4; \
export RUN_TYPE="paddedexome";

### ERRORS:
05 - MarkDuplicates for
-J: No matching job found. Job not submitted.
06 - ValidateSam-III on 11_8_8.226^8008285859^TGI_WES.aln.srt.isec-paddedexome.markd.bam
"11_8_8.226^8008285859^TGI_WES_s05mrkdup: No matching job found. Job not submitted.
07 - Depth of Coverage
"11_8_8.226^8008285859^TGI_WES_s06vldate: No matching job found. Job not submitted.
08 - Base Recalibrator
"11_8_8.226^8008285859^TGI_WES_s06vldate: No matching job found. Job not submitted.
09 - ValidateSam IV on 11_8_8.226^8008285859^TGI_WES.aln.srt.isec-paddedexome.markd.recal.bam
"11_8_8.226^8008285859^TGI_WES_s08bsercl: No matching job found. Job not submitted.
10 - Haplotype Caller
"11_8_8.226^8008285859^TGI_WES_s09vldate: No matching job found. Job not submitted.
11 - VariantEval
"11_8_8.226^8008285859^TGI_WES_s07dcvrge: No matching job found. Job not submitted.
Finished submitting jobs for 11_8_8.226^8008285859^TGI_WES
Finished submitting per-bam jobs for 11_8_8.226^8008285859^TGI_WES
Create directory /gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/TGI_WES/65467^8006305441^TGI_WES
mkdir: cannot create directory ‘/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/TGI_WES/65467^8006305441^TGI_WES’: File exists

### Need to resend 11_8_8.226^8008285859^TGI_WES


##Missing samples & failed steps:
11_15_15.201^8006302878^TGI_WES --> "failed at 06_validate mark duplicates bam"
11_8_8.226^8008285859^TGI_WES --> FINISIHED
65467^8006305441^TGI_WES --> "STOP at markduplicates" --> FINISHED
65533^8005873588^TGI_WES --> "STOP at markduplicates" --> FINISHED
ADNI_0814^8006304109^TGI_WES --> "doesn't go beyond revertbam - in theory it's successfully completed but no bam is generated, but there are some syntax error along the way" --> FINISHED
ADNI_1035^8006303895^TGI_WES --> "doesn't go beyond revertbam - in theory it's successfully completed but no bam is generated, but there are some syntax error along the way" --> FINISHED
ADNI_1263^8006303805^TGI_WES --> "doesn't go beyond revertbam - in theory it's successfully completed but no bam is generated, but there are some syntax error along the way" --> FINISHED
ADNI_1295^8006303449^TGI_WES --> "doesn't go beyond revertbam - in theory it's successfully completed but no bam is generated, but there are some syntax error along the way" --> FINISHED
ADNI_1334^8006303427^TGI_WES --> "doesn't go beyond revertbam - in theory it's successfully completed but no bam is generated, but there are some syntax error along the way" --> FINISHED
ADNI_1373^8006303758^TGI_WES --> "doesn't go beyond revertbam - in theory it's successfully completed but no bam is generated, but there are some syntax error along the way" --> FINISHED
ADNI_1409^8006303784^TGI_WES --> "doesn't go beyond revertbam - in theory it's successfully completed but no bam is generated, but there are some syntax error along the way" --> FINISHED


### 2020-01-21 
## Need to resend the 11_15_15.201^8006302878^TGI_WES
START=1; END=1; 
DELAY=10
EMAIL="fernandezv@wustl.edu"
SHELLDROP=0

export BASE="/gscmnt/gc2645/wgs"; \
export WORKDIR="${BASE}/tmp"; \
export THREADS=16; \
export BWA_PIPE_SORT=1; \
export TIMING=1;\
USER="fernandezv";\
PR="TGI_WES"
WORKLIST="${BASE}/WXS_Aquilla/01-RAW/${PR}-worklist_fullsm.csv-missing"; \
LOOKUP="${BASE}/WXS_Aquilla/01-RAW/${PR}-sm_dna_pr_bamfile.csv"; \
LOOKUP_COL_SM=1; \
LOOKUP_COL_DNA=2; \
LOOKUP_COL_PR=3; \
LOOKUP_COL_BAMFILE=4; \
export RUN_TYPE="paddedexome";

## IT has failed at step06 .. elt's try run it manually to see if we identify the error:
echo -e "06 - ValidateSam-III on ${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.bam"
WORKDIR="/gscmnt/gc2645/wgs/WXS_Aquilla/02-TRANSIT/TGI_WES/11_15_15.201^8006302878^TGI_WES"
BASE="/gscmnt/gc2645/wgs"
	export MEM=6; \
	export BAMFILE="${WORKDIR}/11_15_15.201^8006302878^TGI_WES.aln.srt.isec-paddedexome.markd.bam"; \
	export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
	bsub \
    -J "${FULLSM}_s06vldate" \
    -o "${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}_s06vldate.%J" \
    -u "${EMAIL}" \
    -n1 -W 1440 \
    -R "rusage[mem=8192]" \
    -q research-hpc \
    -a "docker(vifehe/validatesamref)" \
    entrypoint.sh -IGNORE INVALID_TAG_NM -IGNORE MATE_NOT_FOUND; \
Job <1226764> is submitted to queue <research-hpc>

## It has failed but I am not sure why -- see log error at:
fernandezv@virtual-workstation4:/gscmnt/gc2645/wgs/WXS_Aquilla/02-TRANSIT/TGI_WES/11_15_15.201^8006302878^TGI_WES$ cat 11_15_15.201^8006302878^TGI_WES_s06vldate.1226764


# probably I should try to run it interactively to better see what's going on



	echo -e "07 - Depth of Coverage"
	unset MEM; \
  	export RUN_TYPE="paddedexome" ;\
	export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
	export PADDED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Padded.GRCh38.bed" ;\
	export COVERED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Covered.GRCh38.bed" ;\
	export GATK_THREADS="4"
  bsub \
    -w "done(\"${FULLSM}_s06vldate\")" \
    -J "${FULLSM}_s07dcvrge" \
    -o "${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}_s07dcvrge.%J" \
    -u "${EMAIL}" \
    -n1 -W 1440 \
    -q research-hpc \
    -a "docker(vifehe/depthofcoverage)" \
    entrypoint.sh; \
	echo -e "08 - Base Recalibrator"
	export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.bam"; \
	export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
	export DBSNP1="${BASE}/Genome_Ref/GRCh38/20190522_bundle/dbsnp_138.hg38.vcf.gz" ;\
	export ONEKGP1="${BASE}/Genome_Ref/GRCh38/20190522_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz" ;\
	export MILLS_GOLD="${BASE}/Genome_Ref/GRCh38/20190522_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" ;\
  bsub \
    -w "done(\"${FULLSM}_s06vldate\")" \
    -J "${FULLSM}_s08bsercl" \
    -o "${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}_s08bsercl.%J" \
    -u "${EMAIL}" \
    -n1 -W 1440 \
    -q research-hpc \
    -a "docker(vifehe/baserecalv2)" \
    entrypoint.sh; \
	CLEANUP_LIST+=("${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.isec-${RUN_TYPE}.markd.recal.bam");\
	echo -e "09 - ValidateSam IV on ${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.recal.bam"
	export MEM=6; \
	export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.recal.bam"; \
	export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
	bsub \
    -w "done(\"${FULLSM}_s08bsercl\")" \
    -J "${FULLSM}_s09vldate" \
    -o "${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}_s09vldate.%J" \
    -u "${EMAIL}" \
    -n1 -W 360 \
    -R "rusage[mem=8192]" \
    -q research-hpc \
    -a "docker(vifehe/validatesamref)" \
    entrypoint.sh -IGNORE INVALID_VERSION_NUMBER -IGNORE INVALID_TAG_NM -IGNORE MATE_NOT_FOUND; \
	echo -e "10 - Haplotype Caller"
	export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
	export RUN_TYPE="paddedexome" ;\
	export MEM=32; \
	export DBSNP1="${BASE}/Genome_Ref/GRCh38/20190522_bundle/dbsnp_138.hg38.vcf.gz" ;\
	export OUTDIR="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/"
	bsub \
    -w "done(\"${FULLSM}_s09vldate\")" \
    -J "${FULLSM}_s10hcallr" \
    -o "${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}_s10hcallr.%J" \
    -u "${EMAIL}" \
    -n1 -W 2880 \
    -M 46000000 \
    -R "rusage[mem=49152]" \
    -q research-hpc \
    -a "docker(vifehe/haplocallerv2)" \
    entrypoint.sh; \
	echo -e "11 - VariantEval"
	export MEM=2; \
	export GVCF="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.recal.raw.snps.indels.g.vcf.gz"; \
	export CCDS="/gscmnt/gc2645/wgs/Genome_Ref/GRCh38/GRCh38.CCDS.exons.sorted.bed"
	export CLEANUP="${CLEANUP_LIST[@]}"
	bsub \
    -w "done(\"${FULLSM}_s07dcvrge\") && done(\"${FULLSM}_s10hcallr\")" \
    -J "${FULLSM}_s12vnteval" \
    -o "${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}_s12vnteval.%J" \
    -u "${EMAIL}" \
    -n1 -W 360\
    -q research-hpc \
    -a "docker(vifehe/variantevalref:1.2)" \
    entrypoint.sh; \
	echo "Finished submitting jobs for ${FULLSM}"; \
	sleep ${DELAY2:-10}s; \
	echo "Finished submitting per-bam jobs for ${FULLSM}"; \
	done; \
echo "Finished all samples within ${PR} from ${START} to ${END}"