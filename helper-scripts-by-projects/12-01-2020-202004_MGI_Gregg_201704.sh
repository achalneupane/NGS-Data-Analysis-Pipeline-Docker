################################
# 2020-02-03 - MGI-Gregg 
################################

SOURCEDIR="/40/AD/AD_Seq_Data/01.-RawData/MGI_Gregg_201704/bam/"

## This project has 81 bams - but these are coded in non format" Study_FULLSM_RGBASE.bam", instead sth like:
ls MGI_Gregg_201704/bam/*.bam | head
MGI_Gregg_201704/bam/049da3c49e9f47d8991cb7ce06ff23d4.bam
MGI_Gregg_201704/bam/04c9e42871e74cbc99dd669a2e6e6c02.bam
MGI_Gregg_201704/bam/0752f042f11544878b3d02a04febc508.bam
MGI_Gregg_201704/bam/07f9ee31825c4a5480fd6f0888984270.bam

## the file with the recoding info is:
MGI_Gregg_201704/hexhash_mgiid_fullsm_sm.csv
head MGI_Gregg_201704/hexhash_mgiid_fullsm_sm.csv
bamfile,mgiid,fullsm,sm
f0df9157909641358f9eb6558584f03a.bam,H_XV-65491-8015452875,65491^8015452875^MGI_Gregg_201704,65491
30e5e19ff29641729a42489345e872e6.bam,H_XV-64612-0024080575,64612^0024080575^MGI_Gregg_201704,64612
07f9ee31825c4a5480fd6f0888984270.bam,H_XV-66152-8008291106,66152^8008291106^MGI_Gregg_201704,66152
0ff837796c614f9d8204e4af34057cc1.bam,H_XV-61924-8011736582,61924^8011736582^MGI_Gregg_201704,61924

## From this file generate the FULLSM list and the metafile as we usually use it: sm_dna_pr_bamfile
cut -d "," -f3 MGI_Gregg_201704/hexhash_mgiid_fullsm_sm.csv > MGI_Gregg_201704-worklist-fullsm.csv
cut -d "," -f4 MGI_Gregg_201704/hexhash_mgiid_fullsm_sm.csv | tail -n+2 > MGI_Gregg_201704-sm.list
cut -d "," -f3 MGI_Gregg_201704/hexhash_mgiid_fullsm_sm.csv | cut -d "^" -f 2 | tail -n+2 > MGI_Gregg_201704-dna.list
cut -d "," -f3 MGI_Gregg_201704/hexhash_mgiid_fullsm_sm.csv | cut -d "^" -f 3 | tail -n+2 > MGI_Gregg_201704-pr.list
cut -d "," -f1 MGI_Gregg_201704/hexhash_mgiid_fullsm_sm.csv | tail -n+2 > MGI_Gregg_201704-bam.list
paste -d "," MGI_Gregg_201704-sm.list MGI_Gregg_201704-dna.list MGI_Gregg_201704-pr.list MGI_Gregg_201704-bam.list > MGI_Gregg_201704-sm_dna_pr_bamfile.csv

##At MGI:
mkdir MGI_Gregg_201704
SOURCE="/40/AD/AD_Seq_Data/01.-RawData"
PR="MGI_Gregg_201704"
metafile="${PR}-sm_dna_pr_bamfile.csv"
DEST="/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW"
scp -p ${USER}@fenix.psych.wucon.wustl.edu:${SOURCE}/${metafile} ${DEST}/;

## on MGI - Generate the workfile
cut -d, -f1-3 ${DEST}/${metafile} | tr "," "^" | sort -u > ${DEST}/${PR}-worklist_fullsm.csv
WORKLIST="${DEST}/${PR}-worklist_fullsm.csv"


# Check the missing samples on fenix
comm -23 <(sort /home/achal/01.-RawData/missing_check/MGI_Gregg_all_bam.txt| cut -d. -f1) <(sort /home/achal/01.-RawData/missing_check/MGI_Gregg_201704.txt)
## 68317^8012957245^MGI_Gregg_201704
## 72173^8011736646^MGI_Gregg_201704

# achal@virtual-workstation1:/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW$ cat -n $WORKLIST | grep 68317^8012957245^MGI_Gregg_201704
#     82  68317^8012957245^MGI_Gregg_201704
# achal@virtual-workstation1:/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW$ cat -n $WORKLIST | grep 72173^8011736646^MGI_Gregg_201704
#     83  72173^8011736646^MGI_Gregg_201704

## Run these samples using pipeline as defined in  MGI-BAM-UPLOAD-FASTQ-RUN-v1.sh
## Sets to define:
START=82; \
END=82; \
DELAY=10
EMAIL="achal@wustl.edu"
SHELLDROP=0
# export LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active/achal/WXS_Aquilla/01-RAW/all_bam_files:/storage1/fs1/cruchagac/Active/achal/WXS_Aquilla/01-RAW/all_bam_files"; \
export BASE="/gscmnt/gc2645/wgs"; \
export WORKDIR="${BASE}/tmp"; \
export THREADS=16; \
export BWA_PIPE_SORT=1; \
export TIMING=1;\
USER="achal";\
PR="MGI_Gregg_201704"
WORKLIST="${BASE}/WXS_Aquilla/01-RAW/${PR}-worklist_fullsm.csv"; \
LOOKUP="${BASE}/WXS_Aquilla/01-RAW/${PR}-sm_dna_pr_bamfile.csv"; \
LOOKUP_COL_SM=1; \
LOOKUP_COL_DNA=2; \
LOOKUP_COL_PR=3; \
LOOKUP_COL_BAMFILE=4; \
export RUN_TYPE="paddedexome";


for FULLSM in $(sed -n "${START},${END}p" "${WORKLIST}"); do \
	export FULLSM="${FULLSM}"; \
	SOURCE="/40/AD/AD_Seq_Data/01.-RawData/MGI_Gregg_201704/bam";\
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
	BAMFILE="${BASE}/WXS_Aquilla/01-RAW/${PR}/bam/${BAMBASE}"; \
	echo -e "00a - Uploading BAMFILE ${BAMFILE} for sample ${FULLSM}";\
	rsync -avh ${USER}@fenix.psych.wucon.wustl.edu:${SOURCE}/${BAMBASE%.*}.* ${DEST}/;
	echo -e "00b - Starting revertbam for ${FULLSM}\n${BAMFILE}\n${SM}\n${DNA}\n${PR}"; \
	export BAMFILE="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${BAMBASE}"; \
	export MEM=6; \
	export CREATE_RGFILE=1; \
	export DEBUG=1; \
	export REMOVE_INPUT=1; \
	export OUT_DIR="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}" ;\
	echo -e "OUTPUT is in "${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${BAMBASE}_s00rvtbam.%J" \n OUTDIR is ${OUT_DIR}";\
	bsub \
	-J "${FULLSM}_s00rvtbam" \
	-o "${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${BAMBASE}_s00rvtbam.%J" \
	-u "${EMAIL}" \
	-n1 -W 1440 \
	-R "rusage[mem=8192]" \
	-q research-hpc \
	-a 'docker(vifehe/revertbam:1.1)' \
	/bin/bash; \
	echo -e "00c - reading readgroups RG of BAMfile"
	SAMTOOLS="/gscmnt/gc2645/wgs/variant_calling/samtools"; \
	BAMFILE="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${BAMBASE}"; \
	IFS=$'\n' RGS=($(${SAMTOOLS} view -H "${BAMFILE}" | grep "^@RG")); \
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
	export MEM=32; \
	export CLEANUP=1; 
	export REMOVE_INPUT=1; \
	bsub \
	-w "done(\"${FULLSM}_s00rvtbam\")" \
	-J "${RGBASE}_s01alnsrt" \
	-o "${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}_s01alnsrt.%J" \
	-u "${EMAIL}" \
	-n ${THREADS} -W 2160 \
	-M 46000000 \
	-R "rusage[mem=49152]" \
	-q research-hpc \
	-a 'docker(vifehe/bwaref)' \
	/bin/bash; \
	echo -e "02 - Perform ValidateSam-I on newly aligned bamfile ${BAMFILE}"
	export MEM=6; \
	export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.bam"; \
	echo ${BAMFILE}; \
	export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
	bsub \
	-w "done(\"${RGBASE}_s01alnsrt\")" \
	-J "${RGBASE}_s02vldate" \
	-o "${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}_s02vldate.%J" \
	-u "${EMAIL}" \
	-n1 -W 360 \
	-R "rusage[mem=8192]" \
	-q research-hpc \
	-a "docker(vifehe/validatesamref)" \
	entrypoint.sh -IGNORE INVALID_VERSION_NUMBER -IGNORE INVALID_TAG_NM; \
	echo -e "03 - perform intersectbed"
	export RUN_TYPE="paddedexome"; \
	export BEDFILE="${BASE}/Genome_Ref/GRCh38/Capture_Padded.GRCh38.bed"; \
	export COVERED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Covered.GRCh38.bed" ;\
	export PADDED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Padded.GRCh38.bed" ;\
	export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
	export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.bam"; \
echo -e "03 - Starting intersectbed per sample ${FULLSM} and RGBASE ${RGBASE}\nBAMFILE:${BAMFILE}";\
	bsub \
	-w "done(\"${RGBASE}_s02vldate\")" \
	-J "${RGBASE}_s03intsct" \
	-o "${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}_s03intsct.%J" \
	-u "${EMAIL}" \
	-n1 -W 360 \
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
	-R "rusage[mem=8192]" \
	-q research-hpc \
	-a "docker(vifehe/validatesamref)" \
	entrypoint.sh -IGNORE INVALID_VERSION_NUMBER -IGNORE INVALID_TAG_NM -IGNORE MATE_NOT_FOUND; \
	echo -e "04b - Starting bamtocram for FULLSM ${FULLSM} and BAMFILE ${FULLSM}.aln.srt.bam";\
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
    -a "docker(vifehe/variantevalref)" \
    entrypoint.sh; \
	echo "Finished submitting jobs for ${FULLSM}"; \
	sleep ${DELAY2:-10}s; \
	echo "Finished submitting per-bam jobs for ${FULLSM}"; \
	done; \
echo "Finished all samples within ${PR} from ${START} to ${END}"

### 2020-02-03 - Launch first sample to test the pipeline runs ok.
00b - Starting revertbam for 1380^8015455316^MGI_Gregg_201704
/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/MGI_Gregg_201704/bam/3cfd36d731054d96912c810a5473596c.bam
1380
8015455316
MGI_Gregg_201704
OUTPUT is in /gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/MGI_Gregg_201704/1380^8015455316^MGI_Gregg_201704/3cfd36d731054d96912c810a5473596c.bam_s00rvtbam.%J
 OUTDIR is /gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/MGI_Gregg_201704/1380^8015455316^MGI_Gregg_201704
Job <1713570> is submitted to queue <research-hpc>.
00c - reading readgroups RG of BAMfile
01 - Starting bwa per sample 1380^8015455316^MGI_Gregg_201704 and RGBASE 1380^8015455316^MGI_Gregg_201704.2897237469
FQ1:
FQ2:
RGFILE:/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/MGI_Gregg_201704/1380^8015455316^MGI_Gregg_201704/1380^8015455316^MGI_Gregg_201704.2897237469.rgfile
DNA:{8015455316}
PR:MGI_Gregg_201704
Job <1713571> is submitted to queue <research-hpc>.
02 - Perform ValidateSam-I on newly aligned bamfile
/gscmnt/gc2645/wgs/WXS_Aquilla/02-TRANSIT/MGI_Gregg_201704/1380^8015455316^MGI_Gregg_201704/1380^8015455316^MGI_Gregg_201704.2897237469.aln.srt.bam
Job <1713572> is submitted to queue <research-hpc>.
03 - perform intersectbed
Job <1713573> is submitted to queue <research-hpc>.
04 - Perform ValidateSam-II on 1380^8015455316^MGI_Gregg_201704.2897237469.aln.srt.isec-paddedexome.bam
Job <1713574> is submitted to queue <research-hpc>.
04b - Starting bamtocram for FULLSM 1380^8015455316^MGI_Gregg_201704 and BAMFILE 1380^8015455316^MGI_Gregg_201704.aln.srt.bam
Job <1713576> is submitted to queue <research-hpc>.
Finished submitting per-RG jobs for 3cfd36d731054d96912c810a5473596c.bam
05 - MarkDuplicates for /gscmnt/gc2645/wgs/WXS_Aquilla/02-TRANSIT/MGI_Gregg_201704/1380^8015455316^MGI_Gregg_201704/1380^8015455316^MGI_Gregg_201704.2897237469.aln.srt.isec-paddedexome.bam
Job <1713577> is submitted to queue <research-hpc>.
06 - ValidateSam-III on 1380^8015455316^MGI_Gregg_201704.aln.srt.isec-paddedexome.markd.bam
Job <1713578> is submitted to queue <research-hpc>.
07 - Depth of Coverage
Job <1713579> is submitted to queue <research-hpc>.
08 - Base Recalibrator
Job <1713580> is submitted to queue <research-hpc>.
09 - ValidateSam IV on 1380^8015455316^MGI_Gregg_201704.aln.srt.isec-paddedexome.markd.recal.bam
Job <1713581> is submitted to queue <research-hpc>.
10 - Haplotype Caller
Job <1713583> is submitted to queue <research-hpc>.
11 - VariantEval
Job <1713584> is submitted to queue <research-hpc>.
Finished submitting jobs for 1380^8015455316^MGI_Gregg_201704
Finished submitting per-bam jobs for 1380^8015455316^MGI_Gregg_201704
Finished all samples within MGI_Gregg_201704 from 1 to 1


### 2020-02-03 - This sample has processed successfuly so I will launch the remainign ones:
START=2; END=81


