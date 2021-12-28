### 2019-09 - README file for Genentech-WES sequences processing over at MGI

##############
## 2019-09-19
##############

### In Fenix

# generate metadata file for the PROJECT in question
HOMEDIR="/40/AD/AD_Seq_Data/01.-RawData"
PR="201812_MGI_DIANWGS_REDCLOUD"

##List directories within the PR directory
## Macrogen like DIAN have several bams that need to be merged
metafile="${PR}-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv"
vim ${metafile}   ## copy here the metadata for the DIAN project from the excel file
sed -i 's/\t/,/g' ${metafile}

## this file contains the follwoing fields, comma separated
# SM_(MGI),SM,DNA,PR_(project),RGBASE_(MGI sample ID),final.cram,SMDIR,CRAM_location

## Upload this file to MGI
scp ${metafile} fernandezv@virtual-workstation1.gsc.wustl.edu:/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/

## Log onto MGI
#####################################
#### 1 - upload samples using rsync
#####################################

BASE="/gscmnt/gc2645/wgs/"

## I have uplaoded the list with over 547 samples; AN: testing for to samples only
wc -l 201812_MGI_DIANWGS_REDCLOUD-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv
547 201812_MGI_DIANWGS_REDCLOUD-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv

# In a previous test before I ran 34 samples from this list,  
 wc -l 201812_MGI_DIANWGS_REDCLOUD-worklist_fullsm.csv
34 201812_MGI_DIANWGS_REDCLOUD-worklist_fullsm.csv

### I need to make a list of the smaples that I did not process then to process them now.
#cut -d "," -f 7 201812_MGI_DIANWGS_REDCLOUD-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv > LIST1
cat 201812_MGI_DIANWGS_REDCLOUD-worklist_fullsm.csv > 201812_MGI_DIANWGS_REDCLOUD-worklist_fullsm.csv2
#cat LIST1 201812_MGI_DIANWGS_REDCLOUD-worklist_fullsm.csv | sort | uniq -c | awk '{if ($1==1) print $2}' > 201812_MGI_DIANWGS_REDCLOUD-worklist_fullsm.csv2

## My new working list will be 
worklist="201812_MGI_DIANWGS_REDCLOUD-worklist_fullsm.csv2"
## This list contains 513 samples
PR="201812_MGI_DIANWGS_REDCLOUD"
USER="achal"

metafile="${PR}-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv"

# missing sample:
cat $metafile | grep 010-0009-0155433^0010007652^201812_MGI_DIANWGS_REDCLOUD > ${PR}-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc_missing_achal.csv
metafile=${PR}-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc_missing_achal.csv

# a - test that we are reading things correctly
while read -r line; do \
CRAMLOC="$(echo $line | cut -d "," -f 8)"
MGI_ID="$(echo $line | cut -d "," -f 1)"
SMDIR="$(echo $line | cut -d "," -f 7)"
RGBASE="$(echo $line | cut -d "," -f 5)"
echo -e "${CRAMLOC}\n${MGI_ID}\n${SMDIR}\n${RGBASE}"
done < ${metafile} 
"/40/AD/AD_Seq_Data/01.-RawData/201812_MGI_DIANWGS_REDCLOUD/01.-RawData/xfer.genome.wustl.edu/gxfer1/40905687920482/2857909.491.053019/2debfe46db5c48d0a3f7c4e9dc96838c.final.cram
2debfe46db5c48d0a3f7c4e9dc96838c
094-0103-0082230^0010007365^201812_MGI_DIANWGS
094-0103-0082230^0010007365^201812_MGI_DIANWGS.2debfe46db5c48d0a3f7c4e9dc96838c"

# b - upload files: 
vi 201812_MGI_DIANWGS_REDCLOUD-worklist_fullsm_missing_achal.csv
010-0009-0155433^0010007652^201812_MGI_DIANWGS_REDCLOUD

for FULLSM in $(cat 201812_MGI_DIANWGS_REDCLOUD-worklist_fullsm_missing_achal.csv); do \
echo ${FULLSM}; \
line="$(grep ${FULLSM} ${metafile})";\
echo ${line};\
CRAMLOC="$(echo $line | cut -d "," -f 8)";\
MGI_ID="$(echo $line | cut -d "," -f 1)";\
SMDIR="$(echo $line | cut -d "," -f 7)";\
RGBASE="$(echo $line | cut -d "," -f 5)"; \
echo -e "CRAMLOC:${CRAMLOC}\nMGI_ID:${MGI_ID}\nSMDIR:${SMDIR}\nRGBASE:${RGBASE}"; \
DEST="/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/201812_MGI_DIANWGS_REDCLOUD/${SMDIR}"; \
mkdir -p ${DEST}
rsync -avh ${USER}@fenix.psych.wucon.wustl.edu:${CRAMLOC%.*}.* ${DEST}/;\
done 

### Check that all files have uploaded
fernandezv@virtual-workstation1:/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW$ ls 201812_MGI_DIANWGS_REDCLOUD/*/*.final.cram | wc -l
546

#####################################
### 2 - Run vifehe/cramtobam for all samples using metafile data for input
#####################################

# the original metafile is: 201812_MGI_DIANWGS_REDCLOUD-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv
## I will be runnign them on batches of 50 so when I generate the fastq I do not overload the system.
201812_MGI_DIANWGS_REDCLOUD-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv-1-5 --> DONE
201812_MGI_DIANWGS_REDCLOUD-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv-6-50  --> DONE
201812_MGI_DIANWGS_REDCLOUD-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv-51-100 --> DONE
201812_MGI_DIANWGS_REDCLOUD-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv-101-150 --> DONE
201812_MGI_DIANWGS_REDCLOUD-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv-151-200 --> DONE
201812_MGI_DIANWGS_REDCLOUD-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv-201-250 --> DONE 
201812_MGI_DIANWGS_REDCLOUD-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv-251-300 --> DOING 
201812_MGI_DIANWGS_REDCLOUD-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv-301-350  
201812_MGI_DIANWGS_REDCLOUD-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv-351-400
201812_MGI_DIANWGS_REDCLOUD-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv-401-450
201812_MGI_DIANWGS_REDCLOUD-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv-451-500
201812_MGI_DIANWGS_REDCLOUD-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv-501-547


# # kill jobs by pattern:
# bjobs | grep RUN | cut -d\  -f1 | xargs bkill 
# bjobs | grep PEN | cut -d\  -f1 | xargs bkill 

EMAIL="achal@wustl.edu"
# metafile="201812_MGI_DIANWGS_REDCLOUD-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv-251-300";\
# AN: using only two samples
# metafile="201812_MGI_DIANWGS_REDCLOUD-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv";\
metafile=${PR}-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc_missing_achal.csv; \
while read -r line; do \
	SMDIR="$(echo $line | cut -d "," -f 7)"
	PR="201812_MGI_DIANWGS_REDCLOUD"
	MGI_ID="$(echo $line | cut -d "," -f 1)"
	RUN_TYPE="genome"
	SM="$(echo $line | cut -d "," -f 5)"
	CRAMFILE="${BASE}/WXS_Aquilla/01-RAW/${PR}/${SMDIR}/${MGI_ID}.final.cram"
	OUTDIR="${BASE}/WXS_Aquilla/01-RAW/${PR}/${SMDIR}/"
	REF="${BASE}/Genome_Ref/GRCh38/GRCh38_MGI_20190813/all_sequences.fa"
	echo -e "${MGI_ID}\n${SMDIR}"
	export BASE="/gscmnt/gc2645/wgs"
	export WORKDIR="${BASE}/tmp";
	export RAWDIR="${BASE}/WXS_Aquilla/01-RAW/${PR}/${SMDIR}"
	export CRAMFILE="${BASE}/WXS_Aquilla/01-RAW/${PR}/${SMDIR}/${MGI_ID}.final.cram"; \
    export BAMFILE="${BASE}/WXS_Aquilla/01-RAW/${PR}/${SMDIR}/${MGI_ID}.bam"
	export OUTDIR="${BASE}/WXS_Aquilla/01-RAW/${PR}/${SMDIR}"; \
	export RUN_TYPE="genome"; \
	export REF="${BASE}/Genome_Ref/GRCh38/GRCh38_MGI_20190813/all_sequences.fa"
    echo -e "${CRAMFILE}\n"; \
	bsub \
      -J "${MGI_ID}_s00cram2bam" \
      -o "${BASE}/WXS_Aquilla/04-LOGS/${MGI_ID}_s00cram2bam.%J" \
      -u "${EMAIL}" \
      -n1 -W 1360 \
      -q research-hpc \
      -a "docker(vifehe/cram2bam)" \
      entrypoint.sh; 
done < ${metafile}

# job 3962674
### Test if we got Bam files for all samples
ls 201812_MGI_DIANWGS_REDCLOUD/*/*.final.cram | wc -l
546
ls 201812_MGI_DIANWGS_REDCLOUD/*/*.bam | wc -l
46
# This is only after runing bwa- the big chunk below
ls ../02-TRANSIT/201812_MGI_DIANWGS_REDCLOUD/*/*.bam | wc -l
331


#####################################
### 3 - Start pipeline from fastq to variant calling
#####################################

## Need to generate the MGIID-WASHUID lookup file plus the working file for each of the previous 50x sets;
for LOOKUP in $(ls 201812_MGI_DIANWGS_REDCLOUD-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv*); do \
awk -F"," '{print $1","$2","$3","$4","$1".genome.bam"}' ${LOOKUP} > 201812_MGI_DIANWGS_REDCLOUD-MGIID-SMID-DNA-PR-BAMFILE.${LOOKUP##*.}; done

## Start with "complete pipeline following John's script
for LOOKUP in $(ls 201812_MGI_DIANWGS_REDCLOUD-MGIID-SMID-DNA-PR-BAMFILE.csv*); do cut -d, -f2-4 "${LOOKUP}" | tr "," "^" | sort -u > ${PR}-worklist_fullsm.${LOOKUP##*.}; done

FULLSM="${PR}-worklist_fullsm.csv"
# FBASE may vary depending on the project, due to inconsistency in file structure across projects (e.g. FULLSM/ bams/ bam/ bam/FULLSM)
#   Ideally would query BAMFILE from master LOOKUP with all files across all projects

while read -r line; do \
#FULLSM="$(sed -n '1p' "${LOOKUP}" | cut -d, -f2-4 | tr "," "^")"
FULLSM="$(echo ${line} | cut -d, -f2-4 | tr "," "^")"; \
CRAMFILE="$(echo "${line}" | cut -d"," -f1)"; \
SM="$(echo "${FULLSM}" | cut -d^ -f1)"; \
DNA="$(echo "${FULLSM}" | cut -d^ -f2)"; \
PR="$(echo "${FULLSM}" | cut -d^ -f3)"; \
BAMFILE="$(echo ${line} | cut -d, -f5)"; \
FBASE="/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/${PR}/${FULLSM}"; \
FILE_BAM="${FBASE}/${BAMFILE}"; \
echo -e "${FULLSM}\n${CRAMFILE}\n${SM}\n${DNA}\n${PR}\n${BAMFILE}\n${FBASE}\n${FILE_BAM}" ; \
done < ${LOOKUP}

echo -e "${FULLSM}\n${CRAMFILE}\n${SM}\n${DNA}\n${PR}\n${BAMFILE}\n${FBASE}\n${FILE_BAM}" 
010-0004-0085842^0010007564^201812_MGI_DIANWGS_REDCLOUD
3e2b8f9520e144c29debedf0e74d1a8f
010-0004-0085842
0010007564
201812_MGI_DIANWGS_REDCLOUD
3e2b8f9520e144c29debedf0e74d1a8f.genome.bam
/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/201812_MGI_DIANWGS_REDCLOUD/010-0004-0085842^0010007564^201812_MGI_DIANWGS_REDCLOUD
/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/201812_MGI_DIANWGS_REDCLOUD/010-0004-0085842^0010007564^201812_MGI_DIANWGS_REDCLOUD/3e2b8f9520e144c29debedf0e74d1a8f.genome.bam


#################################################################
##### Run pipeline from bam back to fastq and then call variants
#################################################################
PR="201812_MGI_DIANWGS_REDCLOUD"

201812_MGI_DIANWGS_REDCLOUD-MGIID-SMID-DNA-PR-BAMFILE.csv-1-5 --> RUNNING Spet 23rd 
201812_MGI_DIANWGS_REDCLOUD-MGIID-SMID-DNA-PR-BAMFILE.csv-6-50
201812_MGI_DIANWGS_REDCLOUD-MGIID-SMID-DNA-PR-BAMFILE.csv-51-100
201812_MGI_DIANWGS_REDCLOUD-MGIID-SMID-DNA-PR-BAMFILE.csv-101-150
201812_MGI_DIANWGS_REDCLOUD-MGIID-SMID-DNA-PR-BAMFILE.csv-151-200
201812_MGI_DIANWGS_REDCLOUD-MGIID-SMID-DNA-PR-BAMFILE.csv-201-250
201812_MGI_DIANWGS_REDCLOUD-MGIID-SMID-DNA-PR-BAMFILE.csv-251-300
201812_MGI_DIANWGS_REDCLOUD-MGIID-SMID-DNA-PR-BAMFILE.csv-301-350
201812_MGI_DIANWGS_REDCLOUD-MGIID-SMID-DNA-PR-BAMFILE.csv-351-400
201812_MGI_DIANWGS_REDCLOUD-MGIID-SMID-DNA-PR-BAMFILE.csv-401-450
201812_MGI_DIANWGS_REDCLOUD-MGIID-SMID-DNA-PR-BAMFILE.csv-451-500
201812_MGI_DIANWGS_REDCLOUD-MGIID-SMID-DNA-PR-BAMFILE.csv-501-547


# If needed incorporate this at line 159, after declaring ${FULLSM}
#if [ -d ${BASE}/WXS_Aquilla/03-FINAL/${PR}/${FULLSM} ] || [ -d ${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM} ]; then echo "Skipping ${FULLSM}, already run(ning)"; continue; else echo "RUN ${FULLSM}"; fi; \
# missing sample
010-0009-0155433^0010007652^201812_MGI_DIANWGS_REDCLOUD

START=1; \
END=1; \
DELAY=10
EMAIL="achal@wustl.edu"
SHELLDROP=0

export BASE="/gscmnt/gc2645/wgs"; \
export WORKDIR="${BASE}/tmp"; \
export THREADS=16; \
export BWA_PIPE_SORT=1; \
export TIMING=1;\
PR="201812_MGI_DIANWGS_REDCLOUD"; \
WORKLIST="${BASE}/WXS_Aquilla/01-RAW/${PR}-worklist_fullsm_missing_achal.csv"; \
LOOKUP="${BASE}/WXS_Aquilla/01-RAW/${PR}-MGIID-SMID-DNA-PR-BAMFILE_missing_achal.csv"; \
LOOKUP_COL_SM=2; \
LOOKUP_COL_DNA=3; \
LOOKUP_COL_PR=4; \
LOOKUP_COL_BAMFILE=5; \
export RUN_TYPE="paddedexome"; \
for FULLSM in $(sed -n "${START},${END}p" "${WORKLIST}"); do \
	export FULLSM="${FULLSM}"; \
	SM="$(echo "${FULLSM}" | cut -d^ -f1)"; \
	IFS=$'\n' export DNA=($(awk -F, "\$${LOOKUP_COL_SM} == \"${SM}\"" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_DNA} | sort -u)); \
	if [ ${#DNA[@]} -gt 1 ]; then echo "Warning, \${DNA} not unique for ${SM} (n=${#DNA[@]}: ${DNA[@]})"; fi; \
	DNA="$(echo "${FULLSM}" | cut -d^ -f2)"; \
	PR="$(echo "${FULLSM}" | cut -d^ -f3)"; \
	export OUT_DIR="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}"; \
	INPUT_LIST=(); \
	WAIT_LIST=(); \
	for BAMBASE in $(grep "${SM},${DNA},${PR}" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_BAMFILE}); do \
	echo -e "00a - In this round we will analyze ${FULLSM}\n${BAMFILE}\n${SM}\n${DNA}\n${PR}"; \
	export BAMFILE="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${BAMBASE}"; \
	export MEM=16; \
	export CREATE_RGFILE=1; \
	export DEBUG=1; \
	export REMOVE_INPUT=0; \
	export OUT_DIR="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/"; \
	bsub \
	-J "${FULLSM}_s00rvtbam" \
	-o "${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${FULLSM}.${BAMBASE}_s00rvtbam.%J" \
	-u "${EMAIL}" \
	-n1 -W 1440 \
	-R "rusage[mem=18192]" \
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
	export FQ1="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.${FQ1EXT}"; \
	export FQ2="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.${FQ2EXT}"; \
	export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
	export MEM=32; \
	export CLEANUP=1; 
	export REMOVE_INPUT=1; \
	bsub \
	-w "done(\"${FULLSM}_s00rvtbam\")" \
	-J "${RGBASE}_s01alnsrt" \
	-o "${BASE}/WXS_Aquilla/04-LOGS/${RGBASE}_s01alnsrt.%J" \
	-u "${EMAIL}" \
	-n ${THREADS} -W 2160 \
	-M 46000000 \
	-R "rusage[mem=49152]" \
	-q research-hpc \
	-a 'docker(vifehe/bwaref)' \
	/bin/bash; \
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
     -o "${BASE}/WXS_Aquilla/04-LOGS/${RGBASE}_s03intsct.%J"\
     -u "${EMAIL}" \
     -n1 -W 1360 \
     -q research-hpc \
     -a "docker(vifehe/intersectbedref)" \
     entrypoint.sh; \
	export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.isec-${RUN_TYPE}.bam"; \
	export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
	bsub \
	-w "done(\"${RGBASE}_s03intsct\")" \
	-J "${RGBASE}_s04vldate" \
	-o "${BASE}/WXS_Aquilla/04-LOGS/${RGBASE}_s04vldate.%J" \
	-u "${EMAIL}" \
	-n1 -W 1360 \
	-R "rusage[mem=18192]" \
	-q research-hpc \
	-a "docker(vifehe/validatesamref)" \
	entrypoint.sh -IGNORE INVALID_VERSION_NUMBER -IGNORE INVALID_TAG_NM -IGNORE MATE_NOT_FOUND; \
	INPUT_LIST+=("${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.isec-${RUN_TYPE}.bam"); \
	WAIT_LIST+=("&&" "done(\"${RGBASE}_s04vldate\")"); \
	CLEANUP_LIST+=("${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.${FQ1EXT}");\
	CLEANUP_LIST+=("${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${RGBASE}.${FQ2EXT}");\
	sleep ${DELAY1:-1}s; \
	done; \
	echo "Finished submitting per-RG jobs for ${BAMBASE}"; \
	export MEM=32; \
	export OUT_DIR="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}"; \
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
 	unset MEM; \
  	export RUN_TYPE="paddedexome" ;\
	export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
	export PADDED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Padded.GRCh38.bed" ;\
	export COVERED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Covered.GRCh38.bed" ;\
	export GATK_THREADS="4"
  bsub \
    -w "done(\"${FULLSM}_s07vldate\")" \
    -J "${FULLSM}_s08dcvrge" \
    -o "${BASE}/WXS_Aquilla/04-LOGS/${RGBASE}_s08dcvrge.%J" \
    -u "${EMAIL}" \
    -n1 -W 1440 \
    -q research-hpc \
    -a "docker(vifehe/depthofcoverage)" \
    entrypoint.sh; \
	export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.bam"; \
	export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
	export DBSNP1="${BASE}/Genome_Ref/GRCh38/20190522_bundle/dbsnp_138.hg38.vcf.gz" ;\
	export ONEKGP1="${BASE}/Genome_Ref/GRCh38/20190522_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz" ;\
	export MILLS_GOLD="${BASE}/Genome_Ref/GRCh38/20190522_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" ;\
  bsub \
    -w "done(\"${FULLSM}_s07vldate\")" \
    -J "${FULLSM}_s08bsercl" \
    -o "${BASE}/WXS_Aquilla/04-LOGS/${RGBASE}_s08bsercl.%J" \
    -u "${EMAIL}" \
    -n1 -W 1440 \
    -q research-hpc \
    -a "docker(vifehe/baserecalv2)" \
    entrypoint.sh; \
	export MEM=6; \
	export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.recal.bam"; \
	export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
	bsub \
    -w "done(\"${FULLSM}_s08bsercl\")" \
    -J "${FULLSM}_s09vldate" \
    -o "${BASE}/WXS_Aquilla/04-LOGS/${FULLSM}_s09vldate.%J" \
    -u "${EMAIL}" \
    -n1 -W 360 \
    -R "rusage[mem=8192]" \
    -q research-hpc \
    -a "docker(vifehe/validatesamref)" \
    entrypoint.sh -IGNORE INVALID_VERSION_NUMBER -IGNORE INVALID_TAG_NM -IGNORE MATE_NOT_FOUND; \
	export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
	export RUN_TYPE="paddedexome" ;\
	export MEM=32; \
	export BEDFILE="${BASE}/Genome_Ref/GRCh38/Capture_Padded.GRCh38.bed"; \
	export PADDED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Padded.GRCh38.bed" ;\
	export COVERED_BED="${BASE}/Genome_Ref/GRCh38/Capture_Covered.GRCh38.bed" ;\
	export DBSNP1="${BASE}/Genome_Ref/GRCh38/20190522_bundle/dbsnp_138.hg38.vcf.gz" ;\
	export OUTDIR="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/"; \
	bsub \
    -w "done(\"${FULLSM}_s09vldate\")" \
    -J "${FULLSM}_s10hcallr" \
    -o "${BASE}/WXS_Aquilla/04-LOGS/${RGBASE}_s10hcallr.%J" \
    -u "${EMAIL}" \
    -n1 -W 2880 \
    -M 46000000 \
    -R "rusage[mem=49152]" \
    -q research-hpc \
    -a "docker(vifehe/haplocallerv2)" \
    entrypoint.sh; \
	export MEM=2; \
	export GVCF="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${FULLSM}.aln.srt.isec-${RUN_TYPE}.markd.recal.raw.snps.indels.g.vcf.gz"; \
	export CCDS="/gscmnt/gc2645/wgs/Genome_Ref/GRCh38/GRCh38.CCDS.exons.sorted.bed"; \
	bsub \
    -w "done(\"${FULLSM}_s08dcvrge\") && done(\"${FULLSM}_s10hcallr\")" \
    -J "${FULLSM}_s11vnteval" \
    -o "${BASE}/WXS_Aquilla/04-LOGS/${FULLSM}_s11vnteval.%J" \
    -u "${EMAIL}" \
    -n1 -W 1360\
    -q research-hpc \
    -a "docker(vifehe/variantevalref)" \
    entrypoint.sh; \
	echo "Finished submitting jobs for ${FULLSM}"; \
	sleep ${DELAY2:-10}s; \
	echo "Finished submitting per-bam jobs for ${FULLSM}"; \
done; \
echo "Finished all samples within ${PR}"

### CHECK UP if ALL gvcfs have been created:
for x in $(cat ${WORKLIST}); do ls ../02-TRANSIT/201812_MGI_DIANWGS_REDCLOUD/${x}/*.vcf.gz; done

### CLEANING UP FASTQ files for DIAN_REDCLOUD
# 1 - generate a list of gvcf files finished, just keep the FULLSMID info
ls ../02-TRANSIT/201812_MGI_DIANWGS_REDCLOUD/*/*.g.vcf.gz | cut -d "/" -f 4 | sort | uniq > 201812_MGI_DIANWGS_REDCLOUD-fqtodelete.list
ls ../02-TRANSIT/201812_MGI_DIANWGS_REDCLOUD/*/*.g.vcf.gz | cut -d "/" -f 4 | sort | uniq > 201812_MGI_DIANWGS_REDCLOUD-fqtodelete.list2

cat 201812_MGI_DIANWGS_REDCLOUD-fqtodelete.list2 201812_MGI_DIANWGS_REDCLOUD-fqtodelete.list | sort | uniq -c | awk '{if ($1==1) print $2}' > 201812_MGI_DIANWGS_REDCLOUD-fqtodelete.LIST

## delete all fastq & cram within that FULLSMID directory for the samples that gvcf is completed
for x in $(cat 201812_MGI_DIANWGS_REDCLOUD-fqtodelete.list); do rm 201812_MGI_DIANWGS_REDCLOUD/${x}/*.fastq; done
for x in $(cat 201812_MGI_DIANWGS_REDCLOUD-fqtodelete.list); do rm 201812_MGI_DIANWGS_REDCLOUD/${x}/*.final.cram; done
for x in $(cat 201812_MGI_DIANWGS_REDCLOUD-fqtodelete.list); do rm ../02-TRANSIT/201812_MGI_DIANWGS_REDCLOUD/${x}/*.bam; done


### Check poresentce of files
### Test if we got Bam files for all samples
ls 01-RAW/201812_MGI_DIANWGS_REDCLOUD/*/*.final.cram | wc -l
34
ls 01-RAW/201812_MGI_DIANWGS_REDCLOUD/*/*.genome.bam | wc -l
1 --> 01-RAW/201812_MGI_DIANWGS_REDCLOUD/955-0004-0287652^0010007593^201812_MGI_DIANWGS_REDCLOUD/abf5010305ac466eb0d6c4d9c6e1e045.genome.bam
ls 01-RAW/201812_MGI_DIANWGS_REDCLOUD/*/*.r1.fastq | wc -l
538
ls 01-RAW/201812_MGI_DIANWGS_REDCLOUD/*/*.r2.fastq | wc -l
538
ls 02-TRANSIT/201812_MGI_DIANWGS_REDCLOUD/*/*.aln.srt.bam | wc -l
265
ls 02-TRANSIT/201812_MGI_DIANWGS_REDCLOUD/*/*.aln.srt.isec-paddedexome.markd.bam | wc -l
32
ls 02-TRANSIT/201812_MGI_DIANWGS_REDCLOUD/*/*.aln.srt.isec-*.markd.bam | wc -l
32
ls 02-TRANSIT/201812_MGI_DIANWGS_REDCLOUD/*/*.aln.srt.isec-*.markd.recal.table1 | wc -l
32
ls 02-TRANSIT/201812_MGI_DIANWGS_REDCLOUD/*/*.aln.srt.isec-*.markd.exome-coverage.sample_statistics | wc -l
31
ls 02-TRANSIT/201812_MGI_DIANWGS_REDCLOUD/*/*.aln.srt.isec-*.markd.exome-coverage.sample_summary | wc -l
31
ls 02-TRANSIT/201812_MGI_DIANWGS_REDCLOUD/*/*.aln.srt.isec-*.markd.recal.raw.snps.indels.g.vcf.gz | wc -l
82

### CLEAR UP SOME FILES:
fernandezv@virtual-workstation1:/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW$ rm 201812_MGI_DIANWGS_REDCLOUD/*/*.fastq

