### 2019-09 - README file for 2019_MGI_gDNA_LINDSEY sequences processing over at MGI

##############
## 2019-09-19
##############

### In Fenix

# generate metadata file for the PROJECT in question
HOMEDIR="/40/AD/AD_Seq_Data/01.-RawData"
PR="201909_MGI_gDNA_LINDSEY"

##List directories within the PR directory
## Macrogen like DIAN have several bams that need to be merged
metafile="${PR}-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv"
vim ${metafile}   ## copy here the metadata for the DIAN project from the excel file
sed -i 's/\t/,/g' ${metafile}


## this file contains the follwoing fields, comma separated
# SM_(MGI),SM,DNA,PR_(project),RGBASE_(MGI sample ID),cram,SMDIR,CRAM_location

## Upload this file to MGI
scp ${metafile} fernandezv@virtual-workstation1.gsc.wustl.edu:/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/

## Log onto MGI
#####################################
#### 1 - upload samples using rsync
#####################################

BASE="/gscmnt/gc2645/wgs/"

# ################ Skip this if running first batch
# ## I have uplaoded the list with over 547 samples; AN: testing for to samples only
# wc -l 201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv
# # 547 201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv
########################################

#mv 201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv 201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc_full.csv


# extract range of lines
# Doing 6-26

####### here restricts the search upto that line
# 1,20 done
# 21,47 done
# 48,67 done
# 68, 110 done
# 111, 160 done
# 161, 210 done
# 211,260 done
# 261, 308 done
# 309, 355 done; 42nd not done TWGX-MAP_64786-8036056104
# 356, 404 done; 2nd not done, doing in next run
# 405, 450 done
# 451, 500 done; redo 17, 41, 42 in next run 
# 501, 546; done
# 547, 595; done; redo 23
# 596, 629; redo 34,35  MAP_86519^8036056082^201909_MGI_gDNA_LINDSEY and MAP_86619^8036056056^201909_MGI_gDNA_LINDSEY 

# # 629 Lindsey
# sed -n '596,629p' 201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc_full.csv > 201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv

# 8 samples from second batch of Lindsey
cat 201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc_8_FULL.csv  > 201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv




# printf "%s\n" 34p 35p 42p | ed -s  201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv > file1.csv
printf "%s\n" 1p 2p 5p 7p | ed -s  201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv > file1.csv
# printf "%s\n" 33p 34p | ed -s  201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv > file1.csv

# # sed -n -e 23p 201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv > file1.csv


# # cat file1.csv >> 201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv
# cat file1.csv > 201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv


# sed -n -e '17,48p;16483q' 201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv >> file1.csv

# ######## skip this for the first batch#############
# # In a previous test before I ran 34 samples from this list,  
# wc -l 201909_MGI_gDNA_LINDSEY-worklist_fullsm.csv
# # 34 201909_MGI_gDNA_LINDSEY-worklist_fullsm.csv

# ###################################################

# # FULL worklist
# cut -d "," -f 7 201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc_full.csv > 201909_MGI_gDNA_LINDSEY-worklist_fullsm.csv
# FULL worklist
cut -d "," -f 7 201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv > 201909_MGI_gDNA_LINDSEY-worklist_fullsm.csv

### I need to make a list of the smaples that I did not process then to process them now.
cut -d "," -f 7 201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv > LIST1
cat 201909_MGI_gDNA_LINDSEY-worklist_fullsm.csv > 201909_MGI_gDNA_LINDSEY-worklist_fullsm.csv2

# removes the items in the LIST1
#cat LIST1 201909_MGI_gDNA_LINDSEY-worklist_fullsm.csv | sort | uniq -c | awk '{if ($1==1) print $2}' > 201909_MGI_gDNA_LINDSEY-worklist_fullsm.csv2

## My new working list will be 
worklist="201909_MGI_gDNA_LINDSEY-worklist_fullsm.csv2"
## This list contains 513 samples
PR="201909_MGI_gDNA_LINDSEY"
USER="achal"

metafile="${PR}-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv"

# a - test that we are reading things correctly
Snumber=1
while read -r line; do \
CRAMLOC="$(echo $line | cut -d "," -f 8)"
MGI_ID="$(echo $line | cut -d "," -f 1)"
SMDIR="$(echo $line | cut -d "," -f 7)"
RGBASE="$(echo $line | cut -d "," -f 5)"
echo -e "${CRAMLOC}\n${MGI_ID}\n${SMDIR}\n${RGBASE}"
echo "*****************************************************"; \
echo "*****************************************************"; \
echo "Getting sample number**********: " $Snumber; \
echo "*****************************************************"; \
echo "*****************************************************"; \
((Snumber=Snumber+1)); \
done < ${metafile} 
"/40/AD/AD_Seq_Data/01.-RawData/201909_MGI_gDNA_LINDSEY/01.-RawData/xfer.genome.wustl.edu/gxfer1/54456299110049/TWGX-MAP_10064-8036056040.cram
TWGX-MAP_10064-8036056040
MAP_10064^8036056040^201909_MGI_gDNA_LINDSEY
MAP_10064^8036056040^201909_MGI_gDNA_LINDSEY.TWGX-MAP_10064-8036056040
/40/AD/AD_Seq_Data/01.-RawData/201909_MGI_gDNA_LINDSEY/01.-RawData/xfer.genome.wustl.edu/gxfer1/54456299110049/TWGX-MAP_11570-8036056502.cram
TWGX-MAP_11570-8036056502
MAP_11570^8036056502^201909_MGI_gDNA_LINDSEY
MAP_11570^8036056502^201909_MGI_gDNA_LINDSEY.TWGX-MAP_11570-8036056502
"

# I had to add SSH key to transfer files
# b - upload files: 
# $HOME/.ssh/config
 # $HOME/.ssh/mykey

Snumber=1;
for FULLSM in $(cat 201909_MGI_gDNA_LINDSEY-worklist_fullsm.csv2); do \
echo "*****************************************************"; \
echo "*****************************************************"; \
echo "Getting sample number**********: " $Snumber; \
echo "*****************************************************"; \
echo "*****************************************************"; \
((Snumber=Snumber+1)); \
echo ${FULLSM}; \
line="$(grep ${FULLSM} ${metafile})";\
echo ${line};\
CRAMLOC="$(echo $line | cut -d "," -f 8)";\
MGI_ID="$(echo $line | cut -d "," -f 1)";\
SMDIR="$(echo $line | cut -d "," -f 7)";\
RGBASE="$(echo $line | cut -d "," -f 5)"; \
echo -e "CRAMLOC:${CRAMLOC}\nMGI_ID:${MGI_ID}\nSMDIR:${SMDIR}\nRGBASE:${RGBASE}"; \
DEST="/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/${PR}/${SMDIR}"; \
mkdir -p ${DEST}
# echo "Destination: " ${DEST}
# making it specific to add cram and md5 files; no vcf included
# # rsync -e "ssh -i $HOME/.ssh/mykey" -avh ${USER}@fenix.psych.wucon.wustl.edu:${CRAMLOC%.*}.cram* ${DEST}/;\
rsync -avh ${USER}@fenix.psych.wucon.wustl.edu:${CRAMLOC%.*}.cram* ${DEST}/;\
### rsync -avh ${USER}@fenix.psych.wucon.wustl.edu:${CRAMLOC%.*}.* ${DEST}/;\
done

### Check that all files have uploaded
fernandezv@virtual-workstation1:/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW$ ls 201909_MGI_gDNA_LINDSEY/*/*.cram | wc -l
546

#####################################
### 2 - Run vifehe/cramtobam for all samples using metafile data for input
#####################################

# the original metafile is: 201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv
## I will be runnign them on batches of 50 so when I generate the fastq I do not overload the system.
201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv-1-5 --> DONE
201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv-6-50  --> DONE
201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv-51-100 --> DONE
201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv-101-150 --> DONE
201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv-151-200 --> DONE
201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv-201-250 --> DONE 
201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv-251-300 --> DOING 
201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv-301-350  
201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv-351-400
201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv-401-450
201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv-451-500
201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv-501-547


# # # kill jobs by pattern:
# bjobs | grep RUN | cut -d\  -f1 | xargs bkill 
# bjobs | grep PEN | cut -d\  -f1 | xargs bkill 
## to move log files matching LINDSEY:
# mv `ls | grep LINDSEY` JULY9/      

EMAIL="achal@wustl.edu"
# metafile="201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv-251-300";\

Snumber=1
metafile="201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv";\
while read -r line; do \
	echo "*****************************************************"; \
	echo "*****************************************************"; \
	echo "Doing sample number**********: " $Snumber; \
	echo "*****************************************************"; \
	echo "*****************************************************"; \
	((Snumber=Snumber+1)); \
	SMDIR="$(echo $line | cut -d "," -f 7)";
	PR="201909_MGI_gDNA_LINDSEY"
	MGI_ID="$(echo $line | cut -d "," -f 1)"
	RUN_TYPE="genome"
	SM="$(echo $line | cut -d "," -f 5)"
	CRAMFILE="${BASE}/WXS_Aquilla/01-RAW/${PR}/${SMDIR}/${MGI_ID}.cram"
	OUTDIR="${BASE}/WXS_Aquilla/01-RAW/${PR}/${SMDIR}/"
	REF="/gscmnt/production_reference/info/production_reference_GRCh38DH/reference/all_sequences.fa"
	echo -e "${MGI_ID}\n${SMDIR}"
	export BASE="/gscmnt/gc2645/wgs"
	export WORKDIR="${BASE}/tmp";
	export RAWDIR="${BASE}/WXS_Aquilla/01-RAW/${PR}/${SMDIR}"
	export CRAMFILE="${BASE}/WXS_Aquilla/01-RAW/${PR}/${SMDIR}/${MGI_ID}.cram"; \
    export BAMFILE="${BASE}/WXS_Aquilla/01-RAW/${PR}/${SMDIR}/${MGI_ID}.bam"
	export OUTDIR="${BASE}/WXS_Aquilla/01-RAW/${PR}/${SMDIR}"; \
	export RUN_TYPE="genome"; \
	export REF="/gscmnt/production_reference/info/production_reference_GRCh38DH/reference/all_sequences.fa"
    echo -e "${CRAMFILE}\n"; \
	bsub \
      -J "${MGI_ID}_s00cram2bam" \
      -o "${BASE}/WXS_Aquilla/04-LOGS/${MGI_ID}_s00cram2bam.%J" \
      -u "${EMAIL}" \
      -n1 -W 360 \
      -q research-hpc \
      -a "docker(vifehe/cram2bam)" \
      entrypoint.sh; 
done < ${metafile}

#After this step, there will be:
# fernandezv@virtual-workstation1:/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/201909_MGI_gDNA_LINDSEY/MAP_10064^8036056040^201909_MGI_gDNA_LINDSEY$ ls
# TWGX-MAP_10064-8036056040.cram  TWGX-MAP_10064-8036056040.cram.md5  TWGX-MAP_10064-8036056040.genome.bam


### Test if we got Bam files for all samples
ls 201909_MGI_gDNA_LINDSEY/*/*.cram | wc -l
546
ls 201909_MGI_gDNA_LINDSEY/*/*.bam | wc -l
46
# This is only after runing bwa- the big chunk below
ls ../02-TRANSIT/201909_MGI_gDNA_LINDSEY/*/*.bam | wc -l
331


#####################################
### 3 - Start pipeline from fastq to variant calling
#####################################

## Need to generate the MGIID-WASHUID lookup file plus the working file for each of the previous 50x sets;
for LOOKUP in $(ls 201909_MGI_gDNA_LINDSEY-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv*); do \
awk -F"," '{print $1","$2","$3","$4","$1".genome.bam"}' ${LOOKUP} > 201909_MGI_gDNA_LINDSEY-MGIID-SMID-DNA-PR-BAMFILE.${LOOKUP##*.}; done

## Start with "complete pipeline following John's script
for LOOKUP in $(ls 201909_MGI_gDNA_LINDSEY-MGIID-SMID-DNA-PR-BAMFILE.csv*); do cut -d, -f2-4 "${LOOKUP}" | tr "," "^" | sort -u > ${PR}-worklist_fullsm.${LOOKUP##*.}; done

FULLSM="${PR}-worklist_fullsm.csv"
# FBASE may vary depending on the project, due to inconsistency in file structure across projects (e.g. FULLSM/ bams/ bam/ bam/FULLSM)
#   Ideally would query BAMFILE from master LOOKUP with all files across all projects
Snumber=1
while read -r line; do \
echo "*****************************************************"; \
echo "*****************************************************"; \
echo "Doing sample number**********: " $Snumber; \
echo "*****************************************************"; \
echo "*****************************************************"; \
((Snumber=Snumber+1)); \
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
010-0004-0085842^0010007564^201909_MGI_gDNA_LINDSEY
3e2b8f9520e144c29debedf0e74d1a8f
010-0004-0085842
0010007564
201909_MGI_gDNA_LINDSEY
3e2b8f9520e144c29debedf0e74d1a8f.genome.bam
/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/201909_MGI_gDNA_LINDSEY/010-0004-0085842^0010007564^201909_MGI_gDNA_LINDSEY
/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/201909_MGI_gDNA_LINDSEY/010-0004-0085842^0010007564^201909_MGI_gDNA_LINDSEY/3e2b8f9520e144c29debedf0e74d1a8f.genome.bam


#################################################################
##### Run pipeline from bam back to fastq and then call variants
#################################################################
PR="201909_MGI_gDNA_LINDSEY"




# If needed incorporate this at line 159, after declaring ${FULLSM}
#if [ -d ${BASE}/WXS_Aquilla/03-FINAL/${PR}/${FULLSM} ] || [ -d ${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM} ]; then echo "Skipping ${FULLSM}, already run(ning)"; continue; else echo "RUN ${FULLSM}"; fi; \

Snumber=1
START=1; \
END=4; \
DELAY=10
EMAIL="achal@wustl.edu"
SHELLDROP=0

export BASE="/gscmnt/gc2645/wgs"; \
export WORKDIR="${BASE}/tmp"; \
export THREADS=16; \
export BWA_PIPE_SORT=1; \
export TIMING=1;\
PR="201909_MGI_gDNA_LINDSEY"; \
WORKLIST="${BASE}/WXS_Aquilla/01-RAW/${PR}-worklist_fullsm.csv"; \
LOOKUP="${BASE}/WXS_Aquilla/01-RAW/${PR}-MGIID-SMID-DNA-PR-BAMFILE.csv"; \
LOOKUP_COL_SM=2; \
LOOKUP_COL_DNA=3; \
LOOKUP_COL_PR=4; \
LOOKUP_COL_BAMFILE=5; \
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
	export OUT_DIR="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}"; \
	INPUT_LIST=(); \
	WAIT_LIST=(); \
	for BAMBASE in $(grep "${SM},${DNA},${PR}" "${LOOKUP}" | cut -d, -f${LOOKUP_COL_BAMFILE}); do \
	echo -e "00a - In this round we will analyze ${FULLSM}\n${BAMFILE}\n${SM}\n${DNA}\n${PR}"; \
	export BAMFILE="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${BAMBASE}"; \
	export MEM=6; \
	export CREATE_RGFILE=1; \
	export DEBUG=1; \
	export REMOVE_INPUT=1; \
	export OUT_DIR="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/"; \
	bsub \
	-J "${FULLSM}_s00rvtbam" \
	-o "${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}/${FULLSM}.${BAMBASE}_s00rvtbam.%J" \
	-u "${EMAIL}" \
	-n1 -W 1440 \
	-R "rusage[mem=8192]" \
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
	export MEM=6; \
	export BAMFILE="${BASE}/WXS_Aquilla/02-TRANSIT/${PR}/${FULLSM}/${RGBASE}.aln.srt.bam"; \
	echo ${BAMFILE}; \
	export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
	bsub \
	-w "done(\"${RGBASE}_s01alnsrt\")" \
	-J "${RGBASE}_s02vldate" \
	-o "${BASE}/WXS_Aquilla/04-LOGS/${RGBASE}_s02vldate.%J" \
	-u "${EMAIL}" \
	-n1 -W 360 \
	-R "rusage[mem=8192]" \
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
     -n1 -W 360 \
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
	-n1 -W 360 \
	-R "rusage[mem=8192]" \
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
  unset MEM; \
  	export RUN_TYPE="paddedexome" ;\
	export REF="${BASE}/Genome_Ref/GRCh38/Homo_sapiens_assembly38.fasta"; \
	export BEDFILE="${BASE}/Genome_Ref/GRCh38/Capture_Padded.GRCh38.bed"; \
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
    -n1 -W 360\
    -q research-hpc \
    -a "docker(vifehe/variantevalref)" \
    entrypoint.sh; \
	echo "Finished submitting jobs for ${FULLSM}"; \
	sleep ${DELAY2:-10}s; \
	echo "Finished submitting per-bam jobs for ${FULLSM}"; \
done; \
echo "Finished all samples within ${PR}"

### CHECK UP if ALL gvcfs have been created:
WORKLIST="${BASE}/WXS_Aquilla/01-RAW/${PR}-worklist_fullsm.csv"
for x in $(cat ${WORKLIST}); do ls -lh ../02-TRANSIT/201909_MGI_gDNA_LINDSEY/${x}/*.vcf.gz; done

### CLEANING UP FASTQ files for DIAN_REDCLOUD
# 1 - generate a list of gvcf files finished, just keep the FULLSMID info
ls ../02-TRANSIT/201909_MGI_gDNA_LINDSEY/*/*.g.vcf.gz | cut -d "/" -f 4 | sort | uniq > 201909_MGI_gDNA_LINDSEY-fqtodelete.list
ls ../02-TRANSIT/201909_MGI_gDNA_LINDSEY/*/*.g.vcf.gz | cut -d "/" -f 4 | sort | uniq > 201909_MGI_gDNA_LINDSEY-fqtodelete.list2

cat 201909_MGI_gDNA_LINDSEY-fqtodelete.list2 201909_MGI_gDNA_LINDSEY-fqtodelete.list | sort | uniq -c | awk '{if ($1==1) print $2}' > 201909_MGI_gDNA_LINDSEY-fqtodelete.LIST

## delete all fastq & cram within that FULLSMID directory for the samples that gvcf is completed
for x in $(cat 201909_MGI_gDNA_LINDSEY-fqtodelete.list); do rm 201909_MGI_gDNA_LINDSEY/${x}/*.fastq; done
for x in $(cat 201909_MGI_gDNA_LINDSEY-fqtodelete.list); do rm 201909_MGI_gDNA_LINDSEY/${x}/*.cram; done
for x in $(cat 201909_MGI_gDNA_LINDSEY-fqtodelete.list); do rm ../02-TRANSIT/201909_MGI_gDNA_LINDSEY/${x}/*.bam; done


### Check poresentce of files
### Test if we got Bam files for all samples
ls ./201909_MGI_gDNA_LINDSEY/*/*.cram | wc -l
34
ls ./201909_MGI_gDNA_LINDSEY/*/*.genome.bam | wc -l
1 --> 01-RAW/201909_MGI_gDNA_LINDSEY/955-0004-0287652^0010007593^201909_MGI_gDNA_LINDSEY/abf5010305ac466eb0d6c4d9c6e1e045.genome.bam
ls ./201909_MGI_gDNA_LINDSEY/*/*.r1.fastq | wc -l
538
ls ./201909_MGI_gDNA_LINDSEY/*/*.r2.fastq | wc -l
538
ls ../02-TRANSIT/201909_MGI_gDNA_LINDSEY/*/*.aln.srt.bam | wc -l
265
ls ../02-TRANSIT/201909_MGI_gDNA_LINDSEY/*/*.aln.srt.isec-paddedexome.markd.bam | wc -l
32
ls ../02-TRANSIT/201909_MGI_gDNA_LINDSEY/*/*.aln.srt.isec-*.markd.bam | wc -l
32
ls ../02-TRANSIT/201909_MGI_gDNA_LINDSEY/*/*.aln.srt.isec-*.markd.recal.table1 | wc -l
32
ls ../02-TRANSIT/201909_MGI_gDNA_LINDSEY/*/*.aln.srt.isec-*.markd.exome-coverage.sample_statistics | wc -l
31
ls ../02-TRANSIT/201909_MGI_gDNA_LINDSEY/*/*.aln.srt.isec-*.markd.exome-coverage.sample_summary | wc -l
31
ls ../02-TRANSIT/201909_MGI_gDNA_LINDSEY/*/*.aln.srt.isec-*.markd.recal.raw.snps.indels.g.vcf.gz | wc -l
82

### CLEAR UP SOME FILES:
fernandezv@virtual-workstation1:/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW$ rm ./201909_MGI_gDNA_LINDSEY/*/*.fastq



