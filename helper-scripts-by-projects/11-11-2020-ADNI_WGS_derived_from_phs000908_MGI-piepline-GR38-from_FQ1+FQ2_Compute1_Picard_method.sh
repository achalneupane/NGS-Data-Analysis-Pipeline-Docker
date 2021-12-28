### 11-16-2020 - README file for ADNI_WGS on compute1 

##############
## 11-16-2020 
############## 
unset all ENV variables
for c in $(env | cut -d '=' -f 1); do unset $c; done

#######################################
#######################################
#######################################
#######################################
#######################################
#######################################

USER="achal"
PR="ADNI_WGS"
BASE="/storage1/fs1/cruchagac/Active/achal"
DEST="${BASE}/WXS_Aquilla/01-RAW"
# mkdir -p ${PR}
EMAIL="achal@wustl.edu"

# # on Fenix
# cd /80/AD/AD_Seq_Data/02.-BWA_GATK/02-gVCFs/GRCh37/singlerun
# ls -ld /80/AD/AD_Seq_Data/02.-BWA_GATK/02-gVCFs/GRCh37/singlerun/ADNI_*^ADNI_WGS > /home/achal/01.-RawData/ADNI_WGS/ADNI_SM_DNA_PR_BAM_PATH.txt
# cat /home/achal/01.-RawData/ADNI_WGS/ADNI_SM_DNA_PR_BAM_PATH.txt | awk '{print $9}' > /home/achal/01.-RawData/ADNI_WGS/ADNI_SM_DNA_PR_BAM_PATH.csv
# rm /home/achal/01.-RawData/ADNI_WGS/ADNI_SM_DNA_PR_BAM_PATH.txt



## from 809 work list I will only extract those that are already downoaded. 
ADNI_WGS-worklist_fullsm_809.csv

# this is the list of files already downloaded from BYU
# /storage1/fs1/cruchagac/Active/achal/WXS_Aquilla/01-RAW/downloaded_sample_list_11_16_2020_04_50_pm.txt 
grep -w -F -f /storage1/fs1/cruchagac/Active/achal/WXS_Aquilla/01-RAW/downloaded_sample_list_11_16_2020_04_50_pm.txt /storage1/fs1/cruchagac/Active/achal/WXS_Aquilla/01-RAW/ADNI_WGS-worklist_fullsm_809.csv > "${DEST}/${PR}-worklist_fullsm_downloaded.csv"


#################################################################
##### Run pipeline from bam back to fastq and then call variants
#################################################################
bjobs | grep RUN | cut -d\  -f1 | xargs bkill 
bjobs | grep PEN | cut -d\  -f1 | xargs bkill 
# bjobs | grep SSUSP | cut -d\  -f1 | xargs bkill 
#############
ls /storage1/fs1/cruchagac/Active/achal/WXS_Aquilla/01-RAW/ADNI_WGS/all_bam_files | grep .bam$ | cut -d. -f1 > /storage1/fs1/cruchagac/Active/achal/WXS_Aquilla/01-RAW/run_ready.txt
grep -w -F -f /storage1/fs1/cruchagac/Active/achal/WXS_Aquilla/01-RAW/run_ready.txt /storage1/fs1/cruchagac/Active/achal/WXS_Aquilla/01-RAW/ADNI_WGS-worklist_fullsm_809.csv > "${DEST}/${PR}-worklist_fullsm_immediately_available.csv"


USER="achal"
PR="ADNI_WGS"
BASE="/storage1/fs1/cruchagac/Active/achal"
DEST="${BASE}/WXS_Aquilla/01-RAW"
# mkdir -p ${PR}
EMAIL="achal@wustl.edu"

# 1st batch
# WORKLIST="${DEST}/${PR}-worklist_fullsm_immediately_available.csv"
# #2nd batch
# grep -wv -F -f "${DEST}/${PR}-worklist_fullsm_immediately_available.csv" "${DEST}/${PR}-worklist_fullsm_downloaded.csv" > "${DEST}/${PR}-worklist_fullsm_immediately_available_batch2.csv"
# WORKLIST="${DEST}/${PR}-worklist_fullsm_immediately_available_batch2.csv"
## 3rd batch
# WORKLIST="${DEST}/${PR}-worklist_fullsm_immediately_available_batch3.csv"
# 4th batch
WORKLIST="/storage1/fs1/cruchagac/Active/achal/WXS_Aquilla/01-RAW/ADNI_WGS-worklist_fullsm_immediately_available_batch4.csv"
NUM="$(cat $WORKLIST| wc -l)"



## After batch 2
# On MGI get the list of VCF files:
ls ./*/*vcf.gz| cut -d/ -f2
# paste it to :
vi /storage1/fs1/cruchagac/Active/achal/WXS_Aquilla/01-RAW/completed_list_after_batch_2.txt
# Find the samples that needs to be done
grep -wv -F -f /storage1/fs1/cruchagac/Active/achal/WXS_Aquilla/01-RAW/completed_list_after_batch_2.txt "${DEST}/ADNI_WGS-worklist_fullsm_809.csv" > "${DEST}/${PR}-worklist_fullsm_immediately_available_batch3.csv"


#### create LOOKUP and FULLSM using create metafile script for ADNI
## modify number of jobs that can be submitted
# https://confluence.ris.wustl.edu/pages/viewpage.action?pageId=27592450
# bgmod -L 200 /achal

# du --block-size=GB -s /storage1/fs1/cruchagac/Active/achal/WXS_Aquilla/01-RAW/ADNI_WGS
# scp -r ${FILE/Folder} achal@compute1-client-1.ris.wustl.edu:/storage1/fs1/cruchagac/Active/achal/WXS_Aquilla/01-RAW/all_bam_files

126, 127, 128, 129

# mkdir empty_dir
# for d in "$(ls | grep WGS$)"; do
# echo "Deleting dir $d"
# rsync -a --delete empty_dir/ $d/
# done

# From "${DEST}/${PR}-worklist_fullsm_immediately_available.csv"
# 1..25
# 26..50
# 51..75
# 76..101
### from "${DEST}/${PR}-worklist_fullsm_immediately_available_batch2.csv"
# 1..25
# 51..75
# 76..100
# 126..150
# 151..175
# 176..200
# 201..225
# 226..250
# 251..275
# 276..300
# 301..325
# 326..350
# 351..380
# 3rd batch
# 1..50
# 51..100
# 101..150
# 151..200
# 201..250
# 251..300
301..350



### Find samples that are yet to be processed 
# Samples in File 1 but not in File2
# comm -23 <(sort file1) <(sort file2)
# First on MGI
ls ./*/*g.vcf.gz| cut -d/ -f2 > ../../01-RAW/list.txt
# Then on Compute 1
scp fernandezv@virtual-workstation1.gsc.wustl.edu:/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/list.txt .
comm -23 <(sort ADNI_WGS-worklist_fullsm_809.txt) <(sort list.txt) > list1.txt
# Again find samples already processed in current batch
scp fernandezv@virtual-workstation1.gsc.wustl.edu:/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/ADNI_WGS-worklist_fullsm.csv .
comm -23 <(sort list1.txt) <(sort ADNI_WGS-worklist_fullsm.csv) > ADNI_WGS-worklist_fullsm_immediately_available_batch4.csv


START=26
END=27
USER="achal"
BASE="/storage1/fs1/cruchagac/Active/achal"
DEST="${BASE}/WXS_Aquilla/01-RAW"
EMAIL="achal@wustl.edu"
Snumber=1
export PR="ADNI_WGS"
# WORKLIST="${DEST}/${PR}-worklist_fullsm_immediately_available.csv"; 
# WORKLIST="${DEST}/${PR}-worklist_fullsm_immediately_available_batch2.csv"
# WORKLIST="${DEST}/${PR}-worklist_fullsm_immediately_available_batch3.csv"
WORKLIST="${DEST}/${PR}-worklist_fullsm_immediately_available_batch4.csv"
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
## export LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active/achal/WXS_Aquilla/01-RAW/all_bam_files:/storage1/fs1/cruchagac/Active/achal/WXS_Aquilla/01-RAW/all_bam_files"; \
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
# export BAMFILE="/storage1/fs1/cruchagac/Active/achal/WXS_Aquilla/01-RAW/all_bam_files_round2/${BAMBASE}"; \
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







#########################################
# # transfer required files
# scp -r fernandezv@virtual-workstation1.gsc.wustl.edu:/gscmnt/gc2645/wgs/Genome_Ref/GRCh38 .


# Still running



cd /storage1/fs1/cruchagac/Active/achal/WXS_Aquilla/01-RAW
USER="achal"
PR="ADNI_WGS"
BASE="/storage1/fs1/cruchagac/Active/achal"
DEST="${BASE}/WXS_Aquilla/01-RAW"
# mkdir -p ${PR}
EMAIL="achal@wustl.edu"

# WORKLIST="${DEST}/${PR}-worklist_fullsm_immediately_available.csv"
# WORKLIST="${DEST}/${PR}-worklist_fullsm_immediately_available_batch2.csv"
# WORKLIST="${DEST}/${PR}-worklist_fullsm_immediately_available_batch3.csv"
WORKLIST="${DEST}/${PR}-worklist_fullsm_immediately_available_batch4.csv"

# NUM="$(cat "${DEST}/${PR}-worklist_fullsm_need_globus_download.csv"| wc -l)"

#### create LOOKUP and FULLSM using create metafile script for ADNI
## modify number of jobs that can be submitted
# https://confluence.ris.wustl.edu/pages/viewpage.action?pageId=27592450
# bgmod -L 200 /achal

##########################################
# path to save metafile
metadir="${BASE}/WXS_Aquilla/01-RAW/ADNI_WGS"
PR="ADNI_WGS"
cd $metadir

# Number of directories
ls -l ${metadir} | grep "^d" | grep ${PR} | wc -l

# FROM batch # WORKLIST="${DEST}/${PR}-worklist_fullsm_immediately_available.csv"
# 1..25
# 26..50
# 51..75
# 76..101

# FROM batch # WORKLIST="${DEST}/${PR}-worklist_fullsm_immediately_available_batch2.csv"
# 1..25
# 26..50
# 51..75
# 76..100
# 101..125
# 126..150
# 151..175
# 176..200
# 201..225
# 226..250
# 251..275
# 276..300
# 301..325
# 326..350
351..380
# ADNI_4217^LP6005120-DNA_A03^ADNI_WGS
# bjobs 798706

# FROM batch # WORKLIST="${DEST}/${PR}-worklist_fullsm_immediately_available_batch3.csv"
# 1..50
# 51..100
# 101..150
# 151..200
# 201..250
# 251..300
301..350

# FROM batch # WORKLIST="${DEST}/${PR}-worklist_fullsm_immediately_available_batch4.csv"
1..35

Snumber=1
START=1;
END=35; 

#@@@@@@@@ change the number 
#@@@@@@@@ change the number 
#@@@@@@@@ change the number 
# FILE="${metadir}/${PR}_sm_dna-pr_rgbase_FQ1ext_FQ2ext_fullsm_1_to_25_samples.csv" 
# FILE="${metadir}/${PR}_sm_dna-pr_rgbase_FQ1ext_FQ2ext_fullsm_26_to_50_samples.csv" 
# FILE="${metadir}/${PR}_sm_dna-pr_rgbase_FQ1ext_FQ2ext_fullsm_76_to_101_samples.csv" 
# FILE="${metadir}/${PR}_sm_dna-pr_rgbase_FQ1ext_FQ2ext_fullsm_1_to_25_samples_batch2.csv" 
# FILE="${metadir}/${PR}_sm_dna-pr_rgbase_FQ1ext_FQ2ext_fullsm_26_to_50_samples_batch2.csv" 
# FILE="${metadir}/${PR}_sm_dna-pr_rgbase_FQ1ext_FQ2ext_fullsm_51_to_75_samples_batch2.csv" 
# FILE="${metadir}/${PR}_sm_dna-pr_rgbase_FQ1ext_FQ2ext_fullsm_76_to_100_samples_batch2.csv" 
# FILE="${metadir}/${PR}_sm_dna-pr_rgbase_FQ1ext_FQ2ext_fullsm_101_to_125_samples_batch2.csv"
# FILE="${metadir}/${PR}_sm_dna-pr_rgbase_FQ1ext_FQ2ext_fullsm_126_to_150_samples_batch2.csv" 
# FILE="${metadir}/${PR}_sm_dna-pr_rgbase_FQ1ext_FQ2ext_fullsm_151_to_175_samples_batch2.csv" 
# FILE="${metadir}/${PR}_sm_dna-pr_rgbase_FQ1ext_FQ2ext_fullsm_176_to_200_samples_batch2.csv" 
# FILE="${metadir}/${PR}_sm_dna-pr_rgbase_FQ1ext_FQ2ext_fullsm_226_to_250_samples_batch2.csv" 
# FILE="${metadir}/${PR}_sm_dna-pr_rgbase_FQ1ext_FQ2ext_fullsm_276_to_300_samples_batch2.csv" 
# FILE="${metadir}/${PR}_sm_dna-pr_rgbase_FQ1ext_FQ2ext_fullsm_301_to_325_samples_batch2.csv" 
# FILE="${metadir}/${PR}_sm_dna-pr_rgbase_FQ1ext_FQ2ext_fullsm_326_to_350_samples_batch2.csv" 
# FILE="${metadir}/${PR}_sm_dna-pr_rgbase_FQ1ext_FQ2ext_fullsm_351_to_380_samples_batch2.csv" 
# FILE="${metadir}/${PR}_sm_dna-pr_rgbase_FQ1ext_FQ2ext_fullsm_1_to_50_samples_batch3.csv" 
# FILE="${metadir}/${PR}_sm_dna-pr_rgbase_FQ1ext_FQ2ext_fullsm_51_to_100_samples_batch3.csv" 
# FILE="${metadir}/${PR}_sm_dna-pr_rgbase_FQ1ext_FQ2ext_fullsm_151_to_200_samples_batch3.csv" 
# FILE="${metadir}/${PR}_sm_dna-pr_rgbase_FQ1ext_FQ2ext_fullsm_251_to_300_samples_batch3.csv" 
# FILE="${metadir}/${PR}_sm_dna-pr_rgbase_FQ1ext_FQ2ext_fullsm_301_to_350_samples_batch3.csv" 
FILE="${metadir}/${PR}_sm_dna-pr_rgbase_FQ1ext_FQ2ext_fullsm_1_to_35_samples_batch4.csv" 


#@@@@@@@@
#@@@@@@@@
#@@@@@@@@
#@@@@@@@@
rm ${FILE}
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
# echo "${SM},${DNA},${PR},${RG},${R},${FULLSM}" 
done
done




# Failed files:
# LP6005122-DNA_F06 # in 33_50

# transfer file to MGI
scp  ${FILE} fernandezv@virtual-workstation1.gsc.wustl.edu:/gscmnt/gc2645/wgs/WXS_Aquilla/01-RAW/
#@@@@@@@@ change the number 
#@@@@@@@@ change the number 
#@@@@@@@@ change the number 
#@@@@@@@@ change the number 
# also remove failed file from the list
#@@@@@@@@ 
#@@@@@@@@
#@@@@@@@@ 
#@@@@@@@@





#### GLOBUS transfer
## globus
# ./globusconnectpersonal &
# https://docs.globus.org/how-to/globus-connect-personal-linux/

# #grep -w -F -f /home/achal/01.-RawData/ADNI_WGS/globus_transfer/transfer1.txt /home/achal/01.-RawData/ADNI_WGS/ADNI_WGS-worklist_fullsm.csv > worklist.txt
# #grep -w -F -f /home/achal/01.-RawData/ADNI_WGS/globus_transfer/transfer2.txt /home/achal/01.-RawData/ADNI_WGS/ADNI_WGS-worklist_fullsm.csv > /home/achal/01.-RawData/ADNI_WGS/globus_transfer/worklist.txt


# https://app.globus.org/file-manager
# WORKLIST="/home/achal/01.-RawData/ADNI_WGS/globus_transfer/worklist.txt"
# Use Wash U RIS storage1 dtn1
# WORKLIST="/storage1/fs1/cruchagac/Active/achal/WXS_Aquilla/01-RAW/ADNI_WGS-worklist_fullsm_immediately_available.csv"
# 1..101
# WORKLIST="/storage1/fs1/cruchagac/Active/achal/WXS_Aquilla/01-RAW/ADNI_WGS-worklist_fullsm_immediately_available_batch2.csv"

# Samples BATCH4
# /storage1/fs1/cruchagac/Active/achal/WXS_Aquilla/01-RAW/ADNI_WGS-worklist_fullsm_809.csv




# WORKLIST="/storage1/fs1/cruchagac/Active/achal/WXS_Aquilla/01-RAW/ADNI_WGS-worklist_fullsm_immediately_available_batch3.csv"
# 1..25
# 26..50
# 51..75
# 76..100
# 101..125
# 176..200
# 201..225
# 226..250
# 251..275
# 276..300
# 301..325
# 326..350
# 351..375
# 1..50
# 51..100
# 101..150
# 151..200
# 251..300
# 301..350

WORKLIST="/storage1/fs1/cruchagac/Active/achal/WXS_Aquilla/01-RAW/ADNI_WGS-worklist_fullsm_immediately_available_batch4.csv"

Wash U RIS storage1 dtn2
cd /storage1/fs1/cruchagac/Archive/ADNI_WGS

Snumber=1;
START=26; 
END=27;

s1_id="01f0ac4c-9570-11ea-b3c4-0ae144191ee3"
## On July 27, 2021, I changed it to 9888f31d-ef1f-11eb-ab63-d195c983855c


s1_path="storage1/fs1/cruchagac"
for FULLSM in $(sed -n "${START},${END}p" "${WORKLIST}"); do 
echo "Doing Number: ${Snumber}"  
((Snumber=${Snumber}+1))
DNA="$(echo "${FULLSM}" | cut -d^ -f2).bam"
globus transfer "${s1_id}:${s1_path}/Archive/ADNI_WGS/${DNA}" "${s1_id}:${s1_path}/Active/achal/WXS_Aquilla/01-RAW/all_bam_files/${DNA}"
# globus transfer "${s1_id}:${s1_path}/Archive/ADNI_WGS/${DNA}" "${s1_id}:${s1_path}/Active/achal/WXS_Aquilla/01-RAW/all_bam_files_round2/${DNA}"
done

# # to cancel all task:
# https://docs.globus.org/cli/reference/task_cancel/
globus task cancel -a

## to view task
# https://docs.globus.org/cli/reference/task_show/
globus task show

## task monitoring
#https://docs.globus.org/api/transfer/task/
###########################@@@@@@@@@@ 
###########################@@@@@@@@@@ 
###########################@@@@@@@@@@ 
###########################@@@@@@@@@@ 

# check if all bams are there in the directory
bamdir="/storage1/fs1/cruchagac/Active/achal/WXS_Aquilla/01-RAW/ADNI_WGS/all_bam_files"
cd $bamdir
# cat $WORKLIST |cut -d^ -f2| xargs -I {} bash -c 'echo ${bamdir}*{}.bam' | xargs ls -lht
ls -lht $(cat $WORKLIST |cut -d^ -f2| xargs -I {} bash -c 'echo ${bamdir}*{}.bam')
# ls -1ht $(ls | head -17 | xargs -I {} bash -c 'echo {}')



###############################################################
###############################################################
# delete bam files
# # delete bam files
# # delete bam files
# # 1..25
# 26..50

# START=26; \
# END=50; \
# USER="achal"
# PR="ADNI_WGS"
# BASE="/storage1/fs1/cruchagac/Active/achal"
# DEST="${BASE}/WXS_Aquilla/01-RAW"
# EMAIL="achal@wustl.edu"
# Snumber=1
# FULLSM="ADNI_1380^LP6005117-DNA_G04^ADNI_WGS"
# WORKLIST="${DEST}/${PR}-worklist_fullsm_immediately_available.csv"; \
# for FULLSM in $(sed -n "${START},${END}p" "${WORKLIST}"); do \
# PR="ADNI_WGS"; \
# export DEST="${BASE}/WXS_Aquilla/01-RAW"; \
# mkdir -p ${DEST}/${PR}/${FULLSM}; \
# export SM="$(echo "${FULLSM}" | cut -d^ -f1)"; \
# export FULLSM=${FULLSM}; \
# unset DNA; \
# export DNA="$(echo "${FULLSM}" | cut -d^ -f2)"; \
# export PR="$(echo "${FULLSM}" | cut -d^ -f3)"; \
# export OUT_DIR="${BASE}/WXS_Aquilla/01-RAW/${PR}/${FULLSM}"; \
# export BAMBASE=${DNA}.bam; \
# # export BAMBASE="ADNI_1380^LP6005117-DNA_G04^ADNI_WGS_GATKready.bam"; \
# export MEM=32; \
# # rsync -avh /storage1/fs1/cruchagac/Archive/ADNI_WGS/${BAMBASE} ${DEST}/all_bam_files/;
# unset BAMFILE; \
# BAMFILE="/storage1/fs1/cruchagac/Active/achal/WXS_Aquilla/01-RAW/all_bam_files/${BAMBASE}"; \
# rm -f $BAMFILE; \
# done

# -R 'select[mem>69152 && tmp>250] rusage[mem=69152, tmp=250]' \
# -R "rusage[mem=69152]" \
