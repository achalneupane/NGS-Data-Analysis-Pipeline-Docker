cd /40/AD/AD_Seq_Data/01.-RawData/MGI_Gregg_201704/ready
comm -23 <(sort /home/achal/01.-RawData/missing_check/MGI_Gregg_all_bam.txt| cut -d. -f1) <(sort /home/achal/01.-RawData/missing_check/MGI_Gregg_201704.txt)
## 68317^8012957245^MGI_Gregg_201704
## 72173^8011736646^MGI_Gregg_201704
# PR="MGI_Gregg_201704"
# OUTDIR="/home/achal/01.-RawData/${PR}"
# mkdir -p $OUTDIR
# FILE="$PR-SM-DNA-BAM-DIR.csv"
# rm "${OUTDIR}/$FILE"
# for x in $(ls| grep .bam); do 
# i="${x/.bam/}"
# SM="$(echo ${i} | cut -d^ -f1)"
# DNA="$(echo ${i} | cut -d^ -f2)"
# BAM="${x}"
# DIR="/40/AD/AD_Seq_Data/01.-RawData/MGI_Gregg_201704/ready"
# echo "${SM},${DNA},${BAM},${PR},${DIR}" >> "${OUTDIR}/$FILE"
# done
