# To run bgzip to compress VCF: 

for x in {1..22} X Y; do 
echo ${x}
export DIR="/gscmnt/gc2645/wgs/WXS_Aquilla/03-FINAL/VCFs/DIAN_VCF/SNP_VCF"
export FILE="All_DIAN_recal_joint_call_chr${x}.vcf"
export cmd="bgzip ${DIR}/${FILE} && tabix -p vcf ${DIR}/${FILE}.gz"
bsub \
-M 10000000 \
-R "rusage[mem=11152]" \
-n 4 \
-o "/gscmnt/gc2645/wgs/WXS_Aquilla/03-FINAL/VCFs/DIAN_VCF/SNP_VCF/logs/bgzip_All_DIAN_recal_joint_call_chr${x}.%J" \
-q research-hpc \
-a 'docker(achalneupane/htslib)' \
/bin/bash -c "$cmd" 
done




# export DIR="/gscmnt/gc2645/wgs/WXS_Aquilla/02-TRANSIT/DIAN/0141321^8011738217^DIAN"
# export FILE="141321^8011738217^DIAN.raw.snps.indels.g.vcf.gz"
# export cmd="tabix -p vcf ${DIR}/${FILE}"
# bsub \
# -M 10000000 \
# -R "rusage[mem=11152]" \
# -n 4 \
# -o "$DIR/create_tbi_achal.%J" \
# -q research-hpc \
# -a 'docker(achalneupane/htslib)' \
# /bin/bash -c "$cmd" 


