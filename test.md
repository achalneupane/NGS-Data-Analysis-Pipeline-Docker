# ADGC_2021_Analysis

December 30, 2021 

Here, I am performing PCA analysis with all four ethnicities.

WORKING_DIR =/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/01-EOAD-preQC/02-Analysis/ADGC-HapMap-PCA/ADGC

Create symlinks to the plink files:

```r
# NHW
ln -s /40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/ADGC_NHW_Cohort.fam ADGC_NHW_Cohort_without_chr.fam
ln -s /40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/ADGC_NHW_Cohort.bed ADGC_NHW_Cohort_without_chr.bed
sed 's/chr//g' /40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/ADGC_NHW_Cohort.bim >  ADGC_NHW_Cohort_without_chr.bim

# AA
ln -s /40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_AA/ADGC_AA_Cohort.fam ADGC_AA_Cohort_without_chr.fam
ln -s /40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_AA/ADGC_AA_Cohort.bed ADGC_AA_Cohort_without_chr.bed
sed 's/chr//g' /40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_AA/ADGC_AA_Cohort.bim > ADGC_AA_Cohort_without_chr.bim

# Hispanic
ln -s /40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_Hispanic/ADGC_Hispanic_Cohort.fam ADGC_Hispanic_Cohort_without_chr.fam
ln -s /40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_Hispanic/ADGC_Hispanic_Cohort.bed ADGC_Hispanic_Cohort_without_chr.bed
sed 's/chr//g' /40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_Hispanic/ADGC_Hispanic_Cohort.bim > ADGC_Hispanic_Cohort_without_chr.bim

# Asian
ln -s /40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_Asian/ADGC_Asian_Cohort.fam ADGC_Asian_Cohort_without_chr.fam 
ln -s /40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_Asian/ADGC_Asian_Cohort.bed ADGC_Asian_Cohort_without_chr.bed
sed 's/chr//g' /40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_Asian/ADGC_Asian_Cohort.bim > ADGC_Asian_Cohort_without_chr.bim
```

[HapMap PCA](ADGC_2021_Analysis%20cdfc76305a2746a1a8def50f783e624b/HapMap%20PCA%20aa8082bcb6534bd9bbf1614ea988d1fb.md)

[PCA 1000G](ADGC_2021_Analysis%20cdfc76305a2746a1a8def50f783e624b/PCA%201000G%20a964bd5de0a94d4fac05a278c81ea05d.md)

[Selection of samples by Age criteria](ADGC_2021_Analysis%20cdfc76305a2746a1a8def50f783e624b/Selection%20of%20samples%20by%20Age%20criteria%20a366d3e5e3dc4131bbdfe9f37a6d643c.md)

[2022-01-14 Single variant analyses - NHW](ADGC_2021_Analysis%20cdfc76305a2746a1a8def50f783e624b/2022-01-14%20Single%20variant%20analyses%20-%20NHW%20a98cd4b0c5904edf9d326292e75cd59b.md)

[2022-01-14 Single variant analyses - AA](ADGC_2021_Analysis%20cdfc76305a2746a1a8def50f783e624b/2022-01-14%20Single%20variant%20analyses%20-%20AA%206f9b11f484a14545b7534003b928a57a.md)

[2022-01-14 Single variant analyses - Asian](ADGC_2021_Analysis%20cdfc76305a2746a1a8def50f783e624b/2022-01-14%20Single%20variant%20analyses%20-%20Asian%20d10a93624f114f6087e9bcd789100acd.md)

[2022-01-14 Single variant analyses - Hispanics](ADGC_2021_Analysis%20cdfc76305a2746a1a8def50f783e624b/2022-01-14%20Single%20variant%20analyses%20-%20Hispanics%203abc9ef05d3d4da9b8e21f65e723dda9.md)
