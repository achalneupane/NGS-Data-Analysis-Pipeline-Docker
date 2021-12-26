chmod u+x entrypoint.sh

export BASE="/gscmnt/gc2645/wgs"
export FULLSM="MAP_27856^8038741468^202004_USUHS_EOAD-WGS_gDNA_EOLUS"
export BAMDIR="${BASE}/WXS_Aquilla/02-TRANSIT/202004_USUHS_EOAD-WGS_gDNA_EOLUS/${FULLSM}"
export BAMFILE="${BAMDIR}/${FULLSM}.aln.srt.isec-paddedexome.markd.bam"
export INDEX_BAMFILE="${BAMDIR}/${FULLSM}.aln.srt.isec-paddedexome.markd.bai"
export VCFFILE="${BAMDIR}/${FULLSM}.aln.srt.isec-paddedexome.markd.recal.raw.snps.indels.g.vcf.gz"
export WORKDIR="${BASE}/WXS_Aquilla/gvcfTest/verifybamid/${FULLSM}"
export MEM=16; \
bsub \
     -J "${FULLSM}_s13veriBAMID" \
     -o "${BASE}/WXS_Aquilla/gvcfTest/verifybamid/${FULLSM}_s13veriBAMID.%J"\
     -u "${EMAIL}" \
     -n1 -W 1360 \
     -q research-hpc \
     -a "docker(achalneupane/verifybamid:latest)" \
     ${BASE}/WXS_Aquilla/gvcfTest/verifybamid/entrypoint.sh