
#!/bin/bash

# ================================
# QC pipeline for FreeBayes VCF
# ================================

# Input VCF
VCF="results/freebayes/all_samples_filtered.vcf.gz"

# Output directories
PLOTS_DIR="plots"
QC_DIR="qc_metrics"

mkdir -p ${PLOTS_DIR}
mkdir -p ${QC_DIR}

echo "Starting VCF QC pipeline..."

# 1️⃣ General stats using bcftools
echo "Generating bcftools stats..."
bcftools stats ${VCF} > ${QC_DIR}/all_samples.stats

echo "Plotting VCF stats..."
plot-vcfstats -p ${PLOTS_DIR}/ ${QC_DIR}/all_samples.stats

# 2️⃣ Per-sample depth (Mean Depth)
echo "Calculating per-sample depth..."
vcftools --gzvcf ${VCF} --depth --out ${QC_DIR}/depth

# 3️⃣ Missing genotypes per sample
echo "Calculating missing genotypes per sample..."
vcftools --gzvcf ${VCF} --missing-indv --out ${QC_DIR}/missing_indv

# 4️⃣ SNP counts per sample
echo "Counting SNPs per sample..."
vcftools --gzvcf ${VCF} --remove-indels --freq2 --out ${QC_DIR}/snps

# Count total SNPs
echo "Total SNPs:"
awk '{if(NR>1) count++} END {print count}' ${QC_DIR}/snps.frq

# 5️⃣ INDEL counts per sample
echo "Counting INDELs per sample..."
vcftools --gzvcf ${VCF} --keep-only-indels --freq2 --out ${QC_DIR}/indels

# Count total INDELs
echo "Total INDELs:"
awk '{if(NR>1) count++} END {print count}' ${QC_DIR}/indels.frq

# 6️⃣ Allele balance at heterozygous sites
echo "Extracting allele depths (AD) for allele balance..."
bcftools query -f '%CHROM\t%POS[\t%SAMPLE=%AD]\n' ${VCF} > ${QC_DIR}/allele_depth.txt

echo "Allele balance can now be computed in Python or R using:"
echo "AB = AD_ref / (AD_ref + AD_alt)"

# 7️⃣ Transition / Transversion ratio (Ti/Tv)
echo "Calculating Ti/Tv ratio..."
vcftools --gzvcf ${VCF} --TsTv-summary --out ${QC_DIR}/tstv

# 8️⃣ Variant missingness per site
echo "Calculating variant missingness per site..."
vcftools --gzvcf ${VCF} --missing-site --out ${QC_DIR}/missing_site

echo "VCF QC pipeline complete. Outputs saved in ${QC_DIR} and ${PLOTS_DIR}."
