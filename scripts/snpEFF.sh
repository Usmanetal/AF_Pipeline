#!/bin/bash
#=============================================================
# SnpEff setup and annotation pipeline for Aedes funestus
#=============================================================

# Exit on any error
set -e

# Step 0: Setup directories
mkdir -p snpEff/software/snpEff/data/AfunestusFUMOZ
mkdir -p results/freebayes/annotated_vcf

# Step 1: Copy genome and annotation into SnpEff expected location
cp snpEff/data/AfunestusFUMOZ/sequences.fa snpEff/software/snpEff/data/AfunestusFUMOZ/
cp snpEff/data/AfunestusFUMOZ/genes.gff snpEff/software/snpEff/data/AfunestusFUMOZ/

# Step 2: Unzip SnpEff if not already done
mkdir -p snpEff/software
unzip -o snpEff_latest_core.zip -d snpEff/software

# Step 3: Add genome to snpEff.config if missing
CONFIG_FILE=~/.snpEff/snpEff.config
if ! grep -q "AfunestusFUMOZ.genome" "$CONFIG_FILE"; then
    echo "AfunestusFUMOZ.genome : Aedes funestus FUMOZ" >> "$CONFIG_FILE"
fi

# Step 4: Build the custom genome database
java -jar snpEff/software/snpEff/snpEff.jar build -gff3 -v AfunestusFUMOZ

# Step 5: Annotate all filtered VCFs
for f in results/freebayes/*.vcf.gz; do
    base=$(basename "$f" .vcf.gz)
    java -jar snpEff/software/snpEff/snpEff.jar -v AfunestusFUMOZ "$f" \
        > results/freebayes/annotated_vcf/${base}.ann.vcf
done

# Step 6: Generate SnpEff summary HTML and final annotated VCF
for f in results/freebayes/annotated_vcf/*.ann.vcf; do
    base=$(basename "$f" .ann.vcf)
    java -jar snpEff/software/snpEff/snpEff.jar -v -stats results/freebayes/annotated_vcf/${base}.html \
        AfunestusFUMOZ "$f" > results/freebayes/annotated_vcf/${base}.snpEff.vcf
done

# Step 7: Compress and index a representative VCF (all_samples_filtered)
bgzip -c results/freebayes/annotated_vcf/all_samples_filtered.ann.vcf \
    > results/freebayes/annotated_vcf/all_samples_filtered.ann.vcf.gz

tabix -p vcf results/freebayes/annotated_vcf/all_samples_filtered.ann.vcf.gz

# Step 8: Extract basic fields + SnpEff annotation to TSV
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%INFO/DP\t%INFO/AF\t%INFO/ANN\n' \
    results/freebayes/annotated_vcf/all_samples_filtered.ann.vcf.gz \
    > results/freebayes/annotated_vcf/all_samples_filtered.ann.tsv

echo "=============================================================="
echo "Annotation complete. All annotated VCFs and summary HTML files are in:"
echo "results/freebayes/annotated_vcf/"
echo "=============================================================="
