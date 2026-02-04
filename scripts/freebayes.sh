# #!/bin/bash
# set -euo pipefail

# # ======================================================
# # Config
# # ======================================================
# THREADS=8
# REF="data/genome/VectorBase-68_AfunestusFUMOZ_Genome.fasta"
# ALIGN_DIR="results/hisat2_alignments"
# RG_DIR="results/hisat2_alignments/rgfixed"
# VARIANT_DIR="results/freebayes"

# mkdir -p $RG_DIR
# mkdir -p $VARIANT_DIR

# # ======================================================
# # Helper: checkpoint function
# # ======================================================
# checkpoint () {
#     local STEP="$1"
#     local FILE="$2"

#     if [ -f "$FILE" ]; then
#         echo "[SKIP] $STEP — checkpoint exists: $FILE"
#         return 1
#     else
#         echo "[RUN] $STEP..."
#         return 0
#     fi
# }

# # ======================================================
# # 1. Index reference genome
# # ======================================================
# if checkpoint "Index reference genome" "$REF.fai"; then
#     samtools faidx $REF
# fi

# # ======================================================
# # 2. Fix Read Groups + coordinate sort + markdup + index
# # ======================================================
# for bam in $ALIGN_DIR/*sorted.dedup.bam; do
#     base=$(basename "$bam" .sorted.dedup.bam)
#     rg_bam="$RG_DIR/${base}.rg.bam"

#     # 2.1 Fix RG
#     if checkpoint "Fix RG $base" "$rg_bam.bai"; then
#         samtools addreplacerg \
#             -r "ID:$base" \
#             -r "SM:$base" \
#             -r "PL:ILLUMINA" \
#             -r "LB:lib1" \
#             -r "PU:unit1" \
#             -o "$rg_bam" \
#             "$bam"

#         samtools index "$rg_bam"
#     fi
# done

# # ======================================================
# # 3. Generate FreeBayes regions
# # ======================================================
# REGIONS_FILE="regions.txt"
# if checkpoint "Generate FreeBayes regions" "$REGIONS_FILE"; then
#     fasta_generate_regions.py ${REF}.fai 50000 > $REGIONS_FILE
# fi

# # ======================================================
# # 4. Run FreeBayes in parallel
# # ======================================================
# RAW_VCF="$VARIANT_DIR/all_samples_raw.vcf.gz"
# if checkpoint "FreeBayes variant calling" "$RAW_VCF"; then
#     freebayes-parallel $REGIONS_FILE $THREADS \
#         -f $REF \
#         $RG_DIR/*.rg.bam \
#         | bgzip -c > $RAW_VCF

#     tabix -p vcf $RAW_VCF
# fi

# # ======================================================
# # 5. Filter VCF
# # ======================================================
# FILTERED_VCF="$VARIANT_DIR/all_samples_filtered.vcf.gz"
# if checkpoint "Filter VCF" "$FILTERED_VCF"; then
#     bcftools filter \
#         -e 'QUAL<20 || DP<10' \
#         -Oz -o $FILTERED_VCF \
#         $RAW_VCF

#     bcftools index $FILTERED_VCF
# fi

# # ======================================================
# # 6. Extract SNPs only
# # ======================================================
# SNPS_VCF="$VARIANT_DIR/all_samples_snps_only.vcf.gz"
# if checkpoint "Extract SNPs" "$SNPS_VCF"; then
#     bcftools view -m2 -M2 -v snps \
#         -Oz -o $SNPS_VCF \
#         $FILTERED_VCF

#     bcftools index $SNPS_VCF
# fi

# echo "=============================================="
# echo "Pipeline completed successfully!"
# echo "=============================================="


#!/bin/bash
set -euo pipefail

# ======================================================
# Config
# ======================================================
THREADS=8
REF="data/genome/VectorBase-68_AfunestusFUMOZ_Genome.fasta"
RG_DIR="results/hisat2_alignments/rgfixed"
VARIANT_DIR="results/freebayes"

mkdir -p $VARIANT_DIR

# ======================================================
# Helper: checkpoint
# ======================================================
checkpoint () {
    local STEP="$1"
    local FILE="$2"

    if [ -f "$FILE" ]; then
        echo "[SKIP] $STEP — checkpoint exists: $FILE"
        return 1
    else
        echo "[RUN] $STEP..."
        return 0
    fi
}

# ======================================================
# 1. Index reference genome
# ======================================================
if checkpoint "Index reference genome" "$REF.fai"; then
    samtools faidx $REF
fi

# ======================================================
# 2. Generate FreeBayes genome regions
# ======================================================
REGIONS_FILE="$VARIANT_DIR/regions.txt"
if checkpoint "Generate FreeBayes regions" "$REGIONS_FILE"; then
    fasta_generate_regions.py ${REF}.fai 50000 > $REGIONS_FILE
fi

# ======================================================
# 3. Run FreeBayes in parallel on RG-fixed BAMs
# ======================================================
RAW_VCF="$VARIANT_DIR/all_samples_raw.vcf.gz"
if checkpoint "Run FreeBayes" "$RAW_VCF"; then
    echo "[RUN] Running FreeBayes on all RG-fixed BAMs"

    # Select only final RG.bam files (exclude intermediate fixes)
    BAM_LIST=$(ls $RG_DIR/*.rg.bam | grep -v "fixedRG")

    freebayes-parallel $REGIONS_FILE $THREADS \
        -f $REF \
        $BAM_LIST \
        | bgzip -c > $RAW_VCF

    tabix -p vcf $RAW_VCF
fi

# ======================================================
# 4. Filter variants
# ======================================================

FILTERED_VCF="$VARIANT_DIR/all_samples_filtered.vcf.gz"
if checkpoint "Filter VCF" "$FILTERED_VCF"; then
    bcftools filter \
        -e 'QUAL<20 || INFO/DP<10' \
        -Oz -o $FILTERED_VCF \
        $RAW_VCF

    bcftools index $FILTERED_VCF
fi

# ======================================================
# 5. Extract only SNPs
# ======================================================
SNPS_VCF="$VARIANT_DIR/all_samples_snps_only.vcf.gz"
if checkpoint "Extract SNPs" "$SNPS_VCF"; then
    bcftools view -m2 -M2 -v snps \
        -Oz -o $SNPS_VCF \
        $FILTERED_VCF

    bcftools index $SNPS_VCF
fi

echo "=============================================="
echo "FreeBayes pipeline completed successfully!"
echo "=============================================="
