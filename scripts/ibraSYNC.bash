#!/bin/bash
set -e

# ---------------------------
# User-defined variables
# ---------------------------
mpileup2sync_PATH=$1     # Path to mpileup2sync jar
REF=$2                   # Reference FASTA
FAI_FILE=$3              # Reference .fai file
BAM_DIR=$4               # Directory containing BAM files
WORK_DIR=$5              # Output directory
VARSCAN_PATH=$6          # Path to VarScan jar

REGIONS_FILE="$WORK_DIR/regions.txt"

# ---------------------------
# Create output directory
# ---------------------------
mkdir -p "$WORK_DIR"


# ---------------------------
# Generate regions from .fai (if not exists)
# ---------------------------
if [ ! -f "$REGIONS_FILE" ]; then
    echo "Generating regions from $FAI_FILE..."
    cut -f 1,2 "$FAI_FILE" | awk '
    {
        chrom=$1; size=$2; region_size=1000000000;
        for (start=1; start<=size; start+=region_size) {
            end=start+region_size-1;
            if (end>size) end=size;
            print chrom":"start"-"end;
        }
    }' > "$REGIONS_FILE"
fi

# ---------------------------
# Find all BAM files
# ---------------------------
ls "$BAM_DIR"/*.bam > bam_files.txt

# ---------------------------
# Run mpileup, mpileup2sync and VarScan per-region in parallel
# ---------------------------

MAX_JOBS=10  # Adjust based on your system resources

# ---------------------------
# Loop 1: Generate mpileup files
# ---------------------------
echo "Step 1: Generating mpileup files..."
for region in $(cat "$REGIONS_FILE"); do
  (
    TMP_MP="$WORK_DIR/${region//[:\-]/_}.mpileup"
    BAMS=$(awk '{printf " %s", $0}' bam_files.txt)
    samtools mpileup -f "$REF" $BAMS -r "$region" > "$TMP_MP"
  ) &
  
  # Limit parallel jobs
  while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do sleep 5; done
done
wait

# ---------------------------
# Loop 2: mpileup2sync
# ---------------------------
echo "Step 2: Running mpileup2sync..."
for region in $(cat "$REGIONS_FILE"); do
  (
    TMP_MP="$WORK_DIR/${region//[:\-]/_}.mpileup"
    OUT_SYNC="$WORK_DIR/variant_${region//[:\-]/_}.sync"
    if [ -s "$TMP_MP" ]; then
    perl "$mpileup2sync_PATH" \
      --input "$TMP_MP" \
      --output "$OUT_SYNC" \
      --fastq-type illumina \
      --min-qual 30
    fi
  ) &
  
  while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do sleep 5; done
done
wait

# ---------------------------
# Loop 3: VarScan
# ---------------------------
echo "Step 3: Running VarScan..."
for region in $(cat "$REGIONS_FILE"); do
  for mp in AfunF3_*.mpileup; do
  (
    base=$(basename "$mp" .mpileup)
    varscan mpileup2snp "$mp" \
      --min-avg-qual 30 \
      --min-coverage 10 \
      --min-reads2 5 \
      --min-var-freq 0.25 \
      --p-value 0.01 \
      --min-freq-for-hom 0.75 \
      --output-vcf 1 \
      > "variant_${base}.vcf"
  ) &
  
  while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do sleep 5; done
done

wait  # Ensure all parallel jobs are finished

# ---------------------------
# Merge all .sync files
# ---------------------------
echo "Merging sync files..."
find "$WORK_DIR" -name 'variant_*.sync' -print0 | xargs -0 cat > "$WORK_DIR/merged.sync.unsorted"
sort -k1,1 -k2,2n "$WORK_DIR/merged.sync.unsorted" > "$WORK_DIR/merged.sync"

# ---------------------------
# Merge all .vcf files
# ---------------------------
echo "Merging vcf files with header preserved..."

# Extract header from the first VCF
FIRST_VCF=$(find "$WORK_DIR" -name 'variant_*.vcf' | head -n 1)
grep '^#' "$FIRST_VCF" > "$WORK_DIR/merged.vcf"

# Append all non-header lines from all VCFs
find "$WORK_DIR" -name 'variant_*.vcf' -exec grep -v '^#' {} \; >> "$WORK_DIR/merged.vcf.unsorted"

# Sort and append to final VCF
sort -k1,1 -k2,2n "$WORK_DIR/merged.vcf.unsorted" >> "$WORK_DIR/merged.vcf"

# ---------------------------
# Clean up intermediate files
# ---------------------------
echo "Cleaning up intermediate files..."
find "$WORK_DIR" -name 'variant_*.sync' -delete
find "$WORK_DIR" -name 'variant_*.vcf' -delete
rm -f "$WORK_DIR/merged.sync.unsorted"
rm -f "$WORK_DIR/merged.vcf.unsorted"

echo "? All done."
echo "Final sync file: $WORK_DIR/merged.sync"
echo "Final VCF file: $WORK_DIR/merged.vcf"