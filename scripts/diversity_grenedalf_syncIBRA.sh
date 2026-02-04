#!/bin/bash
set -euo pipefail

# Usage:
# ./diversity_grenedalf_sync.sh sync_file output_dir pool_size threads prefix grenedalf_path bed_file file_suffix sample_file

if [ "$#" -ne 9 ]; then
    echo "Usage: $0 sync_file output_dir pool_size threads prefix grenedalf_path bed_file file_suffix sample_file"
    exit 1
fi

SYNC_FILE=$1
OUT_DIR=$2
POOL_SIZE=$3
THREADS=$4
PREFIX=$5
GRENDALF_BIN=$6
BED_FILE=$7
FILE_SUFFIX=$8
SAMPLE_FILE=$9

# -------------------------
# Validation
# -------------------------
for f in "$SYNC_FILE" "$BED_FILE" "$SAMPLE_FILE" "$GRENDALF_BIN"; do
    if [ ! -e "$f" ]; then
        echo "❌ Error: File not found -> $f"
        exit 1
    fi
done

mkdir -p "$OUT_DIR"

LOG_FILE="${OUT_DIR}/${PREFIX}${FILE_SUFFIX}.log"

echo "🧬 Starting Grenedalf Diversity calculation"
echo "SYNC:      $SYNC_FILE"
echo "OUT DIR:   $OUT_DIR"
echo "POOL SIZE: $POOL_SIZE"
echo "THREADS:   $THREADS"
echo "PREFIX:    $PREFIX"
echo "BED:       $BED_FILE"
echo "SUFFIX:    $FILE_SUFFIX"
echo "SAMPLES:   $SAMPLE_FILE"

# -------------------------
# Run Grenedalf
# -------------------------
"$GRENDALF_BIN" diversity \
  --sync-path "$SYNC_FILE" \
  --multi-file-locus-set intersection \
  --filter-sample-min-count 2 \
  --filter-sample-min-read-depth 10 \
  --filter-total-min-read-depth 10 \
  --filter-total-snp-min-frequency 0.01 \
  --window-type regions \
  --window-region-bed "$BED_FILE" \
  --window-average-policy window-length \
  --pool-sizes "$POOL_SIZE" \
  --separator-char tab \
  --out-dir "$OUT_DIR" \
  --file-prefix "$PREFIX" \
  --file-suffix "$FILE_SUFFIX" \
  --allow-file-overwriting \
  --threads "$THREADS" \
  --log-file "$LOG_FILE" \
  # ✅ Ensure theta_pi and Tajima's D are computed
  --compute-theta-pi \
  --compute-tajimas-d

# -------------------------
# Detect output file
# -------------------------
OUTPUT_FILE=""
for ext in txt tsv csv; do
    f="${OUT_DIR}/${PREFIX}${FILE_SUFFIX}.${ext}"
    if [ -f "$f" ]; then
        OUTPUT_FILE="$f"
        break
    fi
done

if [ -z "$OUTPUT_FILE" ]; then
    echo "❌ Error: Grenedalf output not found"
    exit 1
fi

echo "✔ Found output: $OUTPUT_FILE"

# -------------------------
# Clean + reheader
# -------------------------
FINAL_FILE="${OUT_DIR}/${PREFIX}${FILE_SUFFIX}_cleaned.txt"

mapfile -t SAMPLES < "$SAMPLE_FILE"

# Extract base sync name (used by Grenedalf)
SYNC_BASE=$(basename "$SYNC_FILE")
SYNC_BASE=${SYNC_BASE%%.*}

# Identify columns to keep (chrom, start, end, theta_pi, theta_watterson, tajimas_d)
KEEP_REGEX="^(chrom|start|end|.*theta_pi|.*theta_watterson|.*tajimas_d)$"

awk -v keep="$KEEP_REGEX" -v OFS="\t" '
NR==1 {
    for (i=1; i<=NF; i++) {
        if ($i ~ keep) {
            cols[++n] = i
            header[n] = $i
        }
    }
    for (i=1; i<=n; i++) {
        printf "%s%s", header[i], (i<n ? OFS : ORS)
    }
    next
}
{
    out=""
    for (i=1; i<=n; i++) {
        val = $(cols[i])
        if (val=="") val="NA"  # Replace empty with NA
        out = out ? out OFS val : val
    }
    print out
}
' "$OUTPUT_FILE" > "$FINAL_FILE"

# -------------------------
# Replace sample placeholders with actual sample names
# -------------------------
for i in "${!SAMPLES[@]}"; do
    n=$((i+1))
    sed -i "s/${SYNC_BASE}.${n}/${SAMPLES[$i]}/g" "$FINAL_FILE"
done

echo "✔ Cleaned diversity file saved:"
echo "  $FINAL_FILE"
echo "✔ Log file:"
echo "  $LOG_FILE"
echo "🎉 Diversity analysis completed"
