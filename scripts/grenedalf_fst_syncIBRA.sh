#!/bin/bash

# Usage:
# ./fst_grenedalf_sync.sh sync_file output_dir pool_size threads prefix grenedalf_path bed_file file_suffix sample_file
# Example:
# ./fst_grenedalf_sync.sh /path/to/input.sync /path/to/output 40 40 MyPrefix /home/arnaud.tepa/pipeline/grenedalf/bin/grenedalf /path/to/intervals.bed .100kb sample.txt

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

# Validate sample file
if [ ! -f "$SAMPLE_FILE" ]; then
    echo "? Error: Sample file '$SAMPLE_FILE' not found."
    exit 1
fi

# Validate BED file
if [ ! -f "$BED_FILE" ]; then
    echo "? Error: BED file '$BED_FILE' not found."
    exit 1
fi

# Create output directory if not exists
mkdir -p "$OUT_DIR"

# Log file
LOG_FILE="${OUT_DIR}/${PREFIX}${FILE_SUFFIX}.log"

echo "? Starting Grenedalf FST calculation..."
echo "Input SYNC file: $SYNC_FILE"
echo "Output directory: $OUT_DIR"
echo "Pool size: $POOL_SIZE"
echo "Threads: $THREADS"
echo "Prefix: $PREFIX"
echo "Grenedalf path: $GRENDALF_BIN"
echo "BED file: $BED_FILE"
echo "File suffix: $FILE_SUFFIX"
echo "Sample file: $SAMPLE_FILE"

# Run Grenedalf FST using interval windows
$GRENDALF_BIN fst \
  --sync-path "$SYNC_FILE" \
  --multi-file-locus-set intersection \
  --filter-total-min-read-depth 10 \
  --filter-total-snp-min-frequency 0.01 \
  --window-type regions \
  --window-region-bed "$BED_FILE" \
  --window-average-policy window-length \
  --method unbiased-hudson \
  --pool-sizes "$POOL_SIZE" \
  --separator-char tab \
  --out-dir "$OUT_DIR" \
  --file-prefix "$PREFIX" \
  --file-suffix "$FILE_SUFFIX" \
  --allow-file-overwriting \
  --threads "$THREADS" \
  --log-file "$LOG_FILE"

# Detect Grenedalf output file (csv, tsv, txt)
OUTPUT_FILE=""
for ext in csv tsv txt fst; do
    if [ -f "${OUT_DIR}/${PREFIX}${FILE_SUFFIX}.${ext}" ]; then
        OUTPUT_FILE="${OUT_DIR}/${PREFIX}${FILE_SUFFIX}.${ext}"
        break
    fi
done

if [ -z "$OUTPUT_FILE" ]; then
    echo "? Error: No Grenedalf output file found with prefix ${PREFIX}${FILE_SUFFIX} in ${OUT_DIR}"
    exit 1
fi

# Ensure final output has .fst suffix
FINAL_FST_FILE="${OUT_DIR}/${PREFIX}${FILE_SUFFIX}.fst"
if [ "$OUTPUT_FILE" != "$FINAL_FST_FILE" ]; then
    mv "$OUTPUT_FILE" "$FINAL_FST_FILE"
fi

# Reheader function
function reheader_fst_file {
    local input_file="$1"
    local output_file="$2"

    # Read sample names into array
    mapfile -t samples < "$SAMPLE_FILE"

    # Extract header
    header=$(head -1 "$input_file")

    # Replace generic sample names with actual names
    for i in "${!samples[@]}"; do
        n=$((i+1))
        header=$(echo "$header" | sed "s/sample${n}/${samples[$i]}/g")
    done

    # Write new header
    echo -e "$header" > "$output_file"

    # Append data rows
    tail -n +2 "$input_file" >> "$output_file"
}

# Apply reheader
CLEANED_FILE="${OUT_DIR}/${PREFIX}${FILE_SUFFIX}_cleaned.fst"
reheader_fst_file "$FINAL_FST_FILE" "$CLEANED_FILE"
echo "? Reheader completed. Cleaned file saved: $CLEANED_FILE"

echo "? FST calculation and reheader completed!"
echo "Results saved in: $OUT_DIR"
