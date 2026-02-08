mkdir -p results/alignment_qc

for bam in /mnt/d/fastq_katsina_structured_bam/*.sorted.bam; do
  sample=$(basename "$bam" .sorted.bam)
  samtools flagstat "$bam" > results/alignment_qc/${sample}_flagstat.txt
done

echo -e "Sample\tTotal_Reads\tMapped_Reads\tMapped_Percent\tProperly_Paired_Reads\tProperly_Paired_Percent" \
> results/alignment_qc/flagstat_summary_FINAL_CLEAN.txt

for f in results/alignment_qc/*_flagstat.txt; do
  sample=$(basename "$f" _flagstat.txt)

  total=$(grep "in total" "$f" | awk '{print $1}')

  mapped_line=$(grep -m 1 " mapped (" "$f")
  mapped_reads=$(echo "$mapped_line" | awk '{print $1}')
  mapped_pct=$(echo "$mapped_line" | sed -E 's/.*\(([0-9.]+)%.*/\1/')

  proper_line=$(grep -m 1 "properly paired (" "$f")
  proper_reads=$(echo "$proper_line" | awk '{print $1}')
  proper_pct=$(echo "$proper_line" | sed -E 's/.*\(([0-9.]+)%.*/\1/')

  echo -e "${sample}\t${total}\t${mapped_reads}\t${mapped_pct}\t${proper_reads}\t${proper_pct}" \
  >> results/alignment_qc/flagstat_summary_FINAL_CLEAN.txt
done
