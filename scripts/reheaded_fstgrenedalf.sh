awk '
NR==FNR { samples[NR]=$1; next }
{
  line=$0
  for (i=1; i<=length(samples); i++) {
    gsub("merged\\." i, samples[i], line)
  }
  print line
}
' /mnt/d/ibrahimRNAseq/variant_calling/samples.txt \
  /mnt/d/ibrahimRNAseq/fst_results_100kb/Afunestus_100kbFST.100kb.fst \
> /mnt/d/ibrahimRNAseq/fst_results_100kb/Afunestus_100kbFST.100kb_cleaned.fst
