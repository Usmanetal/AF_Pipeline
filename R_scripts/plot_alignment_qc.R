# ---------------------------------------
# RNA-seq Alignment QC (Publication Grade)
# Two-panel figure using ggarrange
# ---------------------------------------

library(tidyverse)
library(ggpubr)

# Load QC table
qc <- read.delim("results/alignment_qc/flagstat_summary_FINAL_CLEAN.txt")

# Preserve sample order
qc$Sample <- factor(qc$Sample, levels = qc$Sample)

# ------------------------
# Plot 1: Mapped reads (%)
# ------------------------
p_mapped <- ggplot(qc, aes(x = Sample, y = Mapped_Percent)) +
  geom_bar(stat = "identity", fill = "#4C72B0", width = 0.7) +
  ylim(0, 100) +
  labs(
    title = "Mapped reads",
    y = "Percentage of reads",
    x = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# --------------------------------
# Plot 2: Properly paired reads (%)
# --------------------------------
p_paired <- ggplot(qc, aes(x = Sample, y = Properly_Paired_Percent)) +
  geom_bar(stat = "identity", fill = "#55A868", width = 0.7) +
  ylim(0, 100) +
  labs(
    title = "Properly paired reads",
    y = "Percentage of reads",
    x = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# ------------------------
# Combine plots (1 row × 2 columns)
# ------------------------
qc_combined <- ggarrange(
  p_mapped, p_paired,
  nrow = 1,
  ncol = 2,
  labels = c("A", "B"),
  font.label = list(size = 14, face = "bold")
)

# ------------------------
# Save high-resolution PNG
# ------------------------
ggsave(
  filename = "results/alignment_qc/alignment_qc_2panel.png",
  plot = qc_combined,
  width = 14,
  height = 5,
  dpi = 300
)
