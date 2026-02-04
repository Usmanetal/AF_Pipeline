import pandas as pd
df= pd.read_csv("results/RNA_AF_with_annotations_and_variant_label.csv")

from collections import defaultdict

rna_cols = [c for c in df.columns if "_RNA_" in c]

groups = defaultdict(list)

for col in rna_cols:
    if "FANG" in col or "FG7" in col:
        groups["FANG"].append(col)
    elif "UNX" in col:
        groups["UNX"].append(col)
    elif "DELTA" in col:
        groups["DELTA"].append(col)

# sort for reproducible replicate numbering
for g in groups:
    groups[g] = sorted(groups[g])
    
rename_map = {}

for group, cols in groups.items():
    for i, col in enumerate(cols, start=1):
        rename_map[col] = f"Nigeria_{group}_{i:03d}"

df = df.rename(columns=rename_map)

df_vaf= df[['CHROM', 'POS', 'REF', 'ALT', 'Gene_Symbol','aa_change','variant_label']+[col for col in df.columns if "Nigeria" in col]]



import pandas as pd
import numpy as np
from scipy.stats import ttest_ind

# -------------------------------
# 0️⃣ Define population replicate groups
# -------------------------------
delta_cols = ['Nigeria_DELTA_001', 'Nigeria_DELTA_002', 'Nigeria_DELTA_003']
fang_cols  = ['Nigeria_FANG_001',  'Nigeria_FANG_002',  'Nigeria_FANG_003']
unx_cols   = ['Nigeria_UNX_001',   'Nigeria_UNX_002',   'Nigeria_UNX_003']

allele_cols = delta_cols + fang_cols + unx_cols

# Ensure numeric
df_vaf[allele_cols] = df_vaf[allele_cols].astype(float)

# -------------------------------
# 1️⃣ Keep only SNPs with gene annotation
# -------------------------------
df_snps = df_vaf.dropna(subset=["Gene_Symbol"]).copy()

# -------------------------------
# 2️⃣ Welch t-test function per SNP
# -------------------------------
def welch_test(row, group1_cols, group2_cols):
    g1 = row[group1_cols].values.astype(float)
    g2 = row[group2_cols].values.astype(float)
    
    # remove NaNs
    g1 = g1[~np.isnan(g1)]
    g2 = g2[~np.isnan(g2)]
    
    if len(g1) < 2 or len(g2) < 2:
        return pd.Series([np.nan, np.nan, np.nan, np.nan, np.nan])
    
    t_stat, p_val = ttest_ind(g1, g2, equal_var=False)  # Welch's t-test
    mean1 = np.mean(g1)
    mean2 = np.mean(g2)
    diff  = mean1 - mean2
    
    return pd.Series([t_stat, p_val, mean1, mean2, diff])

# -------------------------------
# 3️⃣ Apply Welch's t-test for DELTA vs FANG
# -------------------------------
df_snps[['t_stat_delta_fang','pval_delta_fang','mean_delta','mean_fang','diff_delta_fang']] = \
    df_snps.apply(welch_test, axis=1, group1_cols=delta_cols, group2_cols=fang_cols)

# -------------------------------
# 3b️⃣ Optional: DELTA vs UNX
# -------------------------------
df_snps[['t_stat_delta_unx','pval_delta_unx','mean_delta2','mean_unx','diff_delta_unx']] = \
    df_snps.apply(welch_test, axis=1, group1_cols=delta_cols, group2_cols=unx_cols)

# -------------------------------
# 3c️⃣ Optional: FANG vs UNX
# -------------------------------
df_snps[['t_stat_fang_unx','pval_fang_unx','mean_fang2','mean_unx2','diff_fang_unx']] = \
    df_snps.apply(welch_test, axis=1, group1_cols=fang_cols, group2_cols=unx_cols)

# -------------------------------
# 4️⃣ Prepare Manhattan plot data
# Avoid -log10(0) issues by setting a small floor
# -------------------------------
df_snps['pval_delta_fang'] = df_snps['pval_delta_fang'].clip(lower=1e-300)
df_snps['neg_log10_pval_delta_fang'] = -np.log10(df_snps['pval_delta_fang'])

# -------------------------------
# 5️⃣ Top 40 SNPs for heatmap
# -------------------------------
top_snps = df_snps.sort_values('pval_delta_fang').head(40)

# -------------------------------
# 6️⃣ Save results
# -------------------------------
df_snps.to_csv("SNP_welch_results_delta_vs_fang.csv", index=False)
top_snps.to_csv("Top40_SNPs_delta_vs_fang.csv", index=False)

print("✅ Welch t-test per SNP completed!")
print("✅ Columns added: t_stat_delta_fang, pval_delta_fang, mean_delta, mean_fang, diff_delta_fang")
print("✅ Data prepared for Manhattan plot (-log10 p-values) and top 40 SNP heatmap.")


df = pd.read_csv("results/Manhattan_3panel_Welch_plot_data.csv")

import numpy as np

# n_tests = df.shape[0]
# bonf_threshold = -np.log10(0.05 / n_tests)
bonf_threshold = -np.log10(0.05)


filtered_df = (
    df[df['neg_log10_pval_delta_unx'] >= bonf_threshold]
    .sort_values('neg_log10_pval_delta_unx', ascending=False)
)

cols_of_interest = [
    'CHROM', 'POS', 'REF', 'ALT',
    'Gene_Symbol', 'aa_change', 'variant_label',
    'neg_log10_pval_delta_fang'
]

filtered_df = filtered_df[cols_of_interest]

filtered_df.shape

##########################################
# correcting for hochberg
##########################################

from statsmodels.stats.multitest import multipletests
import numpy as np

pval_map = {
    "delta_fang": "pval_delta_fang",
    "delta_unx":  "pval_delta_unx",
    "fang_unx":   "pval_fang_unx"
}

for tag, pcol in pval_map.items():
    qvals = multipletests(
        df[pcol],
        method="fdr_bh"
    )[1]
    df[f"qval_{tag}"] = qvals
    df[f"neg_log10_qval_{tag}"] = -np.log10(qvals)

df.to_csv("sig_df.fdr_bh.csv", index=False)
################################################################
# difining contrasts 0.05 FDR threshold
################################################################

contrast_map = {
    "delta_fang": {
        "qval": "qval_delta_fang",
        "neglogq": "neg_log10_qval_delta_fang"
    },
    "delta_unx": {
        "qval": "qval_delta_unx",
        "neglogq": "neg_log10_qval_delta_unx"
    },
    "fang_unx": {
        "qval": "qval_fang_unx",
        "neglogq": "neg_log10_qval_fang_unx"
    }
}

import pandas as pd

fdr_threshold = 0.05

# Load gene annotation once
df2 = pd.read_csv("data/genome/Afunestus_protein_coding_genes.csv")

def add_gene_info(df):
    """Merge gene info from df2"""
    return df.merge(
        df2[['ID', 'Name', 'description']],
        left_on='Gene_Symbol',
        right_on='ID',
        how='left'
    )

# Loop over all contrasts
for contrast in contrast_map:
    qcol = contrast_map[contrast]["qval"]
    nlogcol = contrast_map[contrast]["neglogq"]

    # Filter for FDR ≤ threshold
    filtered = (
        df[df[qcol] <= fdr_threshold]
        .sort_values(nlogcol, ascending=False)
        .copy()
    )

    # Annotate with gene info
    annotated = add_gene_info(filtered)

    # Save to CSV
    filename = f"{contrast}_FDR{fdr_threshold}.csv"
    annotated.to_csv(filename, index=False)
    print(f"Saved {filename}, {annotated.shape[0]} significant variants")


#########################################################
# fang vs unx (FDR-BH corrected) manhattan plot
#########################################################

contrast = "delta_fang"

qcol = f"qval_{contrast}"
ycol = f"neg_log10_qval_{contrast}"

chrom_order = (
    df['CHROM']
    .drop_duplicates()
    .sort_values()
    .tolist()
)

chrom_order = [c for c in chrom_order if c != "AfunF3_X"] + ["AfunF3_X"]

df = df.copy()

chrom_offsets = {}
current_offset = 0

for chrom in chrom_order:
    chrom_offsets[chrom] = current_offset
    current_offset += df.loc[df['CHROM'] == chrom, 'POS'].max()

df['cumulative_pos'] = (
    df['CHROM'].map(chrom_offsets) + df['POS']
)

xticks = []
xticklabels = []

for chrom in chrom_order:
    chrom_df = df[df['CHROM'] == chrom]
    center = chrom_df['cumulative_pos'].median()
    xticks.append(center)
    xticklabels.append(chrom)

import matplotlib.pyplot as plt
import numpy as np

plt.figure(figsize=(14, 6))

colors = ["#4C72B0", "#DD8452"]  # alternating colors

for i, chrom in enumerate(chrom_order):
    chrom_df = df[df['CHROM'] == chrom]
    plt.scatter(
        chrom_df['cumulative_pos'],
        chrom_df[ycol],
        s=8,
        color=colors[i % 2],
        alpha=0.7,
        edgecolor="none"
    )

# FDR threshold line
fdr_threshold = 0.05
plt.axhline(
    -np.log10(fdr_threshold),
    color="red",
    linestyle="--",
    linewidth=1,
    label="FDR = 0.05"
)

plt.xticks(xticks, xticklabels, rotation=90)
plt.xlabel("Chromosome")
plt.ylabel(r"$-\log_{10}(\mathrm{q\text{-}value})$")
plt.title("Manhattan plot: FANG vs UNX")

plt.legend()
plt.tight_layout()
plt.show()

###########################################################
# Manhattan Plot
#############################################################

fig, ax = plt.subplots(figsize=(14, 10))

# Plot all points chromosome by chromosome
for i, chrom in enumerate(chrom_order):
    chrom_df = df[df['CHROM'] == chrom]

    # All points in light gray
    ax.scatter(
        chrom_df['cumulative_pos'],
        chrom_df[ycol],
        s=6,
        color="lightgrey",
        alpha=0.5
    )

    # Highlight significant points in the chromosome-specific color
    sig_chrom_df = chrom_df[chrom_df[qcol] <= 0.05]
    if not sig_chrom_df.empty:
        ax.scatter(
            sig_chrom_df['cumulative_pos'],
            sig_chrom_df[ycol],
            s=12,
            color=colors[i % 2],
            alpha=0.9
        )

# FDR threshold line
ax.axhline(-np.log10(0.05), color="black", linestyle="--")

# Vertical bars for top 10
ymax = df[ycol].max() * 1.25
ax.set_ylim(0, ymax)
for _, row in top10.iterrows():
    ax.vlines(
        row['cumulative_pos'],
        ymin=0,
        ymax=ymax,
        color="grey",
        alpha=0.15,
        linewidth=6,
        zorder=0
    )

# Annotations for top 10
for _, row in top10.iterrows():
    ax.annotate(
        row['Gene_Symbol'],
        xy=(row['cumulative_pos'], row[ycol]),
        xytext=(0, 6),
        textcoords="offset points",
        ha="center",
        va="bottom",
        fontsize=9,
        fontweight="bold",
        rotation=0,
        clip_on=False,
        zorder=10
    )

# Axis labels
ax.set_xlabel("Chromosome", fontsize=12, fontweight="bold")
ax.set_ylabel(r"$-\log_{10}(\mathrm{q\text{-}value})$", fontsize=12, fontweight="bold")

# X-axis ticks
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels, rotation=0)
# X-axis tick labels
for i, tick in enumerate(ax.get_xticklabels()):
    tick.set_color(colors[i % 2])
    tick.set_fontweight("bold")

# Y-axis tick labels
for tick in ax.get_yticklabels():
    tick.set_fontweight("bold")

# Spines
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.title(f"Manhattan plot: {contrast.replace('_',' ').upper()} (significant highlighted)")
plt.tight_layout()

# Save once
fig.savefig(f"Manhattan_{contrast}_top10.png", dpi=300, bbox_inches="tight")
fig.savefig(f"Manhattan_{contrast}_top10.pdf", bbox_inches="tight")
plt.show()
