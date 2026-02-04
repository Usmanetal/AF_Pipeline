import pandas as pd
snps_df = pd.read_csv("delta_fang_FDR0.05.csv")
snps_fang_unx = pd.read_csv("delta_unx_FDR0.05")

import re

snps_snv = snps_df[
    snps_df["REF"].str.fullmatch(r"[ACGT]", na=False) &
    snps_df["ALT"].str.fullmatch(r"[ACGT]", na=False)
]

snps_snv[["REF", "ALT"]].drop_duplicates()

import pandas as pd

# Read Excel sheet
target = pd.read_excel(
    "gene_interest.xlsx",
    sheet_name="colored"
)

conditions = [
    target["description"].str.contains("copper", case=False, na=False),
    target["description"].str.contains("heme peroxidase", case=False, na=False),
    target["description"].str.contains("Aminomethyltransferase", case=False, na=False),
    target["description"].str.contains("cytochrome b5", case=False, na=False)
]

choices = [
    "SOD",
    "heme peroxidase",
    "Aminomethyltransferase",
    "CYPb5"
]

import numpy as np

target["Name"] = np.select(
    conditions,
    choices,
    default=target["Name"]
)

target["Name"] = target["Name"].replace(
    {"CYP6P9A": "CYP6P9a"}
)
# target = target.rename(columns={"CYP6P9A": "CYP6P9b"})


# Convert columns to numeric (like as.numeric in R)
target["start"] = pd.to_numeric(target["start"], errors="coerce")
target["end"]   = pd.to_numeric(target["end"], errors="coerce")



target_gene = target[["genes", "Name"]].drop_duplicates()
snps_gene_merged = target_gene.merge(
    snps_snv,
    left_on="genes",
    right_on="Gene_Symbol",
    how="inner"   # keeps only SNPs in target genes
)

df = snps_gene_merged.drop(columns=["genes"])


# df["gene_variant"] = (
#     df["Name_x"]
#     .fillna(df["Gene_Symbol"])
#     .str.strip()
#     + " | "
#     + df["variant_label"].str.strip()
# )


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Filter for CYP genes
cyp_snps = df[df["Name_x"].str.contains("cyp", case=False, na=False)].copy()

# Create a new label combining Name and variant_label
cyp_snps['gene_variant'] = cyp_snps['Name_x'] + " | " + cyp_snps['variant_label']

# Select only the allele frequency columns
freq_cols = [c for c in df.columns if c.startswith("Nigeria_")]
freq_data = cyp_snps[freq_cols]

# Set the gene_variant as the index for the heatmap
freq_data.index = cyp_snps['gene_variant']

# Plot the heatmap with frequency values annotated
plt.figure(figsize=(12, max(6, 0.3*len(freq_data))))  # adjust height to number of variants
sns.heatmap(
    freq_data,
    cmap="viridis",
    linewidths=0.5,
    linecolor='gray',
    cbar_kws={'label': 'Allele Frequency'},
    annot=True,        # Show allele frequency values on cells
    fmt=".2f"          # Format to 2 decimal places
)

plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)
plt.title("Allele Frequency Heatmap for CYP Genes", fontsize=14, fontweight='bold')
plt.tight_layout()
plt.show()

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def plot_gene_heatmap(df, gene_keyword, title=None, save_path=None):
    """
    Plot and optionally save an allele frequency heatmap for a specific gene family.

    Parameters:
    - df: DataFrame containing SNP data
    - gene_keyword: str, keyword to filter gene names (e.g., 'cyp', 'gst', 'coe')
    - title: str, optional custom title for the plot
    - save_path: str, optional path to save the figure (e.g., 'heatmap_cyp.png')
    """
    # Filter for the gene of interest
    gene_snps = df[df["Name_x"].str.contains(gene_keyword, case=False, na=False)].copy()

    if gene_snps.empty:
        print(f"No variants found for gene keyword: {gene_keyword}")
        return

    # Create a new label combining Name and variant_label
    gene_snps['gene_variant'] = gene_snps['Name_x'] + " | " + gene_snps['variant_label']

    # Select only the allele frequency columns
    freq_cols = [c for c in df.columns if c.startswith("Nigeria_")]
    freq_data = gene_snps[freq_cols]

    # Set the gene_variant as the index for the heatmap
    freq_data.index = gene_snps['gene_variant']

    # Set title
    if title is None:
        title = f"Allele Frequency Heatmap for {gene_keyword.upper()} Genes"

    # Plot the heatmap
    plt.figure(figsize=(12, max(6, 0.3*len(freq_data))))  # adjust height to number of variants
    sns.heatmap(
        freq_data,
        cmap="viridis",
        linewidths=0.5,
        linecolor='gray',
        cbar_kws={'label': 'Allele Frequency'},
        annot=True,
        fmt=".2f"
    )

    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.title(title, fontsize=14, fontweight='bold')
    plt.tight_layout()

    # Save the figure if save_path is provided
    if save_path is not None:
        plt.savefig(save_path, dpi=300)
        print(f"Heatmap saved to {save_path}")

    plt.show()


# Example usage:
plot_gene_heatmap(df, gene_keyword="cyp", save_path="heatmap_cyp.png")
plot_gene_heatmap(df, gene_keyword="gst", save_path="heatmap_gst.png")
plot_gene_heatmap(df, gene_keyword="obp|nad|chit|heme|sod|aldeh", save_path="heatmap_coe.png")

df[
    ['Name_x', 'CHROM', 'POS', 'REF', 'ALT', 'Gene_Symbol', 'aa_change']
    + [c for c in df.columns if "Nigeria_" in c]
].to_csv("nigeria_variants.csv", index=False)

