import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Create a dataframe with your data
data = {
    "Sample": [
        "9-FANG", "8-FANG", "7-FANG", "15-GGAFUN_UNX3", "14-GGAFUN_UNX2",
        "13-GGAFUN_UNX1", "12-GGAFUN_DELTA_3", "11-GGAFUN_DELTA_2", "10-GGAFUN_DELTA_1"
    ],
    "Mean_Depth": [76.8693, 85.885, 73.9465, 60.9314, 63.0588, 70.7261, 77.1842, 78.2865, 77.9854]
}

df = pd.DataFrame(data)

# Bar plot
plt.figure(figsize=(10,6))
sns.barplot(x="Sample", y="Mean_Depth", data=df, palette="viridis")
plt.title("Mean Sequencing Depth per Sample")
plt.ylabel("Mean Depth")
plt.xlabel("Sample")
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()



# Add missing data
df["F_MISS"] = [0.0483816, 0.0478906, 0.05408, 0.0658617, 0.0626125, 0.0503717, 0.0471871, 0.0458857, 0.0391775]

plt.figure(figsize=(10,6))
sns.barplot(x="Sample", y="F_MISS", data=df, palette="magma")
plt.title("Fraction of Missing Genotypes per Sample")
plt.ylabel("Fraction Missing")
plt.xlabel("Sample")
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Sample data
data = {
    "Sample": [
        "9-FANG", "8-FANG", "7-FANG", "15-GGAFUN_UNX3", "14-GGAFUN_UNX2",
        "13-GGAFUN_UNX1", "12-GGAFUN_DELTA_3", "11-GGAFUN_DELTA_2", "10-GGAFUN_DELTA_1"
    ],
    "Mean_Depth": [76.8693, 85.885, 73.9465, 60.9314, 63.0588, 70.7261, 77.1842, 78.2865, 77.9854],
    "F_MISS": [0.0483816, 0.0478906, 0.05408, 0.0658617, 0.0626125, 0.0503717, 0.0471871, 0.0458857, 0.0391775],
    "SNPs": [100000, 105000, 98000, 90000, 92000, 97000, 102000, 103000, 104000],
    "INDELs": [5000, 4800, 5200, 6000, 5800, 5100, 4900, 4700, 4600],
    # Simulated allele balance values per sample (heterozygous sites)
    "Allele_Balance": [
        np.random.beta(2,2,1000), np.random.beta(2,2,1000), np.random.beta(2,2,1000),
        np.random.beta(2,2,1000), np.random.beta(2,2,1000), np.random.beta(2,2,1000),
        np.random.beta(2,2,1000), np.random.beta(2,2,1000), np.random.beta(2,2,1000)
    ]
}

df = pd.DataFrame(data)


# Set style
sns.set(style="whitegrid")

# Create figure with 4 subplots
fig, axs = plt.subplots(2, 2, figsize=(16,12))

# Subplot 1: Mean depth
sns.barplot(x="Sample", y="Mean_Depth", data=df, palette="viridis", ax=axs[0,0])
axs[0,0].set_title("Mean Sequencing Depth per Sample")
axs[0,0].set_ylabel("Mean Depth")
axs[0,0].set_xlabel("")
axs[0,0].tick_params(axis='x', rotation=45)

# Subplot 2: Fraction missing
sns.barplot(x="Sample", y="F_MISS", data=df, palette="magma", ax=axs[0,1])
axs[0,1].set_title("Fraction of Missing Genotypes per Sample")
axs[0,1].set_ylabel("Fraction Missing")
axs[0,1].set_xlabel("")
axs[0,1].tick_params(axis='x', rotation=45)

# Subplot 3: SNP vs INDEL ratio per sample
snp_indel_ratio = df["SNPs"] / df["INDELs"]

# Create a temporary DataFrame for plotting
ratio_df = pd.DataFrame({
    "Sample": df["Sample"],
    "SNP_INDEL_Ratio": snp_indel_ratio
})

sns.barplot(x="Sample", y="SNP_INDEL_Ratio", data=ratio_df, palette="coolwarm", ax=axs[1, 0])
axs[1, 0].set_title("SNP/INDEL Ratio per Sample")
axs[1, 0].set_ylabel("SNP / INDEL")
axs[1, 0].set_xlabel("")
axs[1, 0].tick_params(axis='x', rotation=45)

# Subplot 4: Allele balance distribution
for i, sample in enumerate(df["Sample"]):
    sns.kdeplot(df["Allele_Balance"][i], label=sample, ax=axs[1,1])
axs[1,1].set_title("Allele Balance Distribution (Heterozygous Sites)")
axs[1,1].set_xlabel("Allele Balance (Ref / Total)")
axs[1,1].set_ylabel("Density")
axs[1,1].legend(title="Sample", bbox_to_anchor=(1.05,1), loc='upper left')

plt.tight_layout()

plt.savefig(
    "allele_balance_distribution.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

