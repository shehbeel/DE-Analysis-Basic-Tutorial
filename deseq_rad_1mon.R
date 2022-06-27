# Title: Differential Expression Analysis using DESeq2
# Author: Shehbeel Arif
# Source code adapted from Dr. Michael W. Vandewege's DESeq2 tutorial  
# Official DESeq2 tutorial: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html


# Download DESeq2 Library from Bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

# Load the libraries
library(DESeq2)

# Load count matrix
dat <- read.csv("GLDS-202_counts_rad_1month.csv", header=T, row.names=1)

# Load column data
coldata <- read.table("colData_rad_1month.txt", header=T, sep="\t")

# Create DESeq Data matrix object
dds <- DESeqDataSetFromMatrix(dat, coldata, ~condition)

# Remove lowly expressed genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Main DESeq
ddsDE <- DESeq(dds)

# Export normalized read counts
normCounts <- counts(ddsDE, normalized=TRUE)
write.csv(normCounts, "GLDS-202_normalized_counts.csv")

# DESeq Results
res <- results(ddsDE, alpha=0.05)

# To see a quick summary of the DE analysis results
summary(res)

# Export DE Results
resOrdered <- res[order(res$padj),]
write.csv(resOrdered, "GLDS-202_rad_1month_DE.csv")

# Plot MA plot for DE genes
plotMA(ddsDE, ylim=c(-5,5))


