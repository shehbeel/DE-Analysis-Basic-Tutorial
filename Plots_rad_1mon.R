# Title: Plotting Differentially Expressed Genes
# Author: Shehbeel Arif
# Source code adapted from Dr. Michael W. Vandewege's DESeq2 tutorial 

# Install libraries
install.packages("ggplot2")
install.packages("pheatmap")

# Load libraries
library(ggplot2)
library(pheatmap)

# Read in normalized gene counts data
normCount <- read.csv("GLDS-202_rad_1month_normalized_counts.csv", row.names=1)

# Read in DE data
deSeqRes <- read.csv("GLDS-202_rad_1month_DE.csv", row.names=1)

# If DE gene has adjusted p-value <= 0.05, then add a label of "yes"
deSeqRes$sig <- ifelse(deSeqRes$padj <= 0.05, "yes", "no")

# Remove "N/A" values from data                       
deSeqRes <- na.omit(deSeqRes)

# MA-plot
ggplot(deSeqRes, aes(x=log2(baseMean), y=log2FoldChange, color=sig)) + 
  geom_point()

# Volcano Plot
ggplot(deSeqRes, aes(x=log2FoldChange, y=-log2(padj), , color=sig)) + 
  geom_point()

# Heatmap
# Plotting all the genes
data <- as.matrix(normCount)
pheatmap(data, scale="row")

# Plotting only significantly expressed genes
signi <- subset(deSeqRes, padj<=0.05)
allSig <- merge(normCount, signi, by=0)

sigCounts <- allSig[,2:13]
row.names(sigCounts) <- allSig$Row.names

pheatmap(log2(sigCounts+1), scale="row")















