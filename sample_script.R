filelist <- list.files("/Users/luqiao/Desktop/604 project/data_cleanup", full.names=TRUE)

# import data
PRAD <- read_delim("data_cleanup/PRAD.tsv", 
                   +     delim = "\t", escape_double = FALSE, 
                   +     trim_ws = TRUE)

merged_sample_sheet <- read_delim("merged_sample_sheet.tsv", 
                                  +     delim = "\t", escape_double = FALSE, 
                                  +     trim_ws = TRUE)

data.matrix(PRAD)

# read data
coldata <- merged_sample_sheet[,c("Project ID","Sample ID", "Sample Type")]

# select project
PRAD_ID <- coldata[coldata$`Project ID` == "TCGA-PRAD",]
head(PRAD_ID)

# select required col
PRAD_ID <- PRAD_ID[,c("Sample ID", "Sample Type")]
PRAD_ID$`Sample Type` <- factor(coldata$`Sample Type`)

# set row name
ID <- PRAD_ID$`Sample ID`
rownames(PRAD_ID) <- ID


# delete gene name col
PRAD <- subset(PRAD, select = -c(gene_name))

PRAD_no_col <- subset(PRAD, select = -c(gene_id))

gene_id <- PRAD$gene_id

rownames(PRAD_no_col) <- gene_id

PRAD_ID_final <- PRAD_ID[,c("Sample Type")]

sample_id <- colnames(PRAD_no_col)

filtered_ID <- PRAD_ID[PRAD_ID$`Sample ID` %in% sample_id,]

rownames(filtered_ID) <- filtered_ID$`Sample ID`

all(rownames(filtered_ID)) %in% colnames(PRAD_no_col)

# sort
PRAD_ordered <- PRAD_no_col[, rownames(filtered_ID)]

rownames(PRAD_ordered) <- gene_id

colnames(filtered_ID) <- sub(" ", "_", colnames(filtered_ID))

data.matrix(filtered_ID)
data.matrix(PRAD_ordered)

filtered_ID$`Sample_Type` <- factor(filtered_ID$`Sample_Type`)

### load packages and create DESeq data structure
library("DESeq2")

test_PRAD <- DESeqDataSetFromMatrix(countData = PRAD_ordered, colData = filtered_ID, design = ~ Sample_Type)

### pre filtering

keep <- rowSums(counts(test_PRAD)) >= 10
pre_filtered_PRAD <- test_PRAD[keep,]


### run DESeq

test_run_PRAD <- DESeq(pre_filtered_PRAD)
res <- results(test_run_PRAD)
res

### test with different conditions, but here we only have one condition to compare so we don't have to do this actually
res_c <- results(test_run_PRAD, contrast=c("Sample_Type","Primary Tumor","Solid Tissue Normal"))
res_c

### Log fold change shrinkage for visualization and ranking
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("apeglm")

library("apeglm")
resultsNames(test_run_PRAD)
resLFC <- lfcShrink(test_run_PRAD, coef = "Sample_Type_Solid.Tissue.Normal_vs_Primary.Tumor", type = "apeglm")
resLFC

### using 'apeglm' for LFC shrinkage. If used in published research, please cite: Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895

resOrdered <- resLFC[order(resLFC$pvalue),]
summary(resOrdered)
# How many adjusted p-values were less than 0.05?
sum(resOrdered$padj < 0.05, na.rm=TRUE)
sum(resOrdered$padj < 0.001, na.rm=TRUE)
sum(resOrdered$padj < 0.0001, na.rm=TRUE)

plotMA(resLFC, ylim=c(-2,2))

plotCounts(test_run_PRAD, gene=which.min(resLFC$padj), intgroup="Sample_Type")

# write to csv
resSig <- subset(resOrdered, padj < 0.0001)
write.csv(as.data.frame(resSig), 
          file="CANCER_vs_NORMAL.csv")



sun