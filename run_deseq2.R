# load packages
library(readr)
library("DESeq2")
# load all data
file_list <- list.files(path = "./data_Cleanup")
num_file <- length(file_list)
merged_sample_sheet <- read_tsv("merged_sample_sheet.tsv")
coldata <- merged_sample_sheet[,c("Project ID","Sample ID", "Sample Type")]
for (x in 1:num_file) {
  file <- paste("data_cleanup/", file_list[x], sep="")
  project <- substr(file_list[x], 1, 4)
  table <- read_tsv(file)
  ID <- coldata[coldata$`Project ID` == paste("TCGA-", project, sep=""),]
  # select required field
  ID <- ID[,c("Sample ID", "Sample Type")]
  # delete gene name col
  gene_id <- table$gene_id
  table <- subset(table, select = -c(gene_name, gene_id))
  rownames(table) <- gene_id
  filtered_ID <- unique(ID[ID$`Sample ID` %in% colnames(table),])
  rownames(filtered_ID) <- filtered_ID$`Sample ID`
  # sort
  table_ordered <- table[, rownames(filtered_ID)]
  rownames(table_ordered) <- gene_id
  colnames(filtered_ID) <- sub(" ", "_", colnames(filtered_ID))
  filtered_ID$`Sample_Type` <- factor(filtered_ID$`Sample_Type`)
  # create DESeq data structure
  deseq_data <- DESeqDataSetFromMatrix(countData = table_ordered, colData = filtered_ID, design = ~ Sample_Type)
  # pre filtering
  keep <- rowSums(counts(deseq_data)) >= 10
  pre_filtered <- deseq_data[keep,]
  # run DESeq
  deseq <- DESeq(pre_filtered)
  res <- results(deseq)
  write.csv(as.data.frame(res), file=paste(project, "_DESeq.csv", sep=""))
}

# all(rownames(filtered_ID)) %in% colnames(PRAD_no_col)

# test_run_PRAD <- DESeq(pre_filtered_PRAD)
# res <- results(test_run_PRAD)
# res
#
# ### test with different conditions, but here we only have one condition to compare so we don't have to do this actually
# res_c <- results(test_run_PRAD, contrast=c("Sample_Type","Primary Tumor","Solid Tissue Normal"))
# res_c
#
# ### Log fold change shrinkage for visualization and ranking
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("apeglm")
#
# library("apeglm")
# resultsNames(test_run_PRAD)
# resLFC <- lfcShrink(test_run_PRAD, coef = "Sample_Type_Solid.Tissue.Normal_vs_Primary.Tumor", type = "apeglm")
# resLFC
#
# ### using 'apeglm' for LFC shrinkage. If used in published research, please cite: Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895
#
# resOrdered <- resLFC[order(resLFC$pvalue),]
# summary(resOrdered)
# # How many adjusted p-values were less than 0.05?
# sum(resOrdered$padj < 0.05, na.rm=TRUE)
# sum(resOrdered$padj < 0.001, na.rm=TRUE)
# sum(resOrdered$padj < 0.0001, na.rm=TRUE)
#
# plotMA(resLFC, ylim=c(-2,2))
#
# plotCounts(test_run_PRAD, gene=which.min(resLFC$padj), intgroup="Sample_Type")
#
# # write to csv
# resSig <- subset(resOrdered, padj < 0.0001)
# write.csv(as.data.frame(resSig),
#           file="CANCER_vs_NORMAL.csv")

