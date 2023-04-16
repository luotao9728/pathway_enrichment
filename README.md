# Project: Analysis of Differential Gene Expression to Identify Common Gene Signatures in Multiple Cancers
## 02-604 Fundamental of Bioinformatics
## Team Name: Gene-ius
## Team Member: Tao Luo, Lu Qiao, I-Shu Wang, Yiyang Zheng 

## General work flow:
1. Data cleanup (use Python)
2. Perform statistics analysis of different genes of different cancer types (DESeq2)
3. Cluster genes based on statistical significance (from ~20,000 to ~4,000 clusters)
4. Cluster of clusters based on statistical significance (from ~4,000 to ~46 clusters)
5. Create Heat Map

## Instruction:
1. Download RNA-seq data in STAR format from TCGA database as well as their sample sheets.
2. Run data_cleanup.py
> python data_cleanup.py <data directory>
3. Run run_deseq2.R
4. Run clustering.py
> python data_cleanup.py <directory to all DESeq files> <directory to all cleanup data> <Enrichr database name>
