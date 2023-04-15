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
