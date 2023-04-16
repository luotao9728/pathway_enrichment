import pandas as pd

## read pvalue
pvalue = pd.read_csv("pvalue.csv")
## read pathway
pathway = pd.read_csv("pathway.csv")
pathway_gene = {}
## Create dict for pathway: [gene1,gene2]
for index, row in pathway.iterrows():
    if row["Term"] not in pathway_gene.keys():
        pathway_gene[row["Term"]] = row["Genes"].split(";")
    else:
        pathway_gene[row["Term"]] += row["Genes"].split(";")
gene_pathway = {}
## Switch dict to gene: [pathway, pathway]
for keys,values in pathway_gene.items():
    for i in values:
        if i not in gene_pathway.keys():
            gene_pathway[i]=[keys]
        else:
            gene_pathway[i].append(keys)
## dict to dataframe
gene_pathway_pd = pd.DataFrame.from_dict(gene_pathway, orient="index").reset_index()
# gene_pathway_pd.to_csv('gene_pathway.csv')
gene_pathway_pd = gene_pathway_pd.rename(columns={"index":"gene_id"})
## merge pvalue and pathway on gene_id
pvalue_pathway = pd.merge(pvalue, gene_pathway_pd, how='left', on = "gene_id")
pvalue_pathway.to_csv("pvalue_pathway.csv")
