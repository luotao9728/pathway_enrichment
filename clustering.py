import pandas as pd
import os, sys
import gseapy as gp

class ClusterDESeq():
    # deseq_filepath: filepath to all DESeq files
    # data_filepath: filepath to all data after cleanup
    # db_list: list of database for clustering genes
    def __init__(self, deseq_filepath, data_filepath, db):
        self.db = db
        df_dict = self.read_deseq(deseq_filepath)
        df_pvalue = self.get_df_pvalue(df_dict)
        id_to_name = self.dict_id_to_name(data_filepath)
        df_pvalue_name = df_pvalue.rename(index=id_to_name)
        df_pvalue_name.to_csv('pvalue.csv')
        gene_name = df_pvalue_name.index.tolist()
        # Potential Database:
        # KEGG_2021_Human 320 8078
        # Human_Gene_Atlas 84 13373
        # WikiPathway_2021_Human 622 7173
        # Tissue_Protein_Expression_from_Human_Proteome_Map 30 6454
        self.query_pathway(gene_name)
        self.merge_pathway()

    def read_deseq(self, filepath):
        files = os.listdir(filepath)
        df_dict = {}
        for file in files:
            project = file[:4]
            df = pd.read_csv(filepath + '\\' + file)
            df_dict[project] = df.rename(columns={'Unnamed: 0': 'gene_id'})
        return df_dict

    def get_df_pvalue(self, df_dict):
        keys = list(df_dict.keys())
        df_pvalue = df_dict[keys[0]][['gene_id', 'pvalue']]
        df_pvalue = df_pvalue.rename(columns={'pvalue': 'BLCA'})
        for project in keys[1:]:
            df_pvalue = pd.merge(df_pvalue, df_dict[project][['gene_id', 'pvalue']], on='gene_id', how='left')
            df_pvalue = df_pvalue.rename(columns={'pvalue': project})
        df_pvalue.index = df_pvalue['gene_id']
        df_pvalue = df_pvalue.drop(df_pvalue.columns[0], axis=1)
        return df_pvalue

    def dict_id_to_name(self, filepath):
        files = os.listdir(filepath)
        id_to_name = {}
        for file in files:
            df = pd.read_csv(filepath + '\\' + file,
                             sep='\t',
                             usecols=['gene_id', 'gene_name'])
            for i in range(df.shape[0]):
                id_to_name[df.iloc[i, 0]] = df.iloc[i, 1]
        return id_to_name

    def query_pathway(self, gene_name):
        length = len(gene_name)
        n = length / 1000
        if length < 1000: n = 1
        n = int(n)
        for i in range(n):
            start = int(i * (length / n))
            end = int((i + 1) * (length / 50))
            gp.enrichr(gene_list=gene_name[start:end],
                       description='pathway',
                       gene_sets=self.db,
                       outdir=self.db + '/' + str(i))
        # combine query
        dir = []
        files = os.listdir(self.db)
        for file in files:
            curr_list = os.listdir(self.db + '/' + file)
            for item in curr_list:
                if item[-3:] == 'txt':
                    dir.append(self.db + '/' + file + '/' + item)
        self.filename = self.db + '.csv'
        with open(self.filename, 'w') as f:
            for file in dir:
                f.write(open(file, 'r').read())

    def merge_pathway(self):
        ## read pvalue
        pvalue = pd.read_csv("pvalue.csv")
        ## read pathway
        pathway = pd.read_csv(self.filename)
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
        pvalue_pathway.to_csv(self.db + "_pathway.csv")

if __name__ == "__main__":
    args = sys.argv
    ClusterDESeq(args[1], args[2], args[3])