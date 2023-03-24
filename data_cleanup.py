# cleanup the original data and output as a new tsv file
# import os, sys
import pandas as pd

def read_file(filepath):
    df = pd.read_csv(filepath, sep='\t', skiprows=[0,2,3,4,5])
    components = filepath.split('/')
    df['sample_type'] = components[1]
    df['tumor_name'] = components[0]
    return df

df = read_file('KICH/CANCER/gdc_download_20230301_212055.202248/0ba21ef5-0829-422e-a674-d3817498c333/4868e8fc-e045-475a-a81d-ef43eabb7066.rna_seq.augmented_star_gene_counts.tsv')
print(df)
