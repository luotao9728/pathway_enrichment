# cleanup the original data and output as a new tsv file
import pandas as pd
import os, sys

def merge_sample_sheet(filepath):
    files = os.listdir(filepath)
    df = pd.DataFrame()
    for file in files:
        df = pd.concat([df, pd.read_csv(filepath + '\\' + file, sep='\t')])
    df.to_csv('merged_sample_sheet.tsv', sep='\t', index=False)
    return df

class data_cleanup():
    def __init__(self):
        # filename should be the same as the project ID (e.g. KICH)
        filename = sys.argv[1]
        file_list = self.get_dir_list(filename)
        df = self.combine_df(file_list)
        df = self.rename_col(df)
        df.to_csv(filename + '.tsv', sep='\t', index=False)

#     def get_dir_list(self, filename):
#         file_list = []
#         for file in ['NORMAL', 'CANCER']:
#             new_filename1 = filename + '/' + file + '/' + os.listdir(filename + '/' + file)[0]
#             dir = os.listdir(new_filename1)
#             for i in range(len(dir)):
#                 new_filename2 = new_filename1 + '/' + dir[i] + '/'
#                 new_filename2 += os.listdir(new_filename2)[0]
#                 if new_filename2[-3:] == 'tsv':
#                     file_list.append(new_filename2)
#         return file_list

#     def read_file(self, file):
#         df = pd.read_csv(file,
#                          sep='\t',
#                          skiprows=[0,2,3,4,5],
#                          usecols=['gene_id', 'gene_name', 'unstranded'])
#         file_name = file.split('/')[-1]
#         df = df.rename(columns={'unstranded': file_name})
#         return df

#     def combine_df(self, file_list):
#         df = self.read_file(file_list[0])
#         for dir in file_list[1:]:
#             df = pd.merge(df, self.read_file(dir), on=['gene_id', 'gene_name'], how='left')
#         return df

    def rename_col(self, df):
        sample_sheet = pd.read_csv('merged_sample_sheet.tsv',
                                   sep='\t',
                                   usecols=['File Name', 'Sample ID'])
        file_names = df.columns[2:]
        rename_dict = {}
        for file_name in file_names:
            match = list(sample_sheet.loc[sample_sheet['File Name'] == file_name]['Sample ID'])
            if len(match) != 0:
                rename_dict[file_name] = match[0]
            else:
                print('File Name:', file_name, 'does not exist!')
        df = df.rename(columns=rename_dict)
        return df

if __name__ == "__main__":
    merge_sample_sheet('sample_sheet')
    data_cleanup()
