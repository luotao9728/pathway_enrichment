# cleanup the original data and output as a new tsv file
import pandas as pd
import os, sys

class data_cleanup():
    def __init__(self):
        # filename should be the same as the project ID (e.g. KICH)
        # filename = sys.argv[0]
        filename = 'KICH'
        file_dict = self.get_dir_list(filename)
        df = self.combine_df(file_dict)
        df.to_csv(filename + '.tsv', sep='\t')

    def get_dir_list(self, filename):
        file_list = []
        for file in ['CANCER', 'NORMAL']:
            new_filename1 = filename + '/' + file + '/' + os.listdir(filename + '/' + file)[0]
            dir = os.listdir(new_filename1)
            for i in range(1, 10):
                new_filename2 = new_filename1 + '/' + dir[i] + '/'
                new_filename2 += os.listdir(new_filename2)[0]
                file_list.append(new_filename2)
        return file_list

    def read_file(self, filepath):
        df = pd.read_csv(filepath, sep='\t', skiprows=[0,2,3,4,5], index_col=False)
        components = filepath.split('/')
        df['sample_type'] = components[1]
        df['tumor_name'] = components[0]
        return df

    def combine_df(self, file_list):
        df = pd.DataFrame()
        for dir in file_list:
            df = pd.concat([df, self.read_file(dir)])
        return df

if __name__ == "__main__":
    data_cleanup()
