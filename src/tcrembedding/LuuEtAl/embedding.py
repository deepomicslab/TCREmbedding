import os
import pandas as pd
import numpy as np
import pickle as pk

class EmbeddingLuuEtAl:

    df_data = None
    TCR_PAD_LENGTH = 30
    EP_PAD_LENGTH = 20
    NCH = 6

    def __init__(self, tcr_pad_length=30, ep_pad_length=20):
        
        self.TCR_PAD_LENGTH = tcr_pad_length
        self.EP_PAD_LENGTH = ep_pad_length

    def load_data(self, file_path):
        '''
        Read cdr3 and epitope data from csv files.

        Csv file example:
            cdr3,antigen.epitope
            CASSSGQLTNTEAFF,GLCTLVAML
            CASSASARPEQFF,GLCTLVAML
            CASSSGLLTADEQFF,GLCTLVAML

        :param file_path: Path of the csv file
        '''
        df = pd.read_csv(file_path)
        # df = df[(df['antigen.epitope'].str.match('^[A-Z]{1,10}$')) &
        #         (~df['antigen.epitope'].str.contains('B|J|O|U|X|Z')) &
        #         (df['cdr3'].str.match('^[A-Z]{1,20}$')) &
        #         (~df['cdr3'].str.contains('B|J|O|U|X|Z'))]
        self.df_data = df

    def pad_seq(self, s, length):
        '''
        Note:
            This function is copied from data_processing.py in the original project.
            See https://github.com/jssong-lab/TCR-Epitope-Binding for further details.
        '''
        return s + ' ' * (length - len(s))

    def encode_seq(self, s, aa_vec):
        '''
        Note:
            This function is copied from data_processing.py in the original project.
            See https://github.com/jssong-lab/TCR-Epitope-Binding for further details.
        '''
        
        s_enc = np.empty((len(s), self.NCH), dtype=np.float32)
        for i, c in enumerate(s):
            s_enc[i] = aa_vec[c]
        return s_enc

    def encode_seq_array(self, arr, aa_vec, pad=True, pad_length=TCR_PAD_LENGTH):
        '''
        Note:
            This function is copied from data_processing.py in the original project.
            See https://github.com/jssong-lab/TCR-Epitope-Binding for source code.
        '''
        if pad:
            arr = arr.map(lambda x: self.pad_seq(x, pad_length))
        enc_arr = arr.map(lambda x: self.encode_seq(x, aa_vec))
        enc_tensor = np.empty((len(arr), pad_length, self.NCH))
        for i, mat in enumerate(enc_arr):
            enc_tensor[i] = mat
        return enc_tensor

    def embed(self):
        '''
        Embed CDR3 and epitope chain.

        :return:
          - X: embedded CDR3 pandas data frame (n, 20, 6) and
          - y: embedded epitope pandas data frame (n, 10, 6)
        '''

        aa_vec = pk.load(open('atchley.pk', 'rb'))

        X = self.encode_seq_array(self.df_data['cdr3'], aa_vec, True, self.TCR_PAD_LENGTH)
        y = self.encode_seq_array(self.df_data['antigen.epitope'], aa_vec, True, self.EP_PAD_LENGTH)

        return X , y


def csv_modifier(input_path='../dataset/merged_data/combined_dataset.csv', output_path='../dataset/testdata_LuuEtAl.csv', rows=200):
    '''
    Modify /dataset/merged_data/combined_dataset.csv to Luu et al data format.

    :param input_path: input csv file path
    :param output_path: output csv file path
    :param rows: first n rows of data from origin csv file
    :return: True if successful
    '''

    df = pd.read_csv(input_path, nrows=rows)

    # Feel free to modify the following statements if you want to use your own data.
    # df = df.drop(columns=['Affinity'])
    # df = df.rename(columns={'Epitope': 'antigen.epitope'})
    df = df.rename(columns={'CDR3b': 'cdr3'})
    df.reindex(columns=['cdr3'])
    # df = df.reindex(columns=['cdr3', 'antigen.epitope'])

    df = df[(df['cdr3'].str.match('^[A-Z]{1,30}$')) &
            (~df['cdr3'].str.contains('B|J|O|U|X|Z'))]

    df.to_csv(path_or_buf=output_path, index=False)

    return True

def get_file_paths(folder_path):
    file_paths_dict = {}
    
    # Iterate over each file in the folder
    for file in os.listdir(folder_path):
        if file.endswith('.csv'):
            # Construct the file path
            file_path = os.path.abspath(os.path.join(folder_path, file))
            # Add to the dictionary: key is file name, value is file path
            file_paths_dict[file] = file_path

    return file_paths_dict

if __name__ == '__main__':
    
    # if you want to use your own datasets and the format is not , you can use csv_modifier() to transefer
    # csv_modifier(file_path, './csv_modified_result/TCRantigenData_detailed1.csv', rows=9033)
    embedding = EmbeddingLuuEtAl('data/testdata_Luu_et_al.csv')
    TCR_encode_result, epitope_encode_result = embedding.embed()
    print(TCR_encode_result.shape)
    print(epitope_encode_result.shape)
