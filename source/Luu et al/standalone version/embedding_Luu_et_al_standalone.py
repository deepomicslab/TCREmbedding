class EmbeddingLuuEtAlStandalone:
    """
    This class takes the embedding operations on cdr3 and epitope sequences from the methods of Luu et al.

    This standalone version does not need to download the source code from GitHub.

    Reference:
        - Article: Predicting TCR-epitope binding specificity using deep metric learning and multimodal learning
        - Authors: Luu, A. M., Leistico, J. R., Miller, T., Kim, S. & Song, J. S.
        - DOI link: https://doi.org/10.3390%2Fgenes12040572
        - GitHub link: https://github.com/jssong-lab/TCR-Epitope-Binding

    """

    df_data = None
    TCR_PAD_LENGTH = 20
    EP_PAD_LENGTH = 10
    NCH = 6

    def __init__(self, file_path, tcr_pad_length=20, ep_pad_length=10):
        '''
        Read cdr3 and epitope data from csv files.

        Csv file example:
            cdr3,antigen.epitope
            CASSSGQLTNTEAFF,GLCTLVAML
            CASSASARPEQFF,GLCTLVAML
            CASSSGLLTADEQFF,GLCTLVAML

        :param file_path: Path of the csv file
        '''
        import pandas as pd
        df = pd.read_csv(file_path)
        df = df[(df['antigen.epitope'].str.match('^[A-Z]{1,10}$')) &
                (~df['antigen.epitope'].str.contains('B|J|O|U|X|Z')) &
                (df['cdr3'].str.match('^[A-Z]{1,20}$')) &
                (~df['cdr3'].str.contains('B|J|O|U|X|Z'))]
        self.df_data = df
        self.TCR_PAD_LENGTH = tcr_pad_length
        self.EP_PAD_LENGTH = ep_pad_length

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
        import numpy as np
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
        import numpy as np
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
        import pickle as pk

        aa_vec = pk.load(open('atchley.pk', 'rb'))

        X = self.encode_seq_array(self.df_data['cdr3'], aa_vec, True, self.TCR_PAD_LENGTH)
        y = self.encode_seq_array(self.df_data['antigen.epitope'], aa_vec, True, self.EP_PAD_LENGTH)

        return X, y


def csv_modifier(input_path='../dataset/merged_data/combined_dataset.csv', output_path='../dataset/testdata_LuuEtAl.csv', rows=200):
    '''
    Modify /dataset/merged_data/combined_dataset.csv to Luu et al data format.

    :param input_path: input csv file path
    :param output_path: output csv file path
    :param rows: first n rows of data from origin csv file
    :return: True if successful
    '''

    import pandas as pd

    df = pd.read_csv(input_path, nrows=rows)

    # Feel free to modify the following statements if you want to use your own data.
    df = df.drop(columns=['Affinity'])
    df = df.rename(columns={'Epitope': 'antigen.epitope'})
    df = df.rename(columns={'CDR3': 'cdr3'})
    df = df.reindex(columns=['cdr3', 'antigen.epitope'])

    df = df[(df['antigen.epitope'].str.match('^[A-Z]{1,10}$')) &
            (~df['antigen.epitope'].str.contains('B|J|O|U|X|Z')) &
            (df['cdr3'].str.match('^[A-Z]{1,20}$')) &
            (~df['cdr3'].str.contains('B|J|O|U|X|Z'))]

    df.to_csv(path_or_buf=output_path, index=False)

    return True


if __name__ == '__main__':
    # if you want to use your own data
    # csv_modifier(rows=200)
    file_path = 'testdata_Luu_et_al.csv'
    embedding = EmbeddingLuuEtAlStandalone(file_path)
    X, y = embedding.embed()
    print(X.shape)
    print(y.shape)