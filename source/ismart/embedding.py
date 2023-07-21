from iSMARTf3 import iSMART_model
import pandas as pd
import numpy as np

class EmbeddingiSMART:
    def __init__(self):
        self.model = None

    def read_csv(self, file_path, use_columns):
        self.data = pd.read_csv(file_path, header=0, sep="\t")[use_columns].head(1000)

    def load_model(self):
        self.model = iSMART_model()

    # run pairwise alignment algorithm to analyze CDR3s
    def encode(self):
        '''
        data: Series. CDR3 sequences.
        -----
        return:
        matrix: ndarray. Each element in the array represents the pairwise distance between every two sequences.
        '''
        model = self.model
        encode_result = model.encode(self.data)

        return encode_result

if __name__ == "__main__":
    encoder = EmbeddingiSMART()
    encoder.read_csv("D:/TCR/TCRantigenData_unique_filt.tsv", use_columns="CDR3b")
    encoder.load_model()
    encode_result = encoder.encode()

    print(encode_result)
    #np.save("iSMART_tcr", encode_result)