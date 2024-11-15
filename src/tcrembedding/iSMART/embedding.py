from iSMARTf3 import iSMART_model
import pandas as pd
import numpy as np
from joblib import Parallel, delayed

class EmbeddingiSMART:
    def __init__(self):
        self.model = None

    def load_data(self, file_path, use_columns, sep=","):
        self.data = pd.read_csv(file_path, header=0, sep=sep)[use_columns]

    def load_model(self):
        self.model = iSMART_model()

    # run pairwise alignment algorithm to analyze CDR3s
    def embed(self):
        '''
        data: Series. CDR3 sequences.
        -----
        return:
        matrix: ndarray. Each element in the array represents the pairwise distance between every two sequences.
        '''
        model = self.model
        def calculate_distance(data_chunk):
            return model.encode(data_chunk)
        num_core = -1
        chunk_size = 1000
        chunks = [self.data[i:i + chunk_size] for i in range(0, len(self.data), chunk_size)]
        encode_result = Parallel(n_jobs=num_core)(delayed(calculate_distance)(chunk) for chunk in chunks if chunk.size != 0)

        return np.vstack(encode_result)

if __name__ == "__main__":

    encoder = EmbeddingiSMART()
    encoder.load_data("data/testdata_iSMART.csv", use_columns="CDR3b")
    encoder.load_model()
    encode_result = encoder.embed()

    print(encode_result.shape)

