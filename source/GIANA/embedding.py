import numpy as np

from giana import GIANA_model
import pandas as pd

class EmbeddingGIANA:
    def __init__(self):
        self.model = None

    def read_csv(self, file_path, use_columns):
        self.data = pd.read_csv(file_path, header=0, sep="\t")[use_columns]

    def load_model(self):
        self.model = GIANA_model()

    # run pairwise alignment algorithm to analyze CDR3s
    def encode(self):
        '''
        data: Series. CDR3 sequences.
        ST: int. Starting position of CDR3 sequence. The first ST letters are omitted. CDR3 sequence length L must be >= ST + 7
        verbose: Bool. if True, it will print intermediate messages.
        -----
        return:
        vectors: list. embedding result, a list of vectors. Each vector in the list corresponds to the embedding result of the respective CDR3 sequence.
        '''
        model = self.model
        vectors = model.EncodeRepertoire(self.data)

        return vectors


if __name__ == "__main__":
    encoder = EmbeddingGIANA()
    encoder.read_csv("D:/TCR/TCRantigenData_unique_filt.tsv", use_columns="CDR3b")
    encoder.load_model()
    vectors = encoder.encode()
    encode_result = np.vstack(vectors)
    print(encode_result)
    print(type(encode_result))
    print(encode_result.shape)

    np.save("D:/TCR/cluster/GIANA_tcr.npy", encode_result)
