import sys
import pandas as pd
import numpy as np

class EmbeddingNetTCR2:

    def __init__(self):
        pass

    def load_data(self, file_path):
        """
        Set path of the input csv file.
        Make sure your data has at least one sequence column consisting CAPITAL alphabets only and no B, J, O, U, X or Z

        :param path: original csv file path
        """

        self.path = file_path

    def enc_list_bl_max_len(self, aa_seqs, blosum, max_seq_len):
        """
        Note:
            This function is copied from utils.py in NetTCR-2.0-main.
            See https://github.com/mnielLab/NetTCR-2.0 for original source code.

        blosum encoding of a list of amino acid sequences with padding
        to a max length

        parameters:
            - aa_seqs : list with AA sequences
            - blosum : dictionnary: key= AA, value= blosum encoding
            - max_seq_len: common length for padding
        returns:
            - enc_aa_seq : list of np.ndarrays containing padded, encoded amino acid sequences
        """

        # encode sequences:
        sequences = []
        for seq in aa_seqs:
            e_seq = np.zeros((len(seq), len(blosum["A"])))
            count = 0
            for aa in seq:
                if aa in blosum:
                    e_seq[count] = blosum[aa]
                    count += 1
                else:
                    sys.stderr.write("Unknown amino acid in peptides: " + aa + ", encoding aborted!\n")
                    sys.exit(2)

            sequences.append(e_seq)

        # pad sequences:
        # max_seq_len = max([len(x) for x in aa_seqs])
        n_seqs = len(aa_seqs)
        n_features = sequences[0].shape[1]

        enc_aa_seq = np.zeros((n_seqs, max_seq_len, n_features))
        for i in range(0, n_seqs):
            enc_aa_seq[i, :sequences[i].shape[0], :n_features] = sequences[i]

        return enc_aa_seq
    
    def embed(self, header, length=None):
        """
        Embed one column from the input csv file.
        If a column has n sequences, after embedding it becomes n * length * 20 matrix.

        :param header: header of the column for embedding in the csv file
        :param length: maximum length of a single sequence with padding
        :return: 3-dimensional matrix of amino acid sequences representation after embedding, in Pandas DataFrame type
        """

        if length is None:
            if 'CDR3b' or 'cdr' in header:
                length = 30
            elif header == 'Epitope':
                length = 20
            else:
                print("Please specific the maximum length!")
                return None

        data = pd.read_csv(self.path)

        # Note:
        #     The following statement is from utils.py in NetTCR-2.0-main.
        #     See https://github.com/mnielLab/NetTCR-2.0 for original source code.
        blosum50_20aa = {
            'A': np.array((5, -2, -1, -2, -1, -1, -1, 0, -2, -1, -2, -1, -1, -3, -1, 1, 0, -3, -2, 0)),
            'R': np.array((-2, 7, -1, -2, -4, 1, 0, -3, 0, -4, -3, 3, -2, -3, -3, -1, -1, -3, -1, -3)),
            'N': np.array((-1, -1, 7, 2, -2, 0, 0, 0, 1, -3, -4, 0, -2, -4, -2, 1, 0, -4, -2, -3)),
            'D': np.array((-2, -2, 2, 8, -4, 0, 2, -1, -1, -4, -4, -1, -4, -5, -1, 0, -1, -5, -3, -4)),
            'C': np.array((-1, -4, -2, -4, 13, -3, -3, -3, -3, -2, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1)),
            'Q': np.array((-1, 1, 0, 0, -3, 7, 2, -2, 1, -3, -2, 2, 0, -4, -1, 0, -1, -1, -1, -3)),
            'E': np.array((-1, 0, 0, 2, -3, 2, 6, -3, 0, -4, -3, 1, -2, -3, -1, -1, -1, -3, -2, -3)),
            'G': np.array((0, -3, 0, -1, -3, -2, -3, 8, -2, -4, -4, -2, -3, -4, -2, 0, -2, -3, -3, -4)),
            'H': np.array((-2, 0, 1, -1, -3, 1, 0, -2, 10, -4, -3, 0, -1, -1, -2, -1, -2, -3, 2, -4)),
            'I': np.array((-1, -4, -3, -4, -2, -3, -4, -4, -4, 5, 2, -3, 2, 0, -3, -3, -1, -3, -1, 4)),
            'L': np.array((-2, -3, -4, -4, -2, -2, -3, -4, -3, 2, 5, -3, 3, 1, -4, -3, -1, -2, -1, 1)),
            'K': np.array((-1, 3, 0, -1, -3, 2, 1, -2, 0, -3, -3, 6, -2, -4, -1, 0, -1, -3, -2, -3)),
            'M': np.array((-1, -2, -2, -4, -2, 0, -2, -3, -1, 2, 3, -2, 7, 0, -3, -2, -1, -1, 0, 1)),
            'F': np.array((-3, -3, -4, -5, -2, -4, -3, -4, -1, 0, 1, -4, 0, 8, -4, -3, -2, 1, 4, -1)),
            'P': np.array((-1, -3, -2, -1, -4, -1, -1, -2, -2, -3, -4, -1, -3, -4, 10, -1, -1, -4, -3, -3)),
            'S': np.array((1, -1, 1, 0, -1, 0, -1, 0, -1, -3, -3, 0, -2, -3, -1, 5, 2, -4, -2, -2)),
            'T': np.array((0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 2, 5, -3, -2, 0)),
            'W': np.array((-3, -3, -4, -5, -5, -1, -3, -3, -3, -3, -2, -3, -1, 1, -4, -4, -3, 15, 2, -3)),
            'Y': np.array((-2, -1, -2, -3, -3, -1, -2, -3, 2, -1, -1, -2, 0, 4, -3, -2, -2, 2, 8, -1)),
            'V': np.array((0, -3, -3, -4, -1, -3, -3, -4, -4, 4, 1, -3, 1, -1, -3, -2, 0, -3, -1, 5))
        }

        column = data[header]

        return self.enc_list_bl_max_len(column, blosum50_20aa, length)


if __name__ == '__main__':

    embedding = EmbeddingNetTCR2('data/testdata_NetTCR-2.0.csv')
    embedding_data = embedding.embed('CDR3b')
    print(embedding_data.shape)

