class EmbeddingNetTCR2:
    """
    This class takes the embedding operations on cdr3 and peptide sequences from NetTCR-2.0.

    reference:
        - Article: NetTCR-2.0 enables accurate prediction of TCR-peptide binding by using paired TCRα and β sequence data
        - Authors: Montemurro, A. et al
        - DOI link: https://doi.org/10.1038%2Fs42003-021-02610-3
        - GitHub link: https://github.com/mnielLab/NetTCR-2.0

    Please download the source code from reference GitHub link and directly put this file into the main directory /NetTCR-2.0-main.

    Files required:
        - /NetTCR-2.0-main/utils.pk
        - /NetTCR-2.0-main/data_processing.py
    """

    def __init__(self, path='data/train_ab_90_alphabeta.csv'):
        """
        Set path of the input csv file.
        Make sure your data has at least one sequence column consisting CAPITAL alphabets only and no B, J, O, U, X or Z

        :param path: original csv file path
        """
        self.path = path

    def embed(self, header, length=None):
        '''
        Embed one column from the input csv file.
        If a column has n sequences, after embedding it becomes n * length * 20 matrix.

        :param header: header of the column for embedding in the csv file
        :param length: maximum length of a single sequence with padding
        :return: 3-dimensional matrix of amino acid sequences representation after embedding, in Pandas DataFrame type
        '''
        import utils
        import pandas as pd

        if length is None:
            if 'CDR' or 'cdr' in header:
                length = 30
            elif header == 'peptide':
                length = 9
            else:
                print("Only CDR, peptide can be embedded.")
                return None

        data = pd.read_csv(self.path)
        encoding = utils.blosum50_20aa

        column = data[header]

        return utils.enc_list_bl_max_len(column, encoding, length)


if __name__ == '__main__':
    embedding = EmbeddingNetTCR2('data/train_ab_90_alphabeta.csv')
    embedding_data = embedding.embed('CDR3b')
    print(embedding_data.shape)
