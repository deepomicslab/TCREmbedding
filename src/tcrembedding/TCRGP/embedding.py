import numpy as np

from tcrgp import subsmatFromAA2, get_pcs
from tcrgp import encode_with_pc, encode_with_pc_2d, get_sequence_lists

class EmbeddingTCRGP:
    """
    This class takes the embedding operations on cdr3 sequences from TCRGP.

    reference:
        - Article: Predicting recognition between T cell receptors and epitopes with TCRGP
        - Authors: Jokinen, E., Huuhtanen, J., Mustjoki, S., Heinonen, M. & Lähdesmäki, H
        - DOI link: https://doi.org/10.1371%2Fjournal.pcbi.1008814
        - GitHub link: https://github.com/emmijokinen/TCRGP
    """

    def __init__(self, datafile='training_data/examples/vdj_human_ATDALMTGY.csv'):
        '''
        Set datafile path.
        :param datafile:
        '''
        
        self.subsmat = subsmatFromAA2('HENS920102')
        self.pc_blo = get_pcs(self.subsmat, d=21)
        self.datafile = datafile
    
    def embed(self, epi, dimension=1):
        '''
        Embed cdr3b sequences.
        :param epi: epitope name in datafile (ignored if balance_controls=False)
        :param dimension: 1 or 2, each sequence will be encoded as 1-d array(if 1) or 2-d array (if 2)
        :return:embedded cdr3b sequences, in numpy ndarray version
        '''
        
        organism = 'human'
        cdr_types = [[], ['cdr3']]
        delim = ','
        clip = [0, 0]

        # Read data file and extract requested CDRs
        epitopes, subjects, cdr_lists, lmaxes, _ = get_sequence_lists(self.datafile, organism, epi, cdr_types, delim, clip, 30,
                                                                            None, None, None, 'CDR3b', None, None,
                                                                            check_v='none', balance_controls=False ,
                                                                            alphabet_db_file_path='data/alphabeta_db.tsv')


        if dimension == 1 :
            # Every sequence is encoded as a 1-d array
            X = encode_with_pc(cdr_lists, lmaxes, self.pc_blo)
        else:
            # Every sequence is encoded as a 2-d array
            X = encode_with_pc_2d(cdr_lists, lmaxes, self.pc_blo)

        return X


if __name__ == '__main__':

    filepath = "data/testdata_TCRGP.csv"
    epitope = 'ATDALMTGY' # epitope name in datafile, ignore if balance control is False
    embedding = EmbeddingTCRGP(filepath)
    embedded_data = embedding.embed(epitope,dimension=1)
    print(embedded_data.shape)
