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
        from tcrgp import subsmatFromAA2, get_pcs
        self.subsmat = subsmatFromAA2('HENS920102')
        self.pc_blo = get_pcs(self.subsmat, d=21)
        self.datafile = datafile

    def embed(self, epi):
        '''
        Embed cdr3b sequences.
        :param epi: epitope name in datafile (ignored if balance_controls=False)
        :return:embedded cdr3b sequences, in numpy ndarray version
        '''
        from tcrgp import encode_with_pc, get_sequence_lists

        organism = 'human'
        cdr_types = [[], ['cdr3']]
        delim = ','
        clip = [0, 0]

        # Read data file and extract requested CDRs
        epitopes, subjects, cdr_lists, lmaxes, _ = get_sequence_lists(self.datafile, organism, epi, cdr_types, delim, clip, None,
                                                                            'va', 'vb', 'cdr3a', 'cdr3b', 'epitope', 'subject',
                                                                            check_v='none', balance_controls=False ,
                                                                            alphabet_db_file_path='data/alphabeta_db.tsv')

        X = encode_with_pc(cdr_lists, lmaxes, self.pc_blo)

        return X


if __name__ == '__main__':
    path = 'training_data/examples/vdj_human_ATDALMTGY.csv'
    epitope = 'ATDALMTGY' # epitope name in datafile, ignore if balance control is False
    embedding = EmbeddingTCRGP(path)
    embedded_data = embedding.embed(epitope)
    print(embedded_data.shape)
