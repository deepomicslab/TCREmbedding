from lib import iSMARTf3

class embeddingismart:

    #run pairwise alignment algorithm to analyze CDR3s
    def embed(data):
        '''
        data: Series. CDR3 sequences.
        -----
        return:
        matrix: ndarray. Each element in the array represents the pairwise distance between every two sequences.
        '''
        matrix = iSMARTf3.embedding(data)

        return matrix
