import giana

class EmbeddingGIANA:
    # Encode CDR3 sequences into 96 dimensional space,a CDR3 sequence is converted into a 96-dimensional vector.
    def __init__():
    
    def embed(data, ST=3, verbose=True):

        '''
        data: Series. CDR3 sequences.
        ST: int. Starting position of CDR3 sequence. The first ST letters are omitted. CDR3 sequence length L must be >= ST + 7
        verbose: Bool. if True, it will print intermediate messages.
        -----
        return:
        vectors: list. embedding result, a list of vectors. Each vector in the list corresponds to the embedding result of the respective CDR3 sequence.
        '''
        vectors = giana.EncodeRepertoire(data, ST, verbose)

        return vectors
