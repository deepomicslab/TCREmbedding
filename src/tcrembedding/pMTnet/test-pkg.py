from TCRembedding import get_embedding_instance

encoder = get_embedding_instance("EmbeddingpMTnet")
TCR_encoded_matrix, antigen_array = encoder.embed("data/testdata_pMTnet.csv")
print(TCR_encoded_matrix.shape)
print(antigen_array.shape)