from TCRembedding import get_embedding_instance

encoder = get_embedding_instance("EmbeddingclusTCR")
encoder.load_data("data/testdata_clusTCR.csv")
encode_result = encoder.embed()
print(encode_result.shape)