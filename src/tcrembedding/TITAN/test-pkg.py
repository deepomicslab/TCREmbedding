from TCRembedding import get_embedding_instance

encoder = get_embedding_instance("EmbeddingTITAN")
encoder.load_data("data/testdata_TITAN.csv", use_columns="CDR3b")
encoder.load_model()
TCR_encode_result = encoder.embed()
epi_encode_result = encoder.embed_epi()
print(TCR_encode_result.shape)
