from TCRembedding import get_embedding_instance

encoder = get_embedding_instance("EmbeddingTCRpeg")
encoder.load_data(file_path="data/testdata_TCRpeg.csv", use_columns="CDR3b")
encode_result = encoder.embed()
print(encode_result.shape)