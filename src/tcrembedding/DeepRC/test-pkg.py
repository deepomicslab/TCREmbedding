from TCRembedding import get_embedding_instance

encoder = get_embedding_instance("EmbeddingDeepRC")
encoder.load_data("/media/lihe/TCR/project/src/TCRembedding/DeepRC/data/testdata_DeepRC.csv", use_columns="CDR3b")
encode_result = encoder.embed()
print(encode_result.shape)