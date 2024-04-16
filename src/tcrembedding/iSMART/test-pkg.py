from TCRembedding import get_embedding_instance

encoder = get_embedding_instance("EmbeddingiSMART")
encoder.load_data("data/testdata_iSMART.csv", use_columns="CDR3b")
encoder.load_model()
encode_result = encoder.embed()

print(encode_result.shape)
