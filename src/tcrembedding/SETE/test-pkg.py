from TCRembedding import get_embedding_instance

encoder = get_embedding_instance("EmbeddingSETE")
encoder.load_data("data/testdata_SETE.csv")
X, y, kmerDict = encoder.embed(k=3) # Only X is encoded.
print(X.shape)