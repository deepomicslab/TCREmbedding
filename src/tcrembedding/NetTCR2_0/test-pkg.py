from TCRembedding import get_embedding_instance

embedding = get_embedding_instance("EmbeddingNetTCR2")
embedding.load_data(file_path="/media/lihe/TCR/project/src/TCRembedding/NetTCR2_0/data/testdata_NetTCR-2.0.csv")
embedding_data = embedding.embed('CDR3b')
print(embedding_data.shape)
