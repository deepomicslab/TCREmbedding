from TCRembedding import get_embedding_instance

EmbeddingATMTCR = get_embedding_instance("EmbeddingATMTCR")
EmbeddingATMTCR.load_data("/media/lihe/TCR/project/src/TCRembedding/ATMTCR/data/testdata_ATM-TCR.csv")
encode_tcr = EmbeddingATMTCR.embed()
encode_pep = EmbeddingATMTCR.embed_epitode()
print(encode_tcr.shape)
print(encode_pep.shape)