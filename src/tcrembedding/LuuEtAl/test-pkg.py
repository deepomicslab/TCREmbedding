from TCRembedding import get_embedding_instance

EmbeddingLuuEtAl = get_embedding_instance("EmbeddingLuuEtAl")
EmbeddingLuuEtAl.load_data('data/testdata_Luu_et_al.csv')
TCR_encode_result, epitope_encode_result = EmbeddingLuuEtAl.embed()
print(TCR_encode_result.shape)
print(epitope_encode_result.shape)