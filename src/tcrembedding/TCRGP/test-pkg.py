from TCRembedding import get_embedding_instance

# Dynamically get the EmbeddingTCRGP class
EmbeddingTCRGP = get_embedding_instance("EmbeddingTCRGP")
filepath = "data/testdata_TCRGP.csv"
epitope = 'ATDALMTGY' # epitope name in datafile, ignore if balance control is False
EmbeddingTCRGP.datafile = filepath
embedded_data = EmbeddingTCRGP.embed(epitope,dimension=1)
print(embedded_data.shape)