from TCRembedding import get_embedding_instance

embedder = get_embedding_instance("EmbeddingcatELMo")
embedder.load_data("data/testdata_catELMo.csv")
embedder.load_model(weights_file='model/weights.hdf5', options_file='model/options.json')
tcr_embeds = embedder.embed()
print(tcr_embeds.shape)

# if you want to embed epitope, you can set the value of the use_columns parameter in read_csv() to the column name of the column where the epitope is located.
# then use embed_epitope() to embed.
embedder.load_data("data/testdata_catELMo.csv", use_columns='Epitope')
epi_embeds = embedder.embed_epitope()
print(epi_embeds.shape)