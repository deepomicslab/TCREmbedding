from TCRembedding import get_embedding_instance

# Load model and data
embedder = get_embedding_instance("EmbeddingTCRanno")
embedder.load_model(model_path = None)  ## set model_path=None to use the default model (provided by TCRanno)
embedder.load_data(file_path='data/testdata_TCRanno.csv', column_name='CDR3b')

# Get embeddings
X = embedder.embed()
print(X.shape)
