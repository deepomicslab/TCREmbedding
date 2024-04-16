from TCRembedding import get_embedding_instance
import numpy as np

encoder = get_embedding_instance("EmbeddingWord2Vec")
encoder.load_data("data/testdata_Word2Vec.csv", use_columns='CDR3b')
encode_result = encoder.embed()
encode_result = np.vstack(encode_result)
print(encode_result.shape)