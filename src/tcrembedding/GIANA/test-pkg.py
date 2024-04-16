from TCRembedding import get_embedding_instance
import numpy as np

encoder = get_embedding_instance("EmbeddingGIANA")
encoder.read_csv("/media/lihe/TCR/project/src/TCRembedding/GIANA/data/testdata_GIANA.csv", use_columns="CDR3b")
encoder.load_model()
vectors = encoder.embed()
encode_result = np.vstack(vectors)
print(encode_result.shape)