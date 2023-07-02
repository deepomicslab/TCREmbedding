from embedding_giana import EmbeddingGIANA
import pandas as pd

row_data = pd.read_csv(f"test_data.txt", delimiter='\t')
data = row_data.iloc[0:10, 0]

vectors = EmbeddingGIANA.embed(data)
print(vectors, len(vectors))