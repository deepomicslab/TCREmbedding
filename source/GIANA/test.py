import embedding
import pandas as pd

row_data = pd.read_csv(f"input.txt", delimiter='\t')
data = row_data.iloc[0:10, 0]

vectors = embedding.embedding(data)
print(vectors, len(vectors))