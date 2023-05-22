from embedding import ClusTCR_embedding
import pandas as pd

row_data =pd.read_csv(f"../data/combined_dataset.csv", header=0)
data = row_data['CDR3'].iloc[:10]
enc = ClusTCR_embedding()
matrix = enc.encode(data)

print(matrix)
