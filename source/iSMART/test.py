import ismart
import pandas as pd

row_data = pd.read_csv(f"lib/input.txt", delimiter='\t')
data = row_data.iloc[0:10, 0]

matrix = ismart.iSMART(data)
matrix_df = pd.DataFrame(matrix)

print(matrix)