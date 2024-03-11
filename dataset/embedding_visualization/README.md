
## File description

This directory stores the processed IEDB dataset files for embedding visualization.

`data/IEDB_uniqueTCR_top10_filt.csv` is the data source for Figure 7 and data split experiments.

`.pkl` files are the *numpy* format of the dataset `.csv` file, recording the embedding of a specific method CDR3 sequences.


## Data source and processing flow
`data/generate_dataset/IEDB_uniqueTCR_test.csv` contains 3710 unique TCR sequences specific to 122 epitope classes from IEDB dataset. 

`data/generate_dataset/IEDB_uniqueTCR_top10_pos.csv` selects a subset of the ten epitopes in the former file that binds the most TCR types (we further perform down sampling to avoid class with too big size). Then we use `samper.py` to generate the corresponding negative dataset named `data/generate_dataset/IEDB_uniqueTCR_top10_neg.csv`. Finally, we merge these two files and do quality control to obtain `data/IEDB_uniqueTCR_top10_filt.csv`.

## Usage
These datasets are used for UMAP 2D visualization of the embeddings in different scenarios and methods.
