# IEDB TCR-ept dataset for reference hitting experiment

## Dataset Description
`IEDB_uniqueTCR_test.csv` contains 3710 unique TCR sequences specific to 122 epitope classes. We filtered classes with too big or too small sizes and selected epitope classes of size between 30 and 500 (i.e.. 30-500 TCR sequences specific to the epitope; resulting in 7 epitope classes)

`IEDB_uniqueTCR_test.pkl` is the dict format of the former csv file, with TCR CDR3 sequences as keys and corresponding epitopes as values.

`IEDB_reference.csv` is the TCR CDR3 sequence library after data filtering. 

`IEDB_reference.pkl` records the epitopes corresponding to these CDR3 sequences (dict format).

(Experimental reference TCRanno project supp.fig s1 and s5)

## Usage

1. The dataset can be used to train models to perform umap/t-SNE 2D visualization of the embeddings.
2. We utilize different embedding methods in the TCRanno framework to find the top k candidate TCRs with the highest similarity scores from the reference database (IEDB_reference.csv). We then determine whether the target epitope of the query TCR sequence exists in the target epitopes of the hit TCRs. The hit rate is used as the evaluation metric. (Reference TCRanno supp.fig s1)