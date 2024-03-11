# Embedding Methods from catELMo



## Description

These python scripts take the embedding operations on cdr3 and epitope sequences from the methods of `catELMo`.

To run `embedding.py`, you have to download the original source code.

In this folder we also provide a csv file `testdata_catELMo.csv` for testing.

In `embedding.py`, within the EmbeddingATMTCR class, `embed()` is used for embedding TCR sequences, and `embed_epitope()` is used for embedding epitope sequences. 

## Data format

An example input csv file is as follows:

```csv
#CDR3b
CASSLGNEQF
CASSLGVATGELF
CASSQEEGGGSWGNTIYF
......
```

If you want to use your own csv data, you need to make sure that your csv file is in the required format.

After embedding, each sequence is embedded into a vector of length 1024.



## Reference:

​    \- Article: Context-Aware Amino Acid Embedding Advances Analysis of TCR-Epitope Interactions

​    \- Authors: Pengfei Zhang, Seojin Bang, Michael Cai, Heewook Lee

​    \- DOI link: https://doi.org/10.7554/eLife.88837.1

​    \- GitHub link: https://github.com/Lee-CBG/catELMo
