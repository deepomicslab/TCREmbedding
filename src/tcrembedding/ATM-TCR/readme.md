# Embedding Methods from ATM-TCR



## Description

These python scripts take the embedding operations on cdr3 and epitope sequences from the methods of `ATM-TCR`.

To run `embedding_ATM-TCR.py`, you have to download the original source code.

In this folder we also provide a csv file `testdata_ATM-TCR.csv` for testing.

In `embedding_ATM-TCR.py`, within the EmbeddingATMTCR class, `encode()` is used for embedding TCR sequences, and `encode_epitope()` is used for embedding epitope sequences. However, modifications are needed for read_pTCR_list() in `data_io_tf.py` when you want to embed epitope sequences.


## Data format

An example input csv file is as follows:

```csv
#antigen.epitope,cdr3,bound
EAAGIGILTV,CASSLGNEQF,1
EAAGIGILTV,CASSLGVATGELF,1
EAAGIGILTV,CASSQEEGGGSWGNTIYF,1
......
```

where the amino acid sequence should consist of upper case letters and not contain B, J, O, U, X and Z. 

If you want to use your own csv data, you need to make sure that your csv file is in the required format.

After embedding, each amino acid sequence of CDR3 will be mapped as a $30\times25$ matrix and epitope as a $20\times25$ matrix. They will be represented as `DataFrame` in pandas.



## Reference:

​    \- Article: ATM-TCR: TCR–epitope binding affinity prediction using a multi-head self-attention model.

​    \- Authors: Cai, M., Bang, S., Zhang, P. & Lee, H.

​    \- DOI link: https://doi.org/10.3389/fimmu.2022.893247

​    \- GitHub link: https://github.com/Lee-CBG/ATM-TCR
