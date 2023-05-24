# Embedding Methods from NetTCR-2.0



## Description

These python scripts take the embedding operations on cdr3 and epitope sequences from the methods of `NetTCR-2.0`.

To run `embedding_NetTCR2.0.py`, you have to download the original source code at first and then put this file into the main folder `/NetTCR-2.0-main`.

We also provide a **standalone** version in the folder `standalone version`. In this folder we also provide a csv data for testing.



## Data format

An example input csv file is as follows:

```csv
CDR3b,peptide,binder
ASSVDGSTYNEQF,GILGFVFTL,1
ASSLTGTPQETQY,GILGFVFTL,1
QQRVRAGISSYEQY,GILGFVFTL,1
......
```

where the amino acid sequence should consist of CAPITAL alphabets and not contain B, J, O, U, X and Z. 

If you want to use your own csv data, you need to make sure that your csv file is in the required format.

After embedding, each amino acid sequence  will be mapped as a $n\times20$ matrix, where $n$ is the maximum length of a single sequence with padding. They will be represented in `Pandas ` `DataFrame`  type.



## Reference:

\- Article: NetTCR-2.0 enables accurate prediction of TCR-peptide binding by using paired TCRα and β sequence data
\- Authors: Montemurro, A. et al
\- DOI link: https://doi.org/10.1038%2Fs42003-021-02610-3
\- GitHub link: https://github.com/mnielLab/NetTCR-2.0