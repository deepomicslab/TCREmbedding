# Embedding Methods from Luu et al



## Description

These python scripts take the embedding operations on cdr3 and epitope sequences from the methods of `Luu et al`.

To run `embedding_Luu_et_at.py`, you have to download the original source code at first and then put this file into the main folder `/TCR-Epitope-Binding-main`.

We also provide a **standalone** version in the folder `standalone version`. In this folder we also provide a csv data for testing.



## Data format

An example input csv file is as follows:

```csv
cdr3,antigen.epitope
CASSSGQLTNTEAFF,GLCTLVAML
CASSASARPEQFF,GLCTLVAML
CASSSGLLTADEQFF,GLCTLVAML
......
```

where the amino acid sequence should consist of upper case letters and not contain B, J, O, U, X and Z. 

If you want to use your own csv data, you need to make sure that your csv file is in the required format.

After embedding, each amino acid sequence of CDR3 will be mapped as a $6\times20$ matrix and epitope as a $6\times10$ matrix. They will be represented as `DataFrame` in pandas.



## Reference:
​    \- Article: Predicting TCR-epitope binding specificity using deep metric learning and multimodal learning
​    \- Authors: Luu, A. M., Leistico, J. R., Miller, T., Kim, S. & Song, J. S.
​    \- DOI link: https://doi.org/10.3390%2Fgenes12040572
​    \- GitHub link: https://github.com/jssong-lab/TCR-Epitope-Binding