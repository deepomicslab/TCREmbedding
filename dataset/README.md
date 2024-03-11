
This directory contains all the datasets supporting the project's Figure 3-7.

## Negative data sampling

For binding prediction tasks and embedding visualization, we chose the best strategy discussed in the *TEINet project* to generate the negative samples. Refer to [sampler script](https://github.com/jiangdada1221/TEINet/blob/master/sampler.py).

## Adaptive and IMGT vgene name conversion 

The [tcrdist3](https://github.com/kmayerb/tcrdist3) project provides a lookup table and script to convert T cell receptor (TCR) V gene names between the Adaptive Biotechnologies and IMGT format. This can help convert V gene information in your dataset to the format required for TCR embedding methods.

For example, TCRBV06-01 (Adaptive) and TRBV6-1*01 (IMGT) refer to the same V gene in the TCR β chain. 

The [adaptive_imgt_mapping.csv](https://github.com/kmayerb/tcrdist3/blob/55d906b19e4c5038f5fdde843eb2edf8293efd88/tcrdist/db/adaptive_imgt_mapping.csv) maps between these two names. The [swap_gene_name.py](https://github.com/kmayerb/tcrdist3/blob/55d906b19e4c5038f5fdde843eb2edf8293efd88/tcrdist/swap_gene_name.py) can help to convert V gene names in your data.
