
This directory contains all the datasets supporting this project. We provide more detailed data sources and processing scripts in the corresponding directory.

Figure 4&5: `binding/combined_dataset_filtered.csv`

Figure 6:  `clustering/TCRantigenData_unique_filt.tsv`

Figure 7:  `embedding_visualization/data/IEDB_uniqueTCR_top10_filt.csv`

### Negative data sampling tool involved in creating the dataset:

For binding prediction tasks and embedding visualization, we chose the best strategy discussed in the *TEINet project* to generate the negative samples. Refer to [sampler script](https://github.com/jiangdada1221/TEINet/blob/master/sampler.py).
