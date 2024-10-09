
## Dataset for Retraining TCR CDR3 Embedding Models
This directory contains the datasets used for retraining the TCR CDR3 embedding models within the TCRembedding project.

### Pretrain Dataset

The `pretrain_dataset_cdr3.zip` file is specifically prepared for retraining data-driven type TCR CDR3 embedding models used in the TCRembedding project. To avoid potential bias and information leakage in downstream experiments, we carefully excluded any CDR3 sequences involved in binding, clustering, and embedding visualization datasets. The curated dataset contains 4,083,694 unique CDR3 amino acid sequences.

#### Source

The dataset was initially sourced from the catELMo project, which provides a comprehensive collection of human repertoires from seven distinct projects. For more detailed information about the original data, please visit the [catELMo project repository](https://github.com/Lee-CBG/catELMo/blob/main/datasets/catELMo.zip) and publication.

- Zhang Pengfei, Bang Seojin, Cai Michael, Lee Heewook (2023). Context-Aware Amino Acid Embedding Advances Analysis of TCR-Epitope Interactions. _eLife_, 12:RP88837. DOI: [10.7554/eLife.88837.2](https://doi.org/10.7554/eLife.88837.2)

### TCRanno Model Retraining Dataset

`TCRanno_retrain_trainset.pkl` and `TCRanno_retrain_validset.pkl` contain epitope-label CDR3 sequences without overlapping the test dataset used in clustering experiments. We retrain the model TCRanno with the hold-out trainset to avoid data leakage in the clustering test.
