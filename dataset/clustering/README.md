## TCRantigenData_unique

`TCRantigenData_unique.tsv`: RAW DATA, TCR cdr3 sequences with antigen labels. Downloaded from GIANA project ([github](https://github.com/s175573/GIANA/blob/master/data/), [paper (see fig. 2)](https://www.nature.com/articles/s41467-021-25006-7)).

`TCRantigenData_unique_filt.tsv`: Processed data, filtered by script  `TCRantigenData_unique_filt.tsv`

columns in the processed data:

-- CDR3b:
The amino acid sequence of the CDR3 region of TCR beta chain recognized by the antigen.

-- Vgene:
V gene encoding the variable region of the TCR beta chain: TRBVx-x*x.

-- Antigen-epitope:
Antigenic epitope sequence recognized by T cells. Format: Antigen:xxx_epitope.

-- Antigen:
The source of antigen. Extracted from the Antigen-epitope column values.

-- Epitope_type_num:
The number of unique epitope sequences recognized by T cells specific for each antigen.

-- TCR_num:
The number of unique TCR beta chain amino acid sequences that recognize each antigen.

Here show top 10 Antigens:

| Antigen | TCR_num |  
|:-:|:-:|  
| CMV      | 22678  |
| Influenza A virus | 9279 | 
| SARS coronavirus 2| 7003|
| EBV| 4659|
|Homo sapiens| 3531|
| HBV| 3190|
|YFV| 2634|
| HIV1| 2008 |
|SARS coronavirus Tor2 |1203|
| HCV| 884|
```
