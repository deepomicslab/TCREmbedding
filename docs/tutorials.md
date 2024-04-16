# Tutorials

## Installation

TCRembedding is a composite of multiple methods for embedding amino acid sequences.It is available on PyPI and can be downloaded and installed via pip:

```
pip install TCR-embedding==1.0.0
```

### Installation Tutorial

Since different methods rely on different runtime environments and there may be version conflicts between the dependent packages, we suggest that you create a virtual environment to use the embedding methods. At the same time, we provide an installation script *env_creator.py*, the script will be based on different embedding methods, create the corresponding virtual environment. The following is an example of how to use it:

(recommended) Based on Linux , python 3.8.

```
python env_creator.py <base_dir> <env_name> [--mirror_url=<url>]
```

base_dir : The base directory where virtual environments will be created.You also need to make sure that the corresponding *requirements.txt* file is in this directory.The *requirements.txt* file for each embedding method is available under src/TCRembedding/method_name/.

env_name : The name of the virtual environment.

mirror_url : The mirror URL for pip installations.

Example:

```
python env_creator.py /media/lihe/TCR/Word2Vec Word2vec --mirror_url=https://pypi.tuna.tsinghua.edu.cn/simple
```

The command to activate the virtual environment is printed at the end of the script run and the user can run the virtual environment according to the instructions.

Example:

```
source /media/lihe/TCR/Word2Vec/Word2vec/Word2vec_venv/bin/activate
```

After entering the virtual environment, use the pip command to install TCRembedding.

```
pip install TCRembedding
```

## Data

All the data used in the paper is publicly available, so we suggest readers refer to the original papers for more details. We also upload the processed data which can be downloaded via *******.

##  Usage Tutorial

### 1.ATM-TCR

1.1 Input file format

|   Epitope   |   CDR3B   |   Affinity   |
| :--: | :--: | :--: |
| EAAGIGILTV | CASSLGNEQF | 1 |
| EAAGIGILTV | CASSLGVATGELF | 1 |
| EAAGIGILTV | CASSQEEGGGSWGNTIYF | 1 |

​	Note : Epitope is placed in the first column, CDR3 in the second column, and Affinity in the third column. The maximum length of the CDR3 sequence does not exceed 30, and the maximum length of the Epitope sequence does not exceed 20.

1.2 All parameters in class EmbeddingATMTCR

​    file_path: The file path specifying the location of the data file to be processed.

​    blosum: BLOSUM (Blocks Substitution Matrix), a scoring matrix used for protein sequence alignment. If provided, it is used to transform sequence data. Default to using the BLOSUM45 matrix.

​    model_name: The model name, specifying the filename to use when saving or loading the model.

​    cuda: A boolean indicating whether to use CUDA acceleration (i.e., run the model on a GPU).

​    seed: The random seed used to ensure repeatability of model training and data splitting.

​    model_type: The type of model, specifying the model architecture used, for example, "attention" indicates a model using an attention mechanism.

​    drop_rate: The dropout rate used in the model's dropout layers to prevent overfitting.

​    lin_size: The size of the linear layer, specifying the number of neurons in the model's linear layers.

​    padding: The type of padding, used for the padding strategy of sequence data, for example, "mid" indicates padding in the middle.

​    heads: The number of heads in the multi-head attention mechanism, applicable only when using attention models.

​    max_len_tcr: The maximum length of TCR sequences, used to determine the truncating or padding length of sequences.

​    max_len_pep: The maximum length of peptide sequences, used to determine the truncating or padding length of sequences.

​    split_type: The type of data splitting, specifying how the data is split into training and test sets, for example, "random" indicates a random split.

1.3  Example Script

```
from TCRembedding import get_embedding_instance

EmbeddingATMTCR = get_embedding_instance("EmbeddingATMTCR")
EmbeddingATMTCR.load_data("data/testdata_ATM-TCR.csv")
encode_tcr = EmbeddingATMTCR.embed()
encode_pep = EmbeddingATMTCR.embed_epitode()
print(encode_tcr.shape)
print(encode_pep.shape)
```

### 2.catELMo

2.1 Input file format

|    CDR3b    |      Epitope       | Affinity |
| :--------: | :----------------: | :---: |
| EAAGIGILTV |     CASSLGNEQF     |   1   |
| EAAGIGILTV |   CASSLGVATGELF    |   1   |
| EAAGIGILTV | CASSQEEGGGSWGNTIYF |   1   |

2.2 All parameters in class EmbeddingcatELMo

​	None

2.3  Example Script

```
from TCRembedding import get_embedding_instance

embedder = get_embedding_instance("EmbeddingcatELMo")
embedder.load_data("data/testdata_catELMo.csv")
embedder.load_model(weights_file='model/weights.hdf5', options_file='model/options.json')
tcr_embeds = embedder.embed()
print(tcr_embeds.shape)

# if you want to embed epitope, you can set the value of the use_columns parameter in load_data() to the column name of the column where the epitope is located.
# then use embed_epitope() to embed.
embedder.load_data("data/testdata_catELMo.csv", use_columns='Epitope')
epi_embeds = embedder.embed_epitope()
print(epi_embeds.shape)
```

### 3.clusTCR

3.1 Inupt file format

|      CDR3b      |
| :-------------: |
| CAISVAGGPGETQYF |
| CASSYGGSPYEQYF  |
|  CATGTQGDQPQHF  |

3.2 All parameters in class EmbeddingclusTCR

​    max_sequence_size (Optional[int]): The maximum sequence length to consider when processing TCR sequences. If specified, it is used to define the maximum length of sequences during the encoding process. Defaults to None, which means the maximum sequence length will be determined dynamically based on the data.        

​    properties (list): A list of properties to use for encoding the TCR sequences. These properties define how TCR sequences are transformed into numerical representations. Defaults to the OPTIMAL list, which is a predefined set of optimal properties.        

​    n_cpus (Union[str, int]): The number of CPUs to use for parallel computing. Can be an integer specifying the exact number of CPUs, or a string for more flexible configurations (e.g., "auto" might denote automatic determination based on available resources). Defaults to 1, meaning processing is not parallelized.   

3.3  Example Script

```
from TCRembedding import get_embedding_instance

encoder = get_embedding_instance("EmbeddingclusTCR")
encoder.load_data("data/testdata_clusTCR.csv")
encode_result = encoder.embed()
print(encode_result.shape)
```

### 4.DeepRC

4.1 Input file format

|   CDR3b    |      Epitope       | Affinity |
| :--------: | :----------------: | :------: |
| EAAGIGILTV |     CASSLGNEQF     |    1     |
| EAAGIGILTV |   CASSLGVATGELF    |    1     |
| EAAGIGILTV | CASSQEEGGGSWGNTIYF |    1     |

4.2 All parameters in class EmbeddingDeepRC

​	None

4.3  Example Script

```
from TCRembedding import get_embedding_instance

encoder = get_embedding_instance("EmbeddingDeepRC")
encoder.load_data("data/testdata_DeepRC.csv", use_columns="CDR3b")
encode_result = encoder.embed()
print(encode_result.shape)
```

### 5.DeepTCR

5.1 Input file format

|   CDR3b    |      Epitope       | Affinity |
| :--------: | :----------------: | :------: |
| EAAGIGILTV |     CASSLGNEQF     |    1     |
| EAAGIGILTV |   CASSLGVATGELF    |    1     |
| EAAGIGILTV | CASSQEEGGGSWGNTIYF |    1     |

5.2 All parameters in class EmbeddingDeepTCR

​	None

5.3 Example Script

```
from TCRembedding import get_embedding_instance

encoder = get_embedding_instance("EmbeddingDeepTCR")
encoder.load_model(model_folder_name="Test_Model", Load_Prev_Data=True)
encoder_result = encoder.embed(encode_data_directory="Data/data")
print(encoder_result.shape)
```

​	Note: The parameter "encode_data_directory" in embed() specifies the path of the input file, but the embed() method will take all the files under the path as input files.

### 6.ERGO-II

6.1 Input file format

|   CDR3b    |      Epitope       | Affinity |
| :--------: | :----------------: | :------: |
| EAAGIGILTV |     CASSLGNEQF     |    1     |
| EAAGIGILTV |   CASSLGVATGELF    |    1     |
| EAAGIGILTV | CASSQEEGGGSWGNTIYF |    1     |

6.2 All parameters in class EmbeddingERGO

​	tcr_encoding_mode : "AE"/"LSTM"

6.3 Example Script

```
from TCRembedding import get_embedding_instance
from TCRembedding.ERGOII.data_loader import Data_Loader

encoder = get_embedding_instance("EmbeddingERGO")
encoder.tcr_encoding_model = "AE"
encoder.load_model(model_path="Models/version_10ve/checkpoints", args_path="Models/10ve/meta_tags.csv")
data_loader = Data_Loader()
tcrb_list, peptide = data_loader.collate("data/testdata_ERGO-II.csv", "AE")
tcrb_batch, pep_batch = encoder.forward(tcrb_list, peptide)
tcrb_encoding, pep_encoding = encoder.embed(tcrb_batch, pep_batch)

tcr_encode_result = tcrb_encoding.detach().numpy()
pep_encode_result = pep_encoding.detach().numpy()
print(tcrb_encoding.shape)
print(pep_encoding.shape)
```

### 7.ESM

7.1 Input file format

​	The input file format of the ESM model is .fasta, so we have to process the sequences that need to be encoded into the .fasta file format.

​	Example：

​	The contents of a csv file are as follows

|   CDR3b    |      Epitope       | Affinity |
| :--------: | :----------------: | :------: |
| EAAGIGILTV |     CASSLGNEQF     |    1     |
| EAAGIGILTV |   CASSLGVATGELF    |    1     |
| EAAGIGILTV | CASSQEEGGGSWGNTIYF |    1     |

​	First,we extract the columns that need to be coded from this csv file, excluding the table headers, and save it as a text file where each row is a tcr sequence.

```
tail -n +2 filename.csv | awk '{print $1}' > filename.tcr
```

​	Then,convert tcr text file to fasta format.

```
index=1;for i in `cat filename.tcr` ;do echo '>'$index && echo $i && let index++ ;done  > filename.fasta
```

7.2 All parameters in class EmbeddingESM

​	toks_per_batch (int): The number of tokens per batch. This parameter controls how many tokens are processed in a single batch during the embedding generation process. Default value is 4096. 

​	repr_layers (list): A list of integers indicating the layers of the model from which representations will be extracted. The layers are indexed starting from 0, with -1 representing the last layer. Default value is [-1], meaning only the last layer's representations are extracted.

​        include (list): A list of string identifiers indicating which type of sequence representations to include. Currently supports ["mean"], which computes the mean of the embeddings across the sequence length. Default is ["mean"].

​        truncation_seq_length (int): The maximum sequence length. Sequences longer than this will be truncated to this length. It's important for managing memory usage and computational efficiency. Default value is 1022.

​        nogpu (bool): A flag indicating whether to force the model to run on CPU even if a GPU is available. Setting this to True can be useful for debugging or environments without a GPU. Default value is False.

7.3 Example Script

```
from TCRembedding import get_embedding_instance

encoder = get_embedding_instance("EmbeddingESM")

model_location = "esm1b_t33_650M_UR50S"
fasta_file = "data/IEDB_uniqueTCR_top10_filter.fasta"
encoder.toks_per_batch = 2048
encoder.repr_layers = [33]
encoder.include = ["mean"]
encoder.truncation_seq_length = 1022
encoder.nogpu = True

encode_result = encoder.run(model_location, fasta_file)
print(encode_result.shape)
```

### 8.GIANA

8.1  Input file format

|   CDR3b    |      Epitope       | Affinity |
| :--------: | :----------------: | :------: |
| EAAGIGILTV |     CASSLGNEQF     |    1     |
| EAAGIGILTV |   CASSLGVATGELF    |    1     |
| EAAGIGILTV | CASSQEEGGGSWGNTIYF |    1     |

8.2 All parameters in class EmbeddingGIANA

​	None

8.3 Example Script

```
from TCRembedding import get_embedding_instance
import numpy as np

encoder = get_embedding_instance("EmbeddingGIANA")
encoder.read_csv("/media/lihe/TCR/project/src/TCRembedding/GIANA/data/testdata_GIANA.csv", use_columns="CDR3b")
encoder.load_model()
vectors = encoder.embed()
encode_result = np.vstack(vectors)
print(encode_result.shape)
```

### 9.ImRex

9.1 Input file format

|   CDR3b    |      Epitope       | Affinity |
| :--------: | :----------------: | :------: |
| EAAGIGILTV |     CASSLGNEQF     |    1     |
| EAAGIGILTV |   CASSLGVATGELF    |    1     |
| EAAGIGILTV | CASSQEEGGGSWGNTIYF |    1     |

9.2 All parameters in class EmbeddingImRex

​    cdr3_range : The minimum and maximum desired cdr3 sequence length.

​    create_neg_dataset : Whether to create negatives by shuffling/sampling, by default True.

​    Note: Should always be set to False when evaluating a dataset that already contains negatives.

9.3 Example Script

```
from TCRembedding import get_embedding_instance
import numpy as np

encoder = get_embedding_instance("EmbeddingImRex")
encoder.load_data("data/testdata_ImRex.csv") 
encode_result = encoder.embed()

iter_tf_dataset = iter(encode_result)
paired_map_list = []

for item in iter_tf_dataset:
    paired_map, affinity = item
    paired_map_list.append(paired_map.numpy())

encode_result = np.stack(paired_map_list)
print(encode_result.shape)

```

### 10.iSMART

10.1 Input file format

|   CDR3b    |      Epitope       | Affinity |
| :--------: | :----------------: | :------: |
| EAAGIGILTV |     CASSLGNEQF     |    1     |
| EAAGIGILTV |   CASSLGVATGELF    |    1     |
| EAAGIGILTV | CASSQEEGGGSWGNTIYF |    1     |

10.2 All parameters in class EmbeddingiSMART

​	None

10.3 Example Script

```
from TCRembedding import get_embedding_instance

encoder = get_embedding_instance("EmbeddingiSMART")
encoder.load_data("data/testdata_iSMART.csv", use_columns="CDR3b") 
encoder.load_model()
encode_result = encoder.embed()

print(encode_result.shape)
```

### 11.Luu et al.

11.1 Input file format

|    CDR3    |  antigen.epitope   |
| :--------: | :----------------: |
| EAAGIGILTV |     CASSLGNEQF     |
| EAAGIGILTV |   CASSLGVATGELF    |
| EAAGIGILTV | CASSQEEGGGSWGNTIYF |

11.2 All parameters in class EmbeddingLuuEtAl

​	tcr_pad_length : The padding length for TCR sequences. This defines the fixed length to which all TCR sequences will be padded or truncated, ensuring consistent input size. Default value is 30.

​	ep_pad_length : The padding length for epitope sequences. This sets the fixed length to which all epitope sequences will be padded or truncated, ensuring consistent input size. Default value is 20.

11.3 Example Script

```
from TCRembedding import get_embedding_instance

EmbeddingLuuEtAl = get_embedding_instance("EmbeddingLuuEtAl")
EmbeddingLuuEtAl.load_data('data/testdata_Luu_et_al.csv')
TCR_encode_result, epitope_encode_result = EmbeddingLuuEtAl.embed()
print(TCR_encode_result.shape)
print(epitope_encode_result.shape)
```

12.NetTCR2.0

12.1 Input file format

|   CDR3b    |      Epitope       | binder |
| :--------: | :----------------: | :----: |
| EAAGIGILTV |     CASSLGNEQF     |   1    |
| EAAGIGILTV |   CASSLGVATGELF    |   1    |
| EAAGIGILTV | CASSQEEGGGSWGNTIYF |   1    |

12.2 All parameters in class EmbeddingNetTCR2

​	None

### 12.3 Example Script

```
from TCRembedding import get_embedding_instance

embedding = get_embedding_instance("EmbeddingNetTCR2")
embedding.load_data(file_path="data/testdata_NetTCR-2.0.csv")
embedding_data = embedding.embed(header='CDR3b')
print(embedding_data.shape)

```

### 13.pMTnet

13.1 Input file format

|   CDR3b    |      Epitope       | affinity |
| :--------: | :----------------: | :------: |
| EAAGIGILTV |     CASSLGNEQF     |    1     |
| EAAGIGILTV |   CASSLGVATGELF    |    1     |
| EAAGIGILTV | CASSQEEGGGSWGNTIYF |    1     |

13.2 All parameters in class EmbeddingpMTnet

​	None

13.3 Example Script

```
from TCRembedding import get_embedding_instance

encoder = get_embedding_instance("EmbeddingpMTnet")
TCR_encoded_matrix, antigen_array = encoder.embed("data/testdata_pMTnet.csv")
print(TCR_encoded_matrix.shape)
print(antigen_array.shape)
```

### 14.SETE

14.1 Input file format

|    CDR3    |      Epitope       | affinity |
| :--------: | :----------------: | :------: |
| EAAGIGILTV |     CASSLGNEQF     |    1     |
| EAAGIGILTV |   CASSLGVATGELF    |    1     |
| EAAGIGILTV | CASSQEEGGGSWGNTIYF |    1     |

14.2 All parameters in class EmbeddingSETE

​	None

14.3 Example Script

```
from TCRembedding import get_embedding_instance

encoder = get_embedding_instance("EmbeddingSETE")
encoder.load_data("data/testdata_SETE.csv")
X, y, kmerDict = encoder.embed(k=3) # Only X is encoded.
print(X.shape)
```

### 15.TCRanno

15.1 Input file format

|   CDR3b    |      Epitope       | affinity |
| :--------: | :----------------: | :------: |
| EAAGIGILTV |     CASSLGNEQF     |    1     |
| EAAGIGILTV |   CASSLGVATGELF    |    1     |
| EAAGIGILTV | CASSQEEGGGSWGNTIYF |    1     |

15.2 All parameters in class EmbeddingTCRanno

​	None

15.3 Example Script

```
from TCRembedding import get_embedding_instance

embedder = get_embedding_instance("EmbeddingTCRanno")
embedder.load_model(model_path = None)  ## set model_path=None to use the default model (provided by TCRanno)
embedder.load_data(file_path='data/testdata_TCRanno.csv', column_name='CDR3b')

X = embedder.embed()
print(X.shape)
```

### 16.TCRGP

16.1 Input file format

|   CDR3b    |      Epitope       | affinity |
| :--------: | :----------------: | :------: |
| EAAGIGILTV |     CASSLGNEQF     |    1     |
| EAAGIGILTV |   CASSLGVATGELF    |    1     |
| EAAGIGILTV | CASSQEEGGGSWGNTIYF |    1     |

16.2 All parameters in class EmbeddingTCRGP

​	None

16.3 Example Script

```
from TCRembedding import get_embedding_instance

EmbeddingTCRGP = get_embedding_instance("EmbeddingTCRGP")
filepath = "data/testdata_TCRGP.csv"
epitope = 'ATDALMTGY' # epitope name in datafile, ignore if balance control is False
EmbeddingTCRGP.datafile = filepath
embedded_data = EmbeddingTCRGP.embed(epitope,dimension=1)
print(embedded_data.shape)
```

### 17.TITAN

17.1 Input file format

|   CDR3b    |      Epitope       | affinity |
| :--------: | :----------------: | :------: |
| EAAGIGILTV |     CASSLGNEQF     |    1     |
| EAAGIGILTV |   CASSLGVATGELF    |    1     |
| EAAGIGILTV | CASSQEEGGGSWGNTIYF |    1     |

17.2 All parameters in class EmbeddingTCRGP

​	model_path (str): The path to the directory where the model is stored. This path is used to load the model                          for further operations such as training, evaluation, or inference. Default value is 'TITAN_model', which assumes there is a folder named 'TITAN_model' in the current directory that contains the model files.

​        params_filepath (str): The path to the JSON file containing model parameters. These parameters are essential                               for initializing the model with the correct configuration before loading its weights. Default value is 'TITAN_model/model_params.json', indicating the model parameters file is located inside the 'TITAN_model' directory.    

17.3 Example Script

```
from TCRembedding import get_embedding_instance

encoder = get_embedding_instance("EmbeddingTITAN")
encoder.load_data("data/testdata_TITAN.csv", use_columns="CDR3b")
encoder.load_model()
TCR_encode_result = encoder.embed()
epi_encode_result = encoder.embed_epi()
print(TCR_encode_result.shape)
```

### 18.Word2Vec

18.1 Input file format

|   CDR3b    |      Epitope       | affinity |
| :--------: | :----------------: | :------: |
| EAAGIGILTV |     CASSLGNEQF     |    1     |
| EAAGIGILTV |   CASSLGVATGELF    |    1     |
| EAAGIGILTV | CASSQEEGGGSWGNTIYF |    1     |

18.2 All parameters in class EmbeddingWord2Vec

​	vector_size (int): The size of the vector to be learnt.

​	model_type :  The context which will be used to infer the representation of the sequence.    If :py:obj:`~immuneML.encodings.word2vec.model_creator.ModelType.ModelType.SEQUENCE` is used, the context of a k-mer is defined by the sequence it occurs in (e.g. if the sequence is CASTTY and k-mer is AST, then its context consists of k-mers CAS, STT, TTY) If :py:obj:`~immuneML.encodings.word2vec.model_creator.ModelType.ModelType.KMER_PAIR` is used, the context for the k-mer is defined as all the k-mers that within one edit distance (e.g. for k-mer CAS, the context includes CAA, CAC, CAD etc.). Valid values for this parameter are names of the ModelType enum.

​	k (int): The length of the k-mers used for the encoding.

​	epochs (int): for how many epochs to train the word2vec model for a given set of sentences (corresponding to epochs parameter in gensim package)

​	window (int): max distance between two k-mers in a sequence (same as window parameter in gensim's word2vec)

18.3 Example Script

```
from TCRembedding import get_embedding_instance
import numpy as np

encoder = get_embedding_instance("EmbeddingWord2Vec")
encoder.load_data("data/testdata_Word2Vec.csv", use_columns='CDR3b')
encode_result = encoder.embed()
encode_result = np.vstack(encode_result)
print(encode_result.shape)
```

### 19.TCRpeg

19.1 Input file format

|   CDR3b    |      Epitope       | affinity |
| :--------: | :----------------: | :------: |
| EAAGIGILTV |     CASSLGNEQF     |    1     |
| EAAGIGILTV |   CASSLGVATGELF    |    1     |
| EAAGIGILTV | CASSQEEGGGSWGNTIYF |    1     |

19.2 All parameters in class EmbeddingTCRpeg

​	hidden_size (int): The number of features in the hidden state of the model.

​	num_layers (int): The number of recurrent layers in the model.

​	embedding_path (str): The path to the embedding file to be used by the model.

​	model_path (str): The path to the pre-trained model file (.pth) to load.

​	device (str): The device to run the model on ('cpu' or 'cuda:index').

​	load_data (bool): Flag to indicate whether to load data upon initialization. Typically false when embedding is the main purpose.

19.3 Example Script

```
from TCRembedding import get_embedding_instance

encoder = get_embedding_instance("EmbeddingTCRpeg")
encoder.load_data(file_path="data/testdata_TCRpeg.csv", use_columns="CDR3b")
encode_result = encoder.embed()
print(encode_result.shape)
```

## Citation

"Citation"

## Contact

"Contact"
