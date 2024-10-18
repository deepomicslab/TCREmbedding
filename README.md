[![Stars](https://img.shields.io/github/stars/deepomicslab/TCREmbedding?logo=GitHub&color=yellow)](https://github.com/deepomicslab/TCREmbedding/stargazers)
[![PyPI](https://img.shields.io/pypi/v/tcrembedding?logo=PyPI)](https://pypi.org/project/tcrembedding/)
[![Downloads](https://static.pepy.tech/badge/tcrembedding)](https://pepy.tech/project/tcrembedding)
[![Docs](https://readthedocs.org/projects/tcrembedding/badge/)](https://tcrembedding.readthedocs.io)

# TCREmbedding - an all-in-one TCR embedding package in Python

TCRembedding is a composite of multiple methods for embedding amino acid sequences. Read the full [documentation] here.

### Installation Tutorial

#### 1.python venv

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
pip install tcrembedding
```

#### 2.conda

In addition to running the *env_creator.py* script to create virtual environments, you can also create and manage virtual environments via *conda*.

Example:

```
conda create --name word2vec python=3.8

conda activate word2vec

pip install -r src/TCRembedding/Word2Vec/requirements.txt
pip install tcrembedding
```

## 

## Citation

If you use `TCREmbedding` in your work, please cite the `TCREmbedding` publication as follows:

> to do


[documentation]: https://tcrembedding.readthedocs.io
