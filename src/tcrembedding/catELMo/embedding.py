import pickle
import pandas as pd
import numpy as np
import torch
from allennlp.commands.elmo import ElmoEmbedder

class EmbeddingcatELMo:
    def __init__(self):
        pass

    def load_model(self, weights_file, options_file):
        self.weights = weights_file
        self.options = options_file
        self.embedder = ElmoEmbedder(self.options, self.weights, cuda_device=2) # cuda_device=-1 for CPU

    def load_data(self, file_path, use_columns="CDR3b"):
        """
        Read TCR sequences from a CSV file.

        Args:
            file_path (str): The path to the CSV file containing TCR sequences.

            column_name (str): Column name of the provided file recording TCRs. Defaults to 'CDR3b'.
        """
        data = pd.read_csv(file_path, sep=",", header=0)[use_columns]
        list_tokenized_sentences = []
        for seq in data.tolist():
            list_tokenized_sentences.append(list(seq))
        self.data = list_tokenized_sentences

    def catELMo_embedding(self, x):
        embed_result = []
        for seq in x:
            embed_result.append(torch.tensor(self.embedder.embed_sentence(seq)).sum(dim=0).mean(dim=0).tolist())
        return embed_result
        #return torch.tensor(self.embedder.embed_sentences(x)).sum(dim=0).mean(dim=0).tolist()

    def embed(self):
        embed_result = self.catELMo_embedding(self.data)
        embed_result = pd.DataFrame(embed_result)
        return embed_result

    def embed_epitope(self):
        embed_result = self.catELMo_embedding(self.data)
        return embed_result

if __name__ == "__main__":

    embedder = EmbeddingcatELMo()
    embedder.read_csv("data/testdata_catELMo.csv")
    tcr_embeds = embedder.embed()
    print(tcr_embeds.shape)

    # if you want to embed epitope, you can set the value of the use_columns parameter in read_csv() to the column name of the column where the epitope is located.
    # then use embed_epitope() to embed.
    embedder.read_csv("data/testdata_catELMo.csv", use_columns='Epitope')
    epi_embeds = embedder.embed_epitope()
    print(epi_embeds.shape)
