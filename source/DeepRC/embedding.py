from typing import Union
import pandas as pd
import torch
import numpy as np

import data_process
from architectures import SequenceEmbeddingCNN, SequenceEmbeddingLSTM, AttentionNetwork
from data_process import DataProcess
from get_feature import GetFeatures

class EmbeddingDeepRC:
    def __init__(self):
        pass

    def read_csv(self, file_path, use_columns):
        # read sequences and filter invalid sequences
        dataprocess = DataProcess()
        seqs = pd.read_csv(file_path, sep="\t", header=0)[use_columns]
        self.data = dataprocess.filter_seqs(seqs)

    def encode(self, type=Union["LSTM", 'CNN']):
        # get sequence feature (position_feature + one-hot_feature)
        getfeature = GetFeatures()
        feature = getfeature.compute_features(sequences=self.data)
        feature = torch.tensor(feature)
        #print(feature)

        # Get sequence embedding h() for single bag (shape: (n_sequences_per_bag, d_v))
        if type == "CNN":
            cnn_embed = SequenceEmbeddingCNN(n_input_features=len(data_process.sequence_characters) + 3)
            embed_result = cnn_embed.forward(feature)
            #print(cnn_embed_result)
        elif type == "LSTM":
            # LSTM embedding
            lstm_embed = SequenceEmbeddingLSTM(n_input_features=len(data_process.sequence_characters) + 3)
            embed_result = lstm_embed.forward(feature, sequence_lengths=18)
            #print(lstm_embed_result)

        # Get attention weights for single bag (shape: (n_sequences_per_bag, 1))
        attention = AttentionNetwork(n_input_features=32)
        attention_weights = attention.forward(embed_result)

        # Calculate attention activations (softmax over n_sequences_per_bag) (shape: (n_sequences_per_bag, 1))
        attention_weights = torch.softmax(attention_weights, dim=0)

        # Apply attention weights to sequence features (shape: (n_sequences_per_bag, d_v))
        emb_seqs_after_attention = embed_result * attention_weights

        return emb_seqs_after_attention


if __name__ == "__main__":
    encoder = EmbeddingDeepRC()
    encoder.read_csv(f"D:/TCR/TCRantigenData_top5.tsv", use_columns="CDR3b")
    encode_result = encoder.encode(type="CNN")
    encode_result = encode_result.detach().numpy()
    print(encode_result)
    print(encode_result.shape)
    np.save("D:/TCR/cluster/DeepRC_tcr_top5.npy", encode_result)