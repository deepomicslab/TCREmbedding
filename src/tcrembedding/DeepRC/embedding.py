import os
import pandas as pd
import torch
import numpy as np
import time
from typing import Union

from . import data_process
from .architectures import SequenceEmbeddingCNN, SequenceEmbeddingLSTM, AttentionNetwork
from .data_process import DataProcess
from .get_feature import GetFeatures

class EmbeddingDeepRC:
    def __init__(self):
        pass

    def load_data(self, file_path, use_columns):
        # read sequences and filter invalid sequences
        dataprocess = DataProcess()
        seqs = pd.read_csv(file_path, sep=",", header=0)[use_columns]
        self.data = dataprocess.filter_seqs(seqs)

    def embed(self):
        # get sequence feature (position_feature + one-hot_feature)
        getfeature = GetFeatures()
        feature = getfeature.compute_features(sequences=self.data)

        return feature

if __name__ == "__main__":

    encoder = EmbeddingDeepRC()
    encoder.read_csv("data/testdata_DeepRC.csv", use_columns="CDR3b")
    encode_result = encoder.embed()
    print(encode_result.shape)
