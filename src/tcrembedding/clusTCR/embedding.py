import os
import pandas as pd
from typing import Union
import numpy as np

from .profiles.properties import OPTIMAL, ALL_PROPS, PROPERTY_COMBINATIONS
from .tools import make_profiles

class EmbeddingclusTCR:

    def __init__(self, max_sequence_size=None, properties: list = OPTIMAL, n_cpus: Union[str, int] = 1):
        self.properties = properties
        self.max_sequence_size = max_sequence_size
        self.n_cpus = n_cpus
        self.profiles = None

    def load_data(self, file_path, use_columns="CDR3b"):
        """Read TCR sequences from a CSV file.

        Args:
            file_path (str): The path to the CSV file containing TCR sequences.
            column_name (str): Column name of the provided file recording TCRs. Defaults to 'full_seq'.
        """
        row_data = pd.read_csv(file_path, sep="\t", header=0)
        self.data = row_data[use_columns]

    def embed(self):
        self.max_sequence_size = max([len(seq) for seq in self.data])
        matrix = make_profiles(self.data, self.properties, self.max_sequence_size, self.n_cpus)
        return matrix


if __name__ == "__main__":
    
    encoder = EmbeddingclusTCR()
    encoder.read_csv("data/testdata_clusTCR.tsv")
    encode_result = encoder.encode()
    print(encode_result.shape)