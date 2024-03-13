import torch
import numpy as np

class GetFeatures():

    def __init__(self, n_input_features=20, add_position_info=True):
        self.n_input_features = n_input_features
        self.add_position_info = add_position_info

    def compute_features(self, sequences):
        # valid amino acids characters
        amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
        n_amino_acids = len(amino_acids)

        # get number of sequences and max length of sequence in data
        n_sequences = len(sequences)
        max_sequence_length = max(len(seq) for seq in sequences)

        # one-hot encoding
        if self.add_position_info == True:
            n_feature = n_amino_acids + 3 * self.add_position_info
        else:
            n_feature = n_amino_acids
        encoded_sequences = np.zeros((n_sequences, max_sequence_length, n_feature), dtype=np.float32)

        for i, sequence in enumerate(sequences):
            for j, amino_acid in enumerate(sequence):
                index = amino_acids.index(amino_acid)
                encoded_sequences[i, j, index] = 1

        # compute position feature
        position_added_features = encoded_sequences
        sequence_lengths = [len(seq) for seq in sequences]
        sequence_lengths = np.array(sequence_lengths)
        half_sequence_lengths = np.asarray(np.ceil(sequence_lengths / 2.), dtype=np.int32)

        for i in range(len(position_added_features)):
            sequence, seq_len, half_seq_len = position_added_features[i], sequence_lengths[i], half_sequence_lengths[i]
            sequence[:seq_len, -1] = np.abs(0.5 - np.linspace(1.0, 0, num=seq_len)) * 2.
            sequence[:half_seq_len, -3] = sequence[:half_seq_len, -1]
            sequence[half_seq_len:seq_len, -2] = sequence[half_seq_len:seq_len, -1]
            sequence[:seq_len, -1] = 1. - sequence[:seq_len, -1]

        # Perform normalization to std=1
        # position_added_features = position_added_features / position_added_features.std()

        return position_added_features