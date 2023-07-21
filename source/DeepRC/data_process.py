from collections import OrderedDict
import numpy as np
import pandas as pd

sequence_characters = (
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')


class DataProcess:
    def __init__(self, file_path: str = "test.csv", sequence_column: str = 'CDR3', verbose: bool = True):
        self.file_path = file_path
        self.sequence_column = sequence_column
        self.verbose = verbose

        # Define AA characters
        self.aas = sequence_characters
        self.aa_ind_dict = OrderedDict(zip(self.aas, range(len(self.aas))))
        self.n_aa = len(self.aas)

        self.seq_lens = None

    def filter_seqs(self, seqs):
        filtered_strings = []
        for seq in seqs:
            valid = True
            for char in seq:
                if char not in sequence_characters:
                    valid = False
                    break
            if valid:
                filtered_strings.append(seq)
        return filtered_strings

    def read_aa_sequence(self, filename):
        """Read sequences of repertoire file and convert to numpy int8 array"""
        try:
            seqs = pd.read_csv(filename, index_col=False, keep_default_na=False, header=0, low_memory=False)[self.sequence_column].values

            # Filter out invalid or excluded/not included sequences
            filted_seqs = self.filter_seqs(seqs)

            # Get max. sequence length
            seq_lens = np.array([len(sequence) for sequence in filted_seqs])
            max_seq_len = seq_lens.max()

            # Convert AA strings to numpy int8 array (padded with -1)
            amino_acid_sequences = np.full(shape=(len(filted_seqs), max_seq_len), dtype=np.int8, fill_value=-1)
            for i, sequence_str in enumerate(filted_seqs):
                amino_acid_sequences[i, :seq_lens[i]] = [self.aa_ind_dict[aa] for aa in sequence_str]
        except Exception as e:
            print(f"\n\n\nFailure in file {filename}\n\n\n")
            raise e

        return amino_acid_sequences
