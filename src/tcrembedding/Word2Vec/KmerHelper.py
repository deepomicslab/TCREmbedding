import itertools
import warnings
from typing import List
from SequenceType import SequenceType

class KmerGenerator:
    @staticmethod
    def generate_kmers(sequence: str, k: int, overlap: bool = True) -> List[str]:
        kmers = []
        step = 1 if overlap else k
        for i in range(0, len(sequence) - k + 1, step):
            kmer = sequence[i:i + k]
            kmers.append(kmer)
        return kmers

    @staticmethod
    def create_all_kmers(k: int, alphabet: list):
        """
        creates all possible k-mers given a k-mer length and an alphabet
        :param k: length of k-mer (int)
        :param alphabet: list of characters from which to make all possible k-mers (list)
        :return: alphabetically sorted list of k-mers
        """
        kmers = [''.join(x) for x in itertools.product(alphabet, repeat=k)]
        kmers.sort()
        return kmers

    @staticmethod
    def create_kmers_from_sequence(sequence: List[str], k: int, overlap: bool = True):
        kmers = []
        step = 1 if overlap else k
        for i in range(0, len(sequence) - k + 1, step):
            kmers.append(sequence[i:i + k])
        return kmers

    @staticmethod
    def get_sequence_alphabet(sequence_type: SequenceType = None):
        """
        :return: alphabetically sorted receptor_sequence alphabet
        """
        seq_type = sequence_type if sequence_type is not None else SequenceType.AMINO_ACID
        if seq_type == SequenceType.AMINO_ACID:
           alphabet = list("ACDEFGHIKLMNPQRSTVWY")
           alphabet.sort()
        elif seq_type == SequenceType.NUCLEOTIDE:
            alphabet = list("ACGT")
            alphabet.sort()
        else:
            raise RuntimeError(
                "EnvironmentSettings: the sequence alphabet cannot be obtained if sequence_type was not set properly. "
                f"Expected AMINO_ACID or NUCLEOTIDE, but got {seq_type} instead.")

        return alphabet

    @staticmethod
    def create_kmers_within_HD(kmer: str, alphabet: list, distance: int = 1):

        assert distance < len(kmer)

        if distance > 1:
            warnings.warn("In create_kmers_within_HD distance larger than 1 is not yet implemented. "
                          "Using default value 1...", Warning)

        pairs = []

        for i in range(len(kmer)):
            for letter in alphabet:
                new_kmer = kmer[0:i] + letter + kmer[i + 1:]
                pairs.append([kmer, new_kmer])

        return pairs
