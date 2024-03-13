from gensim.models import Word2Vec

from SequenceType import SequenceType
from KmerHelper import KmerGenerator


class KmerPairModelCreator:

    def __init__(self, window=5, epochs=10, sequence_type=SequenceType.AMINO_ACID):
        self.epochs = epochs
        self.window = window
        self.sequence_type = sequence_type

    def create_model(self, sequences, k: int, vector_size: int):
        model = Word2Vec(size=vector_size, min_count=1, window=self.window)  # creates an empty model
        all_kmers = KmerGenerator.create_all_kmers(k=k, alphabet=KmerGenerator.get_sequence_alphabet())
        all_kmers = [[kmer] for kmer in all_kmers]
        model.build_vocab(all_kmers)

        for kmer in all_kmers:
            sentences = KmerGenerator.create_kmers_within_HD(kmer=kmer[0],
                                                             alphabet=KmerGenerator.get_sequence_alphabet(),
                                                             distance=1)
            model.train(sentences=sentences, total_words=len(all_kmers), epochs=model.epochs)

        return model
