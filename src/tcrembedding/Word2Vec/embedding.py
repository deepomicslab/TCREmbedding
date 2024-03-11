import os
from KmerHelper import KmerGenerator
from SequenceModelCreator import SequenceModelCreator
from KmerPairModelCreator import KmerPairModelCreator
import numpy as np
import pandas as pd

class ModelType:
    SEQUENCE = "sequence"
    KMER_PAIR = "kmer_pair"

class EmbeddingWord2Vec:
    """
    Word2VecEncoder learns the vector representations of k-mers based on the context (receptor sequence). It works for
    sequence and repertoire datasets. Similar idea was discussed in: Ostrovsky-Berman, M., Frankel, B., Polak, P. & Yaari, G.
    Immune2vec: Embedding B/T Cell Receptor Sequences in ℝN Using Natural Language Processing. Frontiers in Immunology 12, (2021).

    This encoder relies on gensim's implementation of Word2Vec and KmerHelper for k-mer extraction. Currently it works on amino acid level.


    Arguments:

        vector_size (int): The size of the vector to be learnt.

        model_type (:py:obj:`~immuneML.encodings.word2vec.model_creator.ModelType.ModelType`):  The context which will be
        used to infer the representation of the sequence.
        If :py:obj:`~immuneML.encodings.word2vec.model_creator.ModelType.ModelType.SEQUENCE` is used, the context of
        a k-mer is defined by the sequence it occurs in (e.g. if the sequence is CASTTY and k-mer is AST,
        then its context consists of k-mers CAS, STT, TTY)
        If :py:obj:`~immuneML.encodings.word2vec.model_creator.ModelType.ModelType.KMER_PAIR` is used, the context for
        the k-mer is defined as all the k-mers that within one edit distance (e.g. for k-mer CAS, the context
        includes CAA, CAC, CAD etc.).
        Valid values for this parameter are names of the ModelType enum.

        k (int): The length of the k-mers used for the encoding.

        epochs (int): for how many epochs to train the word2vec model for a given set of sentences (corresponding to epochs parameter in gensim package)

        window (int): max distance between two k-mers in a sequence (same as window parameter in gensim's word2vec)
    """
    def __init__(self, vector_size=16, model_type=ModelType.SEQUENCE, k=3, epochs=100, window=8, pretrained_word2vec_model=None):
        self.vector_size = vector_size
        self.model_type = model_type
        self.k = k
        self.epochs = epochs
        self.window = window
        self.pretrained_word2vec_model = pretrained_word2vec_model
        #self.kmer_helper = KmerHelper(k, model_type)

    def read_csv(self, file_path, use_columns):
        self.data = pd.read_csv(file_path, sep=",", header=0)[use_columns]

    def embed(self):
        if self.model_type == ModelType.SEQUENCE:
            model_creator = SequenceModelCreator()
        else:
            model_creator = KmerPairModelCreator()
        model = model_creator.create_model(self.data, self.k, self.vector_size)

        contexts = []
        for sequence in self.data:
            if self.model_type == ModelType.SEQUENCE:
                context = KmerGenerator.create_kmers_from_sequence(sequence, self.k)
            elif self.model_type == ModelType.KMER_PAIR:
                kmers = KmerGenerator.create_kmers_from_sequence(sequence, self.k)
                for kmer in kmers:
                    context = KmerGenerator.create_kmers_within_HD(kmer, KmerGenerator.get_sequence_alphabet())
            else:
                raise ValueError("Invalid model_type")
            contexts.append(context)


        if self.model_type == ModelType.SEQUENCE:
            encoded_sequences = []
            for context in contexts:
                sequence_vector = np.zeros(self.vector_size)
                for kmer in context:
                    if kmer in model.wv:
                        try:
                            word_vector = model.wv[kmer]
                            sequence_vector = np.add(sequence_vector, word_vector)
                        except KeyError:
                            pass
                encoded_sequences.append(sequence_vector)

            return encoded_sequences

        elif self.model_type == ModelType.KMER_PAIR:
            encoded_sequences = []
            for context in contexts:
                sequence_vector = np.zeros(self.vector_size)
                for kmers in context:
                    for kmer in kmers:
                        if kmer in model.wv:
                            try:
                                word_vector = model.wv[kmer]
                                sequence_vector = np.add(sequence_vector, word_vector)
                            except KeyError:
                                pass
                encoded_sequences.append(sequence_vector)

            return encoded_sequences

if __name__ == "__main__":

    encoder = EmbeddingWord2Vec()
    encoder.read_csv(f"D:/TCR/TCRantigenData_detailed1.csv", use_columns='CDR3b')
    encode_result = encoder.embed()
    encode_result = np.vstack(encode_result)
    print(encode_result.shape)

    
