from TCRembedding import get_embedding_instance
import numpy as np
from gensim.models.callbacks import CallbackAny2Vec

class LossLogger(CallbackAny2Vec):
    """Callback to log loss after each epoch, differentiated by kmer batch."""

    def __init__(self, log_file):
        self.epoch = 0
        self.kmer_index = 0
        self.log_file = log_file
        self.loss_previous_step = 0

    def on_epoch_end(self, model):
        loss = model.get_latest_training_loss()
        loss_this_step = loss - self.loss_previous_step
        self.loss_previous_step = loss
        self.epoch += 1
        with open(self.log_file, 'a') as f:
            f.write(f'Kmer {self.kmer_index}, Epoch {self.epoch}: Loss {loss_this_step}\n')

    def reset_epoch(self, kmer_index):
        """Resets the epoch count and previous loss when a new kmer batch starts."""
        self.epoch = 0
        self.kmer_index = kmer_index
        self.loss_previous_step = 0

encoder = get_embedding_instance("EmbeddingWord2Vec")
encoder.pretrained_word2vec_model = '/media/lihe/TCR/Word2Vec/models/sequence_model_4.model' 
encoder.load_data("data/testdata_Word2Vec.csv", use_columns='CDR3b')
encode_result = encoder.embed()
encode_result = np.vstack(encode_result)
print(encode_result.shape)
