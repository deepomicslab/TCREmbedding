import numpy as np
import torch
from tcr2vec.dataset import TCRLabeledDset
from torch.utils.data import DataLoader
from tcr2vec.utils import get_emb
from tcr2vec.model import TCR2vec

class EmbeddingTCR2vec:
    """A class for getting TCR embedding from TCR sequences.

    This class provides an interface to load pre-trained TCR2vec model and get TCR embedding from raw TCR amino acid sequences.

    Attributes:
        tcr2vec_model: The TCR2vec model instance. 
    """

    def __init__(self, model_path, data_path='./data/sample.csv', device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")):

        """Load the pre-trained TCR2vec model and dataset.

        Args:
            model_path (str): The path to the pre-trained TCR2vec model. 
            data_path (str): The path to the CSV file containing TCR sequences.
            
        """
        self.tcr2vec_model = TCR2vec(model_path).to(device)
        self.tcr2vec_model.eval()
        self.data_path = data_path
        

    def encode(self, batch_size = 512, column_name='CDR3.beta'):
        """Get the embeddings for the loaded TCR sequences.
        Args:
            column_name (str): Column name of the provided file recording TCRs.

        Returns:
            np.ndarray: The embeddings for the TCR sequences, with shape [n, embed_size].
        """
        self.dset = TCRLabeledDset(self.data_path, only_tcr=True, use_column=column_name)

        loader = DataLoader(self.dset, batch_size=batch_size, collate_fn=self.dset.collate_fn, shuffle=False)
        emb = get_emb(self.tcr2vec_model, loader)
        return emb


if __name__ == '__main__':
    # Load model and data
    file_path = './data/sample.csv'

    ## pretrained TCR2vec models: 1. CDR3vec_120 (pretrained on CDR3 sequences) 2.TCR2vec_120 (pretrained on TCRb FULL sequences)
    tcr2vec = EmbeddingTCR2vec(model_path = './pretrained_models/CDR3vec_120', data_path='./data/sample.csv' )

    # Get embeddings
    X = tcr2vec.encode(column_name='CDR3.beta')
    print(X.shape)
