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

    def __init__(self):
        self.tcr2vec_model = None

    def load_model(self, model_path, device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")):
        """Load the pre-trained TCR2vec model.

        Args:
            model_path (str): The path to the pre-trained TCR2vec model. 
        """
        self.tcr2vec_model = TCR2vec(model_path).to(device)
        self.tcr2vec_model.eval()

    def read_csv(self, file_path='./data/sample.csv', column_name='full_seq'):
        """Read TCR sequences from a CSV file.

        Args:
            file_path (str): The path to the CSV file containing TCR sequences.
            column_name (str): Column name of the provided file recording TCRs. Defaults to 'full_seq'.
        """
        self.dset = TCRLabeledDset(file_path, only_tcr=True, use_column=column_name)

    def encode(self, batch_size = 512):
        """Get the embeddings for the loaded TCR sequences.

        Returns:
            np.ndarray: The embeddings for the TCR sequences, with shape [n, embed_size].
        """
        loader = DataLoader(self.dset, batch_size=batch_size, collate_fn=self.dset.collate_fn, shuffle=False)
        emb = get_emb(self.tcr2vec_model, loader)
        return emb


if __name__ == '__main__':
    # Load model and data
    file_path = './data/sample.csv'

    embedder = EmbeddingTCR2vec()
    embedder.load_model(model_path = './pretrained_models/TCR2vec_120')
    embedder.read_csv(file_path='./data/sample.csv', column_name='full_seq')

    # Get embeddings
    X = embedder.encode()
    print(X.shape)
