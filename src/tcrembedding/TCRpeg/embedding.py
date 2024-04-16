import pandas as pd
from .tcrpeg.TCRpeg import TCRpeg

class EmbeddingTCRpeg:
    def __init__(self, hidden_size=64, num_layers=3, embedding_path='tcrpeg/data/embedding_32.txt', model_path='tcrpeg/models/tcrpeg.pth', device='cpu', load_data=False):
        """
        Initializes an EmbeddingTCRpeg object with specified TCRpeg model configurations.

        Args:
            hidden_size (int): The number of features in the hidden state of the model.
            num_layers (int): The number of recurrent layers in the model.
            embedding_path (str): The path to the embedding file to be used by the model.
            model_path (str): The path to the pre-trained model file (.pth) to load.
            device (str): The device to run the model on ('cpu' or 'cuda:index').
            load_data (bool): Flag to indicate whether to load data upon initialization. Typically false when embedding is the main purpose.
        """
        self.model = TCRpeg(hidden_size=hidden_size, num_layers=num_layers, embedding_path=embedding_path, load_data=load_data, device=device)
        self.model.create_model(load=True, path=model_path)

    def load_data(self, file_path, use_columns, sep=","):
        self.data = pd.read_csv(file_path, header=0, sep=sep)[use_columns]

    def embed(self):
        """
        Generates embeddings for TCR sequences provided in a CSV file.

        Args:
            file_path (str): The path to the CSV file containing TCR sequences in the 'CDR3' column.

        Returns:
            embs (list): A list of embeddings for the TCR sequences.
        """
        embs = self.model.get_embedding(self.data)
        return embs

# Example usage
if __name__ == "__main__":
    embedder = EmbeddingTCRpeg(hidden_size=128, num_layers=4, device="cpu")
    embeddings = embedder.embed("data/testdata_TCRpeg.csv", use_columns="CDR3b")
    print(embeddings)
