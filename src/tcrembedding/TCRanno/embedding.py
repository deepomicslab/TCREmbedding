import tcranno.model_predict as model_predict
import pandas as pd

class EmbeddingTCRanno:

    def __init__(self):
        self.encoder = None
        self.cdr3s = None

    def load_model(self, model_path=None):
        """
        Args:
            model_path (str): h5 format autoencoder model. Set model_path=None to use the default model (provided by TCRanno)
        """


        self.encoder = model_predict.load_encoder(model_path=model_path)

    def load_data(self, file_path='./data/sample.csv', column_name='aminoAcid'):
        """Read TCR sequences from a CSV file.

        Args:
            file_path (str): The path to the CSV file containing TCR sequences.
            column_name (str): Column name of the provided file recording TCRs. Defaults to 'aminoAcid'.
        """
        self.cdr3s = pd.read_csv(file_path, sep=",")[column_name].tolist()

    def embed(self):
        """
        Get the latent representations of TCRs using the autoencoder model. Each TCR is represented as a 32-dimensional vector.
        """
        embedding = model_predict.get_norm_latent(self.cdr3s, self.encoder)
        return embedding

if __name__ == '__main__':
    # Load model and data
    embedder = EmbeddingTCRanno()
    embedder.load_model(model_path = None)  ## set model_path=None to use the default model (provided by TCRanno)
    embedder.load_data(file_path='data/testdata_TCRanno.csv', column_name='CDR3b')

    # Get embeddings
    X = embedder.embed()
    print(X.shape)
