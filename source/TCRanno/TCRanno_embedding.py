import tcranno.model_predict as model_predict
import pandas as pd

class EmbeddingTCRanno:

    def __init__(self, model_path=None, file_path='./data/sample.csv'):
        """
        Args:
            model_path (str): h5 format autoencoder model. Set model_path=None to use the default model (provided by TCRanno)
            file_path (str): The path to the CSV file containing TCR sequences.
        """
        self.cdr3s = None
        self.data = pd.read_csv(file_path)
        self.encoder = model_predict.load_encoder(model_path=model_path)


    def encode(self, column_name='aminoAcid'):
        """
        Get the latent representations of TCRs using the autoencoder model. Each TCR is represented as a 32-dimensional vector.
        column_name (str): Column name of the provided file recording TCRs. Defaults to 'aminoAcid'.

        """
        self.cdr3s = self.data[column_name].tolist()
        embedding = model_predict.get_norm_latent(self.cdr3s, self.encoder)
        return embedding

if __name__ == '__main__':
    # Load model and data
    file_path = './data/sample.csv'

    TCRanno = EmbeddingTCRanno(model_path=None, file_path='./data/sample.csv') ## set model_path=None to use the default model (provided by TCRanno)

    # Get embeddings
    X = TCRanno.encode(column_name='aminoAcid')
    print(X[:1])
    print(X.shape)