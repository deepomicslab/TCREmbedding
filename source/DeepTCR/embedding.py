from DeepTCR.DeepTCR import DeepTCR_U

class EmbeddingDeepTCR:
    def __init__(self):
        pass

    def load_model(self, train_data_directory, model_folder_name="Test_Model", Load_Prev_Data=False):
        self.model = DeepTCR_U(model_folder_name)
        # Load Data from directories
        self.model.Get_Data(directory=train_data_directory, aggregate_by_aa=True, Load_Prev_Data=Load_Prev_Data,
                       aa_column_beta=1)

        # Train VAE
        self.model.Train_VAE(Load_Prev_Data=Load_Prev_Data)

    def encode(self, encode_data_directory="Data/data", aa_column_beta=0, sep="\t"):
        # Instantiate training object
        DTCRU_rudq = DeepTCR_U('encode')

        # Load Data from directories
        DTCRU_rudq.Get_Data(directory=encode_data_directory, Load_Prev_Data=False, aggregate_by_aa=False,
                            aa_column_beta=aa_column_beta, sep=sep)
        beta_sequences = DTCRU_rudq.beta_sequences

        features = self.model.Sequence_Inference(beta_sequences=beta_sequences)

        return features

if __name__ == "__main__":
    encoder = EmbeddingDeepTCR()
    encoder.load_model(train_data_directory="Test_Model", Load_Prev_Data=True)
    tcr_features = encoder.encode()
    print(tcr_features.shape)
    #np.save("D:/TCR/cluster/DeepTCR_tcr.npy", tcr_features)