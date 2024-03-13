import torch
import torch.nn as nn
import pytorch_lightning as pl
import numpy as np
import pickle

from data_loader import Data_Loader
from Models import AE_Encoder, LSTM_Encoder
from data_loader import get_index_dicts

class EmbeddingERGO(pl.LightningModule):

    def __init__(self, dataset="vdjdb", tcr_encoding_model="LSTM", use_alpha=False, use_vj=False, use_mhc=False,
                 use_t_type=False, cat_encoding="binary", aa_embedding_dim=80, cat_embedding_dim=50, lstm_dim=500,
                 encoding_dim=100, dropout_rate=0.1, lr=1e-3, wd=1e-5):
        super(EmbeddingERGO, self).__init__()
        self.dataset = dataset
        # Model Type
        self.tcr_encoding_model = tcr_encoding_model
        self.use_alpha = use_alpha
        self.use_vj = use_vj
        self.use_mhc = use_mhc
        self.use_t_type = use_t_type
        self.cat_encoding = cat_encoding
        # Dimensions
        self.aa_embedding_dim = aa_embedding_dim
        self.cat_embedding_dim = cat_embedding_dim
        self.lstm_dim = lstm_dim
        self.encoding_dim = encoding_dim
        self.dropout_rate = dropout_rate
        self.lr = lr
        self.wd = wd

        # get train indicies for V,J etc
        if self.cat_encoding == 'embedding':
            with open('Samples/' + self.dataset + '_train_samples.pickle', 'rb') as handle:
                train = pickle.load(handle)
            vatox, vbtox, jatox, jbtox, mhctox = get_index_dicts(train)
            self.v_vocab_size = len(vatox) + len(vbtox)
            self.j_vocab_size = len(jatox) + len(jbtox)
            self.mhc_vocab_size = len(mhctox)

        # TCR Encoder
        if self.tcr_encoding_model == 'AE':
            if self.use_alpha:
                self.tcra_encoder = AE_Encoder(encoding_dim=self.encoding_dim, tcr_type='alpha', max_len=34)
            self.tcrb_encoder = AE_Encoder(encoding_dim=self.encoding_dim, tcr_type='beta')
        elif self.tcr_encoding_model == 'LSTM':
            if self.use_alpha:
                self.tcra_encoder = LSTM_Encoder(self.aa_embedding_dim, self.lstm_dim, self.dropout_rate)
            self.tcrb_encoder = LSTM_Encoder(self.aa_embedding_dim, self.lstm_dim, self.dropout_rate)
            self.encoding_dim = self.lstm_dim
        # Peptide Encoder
        self.pep_encoder = LSTM_Encoder(self.aa_embedding_dim, self.lstm_dim, self.dropout_rate)
        # Categorical
        self.cat_encoding = cat_encoding
        if cat_encoding == 'embedding':
            if self.use_vj:
                self.v_embedding = nn.Embedding(self.v_vocab_size, self.cat_embedding_dim, padding_idx=0)
                self.j_embedding = nn.Embedding(self.j_vocab_size, self.cat_embedding_dim, padding_idx=0)
            if self.use_mhc:
                self.mhc_embedding = nn.Embedding(self.mhc_vocab_size, self.cat_embedding_dim, padding_idx=0)
        # different mlp sizes, depends on model input
        if self.cat_encoding == 'binary':
            self.cat_embedding_dim = 10
        mlp_input_size = self.lstm_dim + self.encoding_dim
        if self.use_vj:
            mlp_input_size += 2 * self.cat_embedding_dim
        if self.use_mhc:
            mlp_input_size += self.cat_embedding_dim
        if self.use_t_type:
            mlp_input_size += 1
        # MLP I (without alpha)
        self.mlp_dim1 = mlp_input_size
        self.hidden_layer1 = nn.Linear(self.mlp_dim1, int(np.sqrt(self.mlp_dim1)))
        self.relu = torch.nn.LeakyReLU()
        self.output_layer1 = nn.Linear(int(np.sqrt(self.mlp_dim1)), 1)
        self.dropout = nn.Dropout(p=self.dropout_rate)
        # MLP II (with alpha)
        if self.use_alpha:
            mlp_input_size += self.encoding_dim
            if self.use_vj:
                mlp_input_size += 2 * self.cat_embedding_dim
            self.mlp_dim2 = mlp_input_size
            self.hidden_layer2 = nn.Linear(self.mlp_dim2, int(np.sqrt(self.mlp_dim2)))
            self.output_layer2 = nn.Linear(int(np.sqrt(self.mlp_dim2)), 1)

    def forward(self, tcrb_list, pep_list):
        # batch output (always)
        tcrb, pep = tcrb_list, pep_list
        if self.tcr_encoding_model == 'LSTM':
            # get lengths for lstm functions
            len_b = torch.sum((tcrb > 0).int(), dim=1)
            #len_b = torch.sum(len_b > 0, dim=1)

        if self.tcr_encoding_model == 'AE':
            pass
        len_p = torch.sum((pep > 0).int(), dim=1)

        if self.tcr_encoding_model == 'LSTM':
            tcrb_batch = (None, (tcrb, len_b))
        elif self.tcr_encoding_model == 'AE':
            tcrb_batch = (None, (tcrb,))
        pep_batch = (pep, len_p)

        return tcrb_batch, pep_batch

    def embed(self, tcrb_batch, pep_batch):
        # PEPTIDE Encoder:
        pep_encoding = self.pep_encoder(*pep_batch)
        # TCR Encoder:
        tcrb = tcrb_batch
        tcrb_encoding = self.tcrb_encoder(*tcrb)

        return tcrb_encoding, pep_encoding

if __name__ == "__main__":

    encoder = EmbeddingERGO(tcr_encoding_model="LSTM", cat_encoding="binary")
    data_loader = Data_Loader()
    tcrb_list, peptide = data_loader.collate("data/testdata_ERGO-II.csv", "LSTM")
    tcr_batch, pep_batch = encoder.forward(tcrb_list, peptide)
    tcra_batch, tcrb_batch = tcr_batch
    tcrb_encoding, pep_encoding = encoder.embed(tcrb_batch, pep_batch)
    tcr_encode_result = tcrb_encoding.detach().numpy()
    pep_encode_result = pep_encoding.detach().numpy()
    print(tcrb_encoding.shape)
    print(pep_encoding.shape)
