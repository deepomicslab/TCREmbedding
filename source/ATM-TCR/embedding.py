import os
import torch
import torch.nn.functional as F
import numpy as np

import data_io_tf
from attention import Net
import data_process
from data_loader import load_embedding

# Constants
PRINT_EVERY_EPOCH = 1

def train(model, device, train_loader, optimizer, epoch):

    model.train()

    for batch in train_loader:

        x_pep, x_tcr, y = batch.X_pep.to(
            device), batch.X_tcr.to(device), batch.y.to(device)

        optimizer.zero_grad()
        yhat = model(x_pep, x_tcr)
        y = y.unsqueeze(-1).expand_as(yhat)
        loss = F.binary_cross_entropy(yhat, y)
        loss.backward()
        optimizer.step()

    if epoch % PRINT_EVERY_EPOCH == 1:
        print('[TRAIN] Epoch {} Loss {:.4f}'.format(epoch, loss.item()))


class EmbeddingATMTCR:

    def __init__(self, file_path, blosum=None, model_name="original.ckpt", cuda=False, seed=1039, model_type="attention", drop_rate=0.25,
                 lin_size=1024, padding="mid", heads=5, max_len_tcr=20, max_len_pep=22, split_type="random"):

        self.infile = file_path
        self.blosum = blosum
        self.model_name = model_name
        self.cuda = cuda
        self.seed = seed
        self.model_type = model_type
        self.drop_rate = drop_rate
        self.lin_size = lin_size
        self.padding = padding
        self.heads = heads
        self.max_len_tcr = max_len_tcr
        self.max_len_pep = max_len_pep
        self.split_type = split_type

        # Set Cuda
        if torch.cuda.is_available() and not self.cuda:
            print("WARNING: You have a CUDA device, so you should probably run with --cuda")
        device = torch.device('cuda' if self.cuda else 'cpu')

        # Set random seed
        torch.manual_seed(self.seed)
        if self.cuda:
            torch.cuda.manual_seed(self.seed)

        # read data
        if self.model_type != 'attention':
            raise ValueError('unknown model name')
        #self.data_tcr, self.data_pep, self.bound = data_io_tf.read_pTCR_list(self.infile)
        # read cluster data
        self.data_tcr = data_io_tf.read_pTCR_list(self.infile)
        # read iedb data
        #self.data_tcr, self.data_epi = data_io_tf.read_pTCR_list(self.infile)
        #read combined_dataset
        #self.data_tcr, self.data_epi = data_io_tf.read_pTCR_list(self.infile)
        # load model
        embedding_matrix = load_embedding(self.blosum)

        model = Net(embedding_matrix, self.blosum, self.heads, self.lin_size, self.max_len_pep, self.max_len_tcr,
                         self.drop_rate).to('cpu')

        # eax1it model
        model_name = self.model_name

        assert model_name in os.listdir('./models')

        model_path = './models/' + model_name
        model.load_state_dict(torch.load(model_path, map_location=torch.device('cpu')))
        self.model = model

    '''
    def embedding(self, peps, tcrs):
        # Load embedding matrix
        embedding_matrix = load_embedding(self.blosum)

        model = Net(embedding_matrix, self.blosum, self.heads, self.lin_size, self.max_len_pep, self.max_len_tcr, self.drop_rate).to('cpu')

        # eax1it model
        model_name = self.model_name

        assert model_name in os.listdir('./models')

        model_name = './models/' + model_name
        model.load_state_dict(torch.load(model_name, map_location=torch.device('cpu')))

        peps = data_process.filter_seqs(peps)
        tcrs = data_process.filter_seqs(tcrs)

        peps_padded = data_process.pad(peps, fix_length=22)
        tcrs_padded = data_process.pad(tcrs, fix_length=20)

        pep_num = data_process.numerialize(peps_padded)
        tcr_num = data_process.numerialize(tcrs_padded)

        peps_tensor = [torch.tensor(seq) for seq in pep_num]
        tcrs_tensor = [torch.tensor(seq) for seq in tcr_num]

        peps_embed = [model.get_embeddings(seq) for seq in peps_tensor]
        tcrs_embed = [model.get_embeddings(seq) for seq in tcrs_tensor]

        peps_embed_ndarray = []
        tcrs_embed_ndarray = []
        for i in peps_embed:
            tmp = i.detach().numpy()
            peps_embed_ndarray.append(tmp)
        for i in tcrs_embed:
            tmp = i.detach().numpy()
            tcrs_embed_ndarray.append(tmp)

        return peps_embed_ndarray, tcrs_embed_ndarray
    '''

    def encode(self):
        tcrs = data_process.filter_seqs(self.data_tcr)
        tcrs_padded = data_process.pad(tcrs, fix_length=20)
        tcr_num = data_process.numerialize(tcrs_padded)
        tcrs_tensor = [torch.tensor(seq).unsqueeze(0) for seq in tcr_num]
        tcrs_encode = [self.model.get_encode(seq) for seq in tcrs_tensor]

        tcrs_encode_ndarray = []
        for i in tcrs_encode:
            tmp = i.detach().numpy()
            tcrs_encode_ndarray.append(tmp)

        encode_result = np.concatenate(tcrs_encode_ndarray, axis=0)
        return encode_result

    def encode_epitode(self):
        peps = data_process.filter_seqs(self.data_epi)
        peps_padded = data_process.pad(peps, fix_length=22)
        pep_num = data_process.numerialize(peps_padded)
        peps_tensor = [torch.tensor(seq).unsqueeze(0) for seq in pep_num]
        peps_encode = [self.model.get_encode(seq) for seq in peps_tensor]

        peps_encode_ndarray = []
        for i in peps_encode:
            tmp = i.detach().numpy()
            peps_encode_ndarray.append(tmp)

        encode_result = np.concatenate(peps_encode_ndarray, axis=0)
        return encode_result

if __name__ == "__main__":
    encoder = EmbeddingATMTCR(file_path=f"D:/TCR/TCRantigenData_top5.tsv")
    encode_tcr = encoder.encode()
    #encode_epi = encoder.encode_epitode()

    print(encode_tcr.shape)
    np.save(f"D:/TCR/cluster/ATM-TCR_tcr_top5.npy", encode_tcr)

    '''print(encode_epi.shape)
    np.save(f"D:/TCR/combined_dataset/ATM-TCR_epi.npy", encode_epi)'''
