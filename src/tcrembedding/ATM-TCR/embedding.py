import os
import numpy as np
import torch

import data_io_tf
from attention import Net
import data_process
from data_loader import load_embedding

class EmbeddingATMTCR:

    def __init__(self, file_path, blosum=None, model_name="original.ckpt", cuda=False, seed=1039, model_type="attention", drop_rate=0.25,
                 lin_size=1024, padding="mid", heads=5, max_len_tcr=20, max_len_pep=22, split_type="random"):
        '''

        file_path: The file path specifying the location of the data file to be processed.

        blosum: BLOSUM (Blocks Substitution Matrix), a scoring matrix used for protein sequence alignment. If provided, it is used to transform sequence data. Default to using the BLOSUM45 matrix.

        model_name: The model name, specifying the filename to use when saving or loading the model.

        cuda: A boolean indicating whether to use CUDA acceleration (i.e., run the model on a GPU).

        seed: The random seed used to ensure repeatability of model training and data splitting.

        model_type: The type of model, specifying the model architecture used, for example, "attention" indicates a model using an attention mechanism.

        drop_rate: The dropout rate used in the model's dropout layers to prevent overfitting.

        lin_size: The size of the linear layer, specifying the number of neurons in the model's linear layers.

        padding: The type of padding, used for the padding strategy of sequence data, for example, "mid" indicates padding in the middle.

        heads: The number of heads in the multi-head attention mechanism, applicable only when using attention models.

        max_len_tcr: The maximum length of TCR sequences, used to determine the truncating or padding length of sequences.

        max_len_pep: The maximum length of peptide sequences, used to determine the truncating or padding length of sequences.

        split_type: The type of data splitting, specifying how the data is split into training and test sets, for example, "random" indicates a random split.

        '''
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
       
        self.data_tcr, self.data_epi, self.bound = data_io_tf.read_pTCR_list(self.infile)

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
    
    def embed(self):
        tcrs = data_process.filter_seqs(self.data_tcr)
        tcrs_padded = data_process.pad(tcrs, fix_length=30)
        tcr_num = data_process.numerialize(tcrs_padded)
        tcrs_tensor = [torch.tensor(seq).unsqueeze(0) for seq in tcr_num]
        tcrs_encode = [self.model.get_encode(seq) for seq in tcrs_tensor]

        tcrs_encode_ndarray = []
        for i in tcrs_encode:
            tmp = i.detach().numpy()
            tcrs_encode_ndarray.append(tmp)

        encode_result = np.concatenate(tcrs_encode_ndarray, axis=0)
        return encode_result

    def embed_epitode(self):
        peps = data_process.filter_seqs(self.data_epi)
        peps_padded = data_process.pad(peps, fix_length=20)
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

    encoder = EmbeddingATMTCR(file_path="data/testdata_ATM-TCR.csv")
    encode_tcr = encoder.embed()
    encode_pep = encoder.embed_epitode()
    print(encode_tcr.shape)
    print(encode_pep.shape)
    
