from utils.utils_aa_embedding import read_blosum, encode_as_blosum, encode_as_onehot
import pandas as pd
import numpy as np
from typing import Union

class EmbeddingTcellMatch:

    def __init__(self):
        self.fn_blosum = None
        self.blosum_encoding = False
        self.chain = "trb"
        self.eos_char = "_"  # End-of-sequence char.
        self.aa_list = [
            'A', 'R', 'N', 'D', 'C',
            'Q', 'E', 'G', 'H', 'I',
            'L', 'K', 'M', 'F', 'P',
            'S', 'T', 'W', 'Y', 'V',
            'B', 'J', 'Z', 'X', '*',
            self.eos_char
        ]
        idx_list = np.arange(0, len(self.aa_list))
        self.dict_aa = dict(zip(self.aa_list, idx_list))
        self.tcr_len = None
        self.pep_len = None

    def _parseinput(self, data, tcr):
        seqs = []
        for index, value in data.items():
            flag = True
            if tcr:
                if value.startswith('C') and value.endswith('F'):
                    for i in list(value):
                        if i not in self.aa_list:
                            flag = False
                            break
                    if flag:
                        seqs.append(value)
                    else:
                        continue
            else:
                for i in list(value):
                    if i not in self.aa_list:
                        flag = False
                        break
                if flag:
                    seqs.append(value)
                else:
                    continue
        return seqs

    def _format_data_aa(
            self,
            x: list,
            fn_blosum: Union[str, None] = None,
            tcr: bool = False,
            chain: str = "trb"
    ) -> np.ndarray:
        """
        Create numeric input data from amino acid sequence data.

        :param x: Input as list (length observations) of lists with strings of amino acid code of each chain.
        :param fn_blosum:
        :return: 4D tensor [observations, chains, amino acid position, amino acid embedding]
            One-hot encoded input.
        """
        if fn_blosum is not None:
            # Blosum encoding
            blosum_embedding = read_blosum(fn_blosum)
            x_encoded = encode_as_blosum(x=x, blosum_embedding=blosum_embedding)
        else:
            # One-hot encoding
            x_encoded = encode_as_onehot(x=x, dict_aa=self.dict_aa, eos_char=self.eos_char)

        if tcr:
            x_encoded = self._format_tcr_chains(x=x_encoded, chain=chain)
        else:
            x_encoded = self._format_antigens(x=x_encoded)
        return x_encoded

    def _format_tcr_chains(
            self,
            x,
            chain="trb"
    ):
        """

        :param x:
        :return:
        """
        if chain == "tra":
            x = np.expand_dims(x[:, 0, :, :], axis=1)
        elif chain == "trb":
            x = np.expand_dims(x[:, 1, :, :], axis=1)
        elif chain == "separate":
            pass
        elif chain == "concat":
            x = np.concatenate([
                np.expand_dims(x[:, 0, :, :], axis=1),
                np.expand_dims(x[:, 1, :, :], axis=1)
            ], axis=2)
        else:
            raise ValueError("self.chains %s not recognized" % chain)

        if self.tcr_len is None:
            self.tcr_len = x.shape[2]
        x = self._pad_tcr(x=x)
        return x

    def _pad_tcr(
            self,
            x: np.ndarray
    ):
        """ Pad TCR to desired length.

        Takes care of chain concatenation: If self.chain is "concat", splits x equally in axis 2 into TRA and TRB
        and pads each to self.tcr_len / 2, then concatenates in axis 2 again. If self.chain is "tra" or "trb", pads
        x to self.tcr_len in axis 2.

        :param x: TCR sequence encoding.
        :return: Padded TCR encoding.
        """
        def pad_block(shape):
            pad_embedding = np.zeros([1, len(self.aa_list)])
            pad_embedding[0, -1] = 1
            return np.zeros(shape) + pad_embedding

        if self.chain.lower() == "concat":
            assert self.tcr_len % 2 == 0, \
                "self.tcr_len (%i) must be divisible by two if 'concat' mode is used for TCR." % self.tcr_len
            assert x.shape[2] % 2 == 0, \
                "dimension 3 of x (%i) must be divisible by two if 'concat' mode is used for TCR." % x.shape[2]
            assert x.shape[2] <= self.tcr_len, \
                "Required tcr length (%i) must at least as high as existing tcr length (%i)." % \
                (self.tcr_len, x.shape[2])
            xa = x[:, :, :int(x.shape[2] / 2), :]
            xb = x[:, :, int(x.shape[2] / 2):, :]
            xa = np.concatenate([
                xa, pad_block([x.shape[0], x.shape[1], int(self.tcr_len / 2) - int(x.shape[2] / 2), x.shape[3]])
            ], axis=2)
            xb = np.concatenate([
                xb, pad_block([x.shape[0], x.shape[1], int(self.tcr_len / 2) - int(x.shape[2] / 2), x.shape[3]])
            ], axis=2)
            x = np.concatenate([xa, xb], axis=2)
        elif self.chain.lower() in ["tra", "trb", "separate"]:
            x = np.concatenate([
                x, pad_block([x.shape[0], x.shape[1], self.tcr_len - x.shape[2], x.shape[3]])
            ], axis=2)
        else:
            raise ValueError("self.chains %s not recognized" % self.chain)
        return x

    def _format_antigens(
            self,
            x
    ):
        """
        :param x:
        :return:
        """
        def pad_block(shape):
            pad_embedding = np.zeros([1, len(self.aa_list)])
            pad_embedding[0, -1] = 1
            return np.zeros(shape) + pad_embedding

        x = np.expand_dims(x[:, 0, :, :], axis=1)
        if self.pep_len is None:
            self.pep_len = x.shape[2]
        else:
            assert self.pep_len >= x.shape[2], \
                "Pre-set antigen length (%i) is smaller than found antigen length %i." % (self.pep_len, x.shape[2])
            x = np.concatenate([
                x, pad_block([x.shape[0], x.shape[1], self.pep_len - x.shape[2], x.shape[3]])
            ], axis=2)
        return x


    def embed(self, file_path, blosum_path=None, chain="trb", blosum_encoding=False):
        row_cdr3_seqs = pd.read_csv(file_path, header=0, delimiter=",")["CDR3"]
        cdr3_seqs = self._parseinput(row_cdr3_seqs, tcr=True)
        row_epi_seqs = pd.read_csv(file_path, header=0, delimiter=",")["Epitope"]
        epi_seqs = self._parseinput(row_epi_seqs, tcr=False)
        cdr3seqslist = list()
        episeqslist = list()

        for seq in cdr3_seqs:
            seq_pair = [None, None]
            if chain == "trb":
                seq_pair[1] = seq
            elif chain == "tra":
                seq_pair[0] = seq
            cdr3seqslist.append(seq_pair)

        for seq in epi_seqs:
            seq_pair = [None]
            seq_pair[0] = seq
            episeqslist.append(seq_pair)

        epis = self._format_data_aa(
            x=episeqslist,
            fn_blosum=blosum_path if blosum_encoding else None,
            tcr=False
        )
        cdr3s = self._format_data_aa(
            x=cdr3seqslist,
            fn_blosum=blosum_path if blosum_encoding else None,
            tcr=True,
            chain=chain
        )
        
        return epis, cdr3s
