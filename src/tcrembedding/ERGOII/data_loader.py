import torch
from torch.utils.data import Dataset
import pandas as pd

class Data_Loader(Dataset):
    def __init__(self):
        self.amino_acids = [letter for letter in 'ARNDCEQGHILKMFPSTWYV']
        self.atox = {amino: index for index, amino in enumerate(['PAD'] + self.amino_acids + ['X'])}
        self.pos_weight_factor = 5

    def aa_convert(self, seq):
        if seq == 'UNK':
            seq = []
        else:
            seq = [self.atox[aa] for aa in seq]
        return seq

    @staticmethod
    def get_max_length(x):
        return len(max(x, key=len))

    def seq_letter_encoding(self, seq):
        def _pad(_it, _max_len):
            return _it + [0] * (_max_len - len(_it))
        return [_pad(it, self.get_max_length(seq)) for it in seq]

    def seq_one_hot_encoding(self, tcr, max_len=28):
        tcr_batch = list(tcr)
        padding = torch.zeros(len(tcr_batch), max_len, 20 + 1)
        # TCR is converted to numbers at this point
        # We need to match the autoencoder atox, therefore -1
        for i in range(len(tcr_batch)):
            # missing alpha
            if tcr_batch[i] == [0]:
                continue
            tcr_batch[i] = tcr_batch[i] + [self.atox['X']]
            for j in range(min(len(tcr_batch[i]), max_len)):
                padding[i, j, tcr_batch[i][j] - 1] = 1
        return padding

    @staticmethod
    def binarize(num):
        l = []
        while num:
            l.append(num % 2)
            num //= 2
        l.reverse()
        # print(l)
        return l

    def filter_seqs(self, seqs):
        filtered_list = []
        for seq in seqs:
            flag = True
            for char in seq:
                if char not in self.amino_acids:
                    flag = False
                    break
            if flag:
                filtered_list.append(seq)
            else:
                continue
        return filtered_list

    def collate(self, file_path, tcr_encoding):
        data = pd.read_csv(file_path, header=0, sep=",")
        tcrb_list, pep_list = data["CDR3b"].tolist(), data["Epitope"].tolist()
        #tcrb_list = data["CDR3"].tolist()
        tcrb_list = self.filter_seqs(tcrb_list)
        pep_list = self.filter_seqs(pep_list)
        # TCRs
        if len(tcrb_list) == 0 or len(pep_list) == 0:
            return None
        tcrb = [self.aa_convert(sample) for sample in tcrb_list]

        if tcr_encoding == 'AE':
            tcrb_lst = torch.FloatTensor(self.seq_one_hot_encoding(tcrb))
        elif tcr_encoding == 'LSTM':
            # we do not send the length, so that ae and lstm batch output be similar
            tcrb_lst = torch.LongTensor(self.seq_letter_encoding(tcrb))
        # Peptide
        peptide = [self.aa_convert(sample) for sample in pep_list]
        pep_lst = torch.LongTensor(self.seq_letter_encoding(peptide))

        return tcrb_lst, pep_lst

def get_index_dicts(train_samples):
    samples = train_samples
    all_va = [sample['va'] for sample in samples if not pd.isna(sample['va'])]
    vatox = {va: index for index, va in enumerate(sorted(set(all_va)), 1)}
    vatox['UNK'] = 0
    all_vb = [sample['vb'] for sample in samples if not pd.isna(sample['vb'])]
    vbtox = {vb: index for index, vb in enumerate(sorted(set(all_vb)), 1)}
    vbtox['UNK'] = 0
    all_ja = [sample['ja'] for sample in samples if not pd.isna(sample['ja'])]
    jatox = {ja: index for index, ja in enumerate(sorted(set(all_ja)), 1)}
    jatox['UNK'] = 0
    all_jb = [sample['jb'] for sample in samples if not pd.isna(sample['jb'])]
    jbtox = {jb: index for index, jb in enumerate(sorted(set(all_jb)), 1)}
    jbtox['UNK'] = 0
    all_mhc = [sample['mhc'] for sample in samples if not pd.isna(sample['mhc'])]
    mhctox = {mhc: index for index, mhc in enumerate(sorted(set(all_mhc)), 1)}
    mhctox['UNK'] = 0
    return [vatox, vbtox, jatox, jbtox, mhctox]
