import os
from os import listdir
from os.path import isfile, join
import numpy as np
import pandas as pd
import pytorch_lightning as pl
from argparse import Namespace
import torch

from Trainer import ERGOLightning
from data_loader import Data_Loader

class EmbeddingERGO():

    def __init__(self, tcr_encoding_model="AE"):
        self.tcr_encoding_model = tcr_encoding_model
        self.model = None

    def get_model(self, hparams, checkpoint_path):
        model = ERGOLightning(hparams)
        checkpoint = torch.load(checkpoint_path, map_location='cpu')
        model.load_state_dict(checkpoint['state_dict'])
        model.eval()
        return model

    def load_model(self, model_path, args_path):
        # get model file from version
        checkpoint_path = model_path
        files = [f for f in listdir(checkpoint_path) if isfile(join(checkpoint_path, f))]
        checkpoint_path = os.path.join(checkpoint_path, files[0])
        # get args from version
        args_path = args_path
        with open(args_path, 'r') as file:
            lines = file.readlines()
            args = {}
            for line in lines[1:]:
                key, value = line.strip().split(',')
                if key in ['dataset', 'tcr_encoding_model', 'cat_encoding']:
                    args[key] = value
                else:
                    args[key] = eval(value)
        hparams = Namespace(**args)
        checkpoint = checkpoint_path
        model = self.get_model(hparams, checkpoint)
        self.model = model
        # return model

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
        model = self.model
        # PEPTIDE Encoder:
        pep_encoding = model.pep_encoder(*pep_batch)
        # TCR Encoder:
        tcra, tcrb = tcrb_batch
        tcrb_encoding = model.tcrb_encoder(*tcrb)

        return tcrb_encoding, pep_encoding


if __name__ == '__main__':
    encoder = EmbeddingERGO(tcr_encoding_model="AE")    
    encoder.load_model(model_path="/media/lihe/TCR/ERGO-II/ERGO-II_paper_logs/paper_models/version_10ve/checkpoints", args_path="/media/lihe/TCR/ERGO-II/ERGO-II_paper_logs/paper_models/10ve/meta_tags.csv")

    data_loader = Data_Loader()
    tcrb_list, peptide = data_loader.collate("/media/lihe/TCR/data/input/TCRantigenData_detailed1.csv", "AE")

    tcrb_batch, pep_batch = encoder.forward(tcrb_list)

    tcrb_encoding, pep_encoding = encoder.embed(tcrb_batch, pep_batch)


    tcr_encode_result = tcrb_encoding.detach().numpy()
    pep_encode_result = pep_encoding.detach().numpy()
    print(tcrb_encoding.shape)
    print(pep_encoding.shape)
    #np.save("ERGOII_cdr3_ae.npy", tcr_encode_result)
