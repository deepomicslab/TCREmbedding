import os
import json
import numpy as np
from pytoda.smiles.smiles_language import SMILESTokenizer
from rdkit import Chem
import pandas as pd

bl62={'A':[4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0],
      'R':[-1,4,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3],
      'N':[-2,0,4,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3],
      'D':[-2,-2,1,4,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3],
      'C':[0,-3,-3,-3,4,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],
      'Q':[-1,1,0,0,-3,4,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2],
      'E':[-1,0,0,2,-4,2,4,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2],
      'G':[0,-2,0,-1,-3,-2,-2,4,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3],
      'H':[-2,0,1,-1,-3,0,0,-2,4,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3],
      'I':[-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3],
      'L':[-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1],
      'K':[-1,2,0,-1,-3,1,1,-2,-1,-3,-2,4,-1,-3,-1,0,-1,-3,-2,-2],
      'M':[-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,4,0,-2,-1,-1,-1,-1,1],
      'F':[-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,4,-4,-2,-2,1,3,-1],
      'P':[-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,4,-1,-1,-4,-3,-2],
      'S':[1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2],
      'T':[0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,4,-2,-2,0],
      'W':[-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,4,2,-3],
      'Y':[-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,4,-1],
      'V':[0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4]}

class EmbeddingTITAN:
    def __init__(self, model_path='TITAN_model', params_filepath='TITAN_model/model_params.json'):
        self.model_path = model_path
        self.params_filepath = params_filepath
        self.model = None

    def load_data(self, file_path, use_columns):
        self.data = pd.read_csv(file_path, header=0, sep=",")[use_columns]

    def load_model(self):
        # Read model params
        params_filepath = self.params_filepath
        params = {}
        with open(params_filepath) as fp:
            params.update(json.load(fp))

        # Load languages
        smiles_language = SMILESTokenizer.from_pretrained(self.params_filepath)
        smiles_language.load_vocabulary(vocab_file="TITAN_model/vocab.json")
        smiles_language.set_encoding_transforms(
            randomize=None,
            add_start_and_stop=params.get('ligand_start_stop_token', False),
            padding=params.get('ligand_padding', True),
            padding_length=params.get('ligand_padding_length', True),
        )
        smiles_language.set_smiles_transforms(
            augment=False,
        )
        self.model = smiles_language
    
    def embed(self):
        num_sequences = len(self.data)
        max_sequence_length = max(len(seq) for seq in self.data)
        encoding_shape = (num_sequences, max_sequence_length, len(bl62['A']))
        encoding_result = np.zeros(encoding_shape)

        for i, sequence in enumerate(self.data):
            for j, residue in enumerate(sequence):
                if residue in bl62:
                    encoding_result[i, j, :] = bl62[residue]
                else:
                    raise ValueError(f"Unknown residue: {residue}")
        return encoding_result

    def embed_epi(self):
        # embedding
        encoding_result = []
        peptide_smiles = []
        for sequence in self.data:
            mol = Chem.MolFromSequence(sequence)
            if mol is None:
                print("%s can not convert to molecule" % sequence)
                continue
            else:
                # get the SMILES representation of sequence
                peptide_smile = Chem.MolToSmiles(mol)
                peptide_smiles.append(peptide_smile)
                # embedding
                embed = self.model.smiles_to_token_indexes(peptide_smile)
                embed = embed.detach().numpy()
                encoding_result.append(embed)

        return peptide_smiles, encoding_result

if __name__ == "__main__":

    encoder = EmbeddingTITAN()
    encoder.load_data("data/testdata_TITAN.csv", use_columns="CDR3b")
    encoder.load_model()
    TCR_encode_result = encoder.embed()
    epi_encode_result = encoder.embed_epi()
    print(TCR_encode_result.shape)
    print(epi_encode_result)
