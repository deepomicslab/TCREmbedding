import os
from io import StringIO
import numpy as np
import pandas as pd
import csv
from tensorflow.keras.models import Model,load_model

########################### Atchley's factors#######################
aa_dict_atchley=dict()

########################### One Hot 'X' is a padding variable ##########################
aa_dict_one_hot = {'A': 0,'C': 1,'D': 2,'E': 3,'F': 4,'G': 5,'H': 6,'I': 7,'K': 8,'L': 9,
           'M': 10,'N': 11,'P': 12,'Q': 13,'R': 14,'S': 15,'T': 16,'V': 17,
           'W': 18,'Y': 19,'X': 20}

########################### Blosum ##########################
BLOSUM50_MATRIX = pd.read_table(StringIO(u"""                                                                                      
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  J  Z  X  *                                                           
A  5 -2 -1 -2 -1 -1 -1  0 -2 -1 -2 -1 -1 -3 -1  1  0 -3 -2  0 -2 -2 -1 -1 -5                                                           
R -2  7 -1 -2 -4  1  0 -3  0 -4 -3  3 -2 -3 -3 -1 -1 -3 -1 -3 -1 -3  0 -1 -5                                                           
N -1 -1  7  2 -2  0  0  0  1 -3 -4  0 -2 -4 -2  1  0 -4 -2 -3  5 -4  0 -1 -5                                                           
D -2 -2  2  8 -4  0  2 -1 -1 -4 -4 -1 -4 -5 -1  0 -1 -5 -3 -4  6 -4  1 -1 -5                                                           
C -1 -4 -2 -4 13 -3 -3 -3 -3 -2 -2 -3 -2 -2 -4 -1 -1 -5 -3 -1 -3 -2 -3 -1 -5                                                           
Q -1  1  0  0 -3  7  2 -2  1 -3 -2  2  0 -4 -1  0 -1 -1 -1 -3  0 -3  4 -1 -5                                                           
E -1  0  0  2 -3  2  6 -3  0 -4 -3  1 -2 -3 -1 -1 -1 -3 -2 -3  1 -3  5 -1 -5                                                           
G  0 -3  0 -1 -3 -2 -3  8 -2 -4 -4 -2 -3 -4 -2  0 -2 -3 -3 -4 -1 -4 -2 -1 -5                                                           
H -2  0  1 -1 -3  1  0 -2 10 -4 -3  0 -1 -1 -2 -1 -2 -3  2 -4  0 -3  0 -1 -5                                                          
I -1 -4 -3 -4 -2 -3 -4 -4 -4  5  2 -3  2  0 -3 -3 -1 -3 -1  4 -4  4 -3 -1 -5                                                           
L -2 -3 -4 -4 -2 -2 -3 -4 -3  2  5 -3  3  1 -4 -3 -1 -2 -1  1 -4  4 -3 -1 -5                                                           
K -1  3  0 -1 -3  2  1 -2  0 -3 -3  6 -2 -4 -1  0 -1 -3 -2 -3  0 -3  1 -1 -5                                                           
M -1 -2 -2 -4 -2  0 -2 -3 -1  2  3 -2  7  0 -3 -2 -1 -1  0  1 -3  2 -1 -1 -5                                                           
F -3 -3 -4 -5 -2 -4 -3 -4 -1  0  1 -4  0  8 -4 -3 -2  1  4 -1 -4  1 -4 -1 -5                                                           
P -1 -3 -2 -1 -4 -1 -1 -2 -2 -3 -4 -1 -3 -4 10 -1 -1 -4 -3 -3 -2 -3 -1 -1 -5                                                           
S  1 -1  1  0 -1  0 -1  0 -1 -3 -3  0 -2 -3 -1  5  2 -4 -2 -2  0 -3  0 -1 -5                                                           
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  2  5 -3 -2  0  0 -1 -1 -1 -5                                                           
W -3 -3 -4 -5 -5 -1 -3 -3 -3 -3 -2 -3 -1  1 -4 -4 -3 15  2 -3 -5 -2 -2 -1 -5                                                           
Y -2 -1 -2 -3 -3 -1 -2 -3  2 -1 -1 -2  0  4 -3 -2 -2  2  8 -1 -3 -1 -2 -1 -5                                                           
V  0 -3 -3 -4 -1 -3 -3 -4 -4  4  1 -3  1 -1 -3 -2  0 -3 -1  5 -3  2 -3 -1 -5                                                           
B -2 -1  5  6 -3  0  1 -1  0 -4 -4  0 -3 -4 -2  0  0 -5 -3 -3  6 -4  1 -1 -5                                                           
J -2 -3 -4 -4 -2 -3 -3 -4 -3  4  4 -3  2  1 -3 -3 -1 -2 -1  2 -4  4 -3 -1 -5                                                           
Z -1  0  0  1 -3  4  5 -2  0 -3 -3  1 -1 -4 -1  0 -1 -2 -2 -3  1 -3  5 -1 -5                                                           
X -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -5                                                           
* -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5  1                                                           
"""), sep='\s+').loc[list(aa_dict_one_hot.keys()), list(aa_dict_one_hot.keys())]
assert (BLOSUM50_MATRIX == BLOSUM50_MATRIX.T).all().all()

ENCODING_DATA_FRAMES = {
    "BLOSUM50": BLOSUM50_MATRIX,
    "one-hot": pd.DataFrame([
        [1 if i == j else 0 for i in range(len(aa_dict_one_hot.keys()))]
        for j in range(len(aa_dict_one_hot.keys()))
    ], index=aa_dict_one_hot.keys(), columns=aa_dict_one_hot.keys())
}

class EmbeddingpMTnet:

    def __init__(self):
        pass

    def preprocess(self, filedir):
        # Preprocess TCR files
        if not os.path.exists(filedir):
            print('Invalid file path: ' + filedir)
            return 0
        dataset = pd.read_csv(filedir, header=0, sep=",")
        dataset = dataset.dropna()
        # Remove antigen that is longer than 15aa
        # dataset = dataset[dataset.Antigen.str.len() < 16]
        TCR_list = dataset['CDR3b'].tolist()
        antigen_list = dataset['Epitope'].tolist()

        return TCR_list, antigen_list

    def aamapping_TCR(self, peptideSeq, aa_dict, encode_dim):
        # Transform aa seqs to Atchley's factors.
        peptideArray = []
        if len(peptideSeq) > encode_dim:
            print('Length: ' + str(len(peptideSeq)) + ' over bound!')
            peptideSeq = peptideSeq[0:encode_dim]
        for aa_single in peptideSeq:
            try:
                peptideArray.append(aa_dict[aa_single])
            except KeyError:
                print('Not proper aaSeqs: ' + peptideSeq)
                peptideArray.append(np.zeros(5, dtype='float64'))
        for i in range(0, encode_dim - len(peptideSeq)):
            peptideArray.append(np.zeros(5, dtype='float64'))
        return np.asarray(peptideArray)

    def peptide_encode_HLA(self, peptide, maxlen, encoding_method):
        '''
        Convert peptide amino acid sequence to one-hot encoding,
        optionally left padded with zeros to maxlen(15).

        The letter 'X' is interpreted as the padding character and
        is assigned a value of zero.

        e.g. encode('SIINFEKL', maxlen=12)
                 := [16,  8,  8, 12,  0,  0,  0,  0,  5,  4,  9, 10]

        Parameters
        ----------
        peptide:string of peptide comprising amino acids
        maxlen : int, default 15
            Pad peptides to this maximum length. If maxlen is None,
            maxlen is set to the length of the first peptide.

        Returns
        -------
        '''
        if len(peptide) > maxlen:
            msg = 'Peptide %s has length %d > maxlen = %d.'
            raise ValueError(msg % (peptide, len(peptide), maxlen))
        peptide = peptide.replace(u'\xa0', u'')  # remove non-breaking space
        o = list(map(lambda x: aa_dict_one_hot[x.upper()] if x.upper() in aa_dict_one_hot.keys() else 20, peptide))
        # if the amino acid is not valid, replace it with padding aa 'X':20
        k = len(o)
        # use 'X'(20) for padding
        o = o[:k // 2] + [20] * (int(maxlen) - k) + o[k // 2:]
        if len(o) != maxlen:
            msg = 'Peptide %s has length %d < maxlen = %d, but pad is "none".'
            raise ValueError(msg % (peptide, len(peptide), maxlen))
        result = ENCODING_DATA_FRAMES[encoding_method].iloc[o]
        return np.asarray(result)

    def antigenMap(self, dataset, maxlen, encoding_method):
        '''Input a list of antigens and get a three dimentional array'''
        m = 0
        for each_antigen in dataset:
            if m == 0:
                antigen_array = self.peptide_encode_HLA(each_antigen, maxlen, encoding_method).reshape(1, maxlen, 21)
            else:
                antigen_array = np.append(antigen_array,
                                          self.peptide_encode_HLA(each_antigen, maxlen, encoding_method).reshape(1, maxlen,
                                                                                                            21), axis=0)
            m = m + 1
        return antigen_array

    def TCRMap(self, dataset, aa_dict, encode_dim):
        # Wrapper of aamapping
        for i in range(0, len(dataset)):
            if i == 0:
                TCR_array = self.aamapping_TCR(dataset[i], aa_dict, encode_dim).reshape(1, encode_dim, 5, 1)
            else:
                TCR_array = np.append(TCR_array,
                                      self.aamapping_TCR(dataset[i], aa_dict, encode_dim).reshape(1, encode_dim, 5, 1),
                                      axis=0)
        return TCR_array
    
    def embed(self, file_dir, encode_dim=80, model_dir="library/h5_file",
               aa_dict_dir="library/Atchley_factors.csv"):
        with open(aa_dict_dir, 'r') as aa:
            aa_reader = csv.reader(aa)
            next(aa_reader, None)
            for rows in aa_reader:
                aa_name = rows[0]
                aa_factor = rows[1:len(rows)]
                aa_dict_atchley[aa_name] = np.asarray(aa_factor, dtype='float')
        # TCR
        TCR_list, antigen_list = self.preprocess(file_dir)
        # read cluster data
        TCR_array = self.TCRMap(TCR_list, aa_dict_atchley, encode_dim)
        TCR_encoder = load_model(model_dir + '/TCR_encoder_30.h5')
        TCR_encoder = Model(TCR_encoder.input, TCR_encoder.layers[-12].output)
        TCR_encoded_result = TCR_encoder.predict(TCR_array)
        TCR_encoded_matrix = pd.DataFrame(data=TCR_encoded_result, index=range(1, len(TCR_list) + 1))

        # Antigen
        antigen_array = self.antigenMap(antigen_list, 15, 'BLOSUM50')

        return TCR_encoded_matrix, antigen_array

if __name__ == "__main__":
    
    encoder = EmbeddingpMTnet()
    TCR_encoded_matrix, antigen_array = encoder.embed("data/testdata_pMTnet.csv")
    print(TCR_encoded_matrix.shape)
    print(antigen_array.shape)

