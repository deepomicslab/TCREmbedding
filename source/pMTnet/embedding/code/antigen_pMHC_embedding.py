import numpy as np
import pandas as pd
import csv
from keras.models import Model, load_model
from keras import backend as K
import os
from io import StringIO

########################### Atchley's factors#######################
aa_dict_atchley=dict()
HLA_seq_lib = {}
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

class antigen_pMHC_embedding:
    def __init__(self):
        pass

    def preprocess(self, filedir, paired=True):
        # Preprocess TCR files
        #print('Processing: ' + filedir)
        if not os.path.exists(filedir):
            print('Invalid file path: ' + filedir)
            return 0
        dataset = pd.read_csv(filedir, header=0)
        # Preprocess HLA_antigen files
        # remove HLA which is not in HLA_seq_lib; if HLA*01:01 not in HLA_seq_lib; then the first HLA startswith input HLA allele will be given
        # Remove antigen that is longer than 15aa
        if paired == 'F':
            HLA_antigen = dataset[['HLA', 'Antigen']].dropna()
            HLA_list = list(HLA_antigen['HLA'])
            antigen_list = list(HLA_antigen['Antigen'])
            ind = 0
            index_list = []
            for i in HLA_list:
                if len([hla_allele for hla_allele in HLA_seq_lib.keys() if hla_allele.startswith(str(i))]) == 0:
                    index_list.append(ind)
                ind = ind + 1
            HLA_antigen = HLA_antigen.drop(HLA_antigen.iloc[index_list].index)
            HLA_antigen = HLA_antigen[HLA_antigen.Antigen.str.len() < 16]
            print(str(max(HLA_antigen.index) - HLA_antigen.shape[0]) + ' antigens longer than 15aa are dropped!')
            TCR_list = list(dataset['CDR3'].dropna())
            antigen_list = list(HLA_antigen['Antigen'])
            HLA_list = list(HLA_antigen['HLA'])
        else:
            dataset = dataset.dropna()
            HLA_list = list(dataset['HLA'])
            ind = 0
            index_list = []
            for i in HLA_list:
                if len([hla_allele for hla_allele in HLA_seq_lib.keys() if hla_allele.startswith(str(i))]) == 0:
                    index_list.append(ind)
                    print('drop ' + i)
                ind = ind + 1
            dataset = dataset.drop(dataset.iloc[index_list].index)
            dataset = dataset[dataset.Antigen.str.len() < 16]
            print(str(max(dataset.index) - dataset.shape[0]) + ' antigens longer than 15aa are dropped!')
            TCR_list = dataset['CDR3'].tolist()
            antigen_list = dataset['Antigen'].tolist()
            HLA_list = dataset['HLA'].tolist()
        return TCR_list, antigen_list, HLA_list

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

    def hla_encode(self, HLA_name, encoding_method):
        '''Convert the HLAs of a sample(s) to a zero-padded (for homozygotes)
        numeric representation.

        Parameters
        ----------
            HLA_name: the name of the HLA
            encoding_method:'BLOSUM50' or 'one-hot'
        '''
        if HLA_name not in HLA_seq_lib.keys():
            if len([hla_allele for hla_allele in HLA_seq_lib.keys() if hla_allele.startswith(str(HLA_name))]) == 0:
                print('cannot find' + HLA_name)
            HLA_name = [hla_allele for hla_allele in HLA_seq_lib.keys() if hla_allele.startswith(str(HLA_name))][0]
        if HLA_name not in HLA_seq_lib.keys():
            print('Not proper HLA allele:' + HLA_name)
        HLA_sequence = HLA_seq_lib[HLA_name]
        HLA_int = [aa_dict_one_hot[char] for char in HLA_sequence]
        while len(HLA_int) != 34:
            # if the pseudo sequence length is not 34, use X for padding
            HLA_int.append(20)
        result = ENCODING_DATA_FRAMES[encoding_method].iloc[HLA_int]
        # Get a numpy array of 34 rows and 21 columns
        return np.asarray(result)

    def peptide_encode_HLA(self, peptide, maxlen, encoding_method):
        '''Convert peptide amino acid sequence to one-hot encoding,
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
                                          self.peptide_encode_HLA(each_antigen, maxlen, encoding_method).reshape(1, maxlen, 21), axis=0)
            m = m + 1
        #print('antigenMap done!')
        return antigen_array

    def HLAMap(self, dataset, encoding_method):
        '''Input a list of HLA and get a three dimentional array'''
        m = 0
        for each_HLA in dataset:
            if m == 0:
                HLA_array = self.hla_encode(each_HLA, encoding_method).reshape(1, 34, 21)
            else:
                HLA_array = np.append(HLA_array, self.hla_encode(each_HLA, encoding_method).reshape(1, 34, 21), axis=0)
            m = m + 1
        #print('HLAMap done!')
        return HLA_array

    def TCRMap(self, dataset, aa_dict, encode_dim):
        # Wrapper of aamapping
        for i in range(0, len(dataset)):
            if i == 0:
                TCR_array = self.aamapping_TCR(dataset[i], aa_dict, encode_dim).reshape(1, encode_dim, 5, 1)
            else:
                TCR_array = np.append(TCR_array,
                                      self.aamapping_TCR(dataset[i], aa_dict, encode_dim).reshape(1, encode_dim, 5, 1),
                                      axis=0)
        #print('TCRMap done!')
        return TCR_array

    def pearson_correlation_f(self, y_true, y_pred):
        fsp = y_pred - K.mean(
            y_pred)  # being K.mean a scalar here, it will be automatically subtracted from all elements in y_pred
        fst = y_true - K.mean(y_true)
        devP = K.std(y_pred)
        devT = K.std(y_true)
        return K.mean(fsp * fst) / (devP * devT)

    def encode(self, file_dir, encode_dim=80, model_dir=r"D:\TCR\program\encode\pMTnet\library\h5_file",
               aa_dict_dir=r"D:\TCR\program\encode\pMTnet\library\Atchley_factors.csv",
               hla_db_dir=r"D:\TCR\program\encode\pMTnet\library\hla_library"):

        # aa_dict_atchley
        with open(aa_dict_dir, 'r') as aa:
            aa_reader = csv.reader(aa)
            next(aa_reader, None)
            for rows in aa_reader:
                aa_name = rows[0]
                aa_factor = rows[1:len(rows)]
                aa_dict_atchley[aa_name] = np.asarray(aa_factor, dtype='float')

        # HLA seq lib
        HLA_ABC = [hla_db_dir + '/A_prot.fasta', hla_db_dir + '/B_prot.fasta', hla_db_dir + '/C_prot.fasta', hla_db_dir + '/E_prot.fasta']
        for one_class in HLA_ABC:
            prot = open(one_class)
            # pseudo_seq from netMHCpan:https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0000796
            pseudo_seq_pos = [7, 9, 24, 45, 59, 62, 63, 66, 67, 79, 70, 73, 74, 76, 77, 80, 81, 84, 95, 97, 99, 114,
                              116, 118, 143, 147, 150, 152, 156, 158, 159, 163, 167, 171]
            # write HLA sequences into a library
            # class I alles
            name = ''
            sequence = ''
            for line in prot:
                if len(name) != 0:
                    if line.startswith('>HLA'):
                        pseudo = ''
                        for i in range(0, 33):
                            if len(sequence) > pseudo_seq_pos[i]:
                                pseudo = pseudo + sequence[pseudo_seq_pos[i]]
                        HLA_seq_lib[name] = pseudo
                        name = line.split(' ')[1]
                        sequence = ''
                    else:
                        sequence = sequence + line.strip()
                else:
                    name = line.split(' ')[1]

        # TCR
        TCR_list, antigen_list, HLA_list=self.preprocess(file_dir)
        TCR_array = self.TCRMap(TCR_list, aa_dict_atchley, encode_dim)

        TCR_encoder = load_model(model_dir + '/TCR_encoder_30.h5')
        TCR_encoder = Model(TCR_encoder.input, TCR_encoder.layers[-12].output)
        TCR_encoded_result = TCR_encoder.predict(TCR_array)
        TCR_encoded_matrix = pd.DataFrame(data=TCR_encoded_result, index=range(1, len(TCR_list) + 1))

        # Antigen
        antigen_array = self.antigenMap(antigen_list, 15, 'BLOSUM50')

        # HLA
        HLA_array = self.HLAMap(HLA_list, 'BLOSUM50')
        HLA_antigen_encoder = load_model(model_dir + '/HLA_antigen_encoder_60.h5',
                                         custom_objects={'pearson_correlation_f': self.pearson_correlation_f})
        HLA_antigen_encoder = Model(HLA_antigen_encoder.input, HLA_antigen_encoder.layers[-2].output)
        HLA_antigen_encoded_result = HLA_antigen_encoder.predict([antigen_array, HLA_array])
        HLA_antigen_encoded_matrix = pd.DataFrame(data=HLA_antigen_encoded_result, index=range(1, len(HLA_list) + 1))

        return TCR_encoded_matrix, antigen_array, HLA_antigen_encoded_matrix

if __name__ == "__main__":
    encoder = antigen_pMHC_embedding()
    TCR_encoded_matrix, antigen_array, HLA_antigen_encoded_matrix = encoder.encode(
        r"D:\TCR\program\encode\pMTnet\embedding\input\test_input.csv")
    print(HLA_antigen_encoded_matrix)
