import numpy as np
from sklearn.decomposition import PCA

class EmbeddingSETE:

    def __init__(self):
        pass
    
    def load_data(self, file_path):
        '''
        Set original csv path.

        :param input_file_path: original input file path.
        '''
        self.input_file_path = file_path

    def splitCDR(self, CDRseq, k=3):
        '''
        Note:
            This function is copied from SETE.py in SETE-master.
            See https://github.com/wonanut/SETE for original source code.
        '''
        retDict = []
        for i in range(len(CDRseq) - k + 1):
            retDict.append(CDRseq[i:i + k])
        return retDict

    def statisticsKmer(self, epiDict, k=3):
        '''
        Note:
            This function is copied from SETE.py in SETE-master.
            See https://github.com/wonanut/SETE for original source code.
        '''
        kmerDict = {}
        for epi in epiDict:
            for i in range(len(epiDict[epi])):
                splitList = self.splitCDR(epiDict[epi][i], k)
                for split in splitList:
                    if split not in kmerDict:
                        kmerDict[split] = 1
                    else:
                        kmerDict[split] += 1
        return kmerDict

    def splitCDR(self, CDRseq, k=3):
        '''
        Note:
            This function is copied from SETE.py in SETE-master.
            See https://github.com/wonanut/SETE for original source code.
        '''
        retDict = []
        for i in range(len(CDRseq) - k + 1):
            retDict.append(CDRseq[i:i + k])
        return retDict

    def buildFeatures(self, epiDict, kmerDict, k=3):
        '''
        Note:
            This function is copied from SETE.py in SETE-master.
            See https://github.com/wonanut/SETE for original source code.
        '''

        counter = 0
        for epi in epiDict:
            counter += len(epiDict[epi])
        retArr = np.zeros((counter, len(kmerDict)))

        kmerList = kmerDict.keys()
        retLabel = []

        iter = 0
        epinum = 0
        for epi in epiDict:
            for cdr in range(len(epiDict[epi])):
                splitlist = self.splitCDR(epiDict[epi][cdr], k)
                retLabel.append(epinum)
                i = 0
                for kmer in kmerList:
                    retArr[iter][i] = splitlist.count(kmer)
                    i += 1
                iter += 1
            epinum += 1
        print()
        return np.array(retArr), np.array(retLabel)

    def pca_analyse(self, X_train, rate=0.9):
        '''
        Perform PCA for the train set.

        Note:
            This function is modified from SETE.py in SETE-master.
            See https://github.com/wonanut/SETE for original source code.
        '''
        pca = PCA(n_components=rate).fit(X_train)
        return pca.transform(X_train)
    
    def embed(self, k=3, remove_duplicate=False):
        '''

        Embed cdr3 sequence into 1d numpy ndarray.

        In SETE, the authors performed PCA downscaling of TCR sequences after k-mer coding as X input. 
        For Epitope, no encoding operation was performed, but it was converted to a label to be predicted.

        If you want to encode a specific column utilizing the SETE form, 
            you need to change the CSV header of that column to CDR3 and 
            you need to make sure that one column has an Epitope header.

        :param k:
        :param remove_duplicate: remove duplicate data if True
        :return: embedded cdr3 sequence and kmer dictionary
        '''
        import pandas as pd

        df = pd.read_csv(self.input_file_path)

        if remove_duplicate:
            head_list = df.columns.values.tolist()
            assert 'epitope' in head_list and 'cdr3b' in head_list and 'vb_gene' in head_list
            subset = ['epitope', 'cdr3b', 'vb_gene']
            if 'vb_gene' in head_list:
                df.drop_duplicates(subset=subset, inplace=True)

        epiDict = {}
        for index, row in df.iterrows():
            if row['Epitope'] not in epiDict:
                epiDict[row['Epitope']] = []
            epiDict[row['Epitope']].append(row['CDR3'])

        statistics_epi = []
        statistics_num = []
        # print('{:22s} {:s}'.format('Epitope', 'Number'))
        for epi in epiDict:
            statistics_epi.append(epi)
            statistics_num.append(len(epiDict[epi]))
            # print('{:22s} {:d}'.format(epi, len(epiDict[epi])))

        kmerDict = self.statisticsKmer(epiDict, k)
        X, y = self.buildFeatures(epiDict, kmerDict, k)
        X = self.pca_analyse( X)
        return X, y, kmerDict


if __name__ == '__main__':
    
    # In SETE, the authors performed PCA downscaling of TCR sequences after k-mer coding as X input. 
    # For Epitope, no encoding operation was performed, but it was converted to a label to be predicted.

    # If you want to encode a specific column utilizing the SETE form, 
    # you need to change the CSV header of that column to CDR3 and 
    # you need to make sure that one column has an Epitope header.

    embedding = EmbeddingSETE('data/testdata_SETE.csv')
    X, y, kmerDict = embedding.embed(k=3) # Only X is encoded.
    print(X.shape)