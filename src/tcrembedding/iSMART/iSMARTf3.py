import numpy as np
from sklearn.manifold import MDS

from substitution_matrices import BLOSUM62

class iSMART_model:
    def __init__(self):
        self.AAstring = 'ACDEFGHIKLMNPQRSTVWY'
        self.AAstringList = list(self.AAstring)

    def InsertGap(self, Seq, n):
        ## Insert n gaps to Seq; n<=2
        if n == 0:
            return [Seq]
        ns = len(Seq)
        SeqList = []
        if n == 1:
            for kk in range(0, ns+1):
                SeqNew = Seq[0:kk]+'-'+Seq[kk:]
                SeqList.append(SeqNew)
        if n == 2:
            for kk in range(0, ns+1):
                SeqNew = Seq[0:kk]+'-'+Seq[kk:]
                for jj in range(0, ns+2):
                    SeqNew0 = SeqNew[0:jj]+'-'+SeqNew[jj:]
                    SeqList.append(SeqNew0)
        return SeqList

    def SeqComparison(self, s1, s2, gap=-6):
        n = len(s1)
        score = 0
        for kk in range(0, n):
            aa = s1[kk]
            bb = s2[kk]
            if aa in ['.','-','*'] or bb in ['.','-','*']:
                if aa != bb:
                    score += gap
                continue
            if aa == bb:
                score += min(4, BLOSUM62[(aa, aa)])
                continue
            KEY = (aa, bb)
            if KEY not in BLOSUM62:
                KEY = (bb, aa)
            if KEY not in BLOSUM62:
                raise "Non-standard amino acid coding!"
            score += BLOSUM62[KEY]
        return score

    def NHLocalAlignment(self, Seq1, Seq2, gap_thr=1, gap=-6):
        n1 = len(Seq1)
        n2 = len(Seq2)
        if n1 < n2:
            Seq = Seq1
            Seq1 = Seq2
            Seq2 = Seq
            nn = n2 - n1
        else:
            nn = n1 - n2
        if nn > gap_thr:
            return -1
        SeqList1 = [Seq1]
        SeqList2 = self.InsertGap(Seq2, nn)
        alns = []
        SCOREList = []
        for s1 in SeqList1:
            for s2 in SeqList2:
                SCOREList.append(self.SeqComparison(s1, s2, gap))
        maxS = max(SCOREList)
        return maxS

    def falign_embed(self, s1, s2, st=3, gapn=1, gap=-6):
        mid1 = s1[st:-2]
        mid2 = s2[st:-2]
        aln = self.NHLocalAlignment(mid1, mid2, gapn, gap)
        score = aln / float(max(len(mid1), len(mid2))) + 4.0
        return score

    def parseinput(self, data):
        CDR3s = []
        for value in data:
            if value.startswith('C') and value.endswith('F'):
                flag = True
                for i in list(value):
                    if i not in ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']:
                        flag = False
                        break
                if flag:
                    CDR3s.append(value)
                else:
                    continue
        unique_set = set(CDR3s)
        unique_list = list(unique_set)
        return unique_list

    def encode(self, data):
        if len(data) == 0:
            print("There is no data")
        else:
            seqs = self.parseinput(data)
            matrix = np.zeros((len(seqs), len(seqs)))
            for i in range(len(seqs)):
                for j in range(i+1, len(seqs)):
                    aln = self.falign_embed(seqs[i], seqs[j])
                    matrix[i, j] = aln
                    matrix[j, i] = aln
        # Isometric embedding
        embedding = MDS(n_components=96, n_init=100, max_iter=1000, eps=0.00001, dissimilarity='precomputed')
        encode_result = embedding.fit_transform(matrix)

        return encode_result
