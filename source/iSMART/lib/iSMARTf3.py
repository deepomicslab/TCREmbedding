#! usr/bin/python3
# -*- coding: utf-8 -*-
## Ultra-fast pairwise alignment algorithm to analyze up to 10^8 CDR3s

import sys, os, re
import numpy as np
from .substitution_matrices import BLOSUM62
from operator import itemgetter
from itertools import chain
import pandas as pd

AAstring='ACDEFGHIKLMNPQRSTVWY'
AAstringList=list(AAstring)
cur_dir=os.path.dirname(os.path.realpath(__file__))+'/'

class CDR3:
    def __init__(self, s, sID, KS, st, ed):
        ## initialize with an input sequence
        ## s: input CDR3 sequence starting from C and ending with the first F in FGXG
        ## sID: unique identifier (increasing integers) given to each CDR3 sequence. Even identical CDR3s should have distinct sIDs
        ## KS: Kmer size
        ## st: the first 0:(st-1) amino acids will not be included in K-merization
        ## ed: the last L-ed amino acids will be skipped
        self.s=s
        self.ID=sID
        L=len(s)
        self.L=L
        sub_s=s[st: (L-ed)]
        Ls=len(sub_s)
        Kmer=[sub_s[x:(x+KS)] for x in range(0,Ls-KS+1)]
        self.Kmer=Kmer

class KmerSet:
    ## Kmer set for fast read searching based on mismatch-allowed Kmer index
    # N 序列数量
    # S
    def __init__(self, Seqs, sIDs, KS, st, ed):
        ## initialize with a list of CDR3s, parse each CDR3 into Kmers, build Kmer-sID dictionary
        ## Seqs and sIDs must have the same length
        if len(Seqs) != len(sIDs):
            raise "Sequence and ID lists have different length. Please check input."
        KmerDict={}
        N=len(Seqs)
        self.N=N
        CDR3Dict={}
        LLs=[]
        for ii in range(0,N):
            s=Seqs[ii]
            sID=sIDs[ii]
            cc=CDR3(s,sID,KS,st,ed)
            CDR3Dict[cc.ID]=cc.Kmer
            KK=cc.Kmer
            LLs.append(cc.L)
            for kk in KK:
                if kk not in KmerDict:
                    KmerDict[kk]=[sID]
                else:
                    KmerDict[kk].append(sID)
        self.KD=KmerDict
        self.KS=KS
        self.CD=CDR3Dict
        self.LL=LLs

    def FindKmerNeighbor(self,kk):
        KS=self.KS
        KS_n1=[]
        for jj in range(KS):
            kk_pre=[kk[0:jj]]*20
            kk_suf=[kk[(jj+1):KS]]*20
            kkn=list(zip(kk_pre,AAstringList,kk_suf))
            KS_n1+=[''.join(list(x)) for x in kkn]
        return KS_n1

    def FindKmerNeighbor2(self,kk):
        ## KS≥6, allowing 2 mismatches. CDR3 length must be ≥ 10
        KS=self.KS
        KS_n1=[]
        for jj in range(KS):
            for ii in range(KS):
                if ii<=jj:
                    continue
                kk_pre=[kk[0:jj]]*20
                kk_mid=[kk[(jj+1):ii]]*20
                kk_suf=[kk[(ii+1):KS]]*400
                kkn=list(zip(kk_pre,AAstringList,kk_mid))
                kkn=[''.join(list(x)) for x in kkn]
                kkn=[[x]*20 for x in kkn]
                kkn=list(chain(*kkn))
                kkn2=list(zip(kkn, AAstringList*20, kk_suf))
                kkn2=[''.join(list(x)) for x in kkn2]
                KS_n1+=kkn2
        return KS_n1

    def KmerIndex(self):
        ## For each K-mer, find its nearest neighbor with 1 character mismatch
        KKs=list(self.KD.keys())
        KS=self.KS
        KKs_set=set(KKs)
        Skk='_'.join(KKs)
        KI_Dict={}
        for kk in KKs:
            KS_n=set(self.FindKmerNeighbor(kk))
            kk_neighbor = KS_n & KKs_set
            KI_Dict[kk]=list(kk_neighbor)
        return KI_Dict

    def updateKD(self, KI):
        ## group sequences sharing motifs with 1-2 mismatches
        KD=self.KD
        KDnew={}
        for kk in KD:
            kkm=KI[kk]
            vvL=itemgetter(*kkm)(KD)
            if isinstance(vvL[0],list):
                vvL=list(chain(*vvL))
            KDnew[kk]=vvL
        return KDnew

def InsertGap(Seq,n):
    ## Insert n gaps to Seq; n<=2
    if n==0:
        return [Seq]
    ns=len(Seq)
    SeqList=[]
    if(n==1):
        for kk in range(0,ns+1):
            SeqNew=Seq[0:kk]+'-'+Seq[kk:]
            SeqList.append(SeqNew)
    if(n==2):
        for kk in range(0,ns+1):
            SeqNew=Seq[0:kk]+'-'+Seq[kk:]
            for jj in range(0,ns+2):
                SeqNew0=SeqNew[0:jj]+'-'+SeqNew[jj:]
                SeqList.append(SeqNew0)
    return SeqList

def SeqComparison(s1,s2,gap=-6):
    n=len(s1)
    CorList=[]
    score=0
    for kk in range(0,n):
        aa=s1[kk]
        bb=s2[kk]
        if aa in ['.','-','*'] or bb in ['.','-','*']:
            if aa!=bb:
                score += gap
            continue
        if aa==bb:
            score += min(4,BLOSUM62[(aa,aa)])
            continue
        KEY=(aa,bb)
        if KEY not in BLOSUM62:
            KEY=(bb,aa)
        if KEY not in BLOSUM62:
            raise "Non-standard amino acid coding!"
        score+=BLOSUM62[KEY]
    return score

def NHLocalAlignment(Seq1,Seq2,gap_thr=1,gap=-6):
    n1=len(Seq1)
    n2=len(Seq2)
    if n1<n2:
        Seq=Seq1
        Seq1=Seq2
        Seq2=Seq
        nn=n2-n1
    else:
        nn=n1-n2
    if nn>gap_thr:
        return -1
    SeqList1=[Seq1]
    SeqList2=InsertGap(Seq2,nn)
    alns=[]
    SCOREList=[]
    for s1 in SeqList1:
        for s2 in SeqList2:
                SCOREList.append(SeqComparison(s1,s2,gap))
    maxS=max(SCOREList)
    return maxS

def falign_embed(s1, s2, st=3, gapn=1, gap=-6):
    mid1=s1[st:-2]
    mid2=s2[st:-2]
    aln=NHLocalAlignment(mid1,mid2,gapn,gap)
    return aln

def parseinput(data):
    CDR3s = []
    for index, value in data.items():
        if value.startswith('C') and value.endswith('F'):
            for i in list(value):
                if i not in ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']:
                    flag = False
                    break
            if flag:
                CDR3s.append(value)
            else:
                continue
    return CDR3s

def embedding(data):
    if data.empty:
        print("There is no data")
    else:
        seqs = parseinput(data)
        matrix = np.zeros((len(seqs), len(seqs)))
        for i in range(len(seqs)):
            for j in range(len(seqs)):
                aln = falign_embed(seqs[i], seqs[j])
                matrix[i, j] = aln
    return matrix
