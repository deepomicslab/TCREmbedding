import numpy as np
import pandas as pd
from collections import defaultdict

class Sampler:
    '''
    The Sampler class implements of the negative sampling methods.
    This class is from TEINet project by Yuepeng Jiang et al. 

    reference:
        - Article: TEINet: a deep learning framework for prediction of TCR–epitope binding specificity
        - Authors: Yuepeng Jiang, Miaozhe Huo, Shuai Cheng Li
        - DOI link: https://doi.org/10.1093/bib/bbad086
        - GitHub link: https://github.com/jiangdada1221/TEINet

    '''
    def __init__(self):
       pass
    
    def construct_neg(self,path,static=False):
        '''
        Mainly to compute useful dicts for sampling. e.g. self.e2c records the interacting TCRs for epitopes
        '''
        d = pd.read_csv(path)
        self.epitopes = d['Epitope'].values #store epitopes
        es = set(self.epitopes)        
        self.tcrs = d['CDR3'].values #store tcrs
        ts = set(self.tcrs)
        self.t2i,self.i2t = dict(),dict()
        if static:
            self.labels = d['Affinity'].values #store the label for static dataset
        for i,t in enumerate(set(self.tcrs)):
            self.t2i[t] = i  
            self.i2t[i] = t
        self.e2i,self.i2e = dict(),dict()
        for i,e in enumerate(set(self.epitopes)):
            self.e2i[e] = i  
            self.i2e[i] = e
        t2e_pos = defaultdict(set) 
        e2c = defaultdict(int)       
        for i in range(len(self.tcrs)):
            t2e_pos[self.tcrs[i]].add(self.epitopes[i])
            e2c[self.epitopes[i]] += 1
        t2e_neg = dict()
        for t in t2e_pos.keys():
            t2e_neg[t] = list(es - t2e_pos[t])
        self.t2e_neg = t2e_neg
        self.e2c = e2c

        e2t = defaultdict(set)
        for i in range(len(self.tcrs)):
            e2t[self.epitopes[i]].add(self.tcrs[i])
        for e in e2t.keys():
            e2t[e] = np.array(list(ts - e2t[e]))
        self.e2t = e2t
                     
    def e2f(self,candidates):
        c = [self.e2c[i] for i in candidates]
        sum_ = np.sum(c)
        return [i / sum_ for i in c]    

    def sample_neg_uniform(self,cdrs,epitopes):
        neg_epitopes = []
        neg_cdrs = []        
        for c in cdrs:
            neg_epitopes.extend(np.random.choice(self.t2e_neg[c],1))
            neg_cdrs.append(c)
        return (list(cdrs), list(epitopes),[1] * len(cdrs)), (list(neg_cdrs),list(neg_epitopes),[0] * len(neg_cdrs))
        # return list(cdrs)+list(cdrs), list(epitopes) + list(neg_epitopes), torch.LongTensor([1]*len(cdrs) + [0]*len(neg_epitopes))
    
    def sample_neg_whole(self,cdrs,epitopes,num=64):
        neg_epitopes = []
        neg_cdrs = [] 
        epis = []       
        for i,c in enumerate(cdrs):
            # epis.append(epitopes[i])
            neg_epitopes.extend(np.random.choice(self.t2e_neg[c],num,replace=False))
            neg_cdrs.extend([c] * num)
        # return (list(cdrs), list(epitopes),[1] * len(cdrs)), (list(neg_cdrs),list(neg_epitopes),[0] * len(neg_cdrs))
        return list(cdrs)+list(neg_cdrs), list(epitopes)+list(neg_epitopes), [1]*len(cdrs) + [0]*len(neg_cdrs)
    
    def sample_neg_fre(self,cdrs,epitopes,num=64):
        neg_epitopes = []
        neg_cdrs = []        
        for c in cdrs:
            neg_epitopes.extend(np.random.choice(self.t2e_neg[c],num,replace=False,p=self.e2f(self.t2e_neg[c])))
            neg_cdrs.extend([c]*num)
        # return (list(cdrs), list(epitopes),[1] * len(cdrs)), (list(neg_cdrs),list(neg_epitopes),[0] * len(neg_cdrs))
        return list(cdrs)+list(neg_cdrs), list(epitopes)+list(neg_epitopes), [1]*len(cdrs) + [0]*len(neg_cdrs)

    def sample_neg_tcr(self,cdrs,epitopes,num=64,reference=False):
        neg_epitopes = []
        neg_cdrs = []           
        for e in epitopes:
            if not reference:                              
                neg_cdrs.extend(np.random.choice(self.e2t[e],num,replace=False))                
            else :                               
                neg_cdrs.extend(np.random.choice(self.reference_tcr,num,replace=False))                                
            neg_epitopes.extend([e]*num)        
        # return (list(cdrs), list(epitopes),[1] * len(cdrs)), (list(neg_cdrs),list(neg_epitopes),[0] * len(neg_cdrs))        
        return list(cdrs)+list(neg_cdrs), list(epitopes)+list(neg_epitopes), [1]*len(cdrs) + [0]*len(neg_cdrs)


if __name__ == '__main__':
    # Create an instance of the Sampler class
    sampler = Sampler()

    # Specify the path to the input CSV file
    input_path = 'combined_dataset_filtered.csv'

    # Initialize necessary data
    sampler.construct_neg(input_path)

    # Generate negative samples
    _, neg_result = sampler.sample_neg_uniform(sampler.tcrs,sampler.epitopes)

    # Unpack the result into separate variables
    neg_cdrs, neg_epitopes, neg_affinity = neg_result
    # Create a dictionary to store the data
    data = {'Epitope': neg_epitopes, 'CDR3': neg_cdrs, 'Affinity': neg_affinity}
    # Convert the data into a DataFrame
    df = pd.DataFrame(data)

    # Save the DataFrame to a CSV file
    df.to_csv('neg_result.csv', index=False)