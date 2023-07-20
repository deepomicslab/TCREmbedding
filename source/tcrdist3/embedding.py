import concurrent.futures
import pandas as pd
from sklearn.manifold import MDS
import numpy as np

from tcrdist.repertoire import TCRrep

class EmbeddingTCRdist:
    def __init__(self):
        self.tcrdist_model = None

    def read_csv(self, file_path):
        """Read TCR sequences from a CSV file.

        Args:
            file_path (str): The path to the CSV file containing TCR sequences.
            column_name (str): Column name of the provided file recording TCRs. Defaults to 'full_seq'.
        """
        self.cell_df = pd.read_csv(file_path, sep="\t", header=0)
        self.cell_df.rename(columns={'CDR3b': 'cdr3_b_aa', 'Vgene': 'v_b_gene', 'TCR_num': 'count'}, inplace=True)

    def load_model(self, cell_df=None,organism='human',chains=['beta'],db_file='alphabeta_gammadelta_db.tsv',clone_df=None,
                    imgt_aligned=True,infer_all_genes=True,infer_cdrs=True,infer_index_cols=True,deduplicate=True,use_defaults=True,
                   store_all_cdr=True,compute_distances=False,index_cols=['index'],cpus=2
                   ):
        '''
        Get an instance object of TCRrep.
        Flexible distance measures for comparing T cell receptors

        The TCRrep Class hold T cell repertoire data, infers CDRs from v-gene name, and computes multi-CDR 'tcrdistance'.

        Attributes
        ----------
        cell_df : pd.DataFrame or None
            Pandas DataFrame containing cell level information
        clone_df : pd.DataFrame or None
            Pandas DataFrame containing clone level information.
            This can be provided directly from a program like MIXCR or
            can be inferred by deduplicating a cell_df.
        organism = str,
            specifies relevant organism for analysis: 'human' or 'mouse'
        chains : list
            specifies relevant chains for single or paried TCR analysis
            ['alpha','beta'], ['alpha'], ['beta'], ['gamma','delta'],  ['gamma'] or ['delta']
        db_file : str
            specifies reference file. The default is 'alphabeta_gammadelta_db.tsv' which
            is preinstalled  with the install python3.7/site-packages/tcrdist/db/
        archive_name : str
            Name for archive file. (only used if archive result is True)
        archive_result: bool
            If True save result to .tar.gz archive
        imgt_aligned : bool
            If True, by default, cdr1, cdr2,and pmhc are inferred aligned to fixed length with gaps.
            If False, cdr1, cdr2,and pmhc are returned as ungapped raw sequences.
        infer_all_genes : bool
            If True, load all_gene reference from 'db_file`.
        infer_cdrs : bool
            If True, infer cdr1, cdr2,and pmhc from the v gene name
        infer_index_cols : bool
            If True, infer index_cols used to deduplicate cell_df.
            If False, index_cols can be specified directly after initialization.
        index_cols : list
            list of index columns used to deduplicate cell_df to clone_df
        deduplicate : bool
            If True, then clone_df will be assigned cell_df grouped_by
            index_cols.
            If False, and clone_df is None, then clone_df will be be
            assigned a copy of cell_df.
        use_defaults : True
            If True, use default metrics, weights, and kargs
        store_all_cdr : True,
            If True, save all distance matrices for each CDR (e.g., pw_cdr3_b_aa).
            If False, only save pw_alpha and pw_beta
        compute_distances : True
            If True, automatically compute distances
        cpus : int,
            Number of cpus to use. In general 1 cpu is sufficient from default Numba metric
            with less than 10^7  pairwise comparisons. However, more cpus will
            result in a speedup for metrics like pw.metrics.nw_hamming_metric for more than 10^6
            pairwise comparisons.
        '''
        if cell_df == None:
            cell_df = self.cell_df

        self.tcrdist_model = TCRrep(organism, chains, db_file, archive_name='tcrdist3.archive', blank=False,
                                    cell_df=cell_df, clone_df=clone_df,imgt_aligned=imgt_aligned,infer_all_genes=infer_all_genes,
                                    infer_cdrs=infer_cdrs,infer_index_cols=infer_index_cols,deduplicate=deduplicate,use_defaults=use_defaults,
                                    store_all_cdr=store_all_cdr,compute_distances=compute_distances,index_cols=index_cols,cpus=cpus)

    def encode(self):
        tr = self.tcrdist_model
        if tr.compute_distances:
            # This is a safety, measure so that a new user doesn't accidently try to compute a pairwise matrix that won't fit in memory
            if tr.clone_df.shape[0] > 10000:
                with concurrent.futures.ThreadPoolExecutor() as executor:
                    # Submit the compute_sparse_rect_distances function to the ThreadPoolExecutor
                    future = executor.submit(tr.compute_sparse_rect_distances, radius=50, chunk_size=100)
                    distance_matrix = future.result()
                distance_matrix = tr.rw_beta
            else:
                tr.compute_distances()
                distance_matrix = tr.pw_beta

        # Isometric embedding
        embedding = MDS(n_components=96, n_init=100, max_iter=1000, eps=0.00001, dissimilarity='precomputed')
        encode_result = embedding.fit_transform(distance_matrix)
        return encode_result

if __name__ == "__main__":
    encoder = EmbeddingTCRdist()
    encoder.read_csv("TCRantigenData_unique_filt.tsv")
    encoder.load_model()
    encode_result = encoder.encode()
    print(encode_result.shape)
    #np.save("TCRdist_tcr.npy", encode_result)

