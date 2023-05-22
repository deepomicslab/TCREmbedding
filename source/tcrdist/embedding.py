import os
from tcrdist.repertoire import TCRrep
from input.vdjdb import parse_vdjdb

def embedding(data=None, chain='beta', sparse=False):
    '''
    Parameters
    ----------
    data : pandas.DataFrame, optional
        Input dataframe containing information about
        CDR3 sequence, V gene and antigen specificity. The default is None.
    chain : str, optional
        TCR chain: alpha, beta or paired. The default is 'beta'.
    sparse : Bool, optional
        Turn on sparse distance computation. The default is False.

    Returns
    -------
    S : numpy.array
        TCRDist distance matrix.
    seq : pandas.Series
        pandas.Series with sequences for which distances have been calculated.
    gt : pandas.DataFrame
        Ground truth. pandas.DataFrame containing information about the
        TCR sequence and its cognate epitope target.

    '''
    if data is None:
        vdjdb = parse_vdjdb(os.path.abspath('input/vdjdb/vdjdb_full.txt'), q=1)
    else:
        vdjdb = data

    if chain == 'beta':
        cdr3 = 'cdr3_b_aa'
        v_name = 'v_b_gene'
        vdjdb = vdjdb.drop(columns=['cdr3.alpha', 'v.alpha'])
        vdjdb = vdjdb.rename(columns={'cdr3.beta': cdr3,
                                      'v.beta': v_name})
    elif chain == 'alpha':
        cdr3 = 'cdr3_a_aa'
        v_name = 'v_a_gene'
        vdjdb = vdjdb.drop(columns=['cdr3.beta', 'v.beta'])
        vdjdb = vdjdb.rename(columns={'cdr3.alpha': cdr3,
                                      'v.alpha': v_name})

    df_epi = vdjdb[[cdr3, v_name, 'antigen.epitope']].dropna().drop_duplicates()
    seq = df_epi.drop(columns=['antigen.epitope']).drop_duplicates().reset_index(drop=True)
    gt = df_epi.rename(columns={cdr3: 'CDR3',
                                v_name: 'V',
                                'antigen.epitope': 'Epitope'})
    if sparse:

        tr = TCRrep(cell_df=seq,
                    organism='human',
                    chains=['beta'],
                    db_file='alphabeta_gammadelta_db.tsv',
                    compute_distances=False)

        tr.cpus = 2
        tr.compute_sparse_rect_distances(radius=200, chunk_size=500)
        S = tr.rw_beta

    else:

        tr = TCRrep(cell_df=seq,
                    organism='human',
                    chains=[chain],
                    db_file='alphabeta_gammadelta_db.tsv',
                    compute_distances=True)

        S = tr.pw_cdr3_b_aa

    return S, seq, gt
