import csv
import pandas as pd
import numpy as np
import sys
from IPython.display import display
from matplotlib import pyplot as plt
pd.set_option('display.max_columns', None)
np.random.seed(621)

mcpas = pd.read_csv(sys.argv[1], encoding='latin', low_memory=False)
def describe_mcpas(df):
    print('Mouse Samples: ', len(df[df['Species'] == 'Mouse']))
    print('Human Samples: ', len(df[df['Species'] == 'Human']))
    print('Unique Antigen Sequences: ', len(df['Epitope.peptide'].unique()))
    print('Cancer Based Samples: ', len(df[df['Category'] == 'Cancer']))
    print('CD8 T-Cell Samples: ', len(df[df['T.Cell.Type'] == 'CD8']))
    print('CD4 T-Cell Samples: ', len(df[df['T.Cell.Type'] == 'CD4']))
    print('MHC I Samples: ', len(df[df['MHC'].str.match('^HLA-[ABC]')==True]))
    print('MHC II Sample: ', len(df[df['MHC'].str.match('^(HLA-D|DR|D)')==True]))
    print('TCR-Epitope Pairs Available: ', len(df))

# McPAS TCR Beta Sequences
mcpas_beta_pairs = mcpas.loc[(mcpas['CDR3.beta.aa'].notnull()) & (mcpas['Epitope.peptide'].notnull())]
print('McPAS TCR Beta Sequences')
print('Unique CDR3 Beta Sequences: ', len(mcpas_beta_pairs['CDR3.beta.aa'].unique()))
describe_mcpas(mcpas_beta_pairs)

# McPAS TCR Alpha Sequences
mcpas_alpha_pairs = mcpas.loc[(mcpas['CDR3.alpha.aa'].notnull()) & (mcpas['Epitope.peptide'].notnull())]
print('\nMcPAS TCR Alpha Sequences')
print('Unique CDR3 Alpha Sequences: ', len(mcpas_alpha_pairs['CDR3.alpha.aa'].unique()))
describe_mcpas(mcpas_alpha_pairs)

# McPAS TCR Alpha AND Beta Sequences
mcpas_both_pairs = mcpas.loc[(mcpas['CDR3.alpha.aa'].notnull()) & (mcpas['Epitope.peptide'].notnull()) & (mcpas['CDR3.beta.aa'].notnull())]
print('\nMcPAS TCR Both Beta and Alpha Available')
describe_mcpas(mcpas_both_pairs)

vdjdb = pd.read_csv(sys.argv[2], encoding='latin', sep='\t')
def describe_vdjdb(df):
    print('Human Samples: ', len(df[df['species'] == 'HomoSapiens']))
    print('Mouse Samples: ', len(df[df['species'] == 'MusMusculus']))
    print('Monkey Samples: ', len(df[df['species'] == 'MacacaMulatta']))
    print('Unique Antigen Sequences: ', len(df['antigen.epitope'].unique()))
    print('MHC Class I: ', len(df[df['mhc.class'] == 'MHCI']))
    print('MHC Class II: ', len(df[df['mhc.class'] == 'MHCII']))
    print('Total Pairs: ', len(df))
    print('Score 0: ', len(df[df['vdjdb.score'] == 0]))
    print('Score 1: ', len(df[df['vdjdb.score'] == 1]))
    print('Score 2: ', len(df[df['vdjdb.score'] == 2]))
    print('Score 3: ', len(df[df['vdjdb.score'] == 3]))
# VDJDB TCR Beta Sequences
vdjdb_beta_pairs = vdjdb.loc[(vdjdb['gene'] == 'TRB') & (vdjdb['antigen.epitope'].notnull())]
print('VDJDB Beta Sequences')
print('Unique CDR3 Beta Sequences: ', len(vdjdb_beta_pairs['cdr3'].unique()))
describe_vdjdb(vdjdb_beta_pairs)
# VDJDB TCR Alpha Sequences
vdjdb_alpha_pairs = vdjdb.loc[(vdjdb['gene'] == 'TRA') & (vdjdb['antigen.epitope'].notnull())]
print('\nVDJDB Alpha Sequences')
print('Unique CDR3 Alpha Sequences: ', len(vdjdb_alpha_pairs['cdr3'].unique()))
describe_vdjdb(vdjdb_alpha_pairs)
# VDJDB TCR Alpha OR Beta Sequences
print('\nVDJDB Alpha or Beta Sequences')
print('Unique CDR3 Alpha or Beta Sequences: ', len(vdjdb['cdr3'].unique()))
describe_vdjdb(vdjdb)

# Filter down the data to relevant columns
def extract_relevant_columns(df):
    df = df[['Epitope - Name', 'Chain 1 - CDR3 Curated', 'Chain 2 - CDR3 Curated']]
    df.columns = ['Peptide', 'Alpha CDR3', 'Beta CDR3']
    return df

iedb_class_i = pd.read_csv(sys.argv[3], low_memory=False, sep='\t')
print('IEDB MHC Class I')
display(iedb_class_i)

iedb_class_ii = pd.read_csv(sys.argv[4], low_memory=False, sep='\t')
print('IEDB MHC Class II')
display(iedb_class_ii)

iedb_class_i = extract_relevant_columns(iedb_class_i)
iedb_class_ii = extract_relevant_columns(iedb_class_ii)

def summary_stats(df):
    print('Unique Peptides: ', len(df['Peptide'].unique()))
    print('Total Pairs: ', len(df))
    
def filter_alpha(df):
    filtered_data = df.loc[df['Alpha CDR3'].notnull()]
    return filtered_data

def filter_beta(df):
    filtered_data = df.loc[df['Beta CDR3'].notnull()]
    return filtered_data

def filter_both_alpha_beta(df):
    filtered_data = df.loc[df['Alpha CDR3'].notnull() & df['Beta CDR3'].notnull()]
    return filtered_data

# Alpha Analysis
print('CDR3 Alpha Pairs')
print('MHC Class I')
summary_stats(filter_alpha(iedb_class_i))
print('MHC Class II')
summary_stats(filter_alpha(iedb_class_ii))

# Beta Analysis
print('\nCDR3 Beta Pairs')
print('MHC Class I')
summary_stats(filter_beta(iedb_class_i))
print('MHC Class II')
summary_stats(filter_beta(iedb_class_ii))

# Alpha and Beta Analysis
print('\nCDR3 Alpha and Beta Pairs')
print('MHC Class I')
summary_stats(filter_both_alpha_beta(iedb_class_i))
print('MHC Class II')
summary_stats(filter_both_alpha_beta(iedb_class_ii))

def final_filter_mcpas(raw_mcpas):
    # Non Null Epitopes
    df = raw_mcpas.loc[(raw_mcpas['Epitope.peptide'].notnull())]
    # Beta Sequences Only
    df = df[['CDR3.beta.aa', 'Epitope.peptide', 'MHC']]
    # Non Null Beta Sequences
    df = df.loc[(df['CDR3.beta.aa'].notnull())]
    # Remove Sequences with *
    df = df[~df['CDR3.beta.aa'].str.contains('\*')]
    # Add MHC Class Label
    df['MHC'] = df['MHC'].replace(regex={r'^HLA-[ABC].*': 'MHCI', r'^(HLA-D|D).*': 'MHCII'})
    df = df.loc[(df['MHC'] == 'MHCI') | (df['MHC'] == 'MHCII')]
    # Rename Columns
    df.columns = ['CDR3', 'Epitope', 'MHC Class']
    return df

def final_filter_vdjdb(raw_vdjdb, non_zero_scores=True):
    # Beta Sequences Only
    df = raw_vdjdb.loc[(raw_vdjdb['gene'] == 'TRB') & (raw_vdjdb['antigen.epitope'].notnull())]
    df = df[['cdr3', 'antigen.epitope', 'mhc.class', 'vdjdb.score']]
    # VDJDB Score > 0
    if non_zero_scores:
        df = df.loc[df['vdjdb.score'] > 0]
    # Rename Columns
    df.columns = ['CDR3', 'Epitope', 'MHC Class', 'VDJDB Score']
    return df

def final_filter_iedb(raw_iedb_class_i, raw_iedb_class_ii):
    # Add MHC Class Label
    raw_iedb_class_i['MHC Class'] = 'MHCI'
    raw_iedb_class_ii['MHC Class'] = 'MHCII'
    df = pd.concat([raw_iedb_class_i, raw_iedb_class_ii])
    # Beta Sequences Only
    df = df[['Beta CDR3', 'Peptide', 'MHC Class']]
    # Non Null Beta Sequences
    df = df.loc[(df['Beta CDR3'].notnull())]
    # Remove Sequences with #, *, or +
    df = df[~df['Beta CDR3'].str.contains('#')]
    df = df[~df['Beta CDR3'].str.contains('\*')]
    df = df[~df['Peptide'].str.contains('\+')]
    df = df[~df['Peptide'].str.contains('beryllium atom')]
    # Remove Whitespace from TCRs
    df['Beta CDR3'] = df['Beta CDR3'].str.replace(" ","")
    # Capitalize All TCRs
    df['Beta CDR3'] = df['Beta CDR3'].str.upper() 
    # Apply Regex
    df = df[df['Beta CDR3'].str.match(r'^[A|R|N|D|B|C|E|Q|Z|G|H|I|L|K|M|F|P|S|T|W|Y|V]+$')]
    # Rename Columns
    df.columns = ['CDR3', 'Epitope', 'MHC Class']
    return df

def final_filter_combined(combined_data):
    # Filter down to MHC I only
    df = combined_data.loc[(combined_data['MHC Class'] == 'MHCI')]
    # Epitope sequences of length 17 or lower
    df = df.loc[(df['Epitope'].str.len() <= 17)]
    # Positively Binding pairs
    df = df.assign(Affinity=1)
    return df

# Perform filtering for training/testing data
final_mcpas = final_filter_mcpas(mcpas)
final_mcpas = final_mcpas.drop_duplicates(['Epitope', 'CDR3'], keep='last')
final_vdjdb = final_filter_vdjdb(vdjdb)
final_vdjdb = final_vdjdb.drop_duplicates(['Epitope', 'CDR3'], keep='last')
final_iedb = final_filter_iedb(iedb_class_i, iedb_class_ii)
final_iedb = final_iedb.drop_duplicates(['Epitope', 'CDR3'], keep='last')

# Grab VDJDB scores of 0 for added summary results
all_vdjdb = final_filter_vdjdb(vdjdb, non_zero_scores=False)
all_vdjdb = all_vdjdb.drop_duplicates(['Epitope', 'CDR3'], keep='last')

# Combine the datasets together and drop duplicates
full_combined_sequences = pd.concat([final_mcpas, final_vdjdb, final_iedb])
combined_sequences = final_filter_combined(full_combined_sequences)
combined_sequences = combined_sequences.drop_duplicates(['Epitope', 'CDR3'], keep='last')
combined_sequences = combined_sequences.reset_index(drop=True)

def determine_overlap(df1, df2):
    """
    Determines the rows which are both in df1 and df2
    """
    df1_sequences_only = df1[['CDR3', 'Epitope']]
    df1_sequences_only = df1_sequences_only.drop_duplicates()
    df2_sequences_only = df2[['CDR3', 'Epitope']]
    df2_sequences_only = df2_sequences_only.drop_duplicates()
    df = pd.merge(df1_sequences_only, df2_sequences_only, how='inner', on=['CDR3', 'Epitope'])
    return df

def print_overlap_statistic(name, overlap_df, original_df):
    overlap_stat = len(overlap_df)/len(original_df) * 100
    print(f'{name} overlap with IEDB:  {len(overlap_df)}/{len(original_df)} ({str(round(overlap_stat, 2))}%)')

overlap_mcpas = determine_overlap(final_mcpas, final_iedb)
overlap_vdjdb = determine_overlap(final_vdjdb, final_iedb)
overlap_all_vdjdb = determine_overlap(all_vdjdb, final_iedb)
overlap_iedb = determine_overlap(final_iedb, final_iedb)
overlap_mcpas_stat = len(overlap_mcpas)/len(final_mcpas) * 100
overlap_vdjdb_stat = len(overlap_vdjdb)/len(final_vdjdb) * 100
overlap_all_vdjdb_stat = len(overlap_all_vdjdb)/len(all_vdjdb) * 100
overlap_iedb_stat = len(overlap_iedb)/len(final_iedb) * 100
print_overlap_statistic('McPAS', overlap_mcpas, final_mcpas)
print_overlap_statistic('VDJDB (Score > 0)', overlap_vdjdb, final_vdjdb)
print_overlap_statistic('VDJDB (Score >= 0)', overlap_all_vdjdb, all_vdjdb)
print_overlap_statistic('IEDB', overlap_iedb, final_iedb)


def save_distribution(distribution, filename):
    with open(filename, 'w') as csv_file:  
        writer = csv.writer(csv_file)
        for key, value in distribution.items():
            writer.writerow([key, value])

unique_epitopes = combined_sequences['Epitope'].unique()
unique_tcrs = combined_sequences['CDR3'].unique()
epitope_distribution = combined_sequences['Epitope'].value_counts()
tcr_distribution = combined_sequences['CDR3'].value_counts()
#save_distribution(epitope_distribution, 'epitope_distribution.csv')
#save_distribution(tcr_distribution, 'tcr_distribution.csv')

# Draw Pie chart of epitope distribution
major_epitope_distribution = epitope_distribution[epitope_distribution >= 3000]
other_epitope_distribution = epitope_distribution[epitope_distribution < 3000]
added_entry = pd.Series([other_epitope_distribution.sum()], index=['Other Epitopes'])
major_epitope_distribution = major_epitope_distribution.append(added_entry)
pcts = [i / major_epitope_distribution.sum() for i in major_epitope_distribution]
pctlabels = [f'{i} ({j* 100:0.1f}%)' for i, j in zip(major_epitope_distribution.index, pcts)]
plt.style.use('Solarize_Light2')
fig = major_epitope_distribution.plot.pie(figsize=(9, 9), labels=None, ylabel="", labeldistance=1.05, rotatelabels=0, legend=True)
plt.legend(pctlabels, bbox_to_anchor=(1.45,0.30), loc='lower right', prop={'size': 14})
fig.set_facecolor(color='white')

plt.savefig('epitope_dist.png', transparent=True, bbox_inches='tight')


epitope_lengths = combined_sequences[['Epitope', 'CDR3']].copy()
epitope_lengths['Epitope Length'] = epitope_lengths['Epitope'].str.len()
print("Mean Epitope Length:", epitope_lengths['Epitope Length'].mean())
print("Mean Unique Epitope Length:", epitope_lengths.drop_duplicates(subset=['Epitope'])['Epitope Length'].mean())
epitope_length_counts = epitope_lengths[['Epitope Length', 'Epitope']].drop_duplicates(subset=['Epitope']).groupby(['Epitope Length']).count()
epitope_lengths = epitope_lengths.groupby(['Epitope', 'Epitope Length']).count()
epitope_lengths = epitope_lengths.groupby(['Epitope Length']).mean()
epitope_lengths.rename(columns = {"CDR3": "Average Cognate CDR3s"}, inplace = True)
epitope_lengths.plot.bar(figsize=(10, 5), title='Average # of Binding TCRs per Epitope Length')
epitope_length_counts.rename(columns = {"Epitope": "Epitopes"}, inplace = True)
epitope_length_counts.plot.bar(figsize=(10, 5), title='Unique Epitopes per Epitope Length')


def combined_summary_stats(df):
    print('Unique Antigen Sequences: ', len(df['Epitope'].unique()))
    print('Unique CDR3 Sequences: ', len(df['CDR3'].unique()))
    print('MHC Class I: ', len(df[df['MHC Class'] == 'MHCI']))
    print('MHC Class II: ', len(df[df['MHC Class'] == 'MHCII']))
    print('VDJDB Pairs: ', len(final_vdjdb))
    print('McPAS Pairs: ', len(final_mcpas))
    print('IEDB Pairs: ', len(final_iedb))
    print('Total Pairs: ', len(df))

print('Combined VDJDB and McPAS Statistics')
combined_summary_stats(combined_sequences)

# Save metadata file
combined_sequences.to_csv('combined_metadata.csv', index=False)
# Get Relevant Columns for training
final_data = combined_sequences[['Epitope', 'CDR3', 'Affinity']]
# Save to File
final_data.to_csv('combined_dataset_positive_only.csv', index=False, header=False)
# Save Reverse Sequences for testing against ERGO
ergo = final_data[['CDR3', 'Epitope', 'Affinity']]
ergo.to_csv('ergo_training_data.csv', index=False, header=False)

print(len(final_iedb['Epitope'].unique()))

def random_recombination(df, epitope_dist, tcr_dist, ratio):
    unique_epitopes = df['Epitope'].unique()
    unique_tcrs = df['CDR3'].unique()
    conversion_df = df[['Epitope', 'CDR3']]
    positive_pairs = set([tuple(x) for x in conversion_df.to_numpy()])

    # We want to weight the tcr choice by frequency in data
    epitope_freq_array = [epitope_dist[peptide] / len(df) for peptide in unique_epitopes]
    tcr_freq_array = [tcr_dist[tcr] / len(df) for tcr in unique_tcrs]
    
    neg_pairs = set()
    for pep in unique_epitopes:
        i = 0
        pairs_to_generate = round(epitope_dist[pep] * ratio)
        while i < pairs_to_generate:
            tcr = np.random.choice(unique_tcrs, p=tcr_freq_array)
            pair = (pep, tcr)
            if pair not in positive_pairs and pair not in neg_pairs:
                neg_pairs.add(pair)
                i += 1
            
    negative_data = pd.DataFrame(neg_pairs, columns = ['Epitope', 'CDR3'])
    negative_data = negative_data.assign(Affinity=0)
    return negative_data

# Ratio of negative to positive data
ratio = 1
negative_data = random_recombination(final_data, epitope_distribution, tcr_distribution, ratio)
negative_data
full_data = pd.concat([final_data, negative_data])
full_data = full_data.reset_index(drop=True)

full_data.to_csv('combined_dataset.csv', index=False, header=False)