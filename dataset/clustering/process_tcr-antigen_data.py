import pandas as pd
import re

    ### these antigen names has alias, we need to merge them
    # 'CMV': ['CMV', 'Human betaherpesvirus 5'],
    # 'EBV': ['EBV', 'Human gammaherpesvirus 4'],
    # 'HIV': ['HIV', 'Human immunodeficiency virus 1'],
    # 'Influenza A': ['Influenza A virus', 'InfluenzaA'],
    # 'HSV': ['HSV', 'Human alphaherpesvirus 2'],
    # 'HCV': ['HCV', 'Hepacivirus C'],
    # 'Homo sapiens': ['HomoSapiens', 'Homo sapiens']
    # 'YFV': ['YFV', 'Yellow fever virus']
    # etc

rename_antigen_dict = {
    'Human betaherpesvirus 5': 'CMV',
    'Human herpesvirus 5 strain AD169': 'CMV AD169',
    'Human gammaherpesvirus 4': 'EBV',
    'Hepatitis B virus': 'HBV',
    'Human immunodeficiency virus 1': 'HIV',
    'InfluenzaA': 'Influenza A virus',
    'Human alphaherpesvirus 2': 'HSV',
    'Hepacivirus C': 'HCV',
    'Yellow fever virus': 'YFV',
    'HomoSapiens': 'Homo sapiens',
    'Severe acute respiratory syndrome coronavirus 2': 'SARS coronavirus 2',
    'Severe acute respiratory syndrome': 'SARS',
    'Human T cell leukemia virus type I': 'HLTV1',
}

# Read data and choose required columns CDR3b, Vgene and Antigen-epitope
df = pd.read_csv('TCRantigenData_unique_raw.tsv', sep='\t', header=None)
df.columns = ['CDR3b', 'Vgene', '-', '-', 'Antigen-epitope']
df = df[['CDR3b', 'Vgene', 'Antigen-epitope']]
df['Antigen-epitope'] = df['Antigen-epitope'].str.replace('Antigen:HIV-1', 'Antigen:HIV1')
df['Antigen-epitope'] = df['Antigen-epitope'].str.replace('Antigen:HTLV-1', 'Antigen:HTLV1')
df['Antigen-epitope'] = df['Antigen-epitope'].str.replace('Antigen:Human T-cell', 'Antigen:Human T cell')
df['Antigen-epitope'] = df['Antigen-epitope'].str.replace('H1N1 subtype', 'Influenza A virus (A/California/07/2009(H1N1))')

# Construct regex pattern to match and replace Antigen-epitope column
pattern = re.compile('Antigen:(.*?)[-_](.*)')

# In Antigen-epitope column, check if each value contains any key of rename_antigen_dict
# If so, use re.sub method to replace with corresponding value, else retain original value
df['Antigen-epitope'] = df['Antigen-epitope'].apply(lambda x: re.sub(pattern, 'Antigen:'+rename_antigen_dict[x.split('Antigen:')[1].split('-')[0].split('_')[0]], x)
                                               if any(key in x for key in rename_antigen_dict.keys())
                                               else x)

# Drop duplicates and reset index
df = df.drop_duplicates().reset_index()

# Extract Antigen from Antigen-epitope column and label unknown values as Unknown
df['Antigen'] = df['Antigen-epitope'].apply(lambda x: x.split('Antigen:')[1].split('-')[0].split('_')[0])
for i in range(len(df)):
    if df.loc[i, 'Antigen-epitope'].startswith('Antigen:_'):
        df.loc[i, 'Antigen'] = 'Unknown'

# Group and choose groups with size ≥ 20
grouped = df.groupby('Antigen')
grouped = grouped.filter(lambda x: len(x) >= 20)  ## 629 from 61106 data is removed here

# Reset index and calculate number of Epitope types for each group, add to df
df = grouped.reset_index(drop=True)
epitope_nums = df.groupby('Antigen')['Antigen-epitope'].nunique()
epitope_nums = epitope_nums.reset_index()
epitope_nums.columns = ['Antigen', 'Epitope_type_num']
epitope_dict = epitope_nums.set_index('Antigen')['Epitope_type_num'].to_dict()
df['Epitope_type_num'] = df['Antigen'].map(epitope_dict)

# Calculate number of TCRs for each group and add to df
group_sizes = df.groupby('Antigen').size()
group_sizes = group_sizes.reset_index()
group_sizes.columns = ['Antigen', 'TCR_num']
tcr_dict = group_sizes.set_index('Antigen')['TCR_num'].to_dict()
df['TCR_num'] = df['Antigen'].map(tcr_dict)

# Output df and group sizes sorted by descending TCR_num
print(df)
group_sizes = group_sizes.sort_values('TCR_num', ascending=False)
print(group_sizes.head(10))

# Output results to file
df.to_csv('TCRantigenData_unique_filt.tsv', sep='\t', index=False)
