import os
import numpy as np
import pandas as pd
import umap
import matplotlib.pyplot as plt
import seaborn as sns
import pickle


ncols = 4  # Define the number of columns in the figure
spot_size = 5
metadata_path = './data/IEDB_uniqueTCR_top10_filt.csv'
embdir = './data/IEDB_uniqueTCR_top10_filt'

def vis_umap(testfiles, selected_idx, labels, ncols):
    testdata_dict = {}
    for file in testfiles:
        model_name = file.split('_')[0]
        ## load pkl
        data = pickle.load(open(os.path.join(embdir, file), 'rb'))
        # print("model_name",model_name)
        data = data[selected_idx]
        if len(data.shape) == 2:
            testdata_dict[model_name] = data

    for variable_name, data in testdata_dict.items():
        print(f'{variable_name}: {data.shape}')

    fig, axs = plt.subplots(nrows=int(np.ceil(len(testdata_dict) / ncols)), ncols=ncols, figsize=(12, 9))
    plt.subplots_adjust(wspace=0.3, hspace=0.3, bottom=0.1)

    for i, (name, data) in enumerate(testdata_dict.items()):

        ax = axs.flatten()[i]
        reducer = umap.UMAP(n_neighbors=20, min_dist=0.3, random_state=42)
        embedding = reducer.fit_transform(data)

        df = pd.DataFrame(embedding, columns=['x', 'y'])
        df['label'] = labels

        if i == 0:
            sns.scatterplot(data=df, x='x', y='y', hue='label', palette="Set3", s=spot_size, ax=ax, alpha=0.8)
            handles, legend_labels = ax.get_legend_handles_labels()
            ax.get_legend().remove()
        else:
            sns.scatterplot(data=df, x='x', y='y', hue='label', palette="Set3", s=spot_size, ax=ax, alpha=0.8,
                            legend=False)
        ax.set_title(name)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.tick_params(axis='x', labelsize=5)
        ax.tick_params(axis='y', labelsize=5)

    for i in range(len(testdata_dict.items()), len(axs.flatten())):
        fig.delaxes(axs.flatten()[i])

    legend = fig.legend(handles, legend_labels, loc='lower center', bbox_to_anchor=(0.5, 0), ncol=5)
    return plt, legend


testfiles = [file for file in os.listdir(embdir) if file.endswith('pkl')]
metadata = pd.read_csv(metadata_path)
## 0 for Negative, 1 for Positive
metadata['Affinity'] = metadata['Affinity'].apply(lambda x: 'Negative' if x == 0 else 'Positive')

## case1: show CDR3s binding with top 10 epitopes (positive)
subdata = metadata[metadata['Affinity']==1]
labels = subdata['Epitope'].values
selected_idx = subdata.index.to_numpy()

plt1, legend = vis_umap(testfiles, selected_idx, labels, ncols)
plt1.savefig('figs/umap_IEDB_top10.svg', dpi=600, bbox_extra_artists=(legend,), bbox_inches='tight')
