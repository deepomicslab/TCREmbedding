import os
import numpy as np
import pandas as pd
import umap
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
from sklearn.model_selection import train_test_split

ncols = 3  # Define the number of columns in the figure
spot_size = 5
metadata_path = './data/IEDB_uniqueTCR_top10_filt.csv'
embdir = './data/binding_task_emb'

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

    fig, axs = plt.subplots(nrows=int(np.ceil(len(testdata_dict) / ncols)), ncols=ncols, figsize=(10, 4))
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

def vis_umap_for_epitope(epitope, testfiles, metadata,  ncols):

    subdata_pos = metadata[(metadata['Epitope'] == epitope) & (metadata['Affinity'] == 'Positive')]
    subdata_neg = metadata[(metadata['Epitope'] == epitope) & (metadata['Affinity'] == 'Negative')]
    subdata = pd.concat([subdata_pos, subdata_neg])
    selected_idx = subdata.index.to_numpy()
    labels = ["Positive" if i in subdata_pos.index else "Negative" for i in selected_idx]

    testdata_dict = {}
    for file in testfiles:
        model_name = file.split('_')[0]
        # load pkl
        data = pickle.load(open(os.path.join(embdir, file), 'rb'))

        if len(data.shape) == 2:
            testdata_dict[model_name] = data

    fig, axs = plt.subplots(nrows=int(np.ceil(len(testdata_dict) / ncols)), ncols=ncols, figsize=(7, 3))
    plt.subplots_adjust(wspace=0.3, hspace=0.3, bottom=0.1)

    for i, (name, data) in enumerate(testdata_dict.items()):
        ax = axs.flatten()[i]
        reducer = umap.UMAP(n_neighbors=20, min_dist=0.3, random_state=42)
        embedding = reducer.fit_transform(data)
        embedding = embedding[selected_idx]

        df = pd.DataFrame(embedding, columns=['x', 'y'])
        df['label'] = labels

        if i == 0:
            sns.scatterplot(data=df, x='x', y='y', hue='label', hue_order=['Positive', 'Negative'], palette=['deeppink', 'grey'], s=spot_size, ax=ax, alpha=0.8)
            handles, legend_labels = ax.get_legend_handles_labels()
            ax.get_legend().remove()
        else:
            sns.scatterplot(data=df, x='x', y='y', hue='label', hue_order=['Positive', 'Negative'],  palette=['deeppink', 'grey'], s=spot_size, ax=ax, alpha=0.8,
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


def vis_umap_for_epitope_subset(epitopes, testfiles, metadata,  ncols=3, figpath='binding_figs/umap_subset.svg'):
    testdata_dict = {}
    for file in testfiles:
        model_name = file.split('_')[0]
        # load pkl
        data = pickle.load(open(os.path.join(embdir, file), 'rb'))
        if len(data.shape) == 2:
            testdata_dict[model_name] = data
    nrows = len(epitopes) + 1

    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=( 1 + 2.5 * ncols, 0.5 + 2.6 * nrows))
    plt.subplots_adjust(wspace=0.35, hspace=0.35, left = 0.1, bottom=0.05, top=0.9)

    ## plot the whole dataset
    for i, (name, data) in enumerate(testdata_dict.items()):
        ax = axs[0, i]
        reducer = umap.UMAP(n_neighbors=20, min_dist=0.3, random_state=42)
        embedding = reducer.fit_transform(data)
        labels = metadata['Affinity'].to_list()
        df = pd.DataFrame(embedding, columns=['x', 'y'])
        df['label'] = labels
        sns.scatterplot(data=df, x='x', y='y', hue='label',  s=spot_size, ax=ax, alpha=0.5,
                        hue_order=['Positive', 'Negative'], palette=['deeppink', 'grey'], legend=False)
        ax.set_title(name)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.tick_params(axis='x', labelsize=5)
        ax.tick_params(axis='y', labelsize=5)
    # axs[0, 0].text(-4, 0, f"All Epitopes", va='center', ha='center', fontsize=14, rotation=90)

    for epitope_idx, epitope in enumerate(epitopes):
        subdata_pos = metadata[(metadata['Epitope'] == epitope) & (metadata['Affinity'] == 'Positive')]
        subdata_neg = metadata[(metadata['Epitope'] == epitope) & (metadata['Affinity'] == 'Negative')]
        subdata = pd.concat([subdata_pos, subdata_neg])
        selected_idx = subdata.index.to_numpy()
        labels = ["Positive" if i in subdata_pos.index else "Negative" for i in selected_idx]

        for i, (name, data) in enumerate(testdata_dict.items()):
            # ax = axs.flatten()[i]
            ax = axs[epitope_idx + 1 , i]
            reducer = umap.UMAP(n_neighbors=20, min_dist=0.3, random_state=42)
            embedding = reducer.fit_transform(data)
            # embedding = np.random.randn(2000, 2)

            embedding = embedding[selected_idx]

            df = pd.DataFrame(embedding, columns=['x', 'y'])
            df['label'] = labels

            if epitope_idx == 0 and i == 0:  # Save the handles to create a legend later
                sns.scatterplot(data=df, x='x', y='y', hue='label', hue_order=['Positive', 'Negative'],
                                palette=['deeppink', 'grey'], s=spot_size, ax=ax, alpha=0.8,
                                legend=True)
                handles, legend_labels = ax.get_legend_handles_labels()
                ax.get_legend().set_visible(False)
            else:
                sns.scatterplot(data=df, x='x', y='y', hue='label', hue_order=['Positive', 'Negative'],
                                palette=['deeppink', 'grey'], s=spot_size, ax=ax, alpha=0.8,
                                legend=False)

            ax.set_title(name)
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.tick_params(axis='x', labelsize=5)
            ax.tick_params(axis='y', labelsize=5)

            ## super title for this epitope
            # axs[epitope_idx + 1, 0].text(-4, 0, f"Epitope: {epitope}", va='center', ha='center', fontsize=14, rotation=90)

    # Add a single legend at the bottom of the figure
    legend = fig.legend(handles, legend_labels, loc='lower center', bbox_to_anchor=(0.5, 0), ncol=2)
    plt.savefig(figpath, dpi=600, bbox_extra_artists=(legend,))

def plot_umap_data_split(testfiles, y_train, y_test, idx_train, idx_test, specific_ept = False, sub_index = None, metadata= None, pos_color='deeppink'):

    testdata_dict = {}
    for file in testfiles:
        model_name = file.split('_')[0]
        # load pkl
        data = pickle.load(open(os.path.join(embdir, file), 'rb'))
        if len(data.shape) == 2:
            testdata_dict[model_name] = data
        print("loading model: ", model_name)
    model_num = len(testfiles)
    fig, axs = plt.subplots(nrows=model_num, ncols=2, figsize=(6, 8))
    plt.subplots_adjust(wspace=0.3, hspace=0.3, bottom=0.05)

    if specific_ept:
        idx_train = np.intersect1d(idx_train, sub_index)
        idx_test = np.intersect1d(idx_test, sub_index)
        y_train = metadata.loc[idx_train, 'Affinity'].to_list()
        y_test = metadata.loc[idx_test, 'Affinity'].to_list()

    for i, (name, data) in enumerate(testdata_dict.items()):
        ax_train = axs.flatten()[i * 2 ]
        ax_val = axs.flatten()[i * 2 + 1]
        reducer = umap.UMAP(n_neighbors=20, min_dist=0.3, random_state=42)
        embedding = reducer.fit_transform(data)

        train_embedding = embedding[idx_train]
        df = pd.DataFrame(train_embedding, columns=['x', 'y'])
        df['label'] = y_train
        sns.scatterplot(data=df, x='x', y='y', hue='label', hue_order=['Positive', 'Negative'],
                        palette=[pos_color, 'grey'], s=spot_size, ax=ax_train, alpha=0.6)
        ax_train.set_title(f'{name}-train set')
        ax_train.set_xlabel('')
        ax_train.set_ylabel('')
        ax_train.tick_params(axis='x', labelsize=5)
        ax_train.tick_params(axis='y', labelsize=5)
        ax_train.legend(loc='lower right',
                        fancybox=True, shadow=False, fontsize='x-small')

        xticks_train = ax_train.get_xticks()
        yticks_train = ax_train.get_yticks()
        xlim_train = ax_train.get_xlim()
        ylim_train = ax_train.get_ylim()

        ## validation
        val_embedding = embedding[idx_test]
        df = pd.DataFrame(val_embedding, columns=['x', 'y'])
        df['label'] = y_test

        sns.scatterplot(data=df, x='x', y='y', hue='label', hue_order=['Positive', 'Negative'],
                        palette=[pos_color, 'grey'], s=spot_size, ax=ax_val, alpha=0.6)
        ax_val.set_title(f'{name}-validation set')
        ax_val.set_xlabel('')
        ax_val.set_ylabel('')
        ax_val.set_xticks(xticks_train)
        ax_val.set_yticks(yticks_train)
        ax_val.set_xlim(xlim_train)
        ax_val.set_ylim(ylim_train)
        ax_val.tick_params(axis='x', labelsize=5)
        ax_val.tick_params(axis='y', labelsize=5)
        ax_val.legend(loc='lower right',
                        fancybox=True, shadow=False, fontsize='x-small')

    return plt


def specific_split(metadata, column='Epitope'):

    arr = metadata[column].to_numpy()
    unique_elements, counts = np.unique(arr, return_counts=True)
    a_size = int(len(arr) * 0.8)
    ## initialize a and b arrays
    a = []
    a_indices = []
    b = []
    b_indices = []
    print('splitting...')

    for element, count in zip(unique_elements, counts):
        element_indices = np.where(arr == element)[0]
        if count < a_size - len(a):
            a.extend([element] * count)
            a_indices.extend(element_indices[:count])
        else:
            b.extend([element] * count)
            b_indices.extend(element_indices[-count:])

    a = np.array(a)
    b = np.array(b)
    a_indices = np.array(a_indices)
    b_indices = np.array(b_indices)

    # check if b and a are disjoint
    result = np.isin(b, a)

    if np.all(~result):
        ## for each element in b, check if it is in a
        print("successfully split the array:")
        a_size = a.size
        b_size = b.size

        print('train set size:', a_size, 'percentage:', a_size / (a_size + b_size) * 100, '%')
        print('test set size:', b_size, 'percentage:', b_size / (a_size + b_size) * 100, '%')

        return a_indices, b_indices
    else:
        print('splitting failed')
        return None, None

testfiles = [file for file in os.listdir(embdir) if file.endswith('pkl')]
metadata = pd.read_csv(metadata_path)
## 0 for Negative, 1 for Positive
metadata['Affinity'] = metadata['Affinity'].apply(lambda x: 'Negative' if x == 0 else 'Positive')


# ## case1: visualize all epitopes and subsets
unique_epitopes = metadata['Epitope'].unique()
print(unique_epitopes)
epitope_show_num = 3 # show n epitopes
## random choose 5 epitopes
np.random.seed(42)
epitopes_to_visualize = np.random.choice(unique_epitopes, epitope_show_num, replace=False)
print("epitopes_to_visualize:")
print(epitopes_to_visualize)
epitopes_to_visualize = unique_epitopes[:epitope_show_num]
vis_umap_for_epitope_subset(epitopes_to_visualize, testfiles, metadata, ncols=3, figpath='binding_figs/umap_subset.svg')



## case2: random split vs epitope split
np.random.seed(42)
X = metadata['CDR3'].values
y = metadata['Affinity'].values
# y = np.array(['Positive' if i == 1 else 'Negative' for i in y])
X_train, X_test, y_train, y_test, idx_train, idx_test = train_test_split(
    X, y, range(len(X)), test_size=0.1, random_state=42
)
plt2 = plot_umap_data_split(testfiles, y_train, y_test, idx_train, idx_test, pos_color='orange')
plt2.savefig('binding_figs/umap_IEDB_random_split.svg', dpi=600)

## case3: split by TCR
np.random.seed(42)
idx_train, idx_test = None, None
while idx_train is None:
    idx_train, idx_test = specific_split(metadata, column='CDR3')
y_train = metadata.loc[idx_train, 'Affinity'].to_list()
y_test = metadata.loc[idx_test, 'Affinity'].to_list()

plt3 = plot_umap_data_split(testfiles, y_train, y_test, idx_train, idx_test, pos_color='limegreen')
plt3.savefig('binding_figs/umap_IEDB_tcr_split.svg', dpi=600)

## case4: split by epitope
np.random.seed(42)
idx_train, idx_test = None, None
while idx_train is None:
    idx_train, idx_test = specific_split(metadata, column='Epitope')
y_train = metadata.loc[idx_train, 'Affinity'].to_list()
y_test = metadata.loc[idx_test, 'Affinity'].to_list()

plt4 = plot_umap_data_split(testfiles, y_train, y_test, idx_train, idx_test, pos_color='lightcoral')
plt4.savefig('binding_figs/umap_IEDB_epitope_split.svg', dpi=600)


