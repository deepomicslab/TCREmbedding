import os
import numpy as np
import pandas as pd
import umap
import matplotlib.pyplot as plt
import seaborn as sns


topk_label = 10  # Define the number of top labels to keep
max_samples_per_class = 100 ## Define the maximum number of samples per class
ncols = 3  # Define the number of columns in the figure
spot_size = 5
metadata = pd.read_csv('./IEDB_uniqueTCR_test.csv')

labels = metadata['Epitope'].values
labels_series = pd.Series(labels)
# Count the occurrences of each label
label_counts = labels_series.value_counts()
# Get the top 10 labels with the highest counts
top_10_labels = label_counts[:topk_label].index

top_10_indices = labels_series[labels_series.isin(top_10_labels)].index
# Subset the metadata using the indices of the top 10 labels
metadata = metadata.loc[top_10_indices]
new_labels_series = pd.Series(metadata['Epitope'].values)


new_metadata = pd.DataFrame()
for label in new_labels_series.unique():
    class_samples = metadata[metadata['Epitope'] == label]
    if len(class_samples) > max_samples_per_class:
        class_samples = class_samples.sample(n=max_samples_per_class, random_state=42)
    new_metadata = pd.concat([new_metadata, class_samples])

new_indices = new_metadata.index
labels = new_metadata['Epitope'].values
print(new_metadata)
new_metadata.to_csv('IEDB_uniqueTCR_top10_filt.csv', index=False)