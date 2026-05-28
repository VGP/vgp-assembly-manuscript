# Import libraries
from matplotlib import pyplot as plt
from matplotlib import ticker
import numpy as np
import pandas as pd
import seaborn as sns
from IPython.display import display
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

# General settings
pd.set_option('display.max_rows', 50)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1000)

# Import data
vgp_assembly_dataset = pd.read_csv('/Users/gformenti/sandbox/vgp-assembly/gfastats.tsv', sep="\t")

# Preprocessing
# Replace 'Class' with 'Lineage' for the dataset
vgp_assembly_dataset['Lineage'] = vgp_assembly_dataset['Tolid'].str[0]
vgp_assembly_dataset.insert(2, 'Lineage', vgp_assembly_dataset.pop('Lineage'))
vgp_assembly_dataset[['A', 'C', 'G', 'T']] = vgp_assembly_dataset['Base composition (A'].str.split(':',
                                                                                                   expand=True).astype(
    int)
vgp_assembly_dataset.drop(['Base composition (A'], axis=1, inplace=True)
vgp_assembly_dataset = vgp_assembly_dataset[~vgp_assembly_dataset['Tolid'].str.contains("alt", case=False)]

# Scaling
scaler = StandardScaler().set_output(transform="pandas")
vgp_assembly_dataset_scaled = scaler.fit_transform(vgp_assembly_dataset.select_dtypes(['number']))
vgp_assembly_dataset_scaled['Lineage'] = vgp_assembly_dataset['Tolid'].str[0]

# Summary
print(f'N genomes: {vgp_assembly_dataset.shape[0]}, N features: {vgp_assembly_dataset.shape[1]}')
print(pd.Index(vgp_assembly_dataset_scaled['Lineage']).value_counts())
print(*vgp_assembly_dataset.columns.tolist(), sep="\n")
with pd.option_context('display.max_rows', None):
    display(vgp_assembly_dataset)
print(vgp_assembly_dataset.describe().apply(lambda s: s.apply('{0:.1f}'.format)))


# All-vs-all plot
def pairplot():
    g = sns.pairplot(vgp_assembly_dataset_scaled, hue='Lineage', height=3, corner=True)
    plt.subplots_adjust(wspace=0.3, hspace=0.3)

    labels = []
    for ax in g.axes[-1, :]:
        label = ax.xaxis.get_label_text()
        labels.append(label)

    for i in range(len(labels)):
        for j in range(i + 1):
            g.axes[i, j].xaxis.set_label_text(labels[j], visible=True)
            g.axes[i, j].yaxis.set_label_text(labels[i], visible=True)

    g.savefig('all_vs_all.png', dpi=300)


# pairplot()

# PCA
def plotPCA():
    n_features = vgp_assembly_dataset_scaled.shape[1] - 1  # drop labels
    print(f'Number of features: {n_features}')
    pca = PCA(n_components=n_features, random_state=42)
    vgp_assembly_dataset_tf = pca.fit_transform(vgp_assembly_dataset_scaled.drop('Lineage', axis=1))
    print(pca.explained_variance_ratio_)
    print(np.sum(pca.explained_variance_ratio_))
    print(np.cumsum(pca.explained_variance_ratio_))

    # EVR
    evr = pca.explained_variance_ratio_
    cum_evr = np.cumsum(pca.explained_variance_ratio_)
    fig, ax1 = plt.subplots(1, 1, figsize=(6, 4))
    ax1.set_xlabel('Principal Components')
    ax1.set_ylabel('Explained variance ratio (EVR)')
    ax1.plot(np.arange(1, n_features + 1, 1), evr, color='tab:red', marker='o', label='EVR')
    ax1.tick_params(axis='y')
    ax1.set_title('Dry Beans: PCA - EVR')
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.set_ylabel('Cumulative EVR')  # we already handled the x-label with ax1
    ax2.plot(np.arange(1, n_features + 1, 1), cum_evr, color='tab:blue', marker='o', label='Cumulative EVR')
    ax2.tick_params(axis='y')
    ax2.axhline(y=0.99, color='tab:purple', linestyle='dashed', label='99% cut-off threshold')
    ax2.axhline(y=0.95, color='tab:orange', linestyle='dashed', label='95% cut-off threshold')
    fig.tight_layout()
    fig.legend(bbox_to_anchor=(0.89, 0.8))
    plt.xticks(range(1, n_features + 1, 2))
    plt.savefig('PCA_EVR.png', dpi=300)

    # PCA 2D Plot
    cmap_brg = plt.get_cmap('brg')
    fig, ax = plt.subplots()
    plt.gca()
    plt.title("vgp-assembly-manuscript PCA")
    plt.xlabel(f"PC1 (% {evr[0] * 100})")
    plt.ylabel(f"PC2 (% {evr[1] * 100})")

    # Use Lineage for coloring
    lineages = pd.unique(vgp_assembly_dataset['Lineage'])
    labels_dict = dict(zip(lineages, range(len(lineages))))
    scatter = ax.scatter(vgp_assembly_dataset_tf[:, 0], vgp_assembly_dataset_tf[:, 1],
                         c=vgp_assembly_dataset['Lineage'].map(labels_dict), cmap=cmap_brg, s=8)

    handles, labels = scatter.legend_elements()
    legend = ax.legend(handles, lineages, loc="upper right", title="Lineage")
    plt.xticks()
    plt.yticks()
    plt.axis('tight')
    plt.savefig('PCA.png', dpi=300)


# Call plotPCA function
plotPCA()


# t-SNE
def plotTSNE():
    vgp_assembly_dataset_embedded = TSNE(random_state=1, n_jobs=-1).fit_transform(
        vgp_assembly_dataset_scaled.drop('Lineage', axis=1))
    df = {'x': vgp_assembly_dataset_embedded[:, 0],
          'y': vgp_assembly_dataset_embedded[:, 1],
          'target': vgp_assembly_dataset_scaled['Lineage']}

    fig, ax = plt.subplots()
    cmap_brg = plt.get_cmap('brg')
    plt.gca()
    plt.title("vgp-assembly-manuscript tSNE")

    # Use Lineage for coloring
    lineages = pd.unique(vgp_assembly_dataset['Lineage'])
    labels_dict = dict(zip(lineages, range(len(lineages))))
    scatter = plt.scatter(df['x'], df['y'], c=df['target'].map(labels_dict), alpha=0.7, s=8, cmap=cmap_brg)

    handles, labels = scatter.legend_elements()
    ax.legend(handles, lineages, loc="upper right", title="Lineage")
    plt.xticks([])
    plt.yticks([])
    plt.grid(None)
    plt.axis('tight')
    plt.savefig('tSNE.png')


# Call plotTSNE function
plotTSNE()