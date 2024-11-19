import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

# Read in the DataFrame
df = pd.read_csv(sys.argv[1], header=None, sep=" ", names=['Accession 1', 'Accession 2', 'D', 'P', 'Hashes'])

# Create a sample dataframe
data = {'A': [1, 2, 3, 4],
        'B': [5, 6, 7, 8],
        'C': [9, 10, 11, 12]}
df = pd.DataFrame(data)

# Create the clustermap

import numpy as np

lis = [('A', 'B', 3),
        ('A', 'D', 4),
        ('B', 'D', 4),
        ('B', 'H', 5),
        ('C', 'L', 2),
        ('D', 'F', 1),
        ('F', 'H', 3),
        ('G', 'H', 2),
        ('G', 'Y', 2),
        ('I', 'J', 6),
        ('I', 'K', 4)]

items = set.union(set([item[0].upper() for item in lis]) , set([item[1].upper() for item in lis]))
value = dict(zip(sorted(items), range(26)))
dist_matrix = np.zeros((len(items) , len(items)))

for i in range(len(lis)):

    # for Upper triangular
    dist_matrix[value[lis[i][0]] , value[lis[i][1]]] = lis[i][2]
    # for lower triangular
    dist_matrix[value[lis[i][1]] , value[lis[i][0]]] = lis[i][2]

    """
    Example:
    [0 3 0]
    [3 0 0]
    [0 0 0]
    """

dist_matrix

sns.clustermap(df)
plt.savefig('mash_heatmap_distances.png')