import pandas as pd
import matplotlib.pyplot as plt
import sys

# Read in the DataFrame
df = pd.read_csv(sys.argv[1], header=None, sep="\t", names=['Accession 1', 'Accession 2', 'Chr2', 'Chr1', 'D', 'P', 'Hashes'])

# creating a histogram
plt.hist(df['D'], bins=100)
plt.savefig('mash_histogram_distances.png')

# creating a histogram log y
plt.hist(df['D'], bins=100)
plt.yscale('log')
plt.savefig('mash_histogram_distances_logy.png')