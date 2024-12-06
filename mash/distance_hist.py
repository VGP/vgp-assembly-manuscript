import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

# Read in the DataFrame
df = pd.read_csv(sys.argv[1], header=None, sep=" ", names=['Accession 1', 'Accession 2', 'D', 'P', 'Hashes'])

# creating a histogram
plt.hist(df['D'], bins=50)
plt.savefig('mash_histogram_distances.png')