import pandas as pd
import sys

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 1000)

# Read in the DataFrame
df = pd.read_csv(sys.argv[1], header=None, sep="\t", names=['Accession', 'Chr 1', 'Chr 2', 'D', 'P', 'Hashes'])
df = df[df['D'] < 0.5]

two_ways = df.groupby(['Chr 2'])['D'].idxmin()

df["is_outlier"] = df.groupby(['Chr 2'], sort=False)['D'].transform(lambda x: -(x - x.mean()) > 4*x.std())
df = df.loc[df["is_outlier"] == True]

two_ways.to_csv('human_min_hits.csv', index=False)
df.to_csv('human_outliers.csv', index=False)
