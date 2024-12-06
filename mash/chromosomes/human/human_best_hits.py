import pandas as pd
import sys

def compute_outliers(x):
    z_score= -(x - x.mean()) / x.std()
    is_outlier = z_score > 3.29 # 1/1000 in normal distribution
    return pd.concat([z_score, is_outlier], ignore_index=True, axis=1)

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 1000)

# read data in
df = pd.read_csv(sys.argv[1], header=None, sep="\t", names=['Accession', 'Chr 1', 'Chr 2', 'D', 'P', 'Hashes'])

# filter very large values
df = df[df['D'] < 0.5]

# detect outliers using Z-score
result = df.groupby(['Chr 2'], sort=False)['D'].apply(compute_outliers)
df[['z_score','is_outlier']]  = result.droplevel(0)
df = df.loc[df["is_outlier"] == True]
df.drop('is_outlier', axis=1, inplace=True)

df.to_csv('human_outliers.csv', index=False)
