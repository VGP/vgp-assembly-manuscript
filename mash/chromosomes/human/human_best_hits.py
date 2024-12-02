import pandas as pd
import sys

def compute_outliers(x):
    D_avg = x.avg()
    D_min = x.min()
    zscore = -(x - x.mean()) / x.std()
    is_outlier = zscore > 3.29 # 1/1000 in normal distribution
    return D_avg, D_min, zscore, is_outlier

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 1000)

# read data in
df = pd.read_csv(sys.argv[1], header=None, sep="\t", names=['Accession', 'Chr 1', 'Chr 2', 'D', 'P', 'Hashes'])

# filter very large values
df = df[df['D'] < 0.5]

# detect outliers using Z-score
df[['D_avg_value', 'D_min_value', 'Z_score', 'is_outlier']] = df.groupby(['Chr 2'], sort=False)['D'].transform(compute_outliers).tolist()

df = df.loc[df["is_outlier"] == True]

df.to_csv('human_outliers.csv', index=False)
