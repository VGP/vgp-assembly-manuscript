import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 1000)

# Read in the DataFrame
df = pd.read_csv(sys.argv[1], header=None, sep=",", names=['Accession 1', 'Tolid 1', 'Accession 2', 'Tolid 2', 'Chr 1', 'Chr 2', 'D'])
df = df[df['D'] < 0.5]
grouped = df.groupby(['Accession 1', 'Accession 2'], sort=False)
result = pd.DataFrame()
outliers = pd.DataFrame()

for name, group in grouped:
    two_ways = group.groupby(['Chr 2'])['D'].idxmin()
    result = pd.concat([result, df.loc[two_ways]], ignore_index=True)

    group["is_outlier"] = group.groupby(['Chr 2'], sort=False)['D'].transform(lambda x: -(x - x.mean()) > 4*x.std())
    group = group.loc[group["is_outlier"] == True]
    outliers = pd.concat([outliers, group], ignore_index=True)

result.to_csv('min_hits.csv', index=False)
outliers.to_csv('outliers.csv', index=False)