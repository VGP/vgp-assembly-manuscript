import pandas as pd
import sys

edges = pd.read_csv(sys.argv[1])
lookup = pd.read_csv(sys.argv[2], sep="\t", header=None)
lookup = lookup.set_axis(["human", "Chr 1"], axis=1, copy=False)

inner_join = pd.merge(edges,
                      lookup,
                      on ='Chr 1',
                      how ='inner')

inner_join.to_csv('human_outliers_with_chr.csv', index=False)