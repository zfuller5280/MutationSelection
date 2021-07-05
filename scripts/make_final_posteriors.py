import pandas as pd
import os, random
import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np
from random import randrange
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
gene_posteriors = pd.read_csv("hs_gene_lookups.raw.csv")

gene_posteriors = gene_posteriors.apply(lambda x: np.where(x.isnull(), x.dropna().sample(len(x), replace=True), x))
gene_posteriors = gene_posteriors.replace(r'\\n','',regex=True)

n_list_path = "neutral_gene_list.4_26.txt"
n_genes = []
with open(n_list_path) as f:
    for line in f:
        line = line.strip("\r\n")
        n_genes.append(line)

for gene in n_genes:
    vals = np.array([0 for val in range(50000)])
    gene_posteriors[gene] = vals

gene_posteriors.to_csv("autosome.posteriors.hs.csv", index=False)
