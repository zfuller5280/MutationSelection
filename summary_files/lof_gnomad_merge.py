import pandas as pd
pd.options.mode.chained_assignment = None
import sys

mu_rates = sys.argv[1]
gnomad = sys.argv[2]

mu_rates = pd.read_csv(mu_rates, sep="\t")
gnomad = pd.read_csv(gnomad, sep="\t")

merged_gnomad = pd.merge(mu_rates, gnomad, on=['SYMBOL'])

chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
                '13', '14', '15', '16', '17', '18', '19', '20', '21', '22' , 'X']

for chrom in chromosomes:
    gnomad_subset = merged_gnomad[merged_gnomad['chromosome']==chrom]
    subset_out = gnomad_subset[['SYMBOL','chromosome','mu_lof','AF_nfe','AN_nfe']]
    subset_out['NFE_k'] = subset_out['AF_nfe'].astype(float) * subset_out['AN_nfe'].astype(float)
    subset_out['NFE_k'] = subset_out['NFE_k'].astype(int)
    subset_out['AN_nfe'] = subset_out['AN_nfe'].astype(int)

    subset_out.to_csv('%s.lof_summary.nfe_canonical.txt'%chrom,sep="\t",index=False)
