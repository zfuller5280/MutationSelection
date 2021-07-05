import pandas as pd
pd.options.mode.chained_assignment = None
import sys
import os
import numpy as np
from hpd import hpd_grid

indir = sys.argv[1]
extension = sys.argv[2]

files = []
for file in os.listdir(indir):
    if file.endswith(".txt"):
        files.append(file)

genes = [x.strip(extension) for x in files]

neut_sims_frq = pd.read_csv("./hg19.exomes.neutral2.sim.input.txt",sep="\t")
neut_sims = pd.read_csv("./hg19.exomes.neutral.sim.input.txt",sep="\t")
canonical_gnomad_stats = pd.read_csv("./canoncial_gnomad.stats.csv")

neut_gnomad_df = neut_sims.merge(canonical_gnomad_stats, left_on='SYMBOL', right_on='gene')
neut_gnomad_df = neut_gnomad_df.merge(neut_sims_frq, on='SYMBOL')
neut_gnomad_df = neut_gnomad_df[neut_gnomad_df['SYMBOL'].isin(genes)]
neut_gnomad_df['nfe_obs_syn_k'] = neut_gnomad_df['AF_nfe']/neut_gnomad_df['AF_nfe_x']
neut_gnomad_df['nfe_avg_syn_AF'] = neut_gnomad_df['AF_nfe']/neut_gnomad_df['possible_syn']

for n, sim_file in enumerate(files):
    gene = genes[n]
    # try:
    sim_df = pd.read_csv("%s/%s"%(indir, sim_file), sep=" ", header = None, error_bad_lines=False)
    sim_df.columns = ['gen','stop','alt_hom','hets','ref_hom','mut']
    sim_df.dropna(inplace=True)
    #print(sim_df)
    sim_df['freq'] = ((sim_df['alt_hom']*2) + sim_df['hets'])/((sim_df['hets'] + sim_df['alt_hom'] + sim_df['ref_hom'])*2)
    gene_sample_size = (neut_gnomad_df[neut_gnomad_df['SYMBOL']==gene]['AN_nfe_y']).astype(int)
    gene_possible_syn = (neut_gnomad_df[neut_gnomad_df['SYMBOL']==gene]['possible_syn']).astype(int)
    try:
        sim_df['samp_count'] = np.random.binomial(gene_sample_size, sim_df['freq'], size=len(sim_df['freq']))
        sim_df['samp_freq'] = sim_df['samp_count'].astype(float).div(int(gene_sample_size))
    except:
        continue
    mean_sim_AF = np.mean(sim_df['samp_freq'])
    median_sim_AF = np.median(sim_df['samp_freq'])
    if sum(sim_df['samp_freq']) > 0:
        hpd_mu, x_mu, y_mu, modes_mu = hpd_grid(sim_df['samp_freq'], roundto=6)
        log_array = np.ma.array(sim_df['samp_freq'], mask=(sim_df['samp_freq']<=0))
        log_hpd_mu, log_x_mu, log_y_mu, log_modes_mu = hpd_grid(np.log10(log_array), roundto=6)
        mean_log_AF = np.mean(np.log10(log_array))
        median_log_AF = np.median(np.log10(log_array))
        std_sim_AF = np.std(sim_df['samp_freq'])
        print(gene, mean_sim_AF, median_sim_AF, std_sim_AF, mean_log_AF, median_log_AF, modes_mu[0], hpd_mu[0][0], hpd_mu[0][1], log_modes_mu[0], log_hpd_mu[0][0], log_hpd_mu[0][1])

    else:
        print(gene, mean_sim_AF, median_sim_AF, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    # except:
    #     print(gene, "NA", "NA")
