import numpy as np
import math
from hpd import hpd_grid
import sys
import pandas as pd
pd.options.mode.chained_assignment = None

infile = sys.argv[1]
name = sys.argv[2]

def try_float(x):
    try:
        float(x)
        return True
    except:
        return False

#col_names = ["s", "h", "hs", "sim_freq", "obs", "prob", "ratio", "mh_ratio", "prop_s","prop_h", "accept", "mutU", "pi_x", "pi_y", "s_xy", "s_yx"]
#raw_samples = pd.read_csv(infile, sep=" ", header=None, names=col_names)
#samples = raw_samples.drop_duplicates(subset=['accept'])
raw_samples = pd.read_csv(infile, sep=",", header=None, names=["s","h"])
samples = raw_samples[raw_samples['s']!='s']
#print(samples)
samples['s'] = samples['s'].astype(float)
samples['h'] = samples['h'].astype(float)
samples['hs'] = samples['s'] * samples['h']
hpd_mu, x_mu, y_mu, modes_mu = hpd_grid(samples['hs'], roundto=10)

log_hpd_mu, log_x_mu, log_y_mu, log_modes_mu = hpd_grid(np.log10(samples['hs']), roundto=8)

modes_mu_max = modes_mu.index(max(modes_mu))
log_modes_mu_max = log_modes_mu.index(max(log_modes_mu))

print(name, hpd_mu[modes_mu_max][0], hpd_mu[modes_mu_max][1], modes_mu[modes_mu_max], log_hpd_mu[log_modes_mu_max][0], log_hpd_mu[log_modes_mu_max][1], log_modes_mu[log_modes_mu_max], np.percentile(samples['hs'],2.5), np.percentile(samples['hs'],97.5), np.mean(samples['hs']), np.median(samples['hs']))
