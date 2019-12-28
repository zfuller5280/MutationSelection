from optparse import OptionParser
import numpy as np
from random import random as rm
import Queue as queue
import multiprocessing as mp
from ctypes import cdll
import subprocess
from subprocess import Popen, PIPE
from operator import add, sub
import random
from scipy.stats import beta
import math
from hpd import hpd_grid


USAGE= """Usage: %prog [options]"""
OPT_DEFAULTS={'infile':'-'}
DESCRIPTION="""Program description: """
EPILOG="""Requirements:"""

def get_options(defaults, usage, description='',epilog=''):
    """Get options, print usage text."""
    parser=OptionParser(usage=usage,description=description,epilog=epilog)
    parser.add_option("-i","--infile",action="store",dest="infile",type="string",
                      default=defaults.get('infile'),
                      help='Name of input gene')
    parser.add_option("-r","--runs",action="store",dest="runs",type="int",
                      default=defaults.get('runs'),
                      help='Number of runs for simulation')
    parser.add_option("-m","--mutation-rate",action="store",dest="mut",type="float",
                      default=defaults.get('mut'),
                      help='Mutation rate (float)')
    parser.add_option("-o","--obs-freq",action="store",dest="obs",type="float",
                      default=defaults.get('obs'),
                      help='Observed mutation frequnecy')

    (options,args)=parser.parse_args()

    return (options, args)

def run_simulation(runs, sel, dom, mutU):
    cmd = ["/Users/zachfuller/mut_population_sims/mut_uncert_sim","%i"%runs,"%s"%sel,"%s"%dom,"%s"%mutU]
    result = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    output = result.stdout.readlines()

    out = [float(x) for x in output[0].split()[5:]]
    freq = ((out[0]*2) + out[1])/(sum(out[:-1])*2)
    out.extend([sel, dom])

    pop_samp = np.random.binomial(113771, freq, 1)
    pop_freq = float(pop_samp)/113771

    mut_u = float(output[0].split()[-1].strip("\n"))
    #print out
    return pop_freq, out, mut_u

def lognuniform(low=1e-6, high=1, size=None, base=np.exp(1)):
    return np.power(base, np.random.uniform(np.log10(low), np.log10(high), size))

def clamp(n, minn=1e-4, maxn=1-(1e-4)):
    return max(min(maxn, n), minn)

def kde_scipy(x, obs, bandwidth=0.3, **kwargs):
    """Kernel Density Estimation with Scipy"""
    # Note that scipy weights its bandwidth by the covariance of the
    # input data.  To make the results comparable to the other methods,
    # we divide the bandwidth by the sample standard deviation here.
    tol = 1e-05
    try:
        kde = gaussian_kde(x, bw_method=bandwidth, **kwargs)
        return kde.integrate_box(obs-tol,obs+tol)
    except:
        return 0.0

def get_beta_params(exp_u, mutU):
    # if exp_u <= 0.5:
    #     skew = min(1.5, (2./3)*(1-2*exp_u)/((exp_u*(1-exp_u))**.5))
    # else:
    #     skew = max(-1.5, (2./3)*(1-2*exp_u)/((exp_u*(1-exp_u))**.5))
    # if skew > 0:
    #     alpha = 2.0*(-skew**2*exp_u**2 + skew**2*exp_u - 4.0*exp_u**2 + 4.0*exp_u + (2.0*exp_u - 1.0)*((skew**2*exp_u**2 - skew**2*exp_u + 4.0*exp_u**2 - 4.0*exp_u + 1.0)**.5) - 1.0)/(skew**2*(exp_u - 1.0))
    #     beta_p = 2.0*(skew**2*exp_u**2 - skew**2*exp_u + 4.0*exp_u**2 - 4.0*exp_u + (-2.0*exp_u + 1.0)*((skew**2*exp_u**2 - skew**2*exp_u + 4.0*exp_u**2 - 4.0*exp_u + 1.0)**.5) + 1.0)/(skew**2*exp_u)
    # else:
    #     alpha = 2.0*(-skew**2*exp_u**2 + skew**2*exp_u - 4.0*exp_u**2 + 4.0*exp_u + (-2.0*exp_u + 1.0)*((skew**2*exp_u**2 - skew**2*exp_u + 4.0*exp_u**2 - 4.0*exp_u + 1.0)**.5) - 1.0)/(skew**2*(exp_u - 1.0))
    #     beta_p = 2.0*(skew**2*exp_u**2 - skew**2*exp_u + 4.0*exp_u**2 - 4.0*exp_u + (2.0*exp_u - 1.0)*((skew**2*exp_u**2 - skew**2*exp_u + 4.0*exp_u**2 - 4.0*exp_u + 1.0)**.5) + 1.0)/(skew**2*exp_u)
    # return alpha, beta_p, skew
    mut_samp = mutU/exp_u
    alpha = mut_samp * 113770
    beta_p = 113770 + 1
    return alpha+1, beta_p

def main():
    (options,args)=get_options(OPT_DEFAULTS, USAGE, DESCRIPTION, EPILOG)
    infile = options.infile
    runs = options.runs
    mutU = options.mut
    obs = options.obs
    if obs < 0:
        print infile, "NA", "NA","NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"
    else:
        #init_s = lognuniform(size=1)[0]
        #init_h = lognuniform(size=1)[0]
        #init_hs = lognuniform(size=1)[0]
        ops = (add, sub)
        init_hs = mutU/obs
        if init_hs > 1: init_hs = .9
        alpha, beta_p = get_beta_params(init_hs, mutU)
        init_hs = clamp(mutU/np.random.beta(a=alpha,b=beta_p))
        init_s, init_h = (init_hs)**.5, (init_hs)**.5
        s, h = init_s, init_h
        hs = init_hs
        #print s, h
        init_freq = run_simulation(1, s, h, mutU)

        #init_freq = run_simulation(1, init_s, init_h, mutU)
        accept = 0
        #hs = s*h
        #diff = 1 - abs(init_freq[0] - obs)
        #pi_x = lognorm.cdf()
        curr_freq = init_freq[0]
        samp_size = 113770

        pi_x = beta.pdf(curr_freq, a=samp_size*obs + 1, b=samp_size-(samp_size*obs)+1)
        #pi_x = min(1-pi_x, pi_x)*2
        hs_L = []

        for run in xrange(runs):

            #prop_s, prop_h = clamp(np.random.lognormal(mean=np.log(s),sigma=.1)), clamp(np.random.lognormal(mean=np.log(h),sigma=.1))
            #prop_param_s, prop_param_h = (float(samp_size)/50)*s, (float(samp_size)/50)*h
            prop_param_s, prop_param_h = 50, 50
            prop_s  = clamp(np.random.beta(a=prop_param_s,b=(prop_param_s-(s*prop_param_s))/s))
            prop_h = clamp(np.random.beta(a=prop_param_h,b=(prop_param_h-(h*prop_param_h))/h))
            #exp_u = hs
            #alpha, beta_p = get_beta_params(exp_u, mutU)


            #prop_s = clamp(mutU/np.random.beta(a=alpha,b=beta_p)/h)
            #prop_hs = clamp(mutU/np.random.beta(a=alpha,b=beta_p))
            #prop_h = prop_hs/prop_s
            #alpha_h, beta_p_h = get_beta_params(h, mutU)
            #prop_h = clamp(mutU/np.random.beta(a=alpha,b=beta_p)/s)
            #print exp_u, alpha, beta_p, prop_s, prop_h
            #print exp_u, alpha, beta_p, mutU/hs, prop_h, prop_s
            #print prop_hs
            sim_result = run_simulation(1, prop_s, prop_h, mutU)
            #sim_result = run_simulation(1, prop_hs, 1, mutU)
            s_xy = beta.pdf(prop_s, a=prop_param_s, b=(prop_param_s-(s*prop_param_s))/s)
            #s_xy = beta.pdf((mutU/(prop_s*h)), a=alpha, b=beta_p)
            #h_xy = beta.pdf((mutU/(prop_h*s)), a=alpha, b=beta_p)
            # s_xy = min(1-s_xy, s_xy)*2
            # #s_xy = min(1-lognorm.cdf(prop_s, s=s, loc=0, scale=s),lognorm.cdf(prop_s, s=s, loc=0, scale=s))
            #s_yx = beta.pdf(s, a=10, b=(10-(prop_s*10))/prop_s)

            # s_yx = min(1-s_yx, s_yx)*2
            # #s_yx = min(1-lognorm.cdf(s, s=prop_s, loc=0, scale=prop_s),lognorm.cdf(s, s=prop_s, loc=0, scale=prop_s))
            #alpha, beta_p = get_beta_params(prop_s*prop_h, mutU)
            #print prop_s*prop_h, alpha, beta_p

            #s_yx = beta.pdf((mutU/(s*prop_h)), a=alpha, b=beta_p)
            #h_yx = beta.pdf((mutU/(h*prop_s)), a=alpha, b=beta_p)
            h_xy = beta.pdf(prop_h, a=prop_param_h, b=(prop_param_h-(h*prop_param_h))/h)
            # h_xy = min(1-h_xy, h_xy)*2
            # #h_xy = min(1-lognorm.cdf(prop_h, s=h, loc=0, scale=h),lognorm.cdf(prop_h, s=h, loc=0, scale=h))
            #prop_param_s, prop_param_h = (float(samp_size)/10)*prop_s, (float(samp_size)/10)*prop_h
            s_yx = beta.pdf(s, a=prop_param_s, b=(prop_param_s-(prop_s*prop_param_s))/prop_s)
            h_yx = beta.pdf(h, a=prop_param_h, b=(prop_param_h-(prop_h*prop_param_h))/prop_h)
            # h_yx = min(1-h_yx, h_yx)*2
            #hs_xy = beta.pdf(mutU/prop_hs, a=alpha, b=beta_p)
            #hs_xy = 1-abs(.5-hs_xy)

            #exp_u = prop_hs

            #alpha, beta_p = get_beta_params(exp_u, mutU)
            #print exp_u, alpha, beta_p, skew
            #print exp_u, alpha, beta_p, mutU/hs
            #hs_yx = beta.pdf(mutU/hs, a=alpha, b=beta_p)
            #hs_yx =1-abs(.5-hs_yx)
            #h_yx = min(1-lognorm.cdf(h, s=prop_h, loc=0, scale=prop_h),lognorm.cdf(h, s=prop_h, loc=0, scale=prop_h))
            pi_x = beta.pdf(curr_freq, a=samp_size*obs + 1, b=samp_size-(samp_size*obs)+1)
            #pi_x = 1-abs(.5-pi_x)
            #pi_x = min(1-lognorm.cdf(curr_freq, s=.1, loc=0, scale=obs), lognorm.cdf(curr_freq, s=.1, loc=0, scale=obs))
            pi_y = beta.pdf(sim_result[0], a=samp_size*obs + 1, b=samp_size-(samp_size*obs)+1)
            #pi_y = 1-abs(.5-pi_y)
            #pi_y = min(1-lognorm.cdf(sim_result[0], s=.1, loc=0, scale=obs),lognorm.cdf(sim_result[0], s=.1, loc=0, scale=obs))

            Q_xy = s_xy * h_xy
            Q_yx = s_yx * h_yx
            #Q_xy = hs_xy
            #Q_yx = hs_yx
            if pi_x*Q_xy > 0:
                mh_r = float(pi_y*Q_yx)/(pi_x*Q_xy)
            else:
                mh_r = 1
            if math.isnan(mh_r) == True: mh_r = 1
            ratio = min(1, mh_r)
            # #print pi_y, Q_yx, pi_x, Q_xy
            prob = rm()
            #print prob, ratio
            if prob <= ratio:
                s, h = prop_s, prop_h
                hs = s*h
                #hs = prop_hs
                accept += 1
                curr_freq = sim_result[0]
                # pi_x = lognorm.cdf(curr_freq, s=.1, loc=0, scale=obs)
                # pi_x = min(1-pi_x, pi_x)
                #diff = prop_diff
            #print s, h, hs, sim_result[0], obs, prob, ratio, mh_r, prop_s, prop_h, accept, Q_xy, Q_yx, pi_x, pi_y, sim_result[-1]
            #print hs, sim_result[0], obs, prob, ratio, mh_r, prop_hs, accept, Q_xy, Q_yx, pi_x, pi_y
            hs_L.append(hs)
        hpd_mu, x_mu, y_mu, modes_mu = hpd_grid(hs_L, roundto=6)
        log_hpd_mu, log_x_mu, log_y_mu, log_modes_mu = hpd_grid(np.log10(hs_L), roundto=6)
        print infile, obs, mutU, np.mean(hs_L), np.percentile(hs_L,2.5), np.percentile(hs_L,97.5), np.median(hs_L), modes_mu[0], hpd_mu[0][0], hpd_mu[0][1], log_modes_mu[0], log_hpd_mu[0][0], log_hpd_mu[0][1]

if __name__ == '__main__':
    main()
