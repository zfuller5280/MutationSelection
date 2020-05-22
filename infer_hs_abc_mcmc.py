#Import libraries
from optparse import OptionParser
import numpy as np
from random import random as rm
import subprocess
from subprocess import Popen, PIPE
from operator import add, sub
import random
from scipy.stats import beta
import math
from hpd import hpd_grid

#Get options from command line
USAGE= """Usage: %prog -i name -r runs -m mutation-rate -o observed-frequency"""
OPT_DEFAULTS={'infile':'-', '-r':1000, '-m':5e-8, '-o':1e-3}
DESCRIPTION="""Program description: This program uses an mcmc approach to simulate
                allele frequencies and infer heterozygous selection coefficients"""
EPILOG="""Requirements: numpy, scipy, hpd"""

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

#Function used to make a call to C++ program to simulate a single run of an allele frequency under a given demographic model
#Default is to use a constant population size (demographic_model=1), change to 0 to run under Schiffels-Durbin model
#Default sample size is 113770 (NFE for gnomad)
def run_simulation(runs, sel, dom, mutU, demographic_model=1, sample_size=113770):
    #Run a single iteration of a forward population genetics simulation
    #Order of arguments is mutation rate, selection, dominance coefficient, runs (1), mutation uncertainty, demographic uncertainty, demographic model, Ne (for contstant size)
    cmd = ["/Users/zachfuller/mut_population_sims/autosome_sim","%s"%mutU,"%s"%sel,"%s"%dom,"%i"%runs,"0","0","%i"%demographic_model,"100000"]
    result = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    #Get the result written to stdout
    output = result.stdout.readlines()

    out = [float(x) for x in output[0].split()[2:6]]
    freq = ((out[0]*2) + out[1])/(sum(out[:-1])*2)
    #Add the selection and dominance coefficient to the output for record keeping
    out.extend([sel, dom])

    #Sample an allele frequency from a binomial distribution matched to sample size
    pop_samp = np.random.binomial(113770, freq, 1)
    pop_freq = float(pop_samp)/113770
    #What was the mutation rate used in the simulation
    mut_u = float(output[0].split()[-1].strip(b"\n"))
    return pop_freq, out, mut_u

#Function to get a random number from a log-uniform distribution
def lognuniform(low=1e-6, high=1, size=None, base=np.exp(1)):
    return np.power(base, np.random.uniform(np.log10(low), np.log10(high), size))

#Set a min and max for a value
def clamp(n, minn=1e-4, maxn=0.99):
    return max(min(maxn, n), minn)

#Function to get kernal density of an observed value in a distribution
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

#Get the alpha and beta parameters for a beta distribution for a given sample size
def get_beta_params(exp_u, mutU, sample_size=113770):

    # return alpha, beta_p, skew
    mut_samp = mutU/exp_u
    alpha = mut_samp * sample_size
    beta_p = sample_size + 1
    return alpha+1, beta_p

#Get the alpha and beta parameters for a beta distribution
def get_proposal_beta_params(hs):

    # u = hs
    # v = (1./500)*hs
    # alpha = ((1.-u)/v-1./u)*(u**2)
    # beta = alpha * (1./u - 1.)
    scale_param = (1./hs) * (-np.log10(hs))
    alpha = hs * scale_param
    beta = (1-hs) * scale_param
    return alpha, beta

#Run the program
def main():
    #Get the options from the command line
    (options,args)=get_options(OPT_DEFAULTS, USAGE, DESCRIPTION, EPILOG)
    infile = options.infile
    runs = options.runs
    mutU = options.mut
    obs = options.obs
    if obs < 0:
        print(infile, "NA", "NA","NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
    else:
        ops = (add, sub)
        #Get a starting point for hs based on deterministic expectation
        if obs > 0:
            init_hs = mutU/obs
        else:
            init_hs = float(mutU)/(float(1)/113770)
        #If observed frequency is less than the mutation rate, set to a maximum value of hs given sample size
        if init_hs > 1: init_hs = float(mutU)/(float(1)/113770)

        #Get the alpha and beta parameters for initial beta distribution of hs and sample size
        alpha, beta_p = get_beta_params(init_hs, mutU)

        #Draw a random hs value from the initial distribution to simulate
        init_hs = clamp(mutU/np.random.beta(a=alpha,b=beta_p))
        #Set h and s so they are equal at initial starting point
        init_s, init_h = (init_hs)**.5, (init_hs)**.5
        s, h, hs = init_s, init_h, init_hs

        #Get the initial frequency from a simulation
        init_freq = run_simulation(1, s, h, mutU)

        #Start keeping track of the number of accepted proposals in the MCMC
        accept = 0

        #What is the current frequency in the MCMC
        curr_freq = init_freq[0]
        samp_size = 113770

        #Calculate likelihood of observing the current simulated allele frequency under binomial sampling given the actual sample size and observed frequency
        pi_x = beta.pdf(curr_freq, a=samp_size*obs + 1, b=samp_size-(samp_size*obs)+1)

        hs_L = []

        #Start the MCMC chain
        for run in range(runs):

            #Get the beta distribution parameters to sample s and h from
            prop_as, prop_bs = get_proposal_beta_params(s)
            prop_ah, prop_bh = get_proposal_beta_params(h)

            #Sample s and h from a beta distribution with parameters from above
            prop_s = clamp(np.random.beta(a=prop_as, b=prop_bs))
            prop_h = clamp(np.random.beta(a=prop_ah, b=prop_bh))

            #Simulate an allele frequency with the sampled s and h
            sim_result = run_simulation(1, prop_s, prop_h, mutU)

            #Get the likelihood of drawing s from the proposed distribution
            s_xy = beta.logpdf(prop_s, a=prop_as, b=prop_bs)
            #Get the likelihood of drawing h from the proposed distribution
            h_xy = beta.logpdf(prop_h, a=prop_ah, b=prop_bh)

            #Get the likelihood of drawing the current s and h values from the proposed distribution for the transition probabilities
            prop_as, prop_bs = get_proposal_beta_params(prop_s)
            prop_ah, prop_bh = get_proposal_beta_params(prop_h)
            s_yx = beta.logpdf(s, a=prop_as, b=prop_bs)
            h_yx = beta.logpdf(h, a=prop_ah, b=prop_bh)

            #Get the likelihood of observing the current allele frequency in the MCMC under binomial sampling given the actual sample size and observed frequency
            pi_x = beta.logpdf(curr_freq, a=samp_size*obs + 1, b=samp_size-(samp_size*obs)+1)

            #Get the likelihood of observing the simulated allele frequency in the MCMC under binomial sampling given the actual sample size and observed frequency
            pi_y = beta.logpdf(sim_result[0], a=samp_size*obs + 1, b=samp_size-(samp_size*obs)+1)

            #Multiply likelihoods of drawing the proposed s and h values to get transition probability of moving from x to y
            Q_xy = s_xy * h_xy
            #Multiply likelihoods of drawing the current s and h values from the proposal distribution to get transition probability of moving from y to x
            Q_yx = s_yx * h_yx

            #Get the metropolis-hastings ratio
            if pi_x*Q_xy > 0:
                mh_r = float(pi_y*Q_yx)/(pi_x*Q_xy)
            else:
                mh_r = 1
            if math.isnan(mh_r) == True: mh_r = 1
            ratio = min(1, mh_r)

            #Get a random uniform number between 0 and 1
            prob = rm()

            #Accept the proposed paramters with proability equal to the metropolis-hastings ratio
            if prob <= ratio:
                s, h = prop_s, prop_h
                hs = s*h

                accept += 1
                curr_freq = sim_result[0]

            #Print the current state of the chain (for debugging and exploring single observed values)
            print(s, h, hs, sim_result[0], obs, prob, ratio, mh_r, prop_s, prop_h, accept, Q_xy, Q_yx, pi_x, pi_y, sim_result[-1])

            hs_L.append(hs)

        #Get the MAP estimate of hs from the posterior distribution and print to stdout (uncomment lines when running for multiple genes or observed values)
        hpd_mu, x_mu, y_mu, modes_mu = hpd_grid(hs_L, roundto=6)
        log_hpd_mu, log_x_mu, log_y_mu, log_modes_mu = hpd_grid(np.log10(hs_L), roundto=6)
        #print(infile, obs, mutU, np.mean(hs_L), np.percentile(hs_L,2.5), np.percentile(hs_L,97.5), np.median(hs_L), modes_mu[0], hpd_mu[0][0], hpd_mu[0][1], log_modes_mu[0], log_hpd_mu[0][0], log_hpd_mu[0][1])

if __name__ == '__main__':
    main()
