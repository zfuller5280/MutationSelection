C++ source code for modules used in ABC implemented in pakman

-simulator.cpp: Main code to run forward simulations to get an allele frequency for an autosomal model. Reads in proposed selection coefficient and epsilon from stdin. Return 'accept' or 'reject' if the proposed parameter is accepted or rejected based on the distance function selected.

-x_simulator.cpp: Main code to run forward simulations to get an allele frequency for an X-chromosome model. Reads in proposed selection coefficient and epsilon from stdin. Return 'accept' or 'reject' if the proposed parameter is accepted or rejected based on the distance function selected.

-perturbation-normal.cpp: Perturbation kernel for ABC-SMC based on a log10-normal distribution. STDEV is included as an argument.

-perturbation-normal-pdf.cpp: Calculate the PDF for samples drawn from the log10-normal distribution perturbation kernel.

-prior_sampler.cpp: Defines a log10-uniform distribution as the prior

-prior-pdf.cpp: Calculates the PDF for a parameter proposed from a log10-uniform distribution 
