C++ source code for modules used in ABC implemented in pakman

-simulator.cpp: Main code to run forward simulations to get an allele frequency for an autosomal model. Reads in proposed selection coefficient and epsilon from stdin. Return 'accept' or 'reject' if the proposed parameter is accepted or rejected based on the distance function selected.

-x_simulator.cpp: Main code to run forward simulations to get an allele frequency for an X-chromosome model. Reads in proposed selection coefficient and epsilon from stdin. Return 'accept' or 'reject' if the proposed parameter is accepted or rejected based on the distance function selected.

-perturber-normal.cpp: Perturbation kernel for ABC-SMC based on a log10-normal distribution. STDEV is included as an argument.

-perturbation-normal-pdf.cpp: Calculate the PDF for samples drawn from the log10-normal distribution perturbation kernel.

-prior_sampler.cpp: Defines a log10-uniform distribution as the prior

-prior-pdf.cpp: Calculates the PDF for a parameter proposed from a log10-uniform distribution 

To compile the simulators: 
##X Chromosome
g++ -O3 -std=c++11 BRand.cpp population_sex_diff.cpp x_simulator.cpp -o x_simulator

##Autosome
g++ -O3 -std=c++11 BRand.cpp population.cpp simulator.cpp -o simulator

All other .cpp files for the modules can be compiled as usual. (e.g., g++ -O3 -std=c++11 prior-pdf.cpp -o prior-pdf)

To run using pakman (all using OpenMPI for parallelization, here illustrated using 5 processors):
1. Standard ABC rejection (here, an exact match with epsilon=0)
mpiexec -n 5 pakman mpi rejection --number-accept 10000 --epsilon 0 --parameter-names=q --simulator="./simulator 1e-6 ${NUM} .5 1 0 0 1 1 11 1 ${SAMP_SIZE}" --prior-sampler="./prior_sampler -6 0" > ${COUNT}.${S}.${H}.${MUT}.mpi.binomial.posterior.txt

2. ABC rejection using a binomial distance kernel
mpiexec -n 5 pakman mpi rejection --number-accept 10000 --epsilon 0 --parameter-names=q --simulator="./simulator 1e-6 ${NUM} .5 1 0 0 1 1 11 2 ${SAMP_SIZE}" --prior-sampler="./prior_sampler -6 0" > ${COUNT}.${S}.${H}.${MUT}.mpi.binomial.posterior.txt

3. SMC with decreasing tolerance (epsilon) based on absolute distance and log10-normal perturbation kernel
mpiexec -n 5 pakman mpi smc --discard-child-stderr --population-size=1000 --epsilons=10,5,2,0 --parameter-names=hs --simulator="./simulator ${MUT} ${NUM} .5 1 0 0 1 1 11 1 ${SAMP_SIZE}" --prior-sampler="./prior_sampler -6 0" --perturber="./perturber-normal 3" --prior-pdf="./prior-pdf 1e-6 1" --perturbation-pdf="./perturbation-pdf-normal 3" > ${NAME}.mpi.posterior.txt


