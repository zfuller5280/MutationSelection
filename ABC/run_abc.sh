#!/bin/bash

## Input order is SYMBOL,chromosome,mu_lof,AF_nfe,AN_nfe,NFE_k

# NAME=$1
# CHROM=$2
# MUT=$3
# FREQ=$4
# SAMP_SIZE=$5
# NUM=$6
#
# NPROC=$7
# NACCEPT=$8

date
echo "Running hs inference with ABC-SMC"
echo "jobId: $AWS_BATCH_JOB_ID"
echo "jobQueue: $AWS_BATCH_JQ_NAME"
echo "computeEnvironment: $AWS_BATCH_CE_NAME"

#echo $NAME
CHECK_FREQ=$(printf "%.14f" $FREQ)
if (( $(echo "$CHECK_FREQ <= 0.0" |bc -l) )); then
	FREQ='1e-6'
fi

MUT_RATE=$(echo "$MUT"|awk '{printf "%.10f \n",$1}')
if (( $(echo "${MUT_RATE} <= 0.0"|bc -l) )); then
        echo "Error: Mutation rate must be positive"
        exit 1
fi

EXP_HS=$(echo "$MUT $FREQ"|awk '{printf "%.7f \n",$1/$2}')

echo "Infer hs for gene: " $NAME
echo "Chromsome: " $CHROM
echo "Mutation rate: " $MUT
echo "Observed LoF Frequency: " $FREQ
echo "Haploid sample size: " $SAMP_SIZE
echo "Observed total LoF number: " $NUM
echo "Number of processors requested for mpi: " $NPROC
echo "Number of acceptances required: " $NACCEPT
#echo "EXPECTED HS: ${EXP_HS}"

if [[ "${CHROM}" == "X" ]]; then
				echo "Running for X Chromosome model..."
        mpiexec -n ${NPROC} pakman mpi smc --discard-child-stderr --population-size=${NACCEPT} \
				 --epsilons=10,5,2,0 --parameter-names=s,h --simulator="./x_simulator ${MUT} ${NUM} 1 0 0 1 1 11 1 ${SAMP_SIZE}" \
				 --prior-sampler="./prior_sampler_x -6 0 0 1" --perturber="./perturber_x 1 0.1" \
				 --prior-pdf="./prior-pdf_x 1e-6 1 0 1" --perturbation-pdf="./perturbation-pdf_x 1 0.1" > outfiles/${NAME}.x.gnomad.canonical.posterior.txt
else
				echo "Running for autosomal model..."
				mpiexec -n ${NPROC} pakman mpi smc --discard-child-stderr --population-size=${NACCEPT} \
				--epsilons=10,5,2,0 --parameter-names=hs --simulator="./simulator ${MUT} ${NUM} .5 1 0 0 1 1 11 1 ${SAMP_SIZE}" \
				--prior-sampler="./prior_sampler -6 0" --perturber="./perturber-ln-normal 1" --prior-pdf="./prior-pdf 1e-6 1" \
				--perturbation-pdf="./perturbation-pdf-ln-normal 1" > outfiles/${NAME}.gnomad.canonical.posterior.txt
fi
