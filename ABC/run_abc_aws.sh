#!/bin/bash

## Input order for infile is SYMBOL,chromosome,mu_lof,AF,AN,k

MAIN_S3_URL="s3://abc-hs"

INPUT_ARRAY_FILE=$1
COMPLETED_GENE_LIST=$4
echo $INPUT_ARRAY_FILE
TMPDIR="$(mktemp -d -t tmp.XXXXXXXXX)" || error_exit "Failed to create temp directory."
TMPFILE="${TMPDIR}/batch-file-temp"
COMPLETEDFILE="${TMPDIR}/completed-file-temp"
#install -m 0600 /dev/null "${TMPFILE}" || error_exit "Failed to create temp file."
LINE=$((AWS_BATCH_JOB_ARRAY_INDEX + 1))
aws s3 cp "${MAIN_S3_URL}/${INPUT_ARRAY_FILE}" - > "${TMPFILE}" || error_exit "Failed to download S3 script."
aws s3 cp "${MAIN_S3_URL}/${COMPLETED_GENE_LIST}" - > "${COMPLETEDFILE}" || error_exit "Failed to download S3 script."
INPUT_LINE=$(sed -n ${LINE}p "${TMPFILE}")
echo $INPUT_LINE
INPUT_LINE_ARRAY=($INPUT_LINE)
NAME=${INPUT_LINE_ARRAY[0]}
CHROM=${INPUT_LINE_ARRAY[1]}
MUT=${INPUT_LINE_ARRAY[2]}
FREQ=${INPUT_LINE_ARRAY[3]}
SAMP_SIZE=${INPUT_LINE_ARRAY[4]}
NUM=${INPUT_LINE_ARRAY[5]}
NPROC=$3
NACCEPT=$2
echo "Setting up compute environment for gene: " $NAME
echo "Chromosome: " $CHROM
echo "Mutation rate: " $MUT
echo "Observed LoF Frequency: " $FREQ
echo "Haploid sample size: " $SAMP_SIZE
echo "Observed total LoF number: " $NUM
echo "The number of acceptances for each SMC run will be: " $NACCEPT


PATH=${PATH}:"pakman/build/src/"
MAIN_S3_URL="s3://abc-hs/expanded_output_data"

exists=$(aws s3 ls ${MAIN_S3_URL}/expanded_output_data/${CHROM}/${NAME}/${NAME}.${CHROM}.mpi.gnomad.canonical.posterior.txt)
if [ -z "$exists" ]; then
	while read -r y z;do
		if [[ "$y" == ${NAME} ]]; then
			echo "The gene ${NAME} already exists at $z, exiting"
			exit 1
		fi
	done < ${COMPLETEDFILE}

  echo "the posterior distriubiton does not exist...proceeding!"
	date
	echo "Running hs inference with ABC-SMC"
	echo "jobId: $AWS_BATCH_JOB_ID"
	echo "jobQueue: $AWS_BATCH_JQ_NAME"
	echo "computeEnvironment: $AWS_BATCH_CE_NAME"
	echo "Array index:" $AWS_BATCH_JOB_ARRAY_INDEX
	echo "default s3 bucket:"

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
	echo "Under deterministic approximation, the expected hs is: " $EXP_HS
	echo "Chromosome: " $CHROM
	echo "Mutation rate: " $MUT
	echo "Observed LoF Frequency: " $FREQ
	echo "Haploid sample size: " $SAMP_SIZE
	echo "Observed total LoF number: " $NUM
	echo "Number of processors requested for mpi: " $NPROC
	echo "Number of acceptances required: " $NACCEPT

	if [[ "${CHROM}" == "X" ]]; then
					echo "Running for X Chromosome model..."
	        mpirun -np ${NPROC} pakman mpi smc --discard-child-stderr --population-size=${NACCEPT} \
					 --epsilons=10,5,2,0 --parameter-names=s,h --simulator="./x_simulator_expand_smallback ${MUT} ${NUM} 1 0 0 1 1 11 1 ${SAMP_SIZE}" \
					 --prior-sampler="./prior_sampler_x -6 0 0 1" --perturber="./perturber_x 1 0.1" \
					 --prior-pdf="./prior-pdf_x 1e-6 1 0 1" --perturbation-pdf="./perturbation-pdf_x 1 0.1" > outfiles/${NAME}.${CHROM}.gnomad.canonical.posterior.txt
	else
					echo "Running for autosomal model..."
					mpirun -np ${NPROC} pakman mpi smc --discard-child-stderr --population-size=${NACCEPT} \
					--epsilons=10,5,2,0 --parameter-names=hs --simulator="./simulator_expand_smallback ${MUT} ${NUM} .5 1 0 0 1 1 11 1 ${SAMP_SIZE}" \
					--prior-sampler="./prior_sampler -6 0" --perturber="./perturber-ln-normal 1" --prior-pdf="./prior-pdf 1e-6 1" \
					--perturbation-pdf="./perturbation-pdf-ln-normal 1" > outfiles/${NAME}.${CHROM}.gnomad.canonical.posterior.txt
	fi

	aws s3 cp --region us-east-2 "outfiles/${NAME}.${CHROM}.gnomad.canonical.posterior.txt" "${MAIN_S3_URL}/expanded_output_data/${CHROM}/${NAME}/${NAME}.${CHROM}.mpi.gnomad.canonical.posterior.txt"

else
  echo "the file already exists, exiting"
	exit 1
fi
