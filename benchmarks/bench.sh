#!/bin/bash



DATA_DIR=/home/yoann/Projects/bioinfo/KFF/tools/data
DATASETS=( "sars.fa" "ecoli.fa" )
# DATASETS=( "ecoli.fa" )
KS=( 23 31 63 127 )
# KS=( 63 )
MIN_M=5
MAX_M=15

set -e

# Get Path to the scripts
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

OUT_DIR=${SCRIPT_DIR}/results
mkdir -p $OUT_DIR

# Compile
cd ${SCRIPT_DIR}/../
cmake . && make -j
cd -

# Prepare env
set +e
PATH_SAVE=$PATH
PATH=${SCRIPT_DIR}/../bin/:${PATH}

# Exec subroutines
for DATASET in ${DATASETS[@]};
do
	mkdir -p ${OUT_DIR}/${DATASET}
	for K in ${KS[@]};
	do
		echo "------------------ K=$K ------------------"
		bash ${SCRIPT_DIR}/generate_kff.sh ${DATA_DIR}/${DATASET} $K ${OUT_DIR}/${DATASET}/${DATASET}_k${K}.kff
		mkdir -p ${OUT_DIR}/${DATASET}/${DATASET}_k${K}
		for (( M=$MIN_M; M<=$MAX_M; M++ ));
		do
			echo $DATASET $K $M
			bash ${SCRIPT_DIR}/sort_kmers.sh ${OUT_DIR}/${DATASET}/${DATASET}_k${K}.kff $M \
				${OUT_DIR}/${DATASET}/${DATASET}_k${K}/${DATASET}_k${K}_m${M}_sorted.kff \
				${OUT_DIR}/${DATASET}/${DATASET}_k${K}/${DATASET}_k${K}_m${M}_compact.kff
		done
	done
done

# Reset env
PATH=$PATH_SAVE
