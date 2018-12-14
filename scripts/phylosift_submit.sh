#!/bin/bash
#
# Using the PBS queue, run Phylosift search and align on a set of fasta sequences.
#
#PBS -l select=1:ncpus=2:mem=8gb
#PBS -q smallq

# Default batch size
DEFAULT_BATCH_SIZE=5

# Phylosift executable
PSBIN=/shared/homes/11025234/phylosift_v1.0.1/bin/phylosift

function printHelp () {
	echo "Usage: [OPTIONS] <input dir> <output dir>"
	echo ""
	echo "Options:"
	echo "  -b N    How many jobs should be batched together (integer)"
	echo "  -h      Help"
	echo ""
	echo "The input directory is expected to contain fasta format sequences"
	echo "whose filenames end in .fasta"
	echo ""
}


#
# Check if invocation was via PBS queue.
# 
# If this is not being run from the queue then
# act as a regular shell script and provide a 
# user interface. Successful invocation will
# then cause the script to submit itself to
# the queue.
#
if [ -z "$PBS_ENVIRONMENT" ]
then
	####################
	# Regular CLI mode #
	####################

	while getopts ":hb:" opt
	do
		case $opt in
			b)
				BATCH_SIZE=$OPTARG
				;;
			h)
				printHelp
				exit 0
				;;
			\?)
				echo "Invalid option: -$OPTARG"
				;;
		esac
	done

	shift $(( OPTIND-1 ))

	if [ $# -ne 2 ]
	then
		printHelp
		exit 1
	fi
	
	if [ ! -d $1 ]
	then
		echo "Input directory [$1] does not exist"
		exit 1
	fi
	
	if [ -d $2 ]
	then
		echo -e "Output directory [$2] already exists\nMove away or delete"
		exit 1
	else
		mkdir -p $2
	fi

	if (( BATCH_SIZE <= 0 ))
	then
		BATCH_SIZE=$DEFAULT_BATCH_SIZE
	fi
	
	# Prepare list of input files to process
	INPUT_FILES="ps.input."`date +%s.txt`
	find $1 -name '*.fasta' > $INPUT_FILES
	TOTAL_JOBS=`wc -l $INPUT_FILES | cut -d ' ' -f 1`
	if (( TOTAL_JOBS <= 0 ))
	then
		echo "No .fasta files were found to process in input folder [$1]"
		exit 1
	fi

	# Submit batch job
	qsub -J 1-$TOTAL_JOBS:$BATCH_SIZE -N PSJOB -v INPUT_FILES=$INPUT_FILES,OUTPUT_TOP=$2,BATCH_SIZE=$BATCH_SIZE $0

else
	##############
	# Queue mode #
	##############

	cd $PBS_O_WORKDIR

	INPUT_FILES=($(<$INPUT_FILES))
	TOTAL_FILES=${#INPUT_FILES[*]}

	n=0
	while ((n < $BATCH_SIZE && n+PBS_ARRAY_INDEX <= TOTAL_FILES))
	do
		# Convert to zero-numbered array index
		IDX=$((n+PBS_ARRAY_INDEX - 1))
		
		OUTPUT_DIR=$OUTPUT_TOP/`basename ${INPUT_FILES[$IDX]}`"."`date +%S.%N`
		mkdir -p $OUTPUT_DIR
		
		# Run phylosift search and align
		$PSBIN search --threads 2 --debug --isolate --besthit --output $OUTPUT_DIR ${INPUT_FILES[$IDX]}
		$PSBIN  align --threads 2 --debug --isolate --besthit --output $OUTPUT_DIR ${INPUT_FILES[$IDX]}

		((n++))
	done

fi
