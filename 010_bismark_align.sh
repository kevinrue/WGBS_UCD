#!/bin/bash

if [ $# -lt 4 ]; then
	echo "Usage: $0 <genome> <rootdir> <outdir> <theads>"
	exit 1
fi

cwd=`pwd`
echo "cwd: $cwd"

genome=$1
echo "genome: $genome"
rootdir=$2
echo "rootdir: $rootdir"
outdir=$3
echo "outdir: $outdir"
threads=$4
echo "threads: $threads"

# Finds folders that contain fq.gz files
folders=`find $rootdir -name "*fq.gz" -exec dirname {} \; | sort | uniq`
echo -e "folders (next line):\n$folders"

if [ ! -e $outdir ]; then
	mkdir -pv $outdir
fi

for folder in `echo $folders`
do
	echo "folder: $folder"
	# Need to preserve batches of sequencing
	# as some fastq files for the same sample have the same name
	batch=$(basename $folder)
	if [ ! -e $outdir/$batch ]; then
		mkdir -pv $outdir/$batch
	fi
	# Guess whether this is a single or paired-end batch
	single=$(ls $folder | grep -c 'trimmed.fq.gz')
	echo "single: $single"
	cmd="bismark -n 1 $genome/ --bowtie1 --unmapped --ambiguous \
		--output_dir $outdir/$batch --multicore $threads"
	if [ $single -gt 0 ]; then
		# Identify all the NOT_BS reads in the folder
		fastqs=$(find $folder -name '*trimmed.fq.gz' | grep -v 'NOT_BS' | \
			tr '\n' ',' | sed -e 's/,$//')
		#echo $fastqs
		cmd_run="$cmd $fastqs"
	else
		# Identify all the forward paired NOT_BS reads in the folder
		paired1s=$(find $folder -name '*_R1_*_val_1*.fq.gz' | \
			grep -v 'NOT_BS' | tr '\n' ',' | sed -e 's/,$//')
		paired2s=$(echo $paired1s | \
			perl -pe 's/_R1_([[:digit:]]{3})_val_1/_R2_\1_val_2/g')
		unpaired=$(find $folder -name '*unpaired*' | grep -v 'NOT_BS' | \
			tr '\n' ',' | sed -e 's/,$//')
		#echo -e "paired1s (next line):\n$paired1s"
		#echo -e "paired2s (next line):\n$paired2s"
		#echo -e "unpaired (next line):\n$unpaired"
		cmd_run="$cmd -1 $paired1s -2 $paired2s; $cmd $unpaired"
	fi
	echo "cmd_run: $cmd_run"
	time(
		$cmd_run
	)
done
