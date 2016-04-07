#!/bin/bash

echo "$0"
echo "$(date -I)"

if [ $# -lt 5 ]; then
	echo "Usage: $0 <genome> <rootdir> <outdir> <threads> <tempdir>"
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
tempdir=$5
echo "tempdir: $tempdir"

# Finds folders that contain fq.gz files
folders=`find $rootdir -name "*fq.gz" -exec dirname {} \; | sort | uniq`
echo -e "folders (next line):\n$folders"

if [ ! -e $outdir ]; then
	mkdir -pv $outdir
fi

if [ ! -e $tempdir ]; then
	mkdir -pv $tempdir
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
	cmd="bismark $genome/ --bowtie2 \
		--output_dir $outdir/$batch -p $threads --temp_dir $tempdir"
	# Identify all the first mates NOT_BS reads in the folder
	# C8_TGACCA_R1_merged_val_1.fq.gz
	paired1s=$(find $folder -name '*_R1_*_val_1*.fq.gz' | \
		grep -v 'NOT_BS' | tr '\n' ',' | sed -e 's/,$//')
	# Deduce the corresponding second mates
	paired2s=$(echo $paired1s | \
		perl -pe 's/_R1_merged_val_1/_R2_merged_val_2/g')
	#echo -e "paired1s (next line):\n$paired1s"
	#echo -e "paired2s (next line):\n$paired2s"
	#echo -e "unpaired (next line):\n$unpaired"
	cmd_run="$cmd -1 $paired1s -2 $paired2s"
	echo "cmd_run: $cmd_run"
	time(
		$cmd_run
	)
done

echo "Completed."
