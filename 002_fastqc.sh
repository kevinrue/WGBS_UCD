#!/bin/bash

echo "$0"
echo "$(date -I)"

if [ $# -lt 4 ]; then
	echo "Usage: $0 <rootdir> <target_file> <outdir> <threads>"
	exit 1
fi

cwd=`pwd`
echo "cwd: $cwd"

rootdir=$1
echo "rootdir: $rootdir"
target_file=$2
echo "target_file: $target_file"
outdir=$3
echo "outdir: $outdir"
threads=$4
echo "threads: $threads"

folders=`find $rootdir -name "$target_file" -exec dirname {} \; | sort | uniq`
echo -e "folders (next line):\n$folders"

if [ ! -e $outdir ]; then
	mkdir -pv $outdir
fi

for folder in `echo $folders`
do
	echo "folder: $folder"
	# Need to separate batches of sequencing
	# as some fastq files for the same sample have the same name
	batch=$(basename $folder)
	if [ ! -e $outdir/$batch ]; then
		mkdir -pv $outdir/$batch
	fi
	fastqs=$(find $folder -name "$target_file" | xargs)
	echo "Processing $(echo $fastqs | wc -w ) files..."
	echo "fastqc --quiet --nogroup --extract --threads $threads --outdir $outdir/$batch $fastqs"
	fastqc --quiet --nogroup --extract --threads $threads --outdir $outdir/$batch $fastqs
done

echo "Completed."
