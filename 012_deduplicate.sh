#!/bin/bash

echo "$0"
echo "$(date -I)"

if [ $# -lt 2 ]; then
	echo "Usage: $0 <rootdir> <threads>"
	exit 1
fi

cwd=`pwd`
echo "cwd: $cwd"

rootdir=$1
echo "rootdir: $rootdir"
threads=$2
echo "threads: $threads"

# Finds folders that contain BAM files
folders=`find $rootdir -name "*bam" -exec dirname {} \; | sort | uniq`
echo -e "folders (next line):\n$folders"

for folder in `echo $folders`
do
	echo "folder: $folder"
	# Identify all the BAM files in the folder
	bams=$(find $folder -name '*bam' | xargs)
	echo -e "BAMs (next lines):\n$bams"
	# Identify whether reverse reads are present in the folder
	paired=$(find $folder -name '*_val_1*.bam' | wc -l)
	echo "paired: $paired"
	cmd="deduplicate_bismark --bam"
	if [ $paired -gt 0 ]
	then
		cmd="$cmd -p"
	else
		cmd="$cmd -s"
	fi
	echo "parallel -j $threads --xapply $cmd ::: $bams"
	time(
		parallel -j $threads --xapply $cmd ::: $bams
	)
done
