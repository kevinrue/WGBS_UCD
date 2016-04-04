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
	bams_paired=$(find $folder -name '*val*bam' | xargs)
	bams_unpaired=$(find $folder -name '*unpaired*bam' | xargs)
	bams_single=$(find $folder -name '*trimmed*bam' | xargs)
	# Count how many elements in each array
	counts_paired=$(echo $bams_paired | wc -w)
	counts_unpaired=$(echo $bams_unpaired | wc -w)
	counts_single=$(echo $bams_single | wc -w)
	echo "paired: $counts_paired"
	echo "unpaired: $counts_unpaired"
	echo "single: $counts_single"
	cmd=""
	if [ $counts_paired -gt 0 ]
	then
		echo "parallel -j $threads --xapply $cmd -p ::: $bams_paired"
		time(
			parallel -j $threads --xapply $cmd -p ::: $bams_paired
		)
	fi
	if [ $counts_unpaired -gt 0 ]
	then
		echo "parallel -j $threads --xapply $cmd -s ::: $bams_unpaired"
		time(
			parallel -j $threads --xapply $cmd -p ::: $bams_unpaired
		)
	fi
	if [ $counts_single -gt 0 ]
	then
		echo "parallel -j $threads --xapply $cmd -s ::: $bams_single"
		time(
			parallel -j $threads --xapply $cmd -p ::: $bams_single
		)
	fi
done
