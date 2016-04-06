#!/bin/bash

echo "$0"
echo "$(date -I)"

if [ $# -lt 2 ]; then
	echo "Usage: $0 <rootdir> <outdir> <threads>"
	exit 1
fi

cwd=`pwd`
echo "cwd: $cwd"

rootdir=$1
echo "rootdir: $rootdir"
outdir=$2
echo "outdir: $outdir"
threads=$3
echo "threads: $threads"

if [ ! -e $outdir ]; then
	mkdir -pv $outdir
fi

# Finds folders that contain BAM files
folders=`find $rootdir -name "*bam" -exec dirname {} \; | sort | uniq`
echo -e "folders (next line):\n$folders"

for folder in `echo $folders`
do
	echo "folder: $folder"
	# Need to preserve batches of sequencing
	# as some fastq files for the same sample have the same name
	batch=$(basename $folder)
	if [ ! -e $outdir/$batch ]; then
		mkdir -pv $outdir/$batch
	fi
	# Identify all the BAM files in the folder
	bams_paired=$(find $folder -name '*val*deduplicated.bam' | xargs)
	bams_unpaired=$(find $folder -name '*unpaired*deduplicated.bam' | xargs)
	bams_single=$(find $folder -name '*trimmed*deduplicated.bam' | xargs)
	# Count how many elements in each array
	counts_paired=$(echo $bams_paired | wc -w)
	counts_unpaired=$(echo $bams_unpaired | wc -w)
	counts_single=$(echo $bams_single | wc -w)
	echo "paired: $counts_paired"
	echo "unpaired: $counts_unpaired"
	echo "single: $counts_single"
	cmd="bismark_methylation_extractor --output $outdir/$batch --bedGraph \
		--buffer_size 10G --scaffolds --cytosine_report \
		--genome_folder bostaurus"
	if [ $counts_paired -gt 0 ]
	then
		cmd="$cmd --no_overlap --paired-end "
		echo "parallel -j $threads --xapply $cmd ::: $bams_paired"
		time(
			parallel -j $threads --xapply $cmd -p ::: $bams_paired
		)
	fi
	if [ $counts_unpaired -gt 0 ]
	then
		cmd="$cmd --single-end"
		echo "parallel -j $threads --xapply $cmd ::: $bams_unpaired"
		time(
			parallel -j $threads --xapply $cmd -p ::: $bams_unpaired
		)
	fi
	if [ $counts_single -gt 0 ]
	then
		cmd="$cmd --single-end"
		echo "parallel -j $threads --xapply $cmd ::: $bams_single"
		time(
			parallel -j $threads --xapply $cmd -p ::: $bams_single
		)
	fi
done

echo "Completed."



