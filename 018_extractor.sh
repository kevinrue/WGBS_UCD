#!/bin/bash

echo "$0"
echo "$(date -I)"

if [ $# -lt 3 ]; then
	echo "Usage: $0 <rootdir> <outdir> <threads> [args]"
	exit 1
fi

cwd=`pwd`
echo "cwd: $cwd"

rootdir=$1
echo "rootdir: $rootdir"; shift
outdir=$1
echo "outdir: $outdir"; shift
threads=$1
echo "threads: $threads"; shift
args=$@
echo $args

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
	# Count how many elements in each array
	counts_paired=$(echo $bams_paired | wc -w)
	echo "paired: $counts_paired"
	cmd="bismark_methylation_extractor --output $outdir/$batch --bedGraph \
		--buffer_size 10G --scaffolds --cytosine_report --gzip \
		--genome_folder bostaurus --no_overlap --paired-end"
		echo "parallel -j $threads --xapply $cmd $args ::: $bams_paired"
		time(
			parallel -j $threads --xapply $cmd $args ::: $bams_paired
		)
done

echo "Completed."



