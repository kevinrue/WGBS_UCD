#!/bin/bash

echo "$0"
echo "$(date -I)"

if [ $# -lt 4 ]; then
	echo "Usage: $0 <db_prefix> <rootdir> <outdir> <threads>"
	exit 1
fi

db_prefix=$1
echo "db_prefix: $db_prefix"
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
	cmd="bwa mem ${db_prefix} -t $threads"
	# Identify all the first mates NOT_BS reads in the folder
	fqFiles=$(find $folder -name '*_R1_*_val_1*.fq.gz' | \
	    grep 'NOT_BS')
	for fqFile1 in $fqFiles
	do
	    echo "fqFile1: $fqFile1"
	    # Deduce the corresponding second mates
	    fqFile2=$(echo $fqFile1 | \
		perl -pe 's/_R1_merged_val_1/_R2_merged_val_2/g')
	    echo "fqFile2: $fqFile2"
	    # Deduce the read group
	    RG=$(basename $fqFile1 | \
		perl -pe 's/([CM][[:digit:]]{1,2}_NOT_BS).*/\1/g')
	    echo "RG: $RG"
	    cmd_bwa="$cmd -M -R @RG\tID:$RG\tSM:$RG $fqFile1 $fqFile2"
	    cmd_samtools="samtools view -uS -"
	    echo "cmd_run: $cmd_bwa | $cmd_samtools > $outdir/$batch/$RG.bam"
	    time(
	    	$cmd_bwa | $cmd_samtools > $outdir/$batch/$RG.bam
	    )
	done
done

echo "Completed."
