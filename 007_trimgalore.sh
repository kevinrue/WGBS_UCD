#!/bin/bash

echo "$0"
echo "$(date -I)"

if [ $# -lt 3 ]; then
	echo "Usage: $0 <outdir> <threads> <rootdirs ...>"
	exit 1
fi

cwd=`pwd`
echo "cwd: $cwd"

md5file=md5.txt

outdir=$1; shift
echo "outdir: $outdir"
threads=$1; shift
echo "threads: $threads"
rootdirs=$@
echo "rootdirs: $rootdirs"

folders=`find $rootdirs -name "*fastq.gz" -exec dirname {} \; | sort | uniq`
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
	# Identify all the forward reads in the folder
	fastq1s=$(find $folder -name '*_R1_*.fastq.gz')
	echo -e "fastq1s:\n$fastq1s"
	# Identify whether reverse reads are present in the folder
	paired=$(find $folder -name '*_R2_*.fastq.gz' | wc -l)
	echo "paired: $paired"
	cmd="trim_galore"
	if [ $paired -gt 0 ]
	then
		cmd="$cmd --paired --retain_unpaired"
	fi
	cmd="$cmd --output_dir $outdir/$batch"
	fastq1s=$(echo $fastq1s | xargs)
	echo "fastq1s count: $(echo $fastq1s | wc -w)"
	if [ $paired -gt 0 ]
	then
		fastq2s=$(echo $fastq1s | perl -p -e 's/_R1_/_R2_/g' | xargs)
		echo "fastq2s count: $(echo $fastq2s | wc -w)"
		echo "parallel -j $threads --xapply $cmd ::: $fastq1s ::: $fastq2s"
		time(
			parallel -j $threads --xapply $cmd ::: $fastq1s ::: $fastq2s
		)
	else
		echo "parallel -j $threads $cmd ::: $fastq1s"
		time(
			parallel -j $threads $cmd ::: $fastq1s
		)
	fi
done


