#!/bin/bash

if [ $# -lt 3 ]; then
	echo "Usage: $0 <rootdir> <outdir> <threads>"
	exit 1
fi

cwd=`pwd`
echo "cwd: $cwd"

md5file=md5.txt

rootdir=$1
echo "rootdir: $rootdir"
outdir=$2
echo "outdir: $outdir"
threads=$3
echo "threads: $threads"

folders=`find $rootdir -name "$md5file" -exec dirname {} \;`
echo -e "folders (next line):\n$folders"

#log=$(basename $0 .sh)_$(date -I)
#echo "log: $log"

#if [ ! -e $rootdir/log ]; then
#	mkdir -pv $rootdir/log
#fi

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
		cmd="$cmd --paired --trim1"
	fi
	cmd="$cmd --output_dir $outdir/$batch --illumina --retain_unpaired"
	fastq1s=$(echo $fastq1s | xargs)
	if [ $paired -gt 0 ]
	then
		fastq2s=$(echo $fastq1s | perl -p -e 's/_R1_/_R2_/' | xargs)
		#echo -e "fastq2s:\n$fastq2s"
		echo "parallel -j $threads --xapply $cmd ::: $fastq1s ::: $fastq2s"
		time(
			parallel -j $threads --xapply $cmd ::: $fastq1s ::: $fastq2s
		)
	else
		echo "parallel -j $threads --xapply $cmd ::: $fastq1s"
		time(
		parallel -j $threads --xapply $cmd ::: $fastq1s
		)
	fi
#	for fastq1 in $fastq1s
#	do
#		echo "$cmd $fastq1 $(echo $fastq1 | perl -p -e 's/_R1_/_R2_/')" | parallel -j 12
#	done
done


