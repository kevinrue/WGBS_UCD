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
	# Need to separate batches of sequencing
	# as some fastq files for the same sample have the same name
	batch=$(basename $folder)
	if [ ! -e $outdir/$batch ]; then
		mkdir -pv $outdir/$batch
	fi
	fastqs=$(find $folder -name '*.fastq.gz' | xargs)
	echo "Processing $(echo $fastqs | wc -w ) files..."
	echo "fastqc --quiet --nogroup --extract --threads $threads --outdir $outdir/$batch $fastqs"
	fastqc --quiet --nogroup --extract --threads $threads --outdir $outdir/$batch $fastqs
done

