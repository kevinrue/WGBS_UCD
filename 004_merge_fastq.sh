#!/bin/bash

echo "$0"
echo "$(date -I)"

if [ $# -lt 3 ]; then
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

# Identify the sample name for all paired-end libraries
samples=$(find $rootdir -name "*_R2_*fastq.gz" -exec basename {} \; | perl -pe 's/([CM][[:digit:]]{1,2}[_NOTBSATGC]*)_.*/\1/g' | sort | uniq)
#echo -e "samples (next lines):\n$samples"

if [ ! -e $outdir ]; then
	mkdir -pv $outdir
fi

for sample in $(echo $samples)
do
	echo "sample: $sample"
	# Define the name of the merged fastq file
	R1outfile="$outdir/${sample}_R1_merged.fastq.gz"
	R2outfile="$outdir/${sample}_R2_merged.fastq.gz"
	# Find the second mates (identify the paired-end runs)
	R2files=$(find $rootdir -name "$sample*_R2_*fastq.gz" | xargs)
#	echo -e "R2files: $(echo $R2files | wc -w)"
	# Deduce the first mates
	R1files=$(echo $R2files | sed -e 's/_R2_/_R1_/g')
#	echo -e "R1files: $(echo $R1files | wc -w)"
	echo "cat $R1files > $R1outfile"
	cat $R1files > $R1outfile
	echo "cat $R2files > $R2outfile"
	cat $R2files > $R2outfile
done

echo "Completed."


