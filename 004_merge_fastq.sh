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

#R2files=$(find $rootdir -name "*_R2_*fastq.gz")
#echo -e "R2files (next lines):\n$R2files"
#echo "R2files: $(echo $R2files | wc -w)"

#R1files=$(echo $R2files | sed -e 's/_R2_/_R1_/g')
#echo -e "R1files (next lines):\n$R1files"
#echo "R1files: $(echo $R1files | wc -w)"

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
	# (Over)Write a script to concatenate the different batches of each sample
	scriptR1=$(echo "$outdir/script_${sample}_R1.sh")
	scriptR2=$(echo "$outdir/script_${sample}_R2.sh")
	echo -n "" > "$scriptR1"
	echo -n "" > "$scriptR2"
	echo "zcat $R1files | gzip -c > $R1outfile" > "$scriptR1"
	echo "zcat $R2files | gzip -c > $R2outfile" > "$scriptR2"
done

scripts=$(find $outdir -name 'script*sh' | xargs)
echo "scripts: $(echo $scripts | wc -w)"

parallel -j $threads bash ::: $scripts

rm -v $scripts


