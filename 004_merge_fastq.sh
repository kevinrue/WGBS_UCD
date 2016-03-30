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
	R2files=$(find $rootdir -name "$sample*_R2_*fastq.gz" | xargs)
	echo -e "R2files: $(echo $R2files | wc -w)"
	R1files=$(echo $R2files | sed -e 's/_R2_/_R1_/g')
	echo -e "R1files: $(echo $R1files | wc -w)"
#	echo "parallel -j $threads --xapply cat {} > \
#		$outdir/${sample}_R1_merged.fastqc.gz ::: $R1files"
	# Merge all R1 files
	echo "cat $R1files > $outdir/${sample}_R1_merged.fastq.gz"
	time(
		zcat $R1files | gzip -c > $outdir/${sample}_R1_merged.fastq.gz
	)
	# Merge all R2 files
	echo "cat $R2files > $outdir/${sample}_R2_merged.fastq.gz"
	time(
		zcat $R2files | gzip -c > $outdir/${sample}_R2_merged.fastq.gz
	)
done




