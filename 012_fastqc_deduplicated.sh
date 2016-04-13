#!/bin/bash

echo "$0"
echo "$(date -I)"

if [ $# -lt 4 ]; then
	echo "Usage: $0 <rootdir> <target_file> <outdir> <threads> [CSVfile]"
	exit 1
fi

cwd=`pwd`
echo "cwd: $cwd"

rootdir=$1
echo "rootdir: $rootdir"
target_file=$2
echo "target_file: $target_file"
outdir=$3
echo "outdir: $outdir"
threads=$4
echo "threads: $threads"

if [ -z $5 ]; then
	outfolder='log'
	CSVfile=$outfolder/"$(basename $0 | sed -e "s/\.sh/_$(date -I).csv/")"
else
	CSVfile=$2
	outfolder=$(dirname $CSVfile)
fi
echo "CSVfile: $CSVfile"

folders=`find $rootdir -name "$target_file" -exec dirname {} \; | sort | uniq`
echo -e "folders (next line):\n$folders"

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
	fastqs=$(find $folder -name "$target_file" | xargs)
	echo "Processing $(echo $fastqs | wc -w ) files..."
	echo "fastqc --quiet --nogroup --extract --threads $threads --outdir $outdir/$batch $fastqs"
	fastqc --quiet --nogroup --extract --threads $threads --outdir $outdir/$batch $fastqs
done

echo "Collating FastQC reports.."

echo "\"Value\",\"QC\",\"File\",\"Sample\",\"Read\",\"Treatment\",\
\"Infection\"" > $CSVfile

for folder in `ls $outdir`
do
	echo "folder: $folder"
	for summaryfile in `find $rootdir/$batch -name summary.txt`
	do
		awk 'BEGIN {FS="\t"; OFS="\",\""} \
		\
		{
			nf=split($3,fs,"_");
			if ($3 ~ /NOT_BS/){
				t="NOT BS"}
			else{
				t="BS"
			}
			if (fs[1] ~ /C/){
				i="Control"}
			else{
				i="M. bovis"
			}
			print "\""$1,$2,$3,fs[1],fs[nf-6],t,i"\""
		}' $summaryfile >> $CSVfile
	done
done

echo "Completed."
