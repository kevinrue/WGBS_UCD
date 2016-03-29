#!/bin/bash

echo "$0"
echo "$(date -I)"

if [ $# -lt 1 ]; then
	echo "Usage: $0 <rootdir> [CSVfile]"
	exit 1
fi

cwd=`pwd`
echo "cwd: $cwd"

rootdir=$1
echo "rootdir: $rootdir"

if [ -z $2 ]; then
	outfolder='log'
	CSVfile=$outfolder/"$(basename $0 | sed -e "s/\.sh/_$(date -I).csv/")"
else
	CSVfile=$2
	outfolder=$(dirname $CSVfile)
fi
echo "CSVfile: $CSVfile"

if [ ! -e $outfolder ]; then
	mkdir -pv $outfolder
fi

echo "\"Value\",\"QC\",\"File\",\"Batch\",\"Lane\",\"Sample\",\
\"Infection\",\"Unpaired\"" > $CSVfile

for folder in `ls $rootdir`
do
	echo "folder: $folder"
	batch=$(basename $folder)
	echo "batch: $batch"
	for summaryfile in `find $rootdir/$batch -name summary.txt`
	do
		awk -v batch=$batch 'BEGIN {FS="\t"; OFS="\",\""} \
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
			if ($3 ~ /unpaired_1/){
				u="R1"}
			else if ($3 ~ /unpaired_2/){
				u="R2"
			}
			else{
				u="Paired"
			}
			print "\""$1,$2,$3,batch,fs[3],fs[1],i,u"\""
		}' $summaryfile >> $CSVfile
	done
done
