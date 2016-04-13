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

echo "\"Value\",\"QC\",\"File\",\"Sample\",\"Read\",\"Treatment\",\
\"Infection\"" > $CSVfile

for folder in `ls $rootdir`
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
