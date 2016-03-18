#!/bin/bash

if [ $# -lt 1 ]; then
	echo "Usage: $0 <rootdir>"
	exit 1
fi

cwd=`pwd`
echo "cwd: $cwd"

rootdir=$1

log=$(basename $0 .sh)_$(date -I).txt
echo "log: $log"

# Overwrite the previous log file
> $cwd/log/$log

if [ ! -e $cwd/log ]; then
	mkdir -p $cwd/log
fi

#for folder in `echo $folders`
for folder in `ls $rootdir`
do
	echo "folder: $folder"
	batch=$(basename $folder)
	echo "batch: $batch"
	for summaryfile in `find $rootdir/$batch -name summary.txt`
	do
		echo "awk -v batch=$batch 'BEGIN {OFS=\"\t\"} {print \$0,batch}' $summaryfile >> $cwd/log/$log"
		awk -v batch=$batch 'BEGIN {OFS="\t"} {print $0,batch}' $summaryfile >> $cwd/log/$log
	done
done
