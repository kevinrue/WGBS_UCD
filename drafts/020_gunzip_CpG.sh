#!/bin/bash

echo "$0"
echo "$(date -I)"

if [ $# -lt 2 ]; then
	echo "Usage: $0 <rootdir> <threads>"
	exit 1
fi

cwd=`pwd`
echo "cwd: $cwd"

rootdir=$1
echo "rootdir: $rootdir"
threads=$2
echo "threads: $threads"


CpGfiles=$(find $rootdir -name "CpG*txt")
echo -e "CpGfiles (next line):\n$CpGfiles"

count=$(echo $CpGfiles | wc -w)
echo "count $count"

CpGfiles=$(echo $CpGfiles | xargs)

uncompressed=$(echo $CpGfiles | perl -pe 's/.gz$//g')

cmd="gunzip -c {1} > {2}"

if [ $count -gt 0 ];
then

	echo "parallel -j $threads --xapply $cmd ::: $CpGfiles ::: $uncompressed"
#	time(
#		parallel -j $threads --xapply $cmd ::: $CpGfiles ::: $uncompressed
#	)

fi

echo "Completed."
