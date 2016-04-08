#!/bin/bash

echo "$0"
echo "$(date -I)"

if [ $# -lt 2 ]; then
	echo "Usage: $0 <rootdir> <threads> "
	exit 1
fi

cwd=`pwd`
echo "cwd: $cwd"

md5file=md5.txt

rootdir=$1
echo "rootdir: $rootdir"
threads=$2
echo "threads: $threads"

BAMfiles=$(find $rootdir -name "*.deduplicated.bam")
echo -e "BAMfiles (next line):\n$BAMfiles"

count=$(echo $BAMfiles | wc -w)

BAMfiles=$(echo $BAMfiles | xargs)
BAMsorted=$(echo $BAMfiles | perl -pe 's/.bam/_sorted.bam/g')
BAMrg=$(echo $BAMfiles | perl -pe 's/.bam/_RG.bam/g')

cmd_sort="samtools sort -o"
cmd_rg="samtools addreplacerg -o"
cmd_index="samtools index"
if [ $count -gt 0 ];
then

	echo "parallel -j $threads --xapply $cmd_sort ::: $BAMsorted ::: $BAMfiles"
#	time(
#		parallel -j $threads --xapply $cmd_sort ::: $BAMsorted ::: $BAMfiles
#	)
	
	echo "parallel -j $threads --xapply $cmd_rg ::: $BAMrg ::: $BAMsorted"
	time(
		parallel -j $threads --xapply $cmd_rg ::: $BAMrg ::: $BAMsorted
	)
	
	echo "parallel -j $threads --xapply $cmd_index ::: $BAMrg"
	time(
		parallel -j $threads --xapply $cmd_index ::: $BAMrg
	)
	
	echo "rm $BAMsorted"
	rm $BAMsorted
	
fi

echo "Completed."
