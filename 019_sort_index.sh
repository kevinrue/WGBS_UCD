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
BAMsortedRG=$(echo $BAMfiles | perl -pe 's/.bam/_picard.bam/g')

RGIDs=$(basename $BAMfiles | perl -pe 's/([CM][[:digit:]]{1,2}).*/\1/g')
echo -e "RGSMs (next line):\n$RGSMs"
RGPUs=$(basename $BAMfiles | perl -pe 's/.*([ATGC]{6}).*/\1/g')
echo -e "RGPUs (next line):\n$RGPUs"

count=$(echo $BAMfiles | wc -w)
echo "count $count"

BAMfiles=$(echo $BAMfiles | xargs)
exit

cmd_sort="samtools sort -o"
cmd_picard="java -jar picard.jar AddOrReplaceReadGroups \
    I={1} \
    O={2} \
    SORT_ORDER=coordinate \
    RGID={3} \
    RGLB={3} \
    RGPL=illumina \
    RGPU={4} \
    RGSM={3}"

if [ $count -gt 0 ];
then
$
	echo "parallel -j $threads --xapply $cmd_picard ::: $BAMfiles ::: $BAMsortedRG"
#	time(
#		parallel -j $threads --xapply $cmd_sort ::: $BAMsorted ::: $BAMfiles
#	)

	echo "parallel -j $threads --xapply $cmd_index ::: $BAMsortedRG"
#	time(
#		parallel -j $threads --xapply $cmd_index ::: $BAMrg
#	)

fi

echo "Completed."
