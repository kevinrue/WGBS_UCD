#!/bin/bash

echo "$0"
echo "$(date -I)"

if [ $# -lt 3 ]; then
	echo "Usage: $0 <rootdir> <outdir> <threads>"
	exit 1
fi

cwd=$(pwd)
echo "cwd: $cwd"

rootdir=$1
echo "rootdir: $rootdir"
outdir=$2
echo "outdir: $outdir"
threads=$3
echo "threads: $threads"

if [ ! -e $outdir ]; then
	mkdir -pv $outdir
fi

alignmentReports=$(find $rootdir -name "*_PE_report.txt" | sort | xargs)
echo "alignmentReports: $(echo $alignmentReports | wc -w)"
#echo $alignmentReports

#dedupReports=$(echo $alignmentReports | sed 's/_PE_report.txt/_pe.deduplication_report.txt/g')
dedupReports=$(find $rootdir -name "*_pe.deduplication_report.txt" | sort | xargs)
#echo $dedupReports

#splittingReports=$(echo $alignmentReports | sed 's/_PE_report.txt/_pe.deduplicated_splitting_report.txt/g')
splittingReports=$(find $rootdir -name "*_pe.deduplicated_splitting_report.txt" | sort | xargs)
#echo $splittingReports

#mBiasReports=$(echo $alignmentReports | sed 's/_PE_report.txt/_pe.deduplicated.M-bias.txt/g')
mBiasReports=$(find $rootdir -name "*_pe.deduplicated.M-bias.txt" | sort | xargs)
#echo $mBiasReports

echo "parallel -j $threads --xargs bismark2report --dir $outdir \
	--alignment_report {1} --dedup_report {2} --splitting_report {3} \
	--mbias_report {4} ::: $alignmentReports ::: $dedupReports ::: \
	$splittingReports ::: $mBiasReports"
	
parallel -j $threads --xargs bismark2report --dir $outdir \
	--alignment_report {1} --dedup_report {2} --splitting_report {3} \
	--mbias_report {4} ::: $alignmentReports ::: $dedupReports ::: \
	$splittingReports ::: $mBiasReports

