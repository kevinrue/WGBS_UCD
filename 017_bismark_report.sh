#!/bin/bash

echo "$0"
echo "$(date -I)"

if [ $# -lt 3 ]; then
	echo "Usage: $0 <rootdir(s)> <outdir> <threads>"
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
echo "dedupReports: $(echo $dedupReports | wc -w)"

#splittingReports=$(echo $alignmentReports | sed 's/_PE_report.txt/_pe.deduplicated_splitting_report.txt/g')
splittingReports=$(find $rootdir -name "*_pe.deduplicated_splitting_report.txt" | sort | xargs)
echo "splittingReports: $(echo $splittingReports | wc -w)"

#mBiasReports=$(echo $alignmentReports | sed 's/_PE_report.txt/_pe.deduplicated.M-bias.txt/g')
mBiasReports=$(find $rootdir -name "*_pe.deduplicated.M-bias.txt" | sort | xargs)
echo "mBiasReports: $(echo $mBiasReports | wc -w)"

#mBiasReports=$(TODO)
nucleotideReports=$(find $rootdir -name "*_pe.nucleotide_stats.txt" | sort | xargs)
echo "nucleotideReports: $(echo $nucleotideReports | wc -w)"

echo "parallel -j $threads --xapply bismark2report --dir $outdir \
	--alignment_report {1} --dedup_report {2} --splitting_report {3} \
	--mbias_report {4} --nucleotide_report {5} ::: $alignmentReports ::: $dedupReports ::: \
	$splittingReports ::: $mBiasReports ::: $nucleotideReports"
	
parallel -j $threads --xapply bismark2report --dir $outdir \
	--alignment_report {1} --dedup_report {2} --splitting_report {3} \
	--mbias_report {4} --nucleotide_report {5} ::: $alignmentReports ::: $dedupReports ::: \
	$splittingReports ::: $mBiasReports ::: $nucleotideReports

