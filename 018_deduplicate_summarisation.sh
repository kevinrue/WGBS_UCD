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

# List unique folders that contain report files
folders=`find $rootdir -name '*deduplication_report.txt' -exec dirname {} \; | \
	sort | uniq`
echo -e "folders (next line):\n$folders"

str_count="([[:digit:]]*)"
str_perc="\(([0-9\.]*)%\)"
str_count_perc="$str_count $str_perc"
str_perc_total="\(([0-9\.]*)% of total\)"
str_count_perc_total="$str_count $str_perc_total"

str_align_in="Total number of alignments analysed in .*:[[:space:]]*"
str_duplicate_align="Total number duplicated alignments removed:[[:space:]]*"
str_duplicate_pos="Duplicated alignments were found at:[[:space:]]*[[:digit:]]* different position(s)"
str_duplicate_pos_perl="Duplicated alignments were found at:[[:space:]]*([[:digit:]]*) different position\(s\)"
str_align_out="Total count of deduplicated leftover sequences:[[:space:]]*"

# Function to extract relevant info from file
extract_info(){
	if [ $# -lt 1 ]; then
		echo "Usage: $0 <trimming_report.txt>"
		exit 1
	fi
	input_align=$(grep "$str_align_in" $1 | 
		perl -pe "s/$str_align_in$str_count/\1/")
	duplicate_align=$(grep "$str_duplicate_align" $1 | \
		perl -pe "s/$str_duplicate_align$str_count_perc/\1\",\"\2/")
	duplicate_pos=$(grep "$str_duplicate_pos" $1 | \
		perl -pe "s/$str_duplicate_pos_perl/\1/")
	output_align=$(grep "$str_align_out" $1 | \
		perl -pe "s/$str_align_out$str_count_perc_total/\1\",\"\2/")
	# return the result
	echo "$input_align\",\"$duplicate_align\",\"$duplicate_pos\",\"\
$output_align"
}

echo "\"Identifier\",\"Filename\",\"Sample\",\"Infection\",\
\"In\",\"Removed\",\"RemovedPct\",\"Positions\",\"Out\",\"OutPct\"" > $CSVfile

for folder in `echo $folders`
do
	echo "folder: $folder"
	# Identify all the forward reads in the folder
	reports=$(find $folder -name '*deduplication_report.txt')
	for report in $reports
	do
		echo "report: $report"
		filename=$(basename $report)
		# Extract sample information from the filename
		identifier=$(echo $filename | perl -pe 's/^(.*)_R1.*/\1/')
		sample=$(echo $filename | perl -pe 's/^([CM]{1}[[:digit:]]{1,2}).*/\1/')
		#echo "sample: $sample"
		infection=$(echo $filename | awk '{
			if ($0 ~ /^C/){i="Control"}
			else{i="M. bovis"}
			print i}' )
#		echo "infection: $infection"
		# Extract information for the forward read
		info=$(extract_info $report)
		echo "\"$identifier\",\"$filename\",\"$sample\",\
\"$infection\",\"$info\"" >> $CSVfile
	done
done

