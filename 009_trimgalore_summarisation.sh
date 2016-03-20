#!/bin/bash

if [ $# -lt 1 ]; then
	echo "Usage: $0 <rootdir>"
	exit 1
fi

cwd=`pwd`
echo "cwd: $cwd"

rootdir=$1
echo "rootdir: $rootdir"

# List unique folders that contain report files
folders=`find trim_galore -name '*trimming_report.txt' -exec dirname {} \; | sort | uniq`
echo -e "folders (next line):\n$folders"

str_input="Total reads processed:[[:space:]]*"
#echo $str_input
str_adapter="Reads with adapters:[[:space:]]*"

# Function to extract relevant info from file
extract_info(){
	if [ $# -lt 1 ]; then
		echo "Usage: $0 <trimming_report.txt>"
		exit 1
	fi
	input_count=$(grep "$str_input" $1 | \
		sed -e "s/$str_input//" | tr -d ',')
	adapter_info=($(grep "$str_adapter" $1 | \
		sed -e "s/$str_adapter//" | tr -d ',' ))
	adapter_count=${adapter_info[0]}
	adapter_perc=${adapter_info[1]//[()%]/}
	echo "$input_count $adapter_count $adapter_perc"
}

for folder in `echo $folders`
do
	echo $folder
	# Identify the batch to annotate the output metrics
	batch=$(basename $folder)
	# Identify all the forward reads in the folder
	report1s=$(find $folder -name '*R1*trimming_report.txt')
	for report1 in $report1s
	do
		echo $report1
		filename=$(basename $report1)
		file_info=(${filename//_/ })
		# Identify the sample to annotate the output metrics
		sample=${file_info[0]}
		echo $sample
		# Extract information for the forward read
		info_forward=$(extract_info $report1)
		echo $info_forward
		exit
	done
	
done

