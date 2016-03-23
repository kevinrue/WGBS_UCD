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

str_perc="\(([0-9\.]*)%\)"
str_count="([[:digit:],]*)"
str_count_perc="$str_count $str_perc"
str_bp="([[:digit:],]*) bp"
str_bp_perc="$str_bp $str_perc"

str_read_in="Total reads processed:[[:space:]]*"
str_adapter="Reads with adapters:[[:space:]]*"
str_read_out="Reads written \(passing filters\):[[:space:]]*"
str_base_in="Total basepairs processed:[[:space:]]*"
str_base_trimmed="Quality-trimmed:[[:space:]]*"
str_base_out="Total written \(filtered\):[[:space:]]*"

str_nts="Bases preceding removed adapters:"

# Function to extract relevant info from file
extract_info(){
	if [ $# -lt 1 ]; then
		echo "Usage: $0 <trimming_report.txt>"
		exit 1
	fi
	input_read=$(grep "$str_read_in" $1 | 
		perl -pe "s/$str_read_in$str_count/\1/" | tr -d ',')
	adapter_info=$(grep "$str_adapter" $1 | \
		perl -pe "s/$str_adapter$str_count_perc/\1 \2/" | tr -d ',')
	output_read=$(grep -E "$str_read_out" $1 | \
		perl -pe "s/$str_read_out$str_count_perc/\1 \2/" | tr -d ',')
	input_base=$(grep "$str_base_in" $1 | \
		perl -pe "s/$str_base_in$str_bp/\1/" | tr -d ',')
	base_trimmed=$(grep "$str_base_trimmed" $1 | \
		perl -pe "s/$str_base_trimmed$str_bp_perc/\1 \2/" | tr -d ',')
	base_out=$(grep -E "$str_base_out" $1 | \
		perl -pe "s/$str_base_out$str_bp_perc/\1 \2/" | tr -d ',')
	A_adapter=$(grep -A 5 "$str_nts" $1 | grep 'A:' | \
		perl -pe "s/ *A: *([[:digit:]\.]*)%/\1/")
	C_adapter=$(grep -A 5 "$str_nts" $1 | grep 'C:' | \
		perl -pe "s/ *C: *([[:digit:]\.]*)%/\1/")
	G_adapter=$(grep -A 5 "$str_nts" $1 | grep 'G:' | \
		perl -pe "s/ *G: *([[:digit:]\.]*)%/\1/")
	T_adapter=$(grep -A 5 "$str_nts" $1 | grep 'T:' | \
		perl -pe "s/ *T: *([[:digit:]\.]*)%/\1/")
	N_adapter=$(grep -A 5 "$str_nts" $1 | grep 'none/other:' | \
		perl -pe "s/ *none\/other: *([[:digit:]\.]*)%/\1/")
	# return the result
	echo "$input_read $adapter_info $output_read $input_base $base_trimmed \
$base_out $A_adapter $C_adapter $G_adapter $T_adapter $N_adapter"
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
		#echo $sample
		# Extract information for the forward read
		info_forward=$(extract_info $report1)
		echo "$batch $sample $info_forward"
		exit
	done
	
done

