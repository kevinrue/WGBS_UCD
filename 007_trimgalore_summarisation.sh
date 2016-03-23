#!/bin/bash

echo "$0"

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
folders=`find $rootdir -name '*trimming_report.txt' -exec dirname {} \; | \
	sort | uniq`
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

str_dec_sep="([[:digit:]]),([[:digit:]])"

# Function to extract relevant info from file
extract_info(){
	if [ $# -lt 1 ]; then
		echo "Usage: $0 <trimming_report.txt>"
		exit 1
	fi
	input_read=$(grep "$str_read_in" $1 | 
		perl -pe "s/$str_read_in$str_count/\1/" | tr -d ',')
	adapter_info=$(grep "$str_adapter" $1 | \
		perl -pe "s/$str_adapter$str_count_perc/\1\",\"\2/" | \
		perl -pe "s/$str_dec_sep/\1\2/g")
	output_read=$(grep -E "$str_read_out" $1 | \
		perl -pe "s/$str_read_out$str_count_perc/\1\",\"\2/" | \
		perl -pe "s/$str_dec_sep/\1\2/g")
	input_base=$(grep "$str_base_in" $1 | \
		perl -pe "s/$str_base_in$str_bp/\1/" | \
		perl -pe "s/$str_dec_sep/\1\2/g")
	base_trimmed=$(grep "$str_base_trimmed" $1 | \
		perl -pe "s/$str_base_trimmed$str_bp_perc/\1\",\"\2/" | \
		perl -pe "s/$str_dec_sep/\1\2/g")
	base_out=$(grep -E "$str_base_out" $1 | \
		perl -pe "s/$str_base_out$str_bp_perc/\1\",\"\2/" | \
		perl -pe "s/$str_dec_sep/\1\2/g")
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
	echo "$input_read\",\"$adapter_info\",\"$output_read\",\"$input_base\",\"\
$base_trimmed\",\"$base_out\",\"$A_adapter\",\"$C_adapter\",\"$G_adapter\",\"\
$T_adapter\",\"$N_adapter"
}

NA_info(){
	echo "NA\",\"NA\",\"NA\",\"NA\",\"NA\",\"NA\",\"NA\",\"NA\",\"NA\",\"\
NA\",\"NA\",\"NA\",\"NA\",\"NA\",\"NA" # 15 NAs
}

echo "\"Batch\",\"Sample\",\"Read1\",\"Adapter1\",\"AdapterPct1\",\
\"ReadOut1\",\"ReadOutPct1\",\"Base1\",\"QualTrim1\",\"QualTrimPct1\",\
\"BaseOut1\",\"BaseOutPct1\",\"A1\",\"C1\",\"G1\",\"T1\",\"N1\",\"Read2\",\
\"Adapter2\",\"AdapterPct2\",\"ReadOut2\",\"ReadOutPct2\",\"Base2\",\
\"QualTrim2\",\"QualTrimPct2\",\"BaseOut2\",\"BaseOutPct2\",\"A2\",\"C2\",\
\"G2\",\"T2\",\"N2\",\"Pass\",\"Pass1\",\"Pass2\"" > $CSVfile

for folder in `echo $folders`
do
	echo "folder: $folder"
	# Identify the batch to annotate the output metrics
	batch=$(basename $folder)
	# Identify all the forward reads in the folder
	report1s=$(find $folder -name '*R1*trimming_report.txt')
	for report1 in $report1s
	do
		echo "report1: $report1"
		filename=$(basename $report1)
		file_info=(${filename//_/ })
		# Identify the sample to annotate the output metrics
		sample=${file_info[0]}
		#echo $sample
		# Extract information for the forward read
		info_forward=$(extract_info $report1)
		# If there is a second mate extract the same info
		report2=$(echo $report1 | sed -e 's/_R1_/_R2_/')
		if [ -e $report2 ]; then
			echo "report2: $report2"
			info_reverse=$(extract_info $report2)
			# Guess the name of the validated forward file
			val1=$(echo $report1 | \
				sed -e 's/\.fastq\.gz_trimming_report\.txt/_val_1.fq.gz/')
#			val2=$(echo $val1 | perl -pe 's/R1_(.*)_val_1/R2_\1_val_2/')
			echo "val1: $val1"
#			echo "val2: $val2"
			# Count the number of validated read pairs (based on forward)
			pass_all=$(($(zcat $val1 | wc -l) / 4 ))
			echo "pass_all: $pass_all"
#			pass_2=$(($(zcat $val2 | wc -l) / 4 ))
#			echo "pass_2: $pass_2"
			unpaired1=$(echo $report1 | \
				sed -e 's/\.fastq\.gz_trimming_report\.txt/_unpaired_1.fq.gz/')
			pass_1=$(($(zcat $unpaired1 | wc -l) / 4 ))
			echo "pass_1: $pass_1"
			unpaired2=$(echo $report2 | \
				sed -e 's/\.fastq\.gz_trimming_report\.txt/_unpaired_2.fq.gz/')
			pass_2=$(($(zcat $unpaired2 | wc -l) / 4 ))
			echo "pass_2: $pass_2"
		# if there is no other mate insert NAs
		else
			info_reverse=$(NA_info)
			trimmed=$(echo $report1 | \
				sed -e 's/\.fastq\.gz_trimming_report\.txt/_trimmed.fq.gz/')
			pass_all=$(($(zcat $trimmed | wc -l) / 4 ))
			pass_1='NA'
			pass_2='NA'
		fi
		echo "\"$batch\",\"$sample\",\"$info_forward\",\"$info_reverse\",\"\
$pass_all\",\"$pass_1\",\"$pass_2\"" >> $CSVfile
	done
done

