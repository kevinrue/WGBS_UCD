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
folders=`find $rootdir -name '*_report.txt' -exec dirname {} \; | \
	sort | uniq`
echo -e "folders (next line):\n$folders"

str_perc="([0-9\.]*)%.*"
str_count="([[:digit:]]*)"
#Number of paired-end alignments with a unique best hit: 7407334
#Number of alignments with a unique best hit from the different alignments:
# Adapt below to seemlessly process single- and paired-end files
str_read_in="Sequence[ pair]*s analysed in total:[[:space:]]*"
str_unique="Number of[[:space:][:alpha:]-]*alignments with a unique best hit[[:space:][:alpha:]]*:[[:space:]]*"
str_efficiency="Mapping efficiency:[[:space:]]*"
str_unaligned="Sequence[ pair]*s with no alignments under any condition:[[:space:]]*"
str_multimap="Sequence[ pair]*s did not map uniquely:[[:space:]]*"
str_seq_issue="Sequence[ pair]*s which were discarded because genomic sequence could not be extracted:[[:space:]]*"
str_OT_start="CT[\/GA]*\/CT:[[:space:]]*"
str_OT_end="[[:space:]]*\(\(converted\) top strand\)"
str_CTOT_start="GA[\/CT]*\/CT:[[:space:]]*"
str_CTOT_end="[[:space:]]*\(complementary to \(converted\) top strand\)"
str_CTOB_start="GA[\/CT]*\/GA:[[:space:]]*"
str_CTOB_end="[[:space:]]*\(complementary to \(converted\) bottom strand\)"
str_OB_start="CT[\/GA]*\/GA:[[:space:]]*"
str_OB_end="[[:space:]]*\(\(converted\) bottom strand\)"
str_m_CpG="Total methylated C's in CpG context:[[:space:]]*"
str_m_CHG="Total methylated C's in CHG context:[[:space:]]*"
str_m_CHH="Total methylated C's in CHH context:[[:space:]]*"
str_m_unk="Total methylated C's in Unknown context:[[:space:]]*"
str_u_CpG="Total unmethylated C's in CpG context:[[:space:]]*"
str_u_CHG="Total unmethylated C's in CpG context:[[:space:]]*"
str_u_CHH="Total unmethylated C's in CpG context:[[:space:]]*"
str_u_unk="Total unmethylated C's in Unknown context:[[:space:]]*"
str_pm_CpG="C methylated in CpG context:[[:space:]]*"
str_pm_CHG="C methylated in CHG context:[[:space:]]*"
str_pm_CHH="C methylated in CHH context:[[:space:]]*"
str_pm_unk="C methylated in [uU]nknown context \(CN or CHN\):[[:space:]]*"

echo "$str_pm_unk$str_perc"

# Function to extract relevant info from file
extract_info(){
	if [ $# -lt 1 ]; then
		echo "Usage: $0 <report.txt>"
		exit 1
	fi
	input_read=$(grep "$str_read_in" $1 | 
		perl -pe "s/$str_read_in$str_count/\1/")
	unique_read=$(grep "$str_unique" $1 | 
		perl -pe "s/$str_unique$str_count/\1/")
	efficiency=$(grep "$str_efficiency" $1 | 
		perl -pe "s/$str_efficiency$str_perc/\1/")
	unaligned=$(grep "$str_unaligned" $1 | 
		perl -pe "s/$str_unaligned$str_count/\1/")
	multimap=$(grep "$str_multimap" $1 | 
		perl -pe "s/$str_multimap$str_count/\1/")
	seq_issue=$(grep "$str_seq_issue" $1 | 
		perl -pe "s/$str_seq_issue$str_count/\1/")
	OT=$(grep -E "$str_OT_end" $1 | 
		perl -pe "s/$str_OT_start$str_count$str_OT_end/\1/")
	CTOT=$(grep -E "$str_CTOT_end" $1 | 
		perl -pe "s/$str_CTOT_start$str_count$str_CTOT_end/\1/")
	CTOB=$(grep -E "$str_CTOB_end" $1 | 
		perl -pe "s/$str_CTOB_start$str_count$str_CTOB_end/\1/")
	OB=$(grep -E "$str_OB_end" $1 | 
		perl -pe "s/$str_OB_start$str_count$str_OB_end/\1/")
	m_CpG=$(grep "$str_m_CpG" $1 | 
		perl -pe "s/$str_m_CpG$str_count/\1/")
	m_CHG=$(grep "$str_m_CHG" $1 | 
		perl -pe "s/$str_m_CHG$str_count/\1/")
	m_CHH=$(grep "$str_m_CHH" $1 | 
		perl -pe "s/$str_m_CHH$str_count/\1/")
	m_unk=$(grep "$str_m_unk" $1 | 
		perl -pe "s/$str_m_unk$str_count/\1/")
	u_CpG=$(grep "$str_u_CpG" $1 | 
		perl -pe "s/$str_u_CpG$str_count/\1/")
	u_CHG=$(grep "$str_u_CHG" $1 | 
		perl -pe "s/$str_u_CHG$str_count/\1/")
	u_CHH=$(grep "$str_u_CHH" $1 | 
		perl -pe "s/$str_u_CHH$str_count/\1/")
	u_unk=$(grep "$str_u_unk" $1 | 
		perl -pe "s/$str_u_unk$str_count/\1/")
	pm_CpG=$(grep "$str_pm_CpG" $1 | 
		perl -pe "s/$str_pm_CpG$str_perc/\1/")
	pm_CHG=$(grep "$str_pm_CHG" $1 | 
		perl -pe "s/$str_pm_CHG$str_perc/\1/")
	pm_CHH=$(grep "$str_pm_CHH" $1 | 
		perl -pe "s/$str_pm_CHH$str_perc/\1/")
	pm_unk=$(grep -E "$str_pm_unk" $1 | 
		perl -pe "s/$str_pm_unk$str_perc/\1/")
	# return the result
	echo "$input_read\",\"$unique_read\",\"$efficiency\",\"$unaligned\",\"\
$multimap\",\"$seq_issue\",\"$OT\",\"$CTOT\",\"$CTOB\",\"$OB\",\"$m_CpG\",\"\
$m_CHG\",\"$m_CHH\",\"$m_unk\",\"$u_CpG\",\"$u_CHG\",\"$u_CHH\",\"$u_unk\",\"\
$pm_CpG\",\"$pm_CHG\",\"$pm_CHH\",\"$pm_unk"
}

echo "\"Batch\",\"Identifier\",\"Filename\",\"Sample\",\"Treatment\",\
\"Infection\",\"Lane\",\"ReadIn\",\"Unique\",\"Efficiency\",\"Unaligned\",\
\"Multimap\",\"SeqIssue\",\"OT\",\"CTOT\",\"CTOB\",\"OB\",\"MethCpG\",\
\"MethCHG\",\"MetCHH\",\"MethUnk\",\"UnmethCpG\",\"UnmethCHG\",\"UnmetCHH\",\
\"UnmethUnk\",\"MethPctCpG\",\"MethPctCHG\",\"MethPctCHH\",\"MethPctUnk\"\
"> $CSVfile

for folder in `echo $folders`
do
	echo "folder: $folder"
	# Identify the batch to annotate the output metrics
	batch=$(basename $folder)
	echo "batch: $batch"
	# Identify all the forward reads in the folder
	report1s=$(find $folder -name '*_report.txt')
	echo -e "folders (next line):\n$report1s"
	for report1 in $report1s
	do
		echo "report1: $report1"
		filename=$(basename $report1)
		# Extract sample information from the filename
		identifier=$(echo $filename | perl -pe 's/^(.*)_R1.*/\1/')
		sample=$(echo $filename | perl -pe 's/^([CM]{1}[[:digit:]]{1,2}).*/\1/')
		#echo "sample: $sample"
		treatment=$(echo $filename | awk '{
			if ($0 ~ /_NOT_BS_/){t="NOT BS"}
			else{t="BS"}
			print t}' )
		echo "treatment: $treatment"
		infection=$(echo $filename | awk '{
			if ($0 ~ /^C/){i="Control"}
			else{i="M. bovis"}
			print i}' )
#		echo "infection: $infection"
		lane=$(echo $filename | perl -pe 's/.*(L[[:digit:]]{3}).*/\1/')
#		echo "lane: $lane"
		# Extract information for the forward read
		info=$(extract_info $report1)
		
		echo "\"$batch\",\"$identifier\",\"$filename\",\"$sample\",\"\
$treatment\",\"$infection\",\"$lane\",\"$info\"" >> $CSVfile
	done
done
