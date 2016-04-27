#!/bin/bash

echo "$0"
echo "$(date -I)"

if [ $# -lt 3 ]; then
	echo "Usage: $0 [OPTIONS] <rootdir> <target_file> <outdir>"
	echo "
Options:
	-s	Skip FastQC step. Only summarise content of outdir in CSV
		file
	
	-t	Number of threads. Only applicable to FastQC step.
"
	exit 1
fi

cwd=`pwd`
echo "cwd: $cwd"

while getopts "st:c:" opt; do
	echo "@: $@"
	echo "opt: $opt"
	case $opt in
		s)
			fastqc=false
			if [ -n "$threads" ]; then
				echo "'-t' not applicable: FastQC step skipped."
			fi
#			echo "OPTIND: $OPTIND"
			shift $((OPTIND-1)); OPTIND=1
			;;
		c)
			CSVfile=$OPTARG
			CSVfolder=$(dirname $CSVfile)
#			echo "OPTIND: $OPTIND"
			shift $((OPTIND-1)); OPTIND=1
			;;
		t)
			threads=$OPTARG
			if [ -n "$fastqc" ]; then
				echo "'-t' not applicable: FastQC step skipped."
			fi
#			echo "OPTIND: $OPTIND"
			shift $((OPTIND-1)); OPTIND=1
			;;
		\?)
			echo "Invalid option: -$OPTARG" >&2
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument." >&2
			exit 1
			;;
	esac
done

echo "rootdir: $rootdir"
echo "target_file: $target_file"
echo "outdir: $outdir"
echo "fastqc: $fastqc"
echo "threads: $threads"
echo "CSVfolder: $CSVfolder"
echo "CSVfile: $CSVfile"

rootdir=$1; shift
target_file=$1; shift
outdir=$1; shift

# Defaults
if [ -z $fastqc ]; then
	fastqc=true
fi
if [ -z $threads ]; then
	threads=1
fi
if [ -z $CSVfolder ]; then
	CSVfolder="log"
fi
if [ -z $CSVfile ]; then
	CSVfile=$CSVfolder/"$(basename $0 | sed -e "s/\.sh/_$(date -I).csv/")"
fi

#echo "rootdir: $rootdir"
#echo "target_file: $target_file"
#echo "outdir: $outdir"
#echo "fastqc: $fastqc"
#echo "threads: $threads"
#echo "CSVfolder: $CSVfolder"
#echo "CSVfile: $CSVfile"

folders=`find $rootdir -name "$target_file" -exec dirname {} \; | sort | uniq`
echo -e "folders (next line):\n$folders"

if [ ! -e $outdir ]; then
	mkdir -pv $outdir
fi
if [ ! -e $CSVfolder ]; then
	mkdir -pv $CSVfolder
fi

if $fastqc; then
	echo "Doing fastqc"
	exit
	for folder in `echo $folders`
	do
		echo "folder: $folder"
		# Need to separate batches of sequencing
		# as some fastq files for the same sample have the same name
		batch=$(basename $folder)
		if [ ! -e $outdir/$batch ]; then
			mkdir -pv $outdir/$batch
		fi
		fastqs=$(find $folder -name "$target_file" | xargs)
		echo "Processing $(echo $fastqs | wc -w ) files..."
		echo "fastqc --quiet --nogroup --extract --threads $threads --outdir $outdir/$batch $fastqs"
		fastqc --quiet --nogroup --extract --threads $threads --outdir $outdir/$batch $fastqs
	done
else
	echo "Skipping FastQC. Jumping to summarisation into CSV file."
fi

echo "Collating FastQC reports.."

echo "\"Value\",\"QC\",\"File\",\"Sample\",\"Treatment\",\
\"Infection\"" > $CSVfile

for folder in `ls $outdir`
do
	echo "folder: $outdir/$folder"
	for summaryfile in `find $outdir/$folder -name summary.txt`
	do
		awk 'BEGIN {FS="\t"; OFS="\",\""} \
		\
		{
			nf=split($3,fs,"_")
			if ($3 ~ /NOT_BS/){
				t="NOT BS"}
			else{
				t="BS"
			}
			if (fs[1] ~ /C/){
				i="Control"}
			else{
				i="M. bovis"
			}
			print "\""$1,$2,$3,fs[1],t,i"\""
		}' $summaryfile >> $CSVfile
	done
done

echo "Completed."
