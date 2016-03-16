#!/bin/bash

if [ $# -lt 1 ]; then
	echo "Usage: $0 <rootdir>"
	exit 1
fi

cwd=`pwd`

md5file=md5.txt

rootdir=$1

folders=`find $rootdir -name "$md5file" -exec dirname {} \;`
echo -e "folders (next line):\n$folders"

log=$(basename $0 .sh)_$(date -I)
echo "log: $log"

if [ ! -e $rootdir/log ]; then
	mkdir -p $rootdir/log
fi

for folder in `echo $folders`
do
	cd $rootdir/$folder
	echo "folder: `pwd`"
	echo "Processing $(ls ./*fastq.gz | wc -l) files..."
	echo "md5sum -c $md5file"
	md5sum -c ./$md5file >> $cwd/log/$log
done

