#!/bin/bash

if [ $# -lt 1 ]; then
	echo "Usage: $0 <rootdir>"
	exit 1
fi

cwd=`pwd`
echo "cwd: $cwd"

md5file=md5.txt

rootdir=$1

folders=`find $rootdir -name "$md5file" -exec dirname {} \;`
echo -e "folders (next line):\n$folders"

log=$(basename $0 .sh)_$(date -I).txt
echo "log: $log"

if [ ! -e $cwd/log ]; then
	mkdir -p $cwd/log
fi

for folder in `echo $folders`
do
	cd $rootdir/$folder
	echo "folder: `pwd`"
	echo "Processing $(wc -l $md5file) files..."
	echo "md5sum -c $md5file"
	md5sum -c ./$md5file >> $cwd/log/$log
	# Restore original state for the loop (critical if rootdir is relative)
	cd $cwd
done

