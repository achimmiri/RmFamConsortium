#!/bin/bash
ProjectName=$1
DirectoryName=$2
#####echo "$DirectoryName"
DirandProjectName="$DirectoryName""/""$ProjectName"

if [[ ! -e $DirandProjectName ]]; then
	mkdir $DirandProjectName
	mkdir $DirandProjectName"/meta data"
	mkdir $DirandProjectName"/converted to msp"
        mkdir $DirandProjectName"/raw data"
	mkdir $DirandProjectName"/raw data/exported as raw msp"
	mkdir $DirandProjectName"/raw data/raw"
	mkdir $DirandProjectName"/RMassBank"
	mkdir $DirandProjectName"/Error-Report"
else
	echo "The Directory name already exists"
fi
