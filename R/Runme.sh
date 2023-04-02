#!/bin/bash
## 1) The arugument 1 is the ProjectName you wish to give for your project
## 2) The arugment 2 is the Directory Location where you want to create the project
ProjectName=$1
DirectoryLocationName=$2
DirandProjectName="$DirectoryLocationName""/""$ProjectName"

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
