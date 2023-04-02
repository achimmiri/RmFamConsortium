#!/bin/bash
## 1) The arugument 1 is the ProjectName you wish to give for your project
## 2) The arugment 2 is the Directory Location where you want to create the project


create_directory(){


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

}


handle_dfiles(){


	ProjectName=$1
        DirectoryLocationName=$2
        DirandProjectName="$DirectoryLocationName""/""$ProjectName"

	dialog --yesno "Use MS-Dail and convert .d files to .msp files and move .msp files to /raw data/exported as raw msp" 10 40

        if [[ $? -eq 1 ]]
	then
		dialog --infobox "rawfiles and metadata files are required fur processig the data" 10 60
		exit -1
	fi

	if [ -n "$(find $DirandProjectName"/raw data/exported as raw msp" -prune -empty 2>/dev/null)" ]
	then
  		dialog --infobox "Empty directory and please check why the files does not exist" 10 60
		exit -1
	fi

  	if [ -n "$(find $DirandProjectName"/meta data" -prune -empty 2>/dev/null)" ]
	then
		dialog --infobox "There us no metadata file there needs to be both raw data and metadata files for processing the data" 20 60
		exit -1
		
	fi
        
	FFind=$(find $DirandProjectName"/meta data" -name "*xlsx")
        FFindCount=$(echo "$FFind"|wc -l)	
        
	if [[ "$FFindCount" -ne 1 ]]
	then
		dialog --infobox "There is more than one metadata file please check and correct it" 20 60
		exit -1
	fi
     
        /usr/bin/env Rscript ValidateMetaData.R "$FFind"  
	
	if [[ $? -eq 1 ]]
	then
		/usr/bin/env Rscript MfAM-Contributions.AllParameetrs.Avilable.R "$FFind" 25 20 0.06
	fi
        ##echo $FFind
        ##/usr/bin/env RScript ValidateMetaData.R 	
	       ###	
	
        
	

}



###create_directory "test-project" "/home/achimmir/temp"

dialog --menu "Select the raw data type" 10 40 20 1 ".d" 2 ".raw" 3 ".mzML" 4 ".mzXML" 2>/tmp/menu
item=$(cat /tmp/menu)

case $item in
	"1") handle_dfiles $1 $2 ;;
	"2") dialog --yesno "Use MS-Dial and convert .raw files to .msp files and move .msp files to /raw data/exported as raw msp" 5 40 ;;
	"3") dialog --yesno "Use MS-dial and convert .mzML to .msp files and move .msp files to  /raw data/exprted as raw msp" 5 40;;
	"4")dialog --yesno "First convert mzXML to mzML using msconvert GUI and move .msp files to /raw data/exprted as raw msp" 5 40;;
esac




