#!/bin/bash
## 1) The arugument 1 is the ProjectName you wish to give for your project
## 2) The arugment 2 is the Directory Location where you want to create the project
## 3) The arugment 3 is the mz ppm
## 4) The argument 4 is the RT ppm
## 5) The argument 5 is the centroiding parameter
## 6) The argument 6 is the API key
## 7) The argument 7 is the Adduct File location 
## 8) The argument 8 is the Database File location


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
	mzPPM=$3
	rtPPM=$4
	rtPPM1=$(echo "scale=1 ; $rtPPM / 100" | bc | sed 's/^\./0./')
	##rtPPM=$((4/100))
	centrD=$5
	## Adding the three new arguments for apikey, adduct ,and database files
	APIK=$6
	AddF=$7
	DBF=$8

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
		####/usr/bin/env Rscript MfAM-Contributions.AllParameetrs.Avilable.R "$FFind" "$mzPPM" "$rtPPM" "$centrD"
		/usr/bin/env Rscript MfAM-Contributions.AllParameetrs.Avilable.R "$FFind" "$mzPPM" "$rtPPM" "$centrD" "$APIK" "$AddF" "$DBF"
		

		CheckFoldernanme=$DirandProjectName/"converted to msp"/mz."$mzPPM"ppm."$rtPPM1"RT
                


		if [ -z "$(ls -A "$CheckFoldernanme")" ]; then
			dialog --infobox "Something went wrong in the script execution check your data" 20 60
   			exit -1
		else
   			echo "$CheckFoldernanme"
			/usr/bin/env Rscript CreateCombinedmsp.R "$CheckFoldernanme"

			###touch "$DirandProjectName"/.combinedmspexists

			
			Fcommsp=$(find "$CheckFoldernanme" -name "*combined.msp" -print0)

			if [ -f "$Fcommsp" ] && [ ! -f "$DirandProjectName"/.combinedmspexists ]; then echo Hello; fi

			##if [ -f "$Fcommsp" ]; then
			if [ -f "$Fcommsp" ] && [ ! -f "$DirandProjectName"/.combinedmspexists ]; then
				
				echo "it is passing the if loop"

				Checkrow=$(/usr/bin/env Rscript CheckExcel.R "$FFind")  
				CheckName=$(grep "NAME:" "$Fcommsp"|wc -l)
				CheckOntolgy=$(grep "Ontology:" "$Fcommsp"|wc -l)
				CheckSmiles=$(grep "SMILES:" "$Fcommsp"|wc -l)
				CheckInchikey=$(grep "INCHIKEY:" "$Fcommsp"|wc -l)
			       
				prefix="[1]"
				Checkrow1=${Checkrow#"$prefix"}
			        
			       
			        if [[ "$Checkrow1" -eq "$CheckName" ]] && [[ "$Checkrow1" -eq "$CheckOntolgy" ]] && [[ "$Checkrow1" -eq "$CheckSmiles" ]] && [[ "$Checkrow1" -eq "$CheckInchikey" ]];then
					
					Fmsp=$(find "$CheckFoldernanme" -name "*.msp" -print0)
					FL=1
					echo "copying the files ...check if it works"
					Fmsp1=${#Fmsp[@]}

					for files in "$CheckFoldernanme"/**msp
					do


						rsync -av "$files" "$DirandProjectName"/RMassBank
					done

					for files in "$FFind"/*xlsx
					do
						 rsync -av "$files" "$DirandProjectName"/RMassBank
					done
					
					echo "files are copied so creating the combinedmspexists"

					touch "$DirandProjectName"/.combinedmspexists

				else
					dialog --infobox "something went wrong in generating the combined msp file check the previous step" 20 60
					
				fi



			else
				dialog --infobox "Something went wrong in combining the msp files" 20 60
				exit -1
			fi
		fi
		###/usr/bin/env Rscript CreateCombinedmsp.R 
	fi


        ##echo $FFind
        ##/usr/bin/env RScript ValidateMetaData.R 	
	       ###	
	
        
	

}



###create_directory "test-project" "/home/achimmir/temp"

dialog --menu "Select the raw data type" 10 40 20 1 ".d" 2 ".raw" 3 ".mzML" 4 ".mzXML" 2>/tmp/menu
item=$(cat /tmp/menu)

if [ -z "$3" ]
then

	mzPPM=25
else
	mzPPM=$3
fi


if [ -z "$4" ]
then
	##rtPPM=20
	rtPPM=0.2

else
	rtPPM=$4
	rtPPM1=$(echo "scale=1 ; $rtPPM / 100" | bc | sed 's/^\./0./')
	
fi

if [ -z "$5" ]
then
	centrD=0.06

else
	centrD=$5
fi

if [ -z "$6" ]
then
	dialog --infobox "The APIKEY is required for processing the chemical transition"
else
	APIK=$6
fi

if [ -z "$7" ]
then
        dialog --infobox "define the location of adduct file"
else
        AddF=$7
fi


if [ -z "$8" ]
then
        dialog --infobox "define the location of Database file"
else
        DBF=$8
fi




case $item in
	####"1") handle_dfiles $1 $2 "$mzPPM" "$rtPPM" "$centD" ;;
	"1") handle_dfiles $1 $2 "$mzPPM" "$rtPPM" "$centrD" "$APIK" "$AddF" "$DBF" ;;
	"2") dialog --yesno "Use MS-Dial and convert .raw files to .msp files and move .msp files to /raw data/exported as raw msp" 5 40 ;;
	"3") dialog --yesno "Use MS-dial and convert .mzML to .msp files and move .msp files to  /raw data/exprted as raw msp" 5 40;;
	"4")dialog --yesno "First convert mzXML to mzML using msconvert GUI and move .msp files to /raw data/exprted as raw msp" 5 40;;
esac




