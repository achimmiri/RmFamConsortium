library(rcdk)
library(rinchi)
library(rromeo)
library(readxl)
library(CAMERA)
library(RMassBank)
library(stringr)
library(stringi)
library(squash)
####################################
####################################
MaKlist<-function(FiN)
{
  #####################################
  #####################################
  gFile=FiN
  #####################################
  ####################################
  if(!is.na(gFile) && length(gFile) >=1){
  ##############################
  lines <- readLines(gFile)
  lst <-split(lines, cumsum(lines==""))
  lst1 <-lapply(lst, function(x) if (x[1] == "") x[-1] else x)
  LL<-sapply(lst, length)
  IR<-which(unname(LL) == 0)
  if(length(IR)>0){LL1<-LL[-(IR)]}else{LL1 <- LL}
  ###### MSP #############################
  lst2 <-lst1[names(LL1)]
  lst2<-purrr::compact(lst2)
  #####################################
  ##print(lst2)
  #####################################
  return(lst2)
  ##############################################
    }## end of test if loop
  ##############################################
}
###########################################################################
############################################################################
args <- commandArgs(TRUE)
DiN<-args[1]
FiN<-args[2]
EXN<-args[3]
INFL=args[4]
########################################
########################################
alis=MaKlist(FiN)
##############################################################################
##############################################################################

File1<-DiN
File2<-EXN
##############################################################################
##############################################################################
## 1) The argument 1 is the RMassBank files path location 
## 2) The argument 2 is the combined msp file location in RMassBank folder
## 3) The excel file location
## 4) ini file location
#############################################################################
##############################################################################
LF<-list.files(path=File1,pattern=".passed.msp")
##############################################################################

RXF<-readxl::read_excel(File2, sheet = 1, col_names = TRUE,skip=1,.name_repair="minimal")
RXF[] <- lapply(RXF, function(x) type.convert(as.character(x)))
RXF1 <- which(is.na(as.character(RXF[["File"]])))
##############################################################################
RXF3=list()
if(length(RXF1) >= 1)
{
 RXF2<- RXF[-RXF1,]
 RXF3=RXF2
}else
{
  RXF3=RXF
}
#####################################
#####################################

#############################################
#############################################
RL<-dim(RXF3)[1]
#####################################
#####################################
readEXCEL<-function(RXF)
{
	RXF<-readxl::read_excel(File2, sheet = 1, col_names = TRUE,skip=1,.name_repair="minimal")
	RXF[] <- lapply(RXF, function(x) type.convert(as.character(x)))
	RXF1 <- which(is.na(as.character(RXF[["File"]])))
	##############################################################################
	RXF3=list()
	if(length(RXF1) >= 1)
	{
  		RXF2<- RXF[-RXF1,]
  		RXF3=RXF2
	}else{
  		RXF3=RXF
	}
	##############################################################################
        return(RXF3)
}

#################################################


RFT<-read.table(paste(DiN,"Filelist.2.csv",sep="/"),sep=",",header=T,fill = TRUE,stringsAsFactors = FALSE)
UFV= unique(RFT$Files)
##########################################################################

##########################################################################
CN1<-c("ID","Name","SMILES","RT","CAS","Formula","Ontology","Adduct","Bpeak","Cenergy","IM","CClass","CONFIDENCE","PUBLICATION\n")
cat(CN1,file =paste(dirname(File2),"Compoundlist.2.csv",sep="/"),sep=",")
############################################################################
dat <- data.frame()
dat1 <- data.frame()
#######################################################
for(i in 1:length(UFV))
{
	######################
        RV=UFV[i]
	#######################
	#######################
	RV1=which(RFT == RV, arr.ind=TRUE)[,1]
	#############################
	#############################
	MV=MaKlist(RV)
	#############################
	##########################
	if(length(MV) == length(RV1))
	{
	####################################	
		print("enter the if loop")
		###################
		for(i in 1:length(MV))
		{
			#####################
			print("enter the if loop ...1")
			######################
			LV=MV[[i]]
                        #######################################################
			Id = RV1[i]
			#######################################################
			Nam = grep("NAME:",LV, value=TRUE)
			Nam1 = stringr::str_trim(sub('NAME: ', '', Nam))
			RT = grep("RETENTIONTIME:",LV, value=TRUE)
                        RT1 = stringr::str_trim(sub('RETENTIONTIME: ', '', RT))
                        ########################################################
			########################################################
                        SM = grep("SMILES:",LV, value=TRUE)
			SM1 = stringr::str_trim(sub('SMILES:', '', SM))
                        CAS = ""
			Form =  grep("FORMULA:",LV, value=TRUE)
                        Form1 = stringr::str_trim(sub('FORMULA: ', '', Form))
                        Ont = grep("Ontology:",LV, value=TRUE)
                        Ont1 = stringr::str_trim(sub('Ontology: ', '', Ont))
			pattern <- "PRECURSORTYPE|ADDUCTIONNAME"
			PTYL = grep(pattern,LV, value=TRUE)
			PTYL1 = stringr::str_trim(sub('PRECURSORTYPE: ', '', PTYL))
			PTYL2 = stringr::str_trim(sub('ADDUCTIONNAME: ', '', PTYL1))
			##NDF<-data.frame(Id,Nam1,SM1,RT1,CAS,Form1,Ont1,PTYL2)
                        Bpeak =  grep("PRECURSORMZ:",LV, value=TRUE)
			Bpeak1 = stringr::str_trim(sub('PRECURSORMZ: ', '', Bpeak))
                        Cenergy =  grep("COLLISIONENERGY:",LV, value=TRUE)
			Cenergy1 = stringr::str_trim(sub('COLLISIONENERGY: ', '', Cenergy))
                        IM = grep("IONMODE:",LV, value=TRUE)
			IM1 = stringr::str_trim(sub('IONMODE: ', '', IM))
			#################################
			INN=grep("INCHI:",LV, value=TRUE)
			INN1=stringr::str_trim(sub("INCHI:",'',INN))
			##################################
			SMV=grep("SMILES:",LV, value=TRUE)
			SMV1=stringr::str_trim(sub("SMILES:",'',SMV))
			##########################################################
			##########################################################
			###FN=which(tools::file_path_sans_ext(as.character(RXF3[["File"]]))==gsub('.passed.msp','',RV))
                        FN=which(tools::file_path_sans_ext(basename(as.character(RXF3[["File"]])))==gsub('.passed.msp','',basename(RV)))
			NAM=which(RXF3[["Name"]]==Nam1)
			FORM=which(RXF3[["Formula"]]==Form1)
			ADUC=which(RXF3[["Adduct"]]==PTYL2)
			CE=which(RXF3[["Collision energy"]]==Cenergy1)
			IM=which(RXF3[["Ionization mode"]]==IM1)
			######RTV=which(RXF3["RT (min)"]==RT1)
			RTV=which(abs(as.numeric(as.character(RXF3[["RT (min)"]]))-as.numeric(RT1)) <= 0.05)
			IKKV=which(as.character(RXF3[["InChI"]])==INN1)
			SMV2=which(as.character(RXF3[["SMILES"]])==SMV1)
			##################################
			InVa=Reduce(intersect, list(FN,NAM))
			###################################			
			###################################
			if(length(InVa) == 1){
				CClass<-as.character(RXF3[InVa,][["Compound class"]])
				CONFIDENCE<-as.character(RXF3[InVa,][["Confidence"]])
				PUBLICATION<-as.character(RXF3[InVa,][["Publication"]])
				NDF<-data.frame(Id,Nam1,SM1,RT1,CAS,Form1,Ont1,PTYL2,Bpeak1,Cenergy1,IM1,CClass,CONFIDENCE,PUBLICATION)
				dat <- rbind(dat,NDF)
				dat1 <- rbind(dat1,Id)
			}else{
				InVa=Reduce(intersect, list(FN,NAM,IM,RTV))
				if(length(InVa) == 1){
					CClass<-as.character(RXF3[InVa,][["Compound class"]])
                                	CONFIDENCE<-as.character(RXF3[InVa,][["Confidence"]])
                                	PUBLICATION<-as.character(RXF3[InVa,][["Publication"]])
                                	NDF<-data.frame(Id,Nam1,SM1,RT1,CAS,Form1,Ont1,PTYL2,Bpeak1,Cenergy1,IM1,CClass,CONFIDENCE,PUBLICATION)
                                	dat <- rbind(dat,NDF)
                                	dat1 <- rbind(dat1,Id)
				}else{
					InVa=Reduce(intersect, list(FN,NAM,IM,FORM,ADUC))
					if(length(InVa) == 1){
						CClass<-as.character(RXF3[InVa,][["Compound class"]])
                                        	CONFIDENCE<-as.character(RXF3[InVa,][["Confidence"]])
                                        	PUBLICATION<-as.character(RXF3[InVa,][["Publication"]])
                                        	NDF<-data.frame(Id,Nam1,SM1,RT1,CAS,Form1,Ont1,PTYL2,Bpeak1,Cenergy1,IM1,CClass,CONFIDENCE,PUBLICATION)
                                        	dat <- rbind(dat,NDF)
                                        	dat1 <- rbind(dat1,Id)
					}else{
						InVa=Reduce(intersect, list(FN,NAM,IM,FORM,ADUC,SMV2))
						if(length(InVa) == 1){
							CClass<-as.character(RXF3[InVa,][["Compound class"]])
                                                	CONFIDENCE<-as.character(RXF3[InVa,][["Confidence"]])
                                                	PUBLICATION<-as.character(RXF3[InVa,][["Publication"]])
                                                	NDF<-data.frame(Id,Nam1,SM1,RT1,CAS,Form1,Ont1,PTYL2,Bpeak1,Cenergy1,IM1,CClass,CONFIDENCE,PUBLICATION)
                                                	dat <- rbind(dat,NDF)
                                                	dat1 <- rbind(dat1,Id)
						}else{
							InVa=Reduce(intersect, list(FN,NAM,IM,FORM,ADUC,IKKV))
							if(length(InVa) == 1){
							CClass<-as.character(RXF3[InVa,][["Compound class"]])
                                                        CONFIDENCE<-as.character(RXF3[InVa,][["Confidence"]])
                                                        PUBLICATION<-as.character(RXF3[InVa,][["Publication"]])
                                                        NDF<-data.frame(Id,Nam1,SM1,RT1,CAS,Form1,Ont1,PTYL2,Bpeak1,Cenergy1,IM1,CClass,CONFIDENCE,PUBLICATION)
                                                        dat <- rbind(dat,NDF)
                                                        dat1 <- rbind(dat1,Id)
							}else{
								InVa=Reduce(intersect, list(FN,NAM,IM,CE))
								if(length(InVa) == 1){
									CClass<-as.character(RXF3[InVa,][["Compound class"]])
                                                        		CONFIDENCE<-as.character(RXF3[InVa,][["Confidence"]])
                                                        		PUBLICATION<-as.character(RXF3[InVa,][["Publication"]])
                                                        		NDF<-data.frame(Id,Nam1,SM1,RT1,CAS,Form1,Ont1,PTYL2,Bpeak1,Cenergy1,IM1,CClass,CONFIDENCE,PUBLICATION)
                                                        		dat <- rbind(dat,NDF)
                                                        		dat1 <- rbind(dat1,Id)
								}else{
									InVa=Reduce(intersect, list(FN,NAM,IM,CE,FORM))
									if(length(InVa) == 1){
										CClass<-as.character(RXF3[InVa,][["Compound class"]])
                                                                        	CONFIDENCE<-as.character(RXF3[InVa,][["Confidence"]])
                                                                        	PUBLICATION<-as.character(RXF3[InVa,][["Publication"]])
                                                                        	NDF<-data.frame(Id,Nam1,SM1,RT1,CAS,Form1,Ont1,PTYL2,Bpeak1,Cenergy1,IM1,CClass,CONFIDENCE,PUBLICATION)
                                                                        	dat <- rbind(dat,NDF)

									}else{
										InVa=Reduce(intersect, list(FN,NAM,IM,CE,FORM,ADUC))
										if(length(InVa) == 1){
											CClass<-as.character(RXF3[InVa,][["Compound class"]])
                                                                                	CONFIDENCE<-as.character(RXF3[InVa,][["Confidence"]])
                                                                                	PUBLICATION<-as.character(RXF3[InVa,][["Publication"]])
                                                                                	NDF<-data.frame(Id,Nam1,SM1,RT1,CAS,Form1,Ont1,PTYL2,Bpeak1,Cenergy1,IM1,CClass,CONFIDENCE,PUBLICATION)
                                                                                	dat <- rbind(dat,NDF)
										}else{
											InVa=Reduce(intersect, list(FN,NAM,IM,CE,FORM,ADUC,SMV2))
											if(length(InVa) == 1){
												CClass<-as.character(RXF3[InVa,][["Compound class"]])
                                                                                        	CONFIDENCE<-as.character(RXF3[InVa,][["Confidence"]])
                                                                                        	PUBLICATION<-as.character(RXF3[InVa,][["Publication"]])
                                                                                        	NDF<-data.frame(Id,Nam1,SM1,RT1,CAS,Form1,Ont1,PTYL2,Bpeak1,Cenergy1,IM1,CClass,CONFIDENCE,PUBLICATION)
                                                                                        	dat <- rbind(dat,NDF)
											}else{
												InVa=Reduce(intersect, list(FN,NAM,IM,CE,FORM,ADUC,IKKV))
												if(length(InVa) == 1){
													CClass<-as.character(RXF3[InVa,][["Compound class"]])
                                                                                                	CONFIDENCE<-as.character(RXF3[InVa,][["Confidence"]])
                                                                                                	PUBLICATION<-as.character(RXF3[InVa,][["Publication"]])
                                                                                                	NDF<-data.frame(Id,Nam1,SM1,RT1,CAS,Form1,Ont1,PTYL2,Bpeak1,Cenergy1,IM1,CClass,CONFIDENCE,PUBLICATION)
                                                                                                	dat <- rbind(dat,NDF)
												}else{
													InVa=Reduce(intersect, list(FN,NAM,RTV))
													if(length(InVa) == 1){
														CClass<-as.character(RXF3[InVa,][["Compound class"]])
                                                                                                        	CONFIDENCE<-as.character(RXF3[InVa,][["Confidence"]])
                                                                                                        	PUBLICATION<-as.character(RXF3[InVa,][["Publication"]])
                                                                                                        	NDF<-data.frame(Id,Nam1,SM1,RT1,CAS,Form1,Ont1,PTYL2,Bpeak1,Cenergy1,IM1,CClass,CONFIDENCE,PUBLICATION)
                                                                                                        	dat <- rbind(dat,NDF)
													}else{
														InVa=Reduce(intersect, list(FN,NAM,IM,SMV2))
														if(length(InVa) == 1){
															CClass<-as.character(RXF3[InVa,][["Compound class"]])
                                                                                                                	CONFIDENCE<-as.character(RXF3[InVa,][["Confidence"]])
                                                                                                                	PUBLICATION<-as.character(RXF3[InVa,][["Publication"]])
                                                                                                                	NDF<-data.frame(Id,Nam1,SM1,RT1,CAS,Form1,Ont1,PTYL2,Bpeak1,Cenergy1,IM1,CClass,CONFIDENCE,PUBLICATION)
                                                                                                                	dat <- rbind(dat,NDF)
														}else{
															InVa=Reduce(intersect, list(FN,NAM,IKKV))

															if(length(InVa) == 1){

															CClass<-as.character(RXF3[InVa,][["Compound class"]])
                                                                                                                        CONFIDENCE<-as.character(RXF3[InVa,][["Confidence"]])
                                                                                                                        PUBLICATION<-as.character(RXF3[InVa,][["Publication"]])
                                                                                                                        NDF<-data.frame(Id,Nam1,SM1,RT1,CAS,Form1,Ont1,PTYL2,Bpeak1,Cenergy1,IM1,CClass,CONFIDENCE,PUBLICATION)
                                                                                                                        dat <- rbind(dat,NDF)
															}else{
																print("there must be some error in generation of msp file check ...msp file again")

															}
														}

													} ##end of else ..where they check ...FN,NAM,SMV2



												} ### end of the else loop ...where they chek ...FN,NAM,RTV

											} ## end of the else loop ..where they check ...FN,NAM,IM,CE,FORM,ADUC,IKKV
										}## end of the else loop ..where they check ...FN,NAM,IM,CE,FORM,ADUC,SMV2


									}### end of the else loop...where they check ...FN,NAM,IM,CE,FORM

								}### end of the else loop...where they check ...FN,NAM,IM,CE...introducing collision energy 


							}### end of else loop ...where they check ...FN,NAM,IM,FORM,ADUC,IKKV

						}##end of else loop ...where they check ..FN,NAM,IM,FORM,ADUC,SMV2 
					}###end of else loop ...where they check...FN,NAM,IM,FORM,ADUC 

				}
			}## end of the else ...where the comparision is filename and compound name
		}## inside for loop
	}###if loop length(MV) == length(RV1)
}### End of the for loop
########################################################################################################################################################
###########################################################################################################
###########################################################################################################
write.table(dat, file = paste(dirname(File2),"Compoundlist.2.csv",sep="/"), sep = ",",quote=TRUE,row.names = F, col.names = F,append=T)
##write.table(dat2, file = paste(dirname(File2),"Compoundlist.2.csv",sep="/"), sep = ",",quote=TRUE,row.names = F, col.names = F,append=T)
#############################################################################################################
############################################################################################################
##RF<-readLines(paste(DiN,"RMB_options.ini",sep="/"))
RF<-readLines(INFL)
################################
################################
AuIN<-grep("authors:",RF)
RF[AuIN]<-paste0("    authors: ",unique(as.character(RXF3[["Authors"]])))
CoIN<-grep("copyright:",RF)
RF[CoIN]<-paste0("    copyright:  Copyright (C) ",unique(as.character(RXF3[["Authors"]])))
PuIN<-grep("publication:",RF)
RF[PuIN]<-paste0("    publication: ",unique(as.character(RXF3[["Publication"]])))
LIIN<-grep("license:",RF)
##RF[LIIN]
inIN<-grep("instrument:",RF)
RF[inIN]<-paste0("    instrument: ",unique(as.character(RXF3[["INSTRUMENT"]])))
intyIN<-grep("instrument_type:",RF)
RF[intyIN]<-paste0("    instrument_type: ",unique(as.character(RXF3[["INSTRUMENT_TYPE"]])))
inIONI<-grep("ionization:",RF)
RF[inIONI]<-paste0("    ionization: ",unique(as.character(RXF3[["IONIZATION"]])))
#################
############################################
concomIN<-grep("confidence_comment",RF)
RF[concomIN]<-paste0("    confidence_comment: ",unique(as.character(RXF3[["Confidence"]])))
##################################################
##################################################
EPIN<-grep("entry_prefix:",RF)
RF[EPIN]<-paste0("    entry_prefix: "," mC")
##################################################
##################################################
inCOCL<-grep("compound_class:",RF)
RF[inCOCL]<-paste0("    compound_class: ",unique(as.character(RXF3[["Compound class"]])))
RF[EPIN]<-paste0("    entry_prefix: "," mC")
msIN<-grep("ms_type:",RF)
RF[msIN]<-paste0("    ms_type: "," MS2")
######################################################
inCE <- grep("COLLISION_ENERGY:",RF)
##RF[inCE]<-paste0("    COLLISION_ENERGY: ",unique(as.character(RXF3[["Collision energy"]])))
######################################################
RF[inCE]<-paste0("    COLLISION_ENERGY: ",gsub("[^0-9]","",unique(as.character(RXF3[["Collision energy"]]))))
#######################################################
#######################################################
##writeLines(RF,file(paste(dirname(File2),"RMB_options.test2.ini",sep="/")),sep="\n")
cat(RF,file=paste(dirname(File2),"RMB_options.test2.ini",sep="/"), append=TRUE, sep = "\n")
