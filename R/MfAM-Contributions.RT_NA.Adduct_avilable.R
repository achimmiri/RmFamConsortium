library(circlize)
library(devtools)
library(gdata)
library(jsonlite)
library(metaMS)
library(plotrix)
library(qdapRegex)
library(rJava)
library(rcdk)
library(rinchi)
library(rromeo)
library(readxl)
library(CAMERA)
library(RMassBank)
library(stringr)
library(stringi)
library(squash)
library(tools)
library(webchem)
library(zeallot)
library(BBmisc)
library(purrr)
library(schoolmath)
library(classyfireR)
library(sen2r)
library(sjmisc)
library(OrgMassSpecR)
library(ChemmineOB)
library("metaMS")
library("rcellminer")
library("R.utils")
library(RCurl)
#########################################################################
##install.packages("sjmisc")
##BiocManager::install("RMassBank")
###BiocManager::install("metaMS")
##BiocManager::install("rcellminerData")
##BiocManager::install("rcellminer")
##BiocManager::install("ChemmineOB")
##BiocManager::install("zlibbioc")
#########################################################################
##devtools::install_github("r-lib/xml2")
##install.packages("RcppEigen",dependencies = TRUE)
##install.packages('Cairo',dependencies = TRUE)
########################## Adding this new package ######################
##BiocManager::install("Rdisop")
#########################################################################
library(Rdisop)
########## Step1: Reading the API key ####################################
##########################################################################
args <- commandArgs(TRUE)
##########################################################################
##########################################################################
### Function to check the number of arguments
CountNumberofArgs<-function(argsV)
{
        if(length(argsV)==8)
        {
                return(length(args))
        }else{
		return(length(args))
                ##message("the number of arguments are less than 8 ..please check and pass the correct arguments")
                ##stop("the number of arguments are less than 8 ..please check and pass the correct arguments")
        }
}
###################################################
argsV1=ifelse(CountNumberofArgs(args) == 8,args,stop("the number of arguments are less than 8 ..please check and pass the correct arguments"))
File1<-argsV1[1]
mz_Tol=as.numeric(argsV1[2])
##RT_Tol=(as.numeric(argsV1[3])/100)
RT_Tol=ifelse(as.numeric(argsV1[3]),(as.numeric(argsV1[3])/100),iflese(is.empty(argsV1[3]),"NA","NA"))
DEFAULT_MZ_TOLERANCE=as.numeric(argsV1[4])
apikey = argsV1[5]
############################################################################
############################################################################
Sys.setenv(CHEMSPIDER_KEY = apikey)
rr_auth(apikey)
## adding two more arguments to read the ADI.4.csv ..which is the adduct file
AIN1 = argsV1[6]
AD = argsV1[7]
##########################################################################
##mz.25ppm.20RT.R
##########################################################################
##File1<-"/mnt/ifs/data/IPB/Projects/2017_005_MS-databases/mFam contributions/maximilian_frey@uni-hohenheim.de/meta data/230181218_Frey_Compound_Spreadsheet_For_MSMS_v17_GB.xlsx"
##########################################################################
##print(File1)
##########################################################################
Fi <- unlist(strsplit(File1, "/"))
Fi1 <-c(paste(Fi[-length(Fi)], collapse = "/"), last(Fi))
Fi2<-Fi1[1]
Fi3<-unlist(strsplit(Fi2, "/"))
Fi4<-c(paste(Fi3[-length(Fi3)], collapse = "/"), last(Fi3))
Fi5<-Fi4[1]
Fi6<-paste(Fi5,"raw data","exported as raw msp",sep="/")
### Input################
Fi7<-paste(Fi6,"/",sep="")
Fi8<-paste(Fi5,"converted to msp",sep="/")
### Output##################
Fi9<-paste(Fi8,"/",sep="")
############################
Fi10<-paste(Fi8,"mz.25ppm.20RT",sep="/")
Fi11<-paste(Fi10,"/",sep="")
############################
Fi12<-paste(Fi8,"mz.40ppm.35RT",sep="/")
Fi13<-paste(Fi12,"/",sep="")
############################
Fi14<-paste(Fi8,"mz.50ppm.40RT",sep="/")
Fi15<-paste(Fi14,"/",sep="")
###############################
NFi5<-unlist(strsplit(Fi5, "/"))
NFi6<-c(paste(NFi5[-length(NFi5)], collapse = "/"), last(NFi5))
N1Fi6<-NFi6[1]
###############################################################
##print(Fi6)
################################################################
####### Step2: Reading the All the required files
##AIN<-read.table("/mnt/ifs/data/IPB/Projects/2017_005_MS-databases/mFam contributions/scripts/ADI.txt",sep="\t",header=F,quote="",stringsAsFactors = FALSE)
##################################################################
AIN1<-read.table("/mnt/ifs/data/IPB/Projects/2017_005_MS-databases/mFam contributions/scripts/ADI.4.csv",sep=",",header=F,quote="",stringsAsFactors = FALSE)
AIN<-AIN1
####################### adding the new database table
AD<-read.table("/mnt/ifs/data/IPB/Projects/2017_005_MS-databases/mFam contributions/scripts/Database_Dec2017.txt",sep="\t",header=T,fill = TRUE,stringsAsFactors = FALSE)
numbers_only <- function(x) !grepl("\\D", x)
####################Step3: READ the meta data file ##############
RXF<-readxl::read_excel(File1, sheet = 1, col_names = TRUE,skip=1,.name_repair="minimal")
RXF[] <- lapply(RXF, function(x) type.convert(as.character(x)))
RXF1 <- which(is.na(as.character(RXF[["File"]])))
#################################################################
RXF3=list()
if(length(RXF1) >= 1)
{
  #print("enter the line ...71")
  RXF2<- RXF[-RXF1,]
  RXF3=RXF2
}else
{
  RXF3=RXF
}
####### THis is test area ################################
##RXF4=RXF3
##RXF3=RXF4[1:25,]
#########################################################

#########################################################
######### Adding this new################################
Lmeda<-tryCatch({base::rle(as.character(RXF3[["File"]]))},warning=function(cond){message("file structure is wrong")})
NLmeda<-tryCatch({Lmeda$lengths},warning=function(cond){message("repeat length calculate mistakes")})
LmeCmu<-tryCatch({abs(cumsum(NLmeda))},warning=function(cond){message("cumulative calculate mistakes happening")})
LmeCmu1<-tryCatch({c(0,LmeCmu)},warning=function(cond){message("cumulative calculate mistakes happening")})
SFileNam<-tryCatch({Lmeda$values},warning=function(cond){message("rle is not able to fetch the information properly")})
#########################################################
##print(LmeCmu1)
#########################################################
PuInKtoSM<-function(getINK)
{
  url<- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"
  out<-tryCatch({jsonlite::fromJSON(paste0(url,getINK, "/JSON"))} ,warning = function(x) {return(NA)})
  prop.names  <-tryCatch({out$PC_Compounds$props[[1]][[1]]},warning = function(x) {return(NA)})
  prop.values <- tryCatch({out$PC_Compounds$props[[1]][[2]]},warning = function(x) {return(NA)})
  sm <-tryCatch({grep("smiles", prop.names[,"label"], ignore.case = TRUE)},warning = function(x) {return(NA)})
  csmiles<-c()
  if(length(sm) >= 1) {
    can <- tryCatch({grep("canonical", prop.names[,"name"], ignore.case = TRUE)},warning= function(x) {return(NA)})
    can1<-tryCatch({prop.values[sm[1],"sval"]},warning= function(x) {return(NA)})
    csmiles<-c(csmiles,can1)
  }else{
    csmiles<-c(csmiles,NA)
  }
  return(csmiles)
  
}
###########################################################
PuInKtoSM1<-function(getINK)
{
  url<- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"
  out<-tryCatch({jsonlite::fromJSON(paste0(url,getINK, "/JSON"))}, warning = function(x) {return(NA)})
  prop.names  <-tryCatch({out$PC_Compounds$props[[1]][[1]]}, warning = function(x) {return(NA)})
  prop.values <- tryCatch({out$PC_Compounds$props[[1]][[2]]}, warning = function(x) {return(NA)})
  sm <-tryCatch({grep("smiles", prop.names[,"label"], ignore.case = TRUE)}, warning = function(x) {return(NA)})
  csmiles<-c()
  if(length(sm) == 2) {
    can <- tryCatch({grep("canonical", prop.names[,"name"], ignore.case = TRUE)}, warning = function(x) {return(NA)})
    can1<-tryCatch({prop.values[sm[2],"sval"]},warning= function(x) {return(NA)})
    csmiles<-c(csmiles,can1)
  }else{
    csmiles<-c(csmiles,NA)
  }
  return(csmiles)
  
}
##########################################################################
##print("coming to area before Ontology")
##########################################################################
#MaKE.ONT.REC<-function(InMEDA)
#{
#  out<-c()
#  ################
#   print("entering the ontology area")
#   ###############
#  if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI="))
#  {
#    ###################################
#    print("enter the line 6")
#    ###################################
#    IN<-as.character(InMEDA[["InChI"]])
#    mol <-tryCatch({rinchi::parse.inchi(IN)},error=function(cond){message("name is empty")})
#    SM<-tryCatch({rcdk::get.smiles(mol[[1]])},error=function(cond){message("name is empty")})
#    IK<-tryCatch({rinchi::get.inchi.key(SM)},error=function(cond){message("name is empty")})
#    IK1<-paste("INCHIKEY:",IK,sep=" ")
#    ############################
#    if(!sjmisc::is_empty(IK)){
#	    ##################
#	    print("enter the line 8")
#      ##########################
#      IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(IK)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
#      if(!sjmisc::is_empty(IKCRV)){
#      ##IKCRV<-tryCatch({classyfireR::get_classification(IK)},warning=function(cond){message("Classifier could not fecth the information")})
#      ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
#      ##########################
#      IK1<-paste("INCHIKEY:",IK,sep=" ")
#      tes2<-paste("Ontology:",ONTV,sep=" ")
#      FINCH<-paste("INCHI:",IN,sep=" ")
#      #################
#      out<-c(out,tes2)
#      out<-c(out,IK1)
#      out<-c(out,FINCH)
#      ########################
#      }else{
#	       IK1<-paste("INCHIKEY:",IK,sep=" ")
#               tes2<-paste("Ontology:","",sep=" ")
#               FINCH<-paste("INCHI:",IN,sep=" ")
#	       ################
#	       out<-c(out,tes2)
#	       out<-c(out,IK1)
#	       out<-c(out,FINCH)
#
#      }
#   ####################################################################
#    }else if(!sjmisc::is_empty(SM))
#    {
#      ##########################################################
#      tes<-tryCatch({classyfireR::submit_query(label = 'query_test', input = SM, type = 'STRUCTURE')},warning=function(cond){message("Classyfire is empty")})
#      tes1<-do.call(paste, c(as.list(tryCatch({tes@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
#      ###########################################################
#      ##IK<-tryCatch({rinchi::get.inchi.key(SM)},warning=function(cond){message("rinchi is not able to fecth")})
#      ##########################################################
#      IK1<-paste("INCHIKEY:",IK,sep=" ")
#      tes2<-paste("Ontology:",tes1,sep=" ")
#      FINCH<-paste("INCHI:",IN,sep=" ")
#      ########################
#      out<-c(out,tes2)
#      out<-c(out,IK1)
#      out<-c(out,FINCH)
#      ########################
#    }else{
#      ############################################
#      F1ONT<-paste("Ontology:","",sep=" ")
#      FINCH<-paste("INCHI:",SM,sep=" ")
#      ##############################################
#      if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){
#        FINCH<-paste("INCHI:",as.character(InMEDA[["InChI"]]),sep=" ")
#        IK<-paste("INCHIKEY:","",sep=" ")
#        ##################
#        out<-c(out,F1ONT)
#        out<-c(out,IK)
#        out<-c(out,FINCH)
#        #################
#      }else{
#	###################################
#        FINCH<-paste("INCHI:","",sep=" ")
#        IK<-paste("INCHIKEY:",as.character(InMEDA[["InChI"]]),sep=" ")
#        ###############
#        out<-c(out,F1ONT)
#        out<-c(out,IK)
#        out<-c(out,FINCH)
#        ####################
#      }
#      #################
#    }
##################################################################################################################
###  ##}else if(!is.na(as.character(RRV[["InChI"]])) & startsWith(as.character(RRV[["InChI"]]),'CAS:')){
###	  cai<-as.character(RRV[["InChI"]])
###	  cai1<-stringr::str_trim(gsub("CAS:","",cai))
###	  tes<-tryCatch({webchem::cir_query(cai1, "smiles")},warning=function(cond){message("some mistake happened in file search files")})
###	  tes1<-tryCatch({tes[[1]]},warning=function(cond){message("some mistake happened in file search files")})
###	  IK<-tryCatch({rinchi::get.inchi.key(tes1)},warning=function(cond){message("rinchi could not fetch missing")})
###	  IN<-tryCatch({rinchi::get.inchi(tes1)},warning=function(cond){message("rinchi could not fetch missing")})
### ## }
######################################################################################################################
#  }else if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]]))){
#    ############################
#    print("enter the line ...53")
#    print(as.character(InMEDA[["InChI"]]))
#    ############################
#    if(tryCatch({webchem::is.inchikey(as.character(InMEDA[["InChI"]]))},warning=function(cond){message("inchikey validation failed")}))
#    {
#      print("enter the line 55")
#      #######################################
#      tes<-tryCatch({webchem::get_cid(stringr::str_trim(as.character(InMEDA[["InChI"]])), from = "inchikey")},error=function(cond){message("webchem not able to get cid from Inchi")})
#      tes1<-tryCatch({tes$cid},error=function(cond){message("Inchi to CID did not convert")})
#      tes2<-tryCatch({webchem::pc_prop(as.numeric(tes1), properties = c("MolecularFormula", "MolecularWeight","CanonicalSMILES","InChI","InChIKey"))},error=function(cond){message("Inchi to CID did not convert so did not get properties")})
#      IN<-tryCatch({tes2$InChI},error=function(cond){message("Inchi to Inchikey failed because of CID not converting")})
#      #######################################
#      #######################################
#      FINCH<-paste("INCHI:",IN,sep=" ")
#      #####################################
#      teIK1<-as.character(InMEDA[["InChI"]])
#      IK1<-paste("INCHIKEY:",teIK1,sep=" ")
#      ###################################
#      if(!sjmisc::is_empty(teIK1)){
#        #####################################
#        print("enter the line ...57")		
#        ########################################
#	IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(teIK1)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
#        if(!sjmisc::is_empty(IKCRV)){
#	##############################
#        ##IKCRV<-tryCatch({classyfireR::get_classification(teIK1)},warning=function(cond){message("Classifier could not fecth the information")})
#        ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
#        ####################################
#        ####################################
#        tes2<-paste("Ontology:",ONTV,sep=" ")
#        #####################################
#        out<-c(out,tes2)
#        out<-c(out,IK1)
#        out<-c(out,FINCH)
#        ####################
#	}else{
#		##############################
#		IK1<-paste("INCHIKEY:",teIK1,sep=" ")
#		tes2<-paste("Ontology:","",sep=" ")
#		FINCH<-paste("INCHI:",IN,sep=" ")
#		#############################
#		out<-c(out,tes2)
#		out<-c(out,IK1)
#		out<-c(out,FINCH)
#		##############################
#	}
#        ##print(out)
#        ####################################
#      }else{
#        ####################################
#        print("enter the line ...89")
#        ################################	
#        SMV<-tryCatch({PuInKtoSM(teIK1)},error=function(cond){message("Pubchem fetch is empty")})
#        IK<-tryCatch({rinchi::get.inchi.key(SMV)},error=function(cond){message("rinchi did not get inchi key")})
#        IV<-tryCatch({rinchi::get.inchi(SMV)},error=function(cond){message("rinchi did not get inchi")})
#        FINCH<-paste("INCHI:",IV,sep=" ")
#        IK1<-paste("INCHIKEY:",IK,sep=" ")
#        ###################################
#        if(!sjmisc::is_empty(IK)){
#          #################################
#	  IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(IK)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
#	  if(!sjmisc::is_empty(IKCRV)){
#		  #############################
#		  print("enter the line ...91")
#		  ##############################
#          ##IKCRV<-tryCatch({classyfireR::get_classification(IK)},warning=function(cond){message("Classifier could not fecth the information")})
#          ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
#          ####################################
#          ##IK1<-paste("INCHIKEY:",IK,sep=" ")
#          tes2<-paste("Ontology:",ONTV,sep=" ")
#          ####################################
#          out<-c(out,tes2)
#          out<-c(out,IK1)
#          out<-c(out,FINCH)
#	  #####################################
#	  }else{
#		  ####################################
#		  IK1<-paste("INCHIKEY:",IK,sep=" ")
#		  tes2<-paste("Ontology:","",sep=" ")
#		  FINCH<-paste("INCHI:",IV,sep=" ")
#		  ######################################
#		  out<-c(out,tes2)
#		  out<-c(out,IK1)
#		  out<-c(out,FINCH)
#		  #####################################
#	  }
#        #####################################################
#        }else if(!sjmisc::is_empty(SMV)){
#		print("enter the line 95")
#	  ##########################################################	
#          tes<-tryCatch({classyfireR::submit_query(label = 'query_test', input = SMV, type = 'STRUCTURE')},warning=function(cond){message("Classyfire is empty")})
#          tes1<-do.call(paste, c(as.list(tryCatch({tes@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
#          ###########################################################
#          IK<-tryCatch({rinchi::get.inchi.key(SMV)},error=function(cond){message("rinchi is not able to fecth")})
#          IK1<-paste("INCHIKEY:",IK,sep=" ")
#          tes2<-paste("Ontology:",tes1,sep=" ")
#          ########################
#          out<-c(out,tes2)
#          out<-c(out,IK1)
#          out<-c(out,FINCH)
#	  #######################
#        }else{
#          F1ONT<-paste("Ontology:","",sep=" ")
#          if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]]))& startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){
#	    ####################################
#            FINCH<-paste("INCHI:",as.character(InMEDA[["InChI"]]),sep=" ")
#            IK<-paste("INCHIKEY:","",sep=" ")
#            ###############
#            out<-c(out,F1ONT)
#            out<-c(out,IK)
#            out<-c(out,FINCH)
#            ###############
#          }else{
#            ##################################
#            FINCH<-paste("INCHI:","",sep=" ")
#            IK<-paste("INCHIKEY:",as.character(InMEDA[["InChI"]]),sep=" ")
#            ###############
#            out<-c(out,F1ONT)
#            out<-c(out,IK)
#            out<-c(out,FINCH)
#            ####################
#          }
#        }
#      }
#    }else{
#      ###################################################
#      F1ONT<-paste("Ontology:","",sep=" ")
#      if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){
#	################################
#        FINCH<-paste("INCHI:",as.character(InMEDA[["InChI"]]),sep=" ")
#        IK<-paste("INCHIKEY:","",sep=" ")
#        ###############
#        out<-c(out,F1ONT)
#        out<-c(out,IK)
#        out<-c(out,FINCH)
#        ###############
#      }else{
#	####################################
#        FINCH<-paste("INCHI:","",sep=" ")
#        IK<-paste("INCHIKEY:",as.character(InMEDA[["InChI"]]),sep=" ")
#        ###############
#        out<-c(out,F1ONT)
#        out<-c(out,IK)
#        out<-c(out,FINCH)
#        ####################
#      }
#    }
#  ################################################  
#    
#
#
#  }else{
#    #######################################
#    print("entering this line...174")
#  ##########################################
#    PCID<- as.character(InMEDA[["PubChem CID"]])
#    PCSM<- as.character(InMEDA[["SMILES"]])
#    if(!sjmisc::is_empty(PCID))
#    {
#      #######################################
#      PCID1<-as.numeric(PCID)
#      tes<-tryCatch({webchem::pc_prop(PCID1, properties = c("MolecularFormula", "MolecularWeight","CanonicalSMILES","InChI","InChIKey"))},warning=function(cond){message("Did not get properties from Pubchem CID")})
#      gSMI<-tryCatch({tes$CanonicalSMILES},error=function(cond){message("smiles is not found")})
#      IK<-tryCatch({tes$InChIKey},error=function(cond){message("some mistake happened in file search files")})
#      IN<-tryCatch({tes$InChI},error=function(cond){message("some mistake happened in file search files")})
#      IK1<-paste("INCHIKEY:",IK,sep=" ")
#      #########################################################################
#      #########################################################################
#      if(!sjmisc::is_empty(IK)){
#      ###############################
#	IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(IK)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
#        if(!sjmisc::is_empty(IKCRV)){
#        ##IKCRV<-tryCatch({classyfireR::get_classification(IK)},warning=function(cond){message("Classifier could not fecth the information")})
#        ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
#        ##################################
#        IK1<-paste("INCHIKEY:",IK,sep=" ")
#        tes2<-paste("Ontology:",ONTV,sep=" ")
#        ##teIK2<-tryCatch({webchem::cs_convert(IK,from="inchikey",to="inchi")},error=function(cond){message("webchecm could not fetch the info")})
#        FINCH<-paste("INCHI:",IN,sep=" ")
#        ###################################
#        out<-c(out,tes2)
#        out<-c(out,IK1)
#        out<-c(out,FINCH)
#	################################
#	}else{
#		IK1<-paste("INCHIKEY:",IK,sep=" ")
#		tes2<-paste("Ontology:","",sep=" ")
#		FINCH<-paste("INCHI:",IN,sep=" ")
#		out<-c(out,tes2)
#		out<-c(out,IK1)
#		out<-c(out,FINCH)
#	}
#
#        #################################
#      }else if(!sjmisc::is_empty(gSMI)){
#	      ##################################
#	      print("enter the line 184")
#        ##################################################
#        tes<-tryCatch({classyfireR::submit_query(label = 'query_test', input = gSMI, type = 'STRUCTURE')},warning=function(cond){message("adduct value is missing")})
#        tes1<-do.call(paste, c(as.list(tryCatch({tes@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
#        ################################################
#        tes2<-paste("Ontology:",tes1,sep=" ")
#        ##teIK2<-tryCatch({webchem::cs_convert(IK,from="inchikey",to="inchi")},error=function(cond){message("webchecm could not fetch the info")})
#        FINCH<-paste("INCHI:",IN,sep=" ")
#        ##################
#        out<-c(out,tes2)
#        out<-c(out,IK1)
#        out<-c(out,FINCH)
#        #################
#      }else{
#        F1ONT<-paste("Ontology:","",sep=" ")
#        if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){
#	  #########################################
#          FINCH<-paste("INCHI:",as.character(InMEDA[["InChI"]]),sep=" ")
#          IK<-paste("INCHIKEY:","",sep=" ")
#          #################
#          out<-c(out,F1ONT)
#          out<-c(out,IK)
#          out<-c(out,FINCH)
#          ################
#        }else{
#          FINCH<-paste("INCHI:","",sep=" ")
#          IK<-paste("INCHIKEY:",as.character(InMEDA[["InChI"]]),sep=" ")
#          #################
#          out<-c(out,F1ONT)
#          out<-c(out,IK)
#          out<-c(out,FINCH)
#          ####################
#        }
#      }
#    }else{
#      #############################################
#	    print("enter the line 194")
#	    ##########################################
#      SMV<-as.character(InMEDA[["SMILES"]])
#      #############################################
#      if(!sjmisc::is_empty(SMV)){
#	############################################      
#        IK<-tryCatch({rinchi::get.inchi.key(SMV)},error=function(cond){message("rinchi could not fetch inchikey missing")})
#        SMV1<-tryCatch({rinchi::get.inchi(SMV)},error=function(cond){message("rinchi could not fetch inchi missing")})
#        IK1<-paste("INCHIKEY:",IK,sep=" ")
#        FINCH<-paste("INCHI:",SMV1,sep=" ")
#        ###################################
#        if(!sjmisc::is_empty(IK)){
#          ######################################
#	  IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(IK)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
#          if(!sjmisc::is_empty(IKCRV)){
#          ##IKCRV<-tryCatch({classyfireR::get_classification(IK)},warning=function(cond){message("Classifier could not fecth the information")})
#          ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
#          #############################
#          IK1<-paste("INCHIKEY:",IK,sep=" ")
#          tes2<-paste("Ontology:",ONTV,sep=" ")
#          ############################
#          out<-c(out,tes2)
#          out<-c(out,IK1)
#          out<-c(out,FINCH)
#	  #########################
#	  }else{
#		  IK1<-paste("INCHIKEY:",IK,sep=" ")
#		  tes2<-paste("Ontology:","",sep=" ")
#		  ##FINCH<-paste("INCHI:",IN,sep=" ")
#		  FINCH<-paste("INCHI:",SMV1,sep=" ")
#		  ################################
#		  out<-c(out,tes2)
#		  out<-c(out,IK1)
#		  out<-c(out,FINCH)
#	  }
#          ############################
#        }else if(!sjmisc::is_empty(SMV)){
#          ###############################
#          ##F1ONT<-paste("Ontology:","",sep=" ")
#          ###############################
#          tes<-tryCatch({classyfireR::submit_query(label = 'query_test', input = SMV, type = 'STRUCTURE')},warning=function(cond){message("Classfire not able to fetch empty")})
#          tes1<-do.call(paste, c(as.list(tryCatch({tes@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
#          tes2<-paste("Ontology:",tes1,sep=" ")
#          teIK2<-tryCatch({webchem::cs_convert(IK,from="inchikey",to="inchi")},error=function(cond){message("webchecm could not fetch the info")})
#          FINCH<-paste("INCHI:",teIK2,sep=" ")
#          ########################
#          out<-c(out,tes2)
#          out<-c(out,IK1)
#          out<-c(out,FINCH)
#          ##################
#        }
#      }else{
#        ##############################################
#        F1ONT<-paste("Ontology:","",sep=" ")
#        if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){
#          FINCH<-paste("INCHI:",as.character(InMEDA[["InChI"]]),sep=" ")
#          IK<-paste("INCHIKEY:","",sep=" ")
#          ###############
#          out<-c(out,F1ONT)
#          out<-c(out,IK)
#          out<-c(out,FINCH)
#          ###############
#        }else{
#          FINCH<-paste("INCHI:","",sep=" ")
#          IK<-paste("INCHIKEY:",as.character(InMEDA[["InChI"]]),sep=" ")
#          ###############
#          out<-c(out,F1ONT)
#          out<-c(out,IK)
#          out<-c(out,FINCH)
#          ####################
#        }
#      }
#    }
#  }
#  return(out)
#}
#
##########################################################################################
##########################################################################################
MaKE.ONT.REC<-function(InMEDA)
{
  out<-c()
  ################
  print("entering the ontology area")
  ###############
  if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI="))
  {
    ###################################
    print("enter the line 6")
    ###################################
    IN<-as.character(InMEDA[["InChI"]])
    mol <-tryCatch({rinchi::parse.inchi(IN)},error=function(cond){message("name is empty")})
    SM<-tryCatch({rcdk::get.smiles(mol[[1]])},error=function(cond){message("name is empty")})
    IK<-tryCatch({rinchi::get.inchi.key(SM)},error=function(cond){message("name is empty")})
    IK1<-paste("INCHIKEY:",IK,sep=" ")
    ############################
    if(!sjmisc::is_empty(IK)){
      ########################
      print("enter the line 8")
      ##########################
      IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(IK)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
      if(!sjmisc::is_empty(IKCRV)){
        #############################
        ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
        ##########################
        IK1<-paste("INCHIKEY:",IK,sep=" ")
        tes2<-paste("Ontology:",ONTV,sep=" ")
        FINCH<-paste("INCHI:",IN,sep=" ")
        #########################
        out<-c(out,tes2)
        out<-c(out,IK1)
        out<-c(out,FINCH)
        ########################
      }else{
        IK1<-paste("INCHIKEY:",IK,sep=" ")
        tes2<-paste("Ontology:","",sep=" ")
        FINCH<-paste("INCHI:",IN,sep=" ")
        ################
        out<-c(out,tes2)
        out<-c(out,IK1)
        out<-c(out,FINCH)
        
      }
    ####################################################################
    }else if(!sjmisc::is_empty(SM)){
      ##############################################################
      tes<-tryCatch({classyfireR::submit_query(label = 'query_test', input = SM, type = 'STRUCTURE')},error=function(cond){message("Classyfire is empty")})
      tes1<-do.call(paste, c(as.list(tryCatch({tes@classification$Classification},error=function(cond){message("Classifier could not fecth the information")})), sep = ","))
      ###########################################################
      ###########################################################
      IK1<-paste("INCHIKEY:",IK,sep=" ")
      tes2<-paste("Ontology:",tes1,sep=" ")
      FINCH<-paste("INCHI:",IN,sep=" ")
      #########################
      out<-c(out,tes2)
      out<-c(out,IK1)
      out<-c(out,FINCH)
      ########################
    }else{
      #############################################
      F1ONT<-paste("Ontology:","",sep=" ")
      ##FINCH<-paste("INCHI:",SM,sep=" ")
      ##############################################
      if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){
        FINCH<-paste("INCHI:",as.character(InMEDA[["InChI"]]),sep=" ")
        IK<-paste("INCHIKEY:","",sep=" ")
        ##################
        out<-c(out,F1ONT)
        out<-c(out,IK)
        out<-c(out,FINCH)
        #################
      }else{
        ###################################
        FINCH<-paste("INCHI:","",sep=" ")
        IK<-paste("INCHIKEY:",as.character(InMEDA[["InChI"]]),sep=" ")
        ###############
        out<-c(out,F1ONT)
        out<-c(out,IK)
        out<-c(out,FINCH)
        ####################
      }
      #################
    }## end of else..else if ...if
  }else if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]]))){
    ###########################
    print("enter the line ...53")
    print(as.character(InMEDA[["InChI"]]))
    ############################
    if(tryCatch({webchem::is.inchikey(as.character(InMEDA[["InChI"]]))},error=function(cond){message("inchikey validation failed")}))
    {
      #########################################
      print("entering the inchikey area")
      ########################################
      tes<-tryCatch({webchem::get_cid(stringr::str_trim(as.character(InMEDA[["InChI"]])), from = "inchikey")},error=function(cond){message("webchem not able to get cid from Inchikey")})
      tes1<-tryCatch({tes$cid},error=function(cond){message("Inchikey to CID did not convert")})
      tes2<-tryCatch({webchem::pc_prop(as.numeric(tes1[1]), properties = c("MolecularFormula", "MolecularWeight","CanonicalSMILES","InChI","InChIKey"))},error=function(cond){message("Inchikey to CID did not convert so did not get properties")})
      IN<-tryCatch({tes2$InChI},error=function(cond){message("Inchi to Inchikey failed because of CID not converting")})
      SM<-tryCatch({tes2$CanonicalSMILES},error=function(cond){message("Inchi to Inchikey failed because of CID not converting")})
      #######################################
      print("check if this is the error area")
      print(IN)
      #######################################
      FINCH<-paste("INCHI:",IN,sep=" ")
      IK<-as.character(InMEDA[["InChI"]])
      IK1<-paste("INCHIKEY:",IK,sep=" ")
      ############################
      if(!sjmisc::is_empty(IK)){
        ########################
        print("enter the line 8")
        ##########################
        IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(IK)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
        if(!sjmisc::is_empty(IKCRV)){
          
          ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
          ##########################
          IK1<-paste("INCHIKEY:",IK,sep=" ")
          tes2<-paste("Ontology:",ONTV,sep=" ")
          FINCH<-paste("INCHI:",IN,sep=" ")
          #########################
          out<-c(out,tes2)
          out<-c(out,IK1)
          out<-c(out,FINCH)
          ########################
        }else{
          IK1<-paste("INCHIKEY:",IK,sep=" ")
          tes2<-paste("Ontology:","",sep=" ")
          FINCH<-paste("INCHI:",IN,sep=" ")
          ################
          out<-c(out,tes2)
          out<-c(out,IK1)
          out<-c(out,FINCH)
          
        }## end of else
      }else if(!sjmisc::is_empty(SM)){
        tes<-tryCatch({classyfireR::submit_query(label = 'query_test', input = SM, type = 'STRUCTURE')},error=function(cond){message("Classyfire is empty")})
        tes1<-do.call(paste, c(as.list(tryCatch({tes@classification$Classification},error=function(cond){message("Classifier could not fecth the information")})), sep = ","))
        IK1<-paste("INCHIKEY:",IK,sep=" ")
        tes2<-paste("Ontology:",tes1,sep=" ")
        FINCH<-paste("INCHI:",IN,sep=" ")
        ########################
        out<-c(out,tes2)
        out<-c(out,IK1)
        out<-c(out,FINCH)
        
        
      }else{
        F1ONT<-paste("Ontology:","",sep=" ")
        ##FINCH<-paste("INCHI:",IN,sep=" ")
        ##############################################
        if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){
          FINCH<-paste("INCHI:",as.character(InMEDA[["InChI"]]),sep=" ")
          IK<-paste("INCHIKEY:","",sep=" ")
          ##################
          out<-c(out,F1ONT)
          out<-c(out,IK)
          out<-c(out,FINCH)
          #################
        }else{
          ###################################
          FINCH<-paste("INCHI:","",sep=" ")
          IK<-paste("INCHIKEY:",as.character(InMEDA[["InChI"]]),sep=" ")
          ###############
          out<-c(out,F1ONT)
          out<-c(out,IK)
          out<-c(out,FINCH)
          ####################
        }
        
      }## end of else ...else if ..else
    }else{
      ##inchikey validation failed ..so must be CAS or try to get ontology from smiles ...
      if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["InChI"]]))) & startsWith(as.character(InMEDA[["InChI"]]),'CAS:')){
        ###############################
        print("enter the CAS area ONtology")
        ##############################
        CV<-stringr::str_trim(as.character(InMEDA[["InChI"]]))
        CV1<-stringr::str_replace(CV,pattern='CAS:',replacement ="")
        CV2<-stringr::str_trim(as.character(CV1))
        ##############################
        PCID<-tryCatch({webchem::get_cid(CV2, from = "xref/rn",match="first")},error=function(cond){message("Pubchem Id is empty")})
        PCID1<-tryCatch({webchem::pc_prop(as.numeric(PCID[[2]][1]), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES","InChI","InChIKey"))},error=function(cond){message("Pubchem Id is empty")})
        ###############################
        gSMI<-tryCatch({PCID1$CanonicalSMILES},error=function(cond){message("smiles is not found")})
        IK<-tryCatch({PCID1$InChIKey},error=function(cond){message("some mistake happened in file search files")})
        IN<-tryCatch({PCID1$InChI},error=function(cond){message("some mistake happened in file search files")})
        IK1<-paste("INCHIKEY:",IK,sep=" ")
        #############################
        if(!sjmisc::is_empty(IK)){
          ##########################
          IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(IK)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
          if(!sjmisc::is_empty(IKCRV)){
            
            ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
            ##########################
            IK1<-paste("INCHIKEY:",IK,sep=" ")
            tes2<-paste("Ontology:",ONTV,sep=" ")
            FINCH<-paste("INCHI:",IN,sep=" ")
            #################
            out<-c(out,tes2)
            out<-c(out,IK1)
            out<-c(out,FINCH)
            ########################
          }## ikcrv END
        }else{
          if(!sjmisc::is_empty(as.character(InMEDA[["SMILES"]])) & !startsWith(as.character(InMEDA[["SMILES"]]),'not available')){
            ##IK<-as.character(InMEDA[["SMILES"]])
            F1SM<-stringr::str_trim(as.character(InMEDA[["SMILES"]]))
            F1SM1<-tryCatch({rinchi::get.inchi.key(F1SM)},error=function(cond){message("webchecm could not fetch the info")})
            #########
            if(!sjmisc::is_empty(F1SM1)){
              IK<-F1SM1
              IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(IK)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
              if(!sjmisc::is_empty(IKCRV)){
                ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
                ##########################
                IK1<-paste("INCHIKEY:",IK,sep=" ")
                tes2<-paste("Ontology:",ONTV,sep=" ")
                FINCH<-paste("INCHI:",IN,sep=" ")
                #################
                out<-c(out,tes2)
                out<-c(out,IK1)
                out<-c(out,FINCH)
              }##IKCRV
            }else{
              ## not able to get inchikey from smiles too check pubchemID
              if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(InMEDA[["PubChem CID"]]),'not available')){
                FPUCID<-stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))
                FPUCID1<-as.numeric(FPUCID)
                FINSM<-tryCatch({webchem::pc_prop(FPUCID1)},error=function(cond){message("webchecm could not fetch the info")})
                FIINK<-tryCatch({FINSM$InChIKey},error=function(cond){message("webchecm could not fetch the info")})
                #############################
                IK<-FIINK
                #############################
                if(!sjmisc::is_empty(IK)){
                  ##########################
                  IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(IK)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
                  if(!sjmisc::is_empty(IKCRV)){
                    
                    ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
                    ##########################
                    IK1<-paste("INCHIKEY:",IK,sep=" ")
                    tes2<-paste("Ontology:",ONTV,sep=" ")
                    FINCH<-paste("INCHI:",IN,sep=" ")
                    #################
                    out<-c(out,tes2)
                    out<-c(out,IK1)
                    out<-c(out,FINCH)
                    ########################
                  }## ikcrv END
                }else{
                  F1ONT<-paste("Ontology:","",sep=" ")
                  ##FINCH<-paste("INCHI:",SM,sep=" ")
                  ##############################################
                  if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){
                    FINCH<-paste("INCHI:",as.character(InMEDA[["InChI"]]),sep=" ")
                    IK<-paste("INCHIKEY:","",sep=" ")
                    ##################
                    out<-c(out,F1ONT)
                    out<-c(out,IK)
                    out<-c(out,FINCH)
                    #################
                  }else{
                    ###################################
                    FINCH<-paste("INCHI:","",sep=" ")
                    IK<-paste("INCHIKEY:",as.character(InMEDA[["InChI"]]),sep=" ")
                    ###############
                    out<-c(out,F1ONT)
                    out<-c(out,IK)
                    out<-c(out,FINCH)
                    ####################
                  }
                  #################
                }## end of else
                
                
              }else{
                F1ONT<-paste("Ontology:","",sep=" ")
                ##FINCH<-paste("INCHI:",SM,sep=" ")
                ##############################################
                if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){
                  FINCH<-paste("INCHI:",as.character(InMEDA[["InChI"]]),sep=" ")
                  IK<-paste("INCHIKEY:","",sep=" ")
                  ##################
                  out<-c(out,F1ONT)
                  out<-c(out,IK)
                  out<-c(out,FINCH)
                  #################
                }else{
                  ###################################
                  FINCH<-paste("INCHI:","",sep=" ")
                  IK<-paste("INCHIKEY:",as.character(InMEDA[["InChI"]]),sep=" ")
                  ###############
                  out<-c(out,F1ONT)
                  out<-c(out,IK)
                  out<-c(out,FINCH)
                  ####################
                }
              }## end of else## CID
              
              
            }## check else..smiles
          }else{
            #############################################
            F1ONT<-paste("Ontology:","",sep=" ")
            ##FINCH<-paste("INCHI:",SM,sep=" ")
            ##############################################
            if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){
              FINCH<-paste("INCHI:",as.character(InMEDA[["InChI"]]),sep=" ")
              IK<-paste("INCHIKEY:","",sep=" ")
              ##################
              out<-c(out,F1ONT)
              out<-c(out,IK)
              out<-c(out,FINCH)
              #################
            }else{
              ###################################
              FINCH<-paste("INCHI:","",sep=" ")
              IK<-paste("INCHIKEY:",as.character(InMEDA[["InChI"]]),sep=" ")
              ###############
              out<-c(out,F1ONT)
              out<-c(out,IK)
              out<-c(out,FINCH)
              ####################
            }
            
          }## inner smiles ..end ..else
        }## end of smiles## CAS
        
      }else{
        #############################################
        F1ONT<-paste("Ontology:","",sep=" ")
        ##FINCH<-paste("INCHI:",SM,sep=" ")
        ##############################################
        if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){
          FINCH<-paste("INCHI:",as.character(InMEDA[["InChI"]]),sep=" ")
          IK<-paste("INCHIKEY:","",sep=" ")
          ##################
          out<-c(out,F1ONT)
          out<-c(out,IK)
          out<-c(out,FINCH)
          #################
        }else{
          ###################################
          FINCH<-paste("INCHI:","",sep=" ")
          IK<-paste("INCHIKEY:",as.character(InMEDA[["InChI"]]),sep=" ")
          ###############
          out<-c(out,F1ONT)
          out<-c(out,IK)
          out<-c(out,FINCH)
          ####################
        }
        
      }####
      ##################################    
    }### end of else ..so starting checking CAS..smiles ..so ..on 
  }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["InChI"]]))) & startsWith(as.character(InMEDA[["InChI"]]),'CAS:')){
    ###############################
    print("enter the CAS area ONtology")
    ##############################
    CV<-stringr::str_trim(as.character(InMEDA[["InChI"]]))
    CV1<-stringr::str_replace(CV,pattern='CAS:',replacement ="")
    CV2<-stringr::str_trim(as.character(CV1))
    ##############################
    PCID<-tryCatch({webchem::get_cid(CV2, from = "xref/rn",match="first")},error=function(cond){message("Pubchem Id is empty")})
    PCID1<-tryCatch({webchem::pc_prop(as.numeric(PCID[[2]][1]), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES","InChI","InChIKey"))},error=function(cond){message("Pubchem Id is empty")})
    ###############################
    gSMI<-tryCatch({PCID1$CanonicalSMILES},error=function(cond){message("smiles is not found")})
    IK<-tryCatch({PCID1$InChIKey},error=function(cond){message("some mistake happened in file search files")})
    IN<-tryCatch({PCID1$InChI},error=function(cond){message("some mistake happened in file search files")})
    IK1<-paste("INCHIKEY:",IK,sep=" ")
    #############################
    if(!sjmisc::is_empty(IK)){
      ##########################
      IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(IK)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
      if(!sjmisc::is_empty(IKCRV)){
        
        ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
        ##########################
        IK1<-paste("INCHIKEY:",IK,sep=" ")
        tes2<-paste("Ontology:",ONTV,sep=" ")
        FINCH<-paste("INCHI:",IN,sep=" ")
        #################
        out<-c(out,tes2)
        out<-c(out,IK1)
        out<-c(out,FINCH)
        ########################
      }## ikcrv END
    }else{
      if(!sjmisc::is_empty(as.character(InMEDA[["SMILES"]])) & !startsWith(as.character(InMEDA[["SMILES"]]),'not available')){
        IK<-as.character(InMEDA[["SMILES"]])
        F1SM<-stringr::str_trim(as.character(InMEDA[["SMILES"]]))
        F1SM1<-tryCatch({rinchi::get.inchi.key(F1SM)},error=function(cond){message("webchecm could not fetch the info")})
        #########
        if(!sjmisc::is_empty(F1SM1)){
          IK<-F1SM1
          IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(IK)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
          if(!sjmisc::is_empty(IKCRV)){
            ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
            ##########################
            IK1<-paste("INCHIKEY:",IK,sep=" ")
            tes2<-paste("Ontology:",ONTV,sep=" ")
            FINCH<-paste("INCHI:",IN,sep=" ")
            #################
            out<-c(out,tes2)
            out<-c(out,IK1)
            out<-c(out,FINCH)
	    ####################
          }##IKCRV
        }else{
          ## not able to get inchikey from smiles too check pubchemID
          if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(InMEDA[["PubChem CID"]]),'not available')){
            FPUCID<-stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))
            FPUCID1<-as.numeric(FPUCID)
            FINSM<-tryCatch({webchem::pc_prop(FPUCID1)},error=function(cond){message("webchecm could not fetch the info")})
            FIINK<-tryCatch({FINSM$InChIKey},error=function(cond){message("webchecm could not fetch the info")})
            #############################
            IK<-FIINK
            #############################
            if(!sjmisc::is_empty(IK)){
              ##########################
              IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(IK)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
              if(!sjmisc::is_empty(IKCRV)){
                
                ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
                ##########################
                IK1<-paste("INCHIKEY:",IK,sep=" ")
                tes2<-paste("Ontology:",ONTV,sep=" ")
                FINCH<-paste("INCHI:",IN,sep=" ")
                #################
                out<-c(out,tes2)
                out<-c(out,IK1)
                out<-c(out,FINCH)
                ########################
              }## ikcrv END
            }else{
              F1ONT<-paste("Ontology:","",sep=" ")
              ##FINCH<-paste("INCHI:",SM,sep=" ")
              ##############################################
              if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){
                FINCH<-paste("INCHI:",as.character(InMEDA[["InChI"]]),sep=" ")
                IK<-paste("INCHIKEY:","",sep=" ")
                ##################
                out<-c(out,F1ONT)
                out<-c(out,IK)
                out<-c(out,FINCH)
                #################
              }else{
                ###################################
                FINCH<-paste("INCHI:","",sep=" ")
                IK<-paste("INCHIKEY:",as.character(InMEDA[["InChI"]]),sep=" ")
                ###############
                out<-c(out,F1ONT)
                out<-c(out,IK)
                out<-c(out,FINCH)
                ####################
              }
              #################
            }## end of else
            
            
          }else{
            F1ONT<-paste("Ontology:","",sep=" ")
            ##FINCH<-paste("INCHI:",SM,sep=" ")
            ##############################################
            if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){
              FINCH<-paste("INCHI:",as.character(InMEDA[["InChI"]]),sep=" ")
              IK<-paste("INCHIKEY:","",sep=" ")
              ##################
              out<-c(out,F1ONT)
              out<-c(out,IK)
              out<-c(out,FINCH)
              #################
            }else{
              ###################################
              FINCH<-paste("INCHI:","",sep=" ")
              IK<-paste("INCHIKEY:",as.character(InMEDA[["InChI"]]),sep=" ")
              ###############
              out<-c(out,F1ONT)
              out<-c(out,IK)
              out<-c(out,FINCH)
              ####################
            }
          }## end of else## CID
          
          
        }## check else..smiles
      }else{
        #############################################
        F1ONT<-paste("Ontology:","",sep=" ")
        ##FINCH<-paste("INCHI:",SM,sep=" ")
        ##############################################
        if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){
          FINCH<-paste("INCHI:",as.character(InMEDA[["InChI"]]),sep=" ")
          IK<-paste("INCHIKEY:","",sep=" ")
          ##################
          out<-c(out,F1ONT)
          out<-c(out,IK)
          out<-c(out,FINCH)
          #################
        }else{
          ###################################
          FINCH<-paste("INCHI:","",sep=" ")
          IK<-paste("INCHIKEY:",as.character(InMEDA[["InChI"]]),sep=" ")
          ###############
          out<-c(out,F1ONT)
          out<-c(out,IK)
          out<-c(out,FINCH)
          ####################
        }
        
      }## inner smiles ..end ..else
    }## end of smiles## CAS
    
  }else{
    ##########################################
    print("entering this line...174")
    #################################
    PCID<- as.character(InMEDA[["PubChem CID"]])
    PCSM<- as.character(InMEDA[["SMILES"]])
    #############################
    if(!sjmisc::is_empty(PCID))
    {
      #######################################
      PCID1<-as.numeric(PCID)
      ######################################
      tes<-tryCatch({webchem::pc_prop(PCID1, properties = c("MolecularFormula", "MolecularWeight","CanonicalSMILES","InChI","InChIKey"))},error=function(cond){message("Did not get properties from Pubchem CID")})
      gSMI<-tryCatch({tes$CanonicalSMILES},error=function(cond){message("smiles is not found")})
      IK<-tryCatch({tes$InChIKey},error=function(cond){message("some mistake happened in file search files")})
      IN<-tryCatch({tes$InChI},error=function(cond){message("some mistake happened in file search files")})
      IK1<-paste("INCHIKEY:",IK,sep=" ")
      #############################
      if(!sjmisc::is_empty(IK)){
        ###############################
        IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(IK)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
        if(!sjmisc::is_empty(IKCRV)){
          ####################
          ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
          ##################################
          IK1<-paste("INCHIKEY:",IK,sep=" ")
          tes2<-paste("Ontology:",ONTV,sep=" ")
          FINCH<-paste("INCHI:",IN,sep=" ")
          #################################
          out<-c(out,tes2)
          out<-c(out,IK1)
          out<-c(out,FINCH)
        }else{
          IK1<-paste("INCHIKEY:",IK,sep=" ")
          tes2<-paste("Ontology:","",sep=" ")
          FINCH<-paste("INCHI:",IN,sep=" ")
          out<-c(out,tes2)
          out<-c(out,IK1)
          out<-c(out,FINCH)
        }###else
      }else if(!sjmisc::is_empty(gSMI)){
        ##################################################
        print("enter the line 184")
        ##################################################
        tes<-tryCatch({classyfireR::submit_query(label = 'query_test', input = gSMI, type = 'STRUCTURE')},error=function(cond){message("adduct value is missing")})
        tes1<-do.call(paste, c(as.list(tryCatch({tes@classification$Classification},error=function(cond){message("Classifier could not fecth the information")})), sep = ","))
        ###############################################
        tes2<-paste("Ontology:",tes1,sep=" ")
        FINCH<-paste("INCHI:",IN,sep=" ")
        ##################
        out<-c(out,tes2)
        out<-c(out,IK1)
        out<-c(out,FINCH)
        #################
      }else{
        F1ONT<-paste("Ontology:","",sep=" ")
        if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){
          #########################################
          FINCH<-paste("INCHI:",as.character(InMEDA[["InChI"]]),sep=" ")
          IK<-paste("INCHIKEY:","",sep=" ")
          #################
          out<-c(out,F1ONT)
          out<-c(out,IK)
          out<-c(out,FINCH)
          ################
        }else{
          FINCH<-paste("INCHI:","",sep=" ")
          IK<-paste("INCHIKEY:",as.character(InMEDA[["InChI"]]),sep=" ")
          #################
          out<-c(out,F1ONT)
          out<-c(out,IK)
          out<-c(out,FINCH)
          ####################
        }##else..if
        
      }###else...else if ..if
    }else{
      ######check the smiles is not empty
      ###PCSM<-SMV
      ################################	    
      SMV<-PCSM
      if(!sjmisc::is_empty(SMV)){
        IK<-tryCatch({rinchi::get.inchi.key(SMV)},error=function(cond){message("rinchi could not fetch inchikey missing")})
        SMV1<-tryCatch({rinchi::get.inchi(SMV)},error=function(cond){message("rinchi could not fetch inchi missing")})
        IK1<-paste("INCHIKEY:",IK,sep=" ")
        FINCH<-paste("INCHI:",SMV1,sep=" ")
        #########################
        if(!sjmisc::is_empty(IK)){
          IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(IK)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
          if(!sjmisc::is_empty(IKCRV)){
            ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
            ###########################
            IK1<-paste("INCHIKEY:",IK,sep=" ")
            tes2<-paste("Ontology:",ONTV,sep=" ")
            ###########################
            out<-c(out,tes2)
            out<-c(out,IK1)
            out<-c(out,FINCH)
            
          }else{
            ### IKCRV is empty
            IK1<-paste("INCHIKEY:",IK,sep=" ")
            tes2<-paste("Ontology:","",sep=" ")
            FINCH<-paste("INCHI:",SMV1,sep=" ")
            ################
            out<-c(out,tes2)
            out<-c(out,IK1)
            out<-c(out,FINCH)
	    #################
          }## end of IKCRV
          
        }else if(!sjmisc::is_empty(SMV)){
	  ################################################	
          tes<-tryCatch({classyfireR::submit_query(label = 'query_test', input = SMV, type = 'STRUCTURE')},error=function(cond){message("Classfire not able to fetch empty")})
          tes1<-do.call(paste, c(as.list(tryCatch({tes@classification$Classification},error=function(cond){message("Classifier could not fecth the information")})), sep = ","))
          tes2<-paste("Ontology:",tes1,sep=" ")
          teIK2<-SMV1
          ##teIK2<-tryCatch({PuInKtoSM(IK)},error=function(cond){message("PUBCHEM not able to convert inchikey to smile")})
          ##teIK2<-tryCatch({webchem::cs_convert(IK,from="inchikey",to="inchi")},error=function(cond){message("webchecm could not fetch the info")})
          FINCH<-paste("INCHI:",teIK2,sep=" ")
          ########################
          out<-c(out,tes2)
          out<-c(out,IK1)
          out<-c(out,FINCH)
	  ##########################
        }else{
          ##############################################
          F1ONT<-paste("Ontology:","",sep=" ")
          if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){
            FINCH<-paste("INCHI:",as.character(InMEDA[["InChI"]]),sep=" ")
            IK<-paste("INCHIKEY:","",sep=" ")
            ###############
            out<-c(out,F1ONT)
            out<-c(out,IK)
            out<-c(out,FINCH)
            ###############
          }else{
            FINCH<-paste("INCHI:","",sep=" ")
            IK<-paste("INCHIKEY:",as.character(InMEDA[["InChI"]]),sep=" ")
            ###############
            out<-c(out,F1ONT)
            out<-c(out,IK)
            out<-c(out,FINCH)
            #################################
          }
          
        }## end of else if ...else ...if
      #########
      }
      ##########
    }## entering the else
    
  }## end of else
##########################
  return(out)
###########################

}## end of ontology function


##########################################################################################
##########################################################################################
finADWC<-function(gs,dIV)
{
  gs1<-unlist(strsplit(gs, "(?=[+-])", perl = TRUE))
  tes1<-c()
  for(i in c(1:length(gs1)))
  {

    if(gs1[i]=="+"){
      tes1<-c(tes1,gs1[i])
    }else if(gs1[i]=="-"){
      tes1<-c(tes1,gs1[i])
    }else{
      cs<-gsub("\\[", "", gs1[i])
      cs1<-gsub("\\]", "",cs)
      ####################################
      ##nt<-paste(cs1,"/",dIV,sep="")
      ##if(grepl("^[[:digit:]]+",cs1)){
      ####################################
      if(grepl("^[[:digit:]]+",cs1) & !grepl("^[0-9]*[.])?[0-9]+",cs1)){
        VL<-comwithNu(cs1)
        nt<-paste(VL,"/",dIV,sep="")
        tes1<-c(tes1,nt)
      }else{
        nt<-paste(cs1,"/",dIV,sep="")
        tes1<-c(tes1,nt)
      }
      ###################################
      ##nt<-paste(gs1[i],"/",dIV,sep="")
      ##tes1<-c(tes1,nt)
      ###################################
    }
  }
  return(paste(tes1,collapse = ""))
}


##########I###########################################################
##comwithNu<-function(tes)
##{
  ##tes1<-c()
  ##for(elem in tes){
    ##if(numbers_only(elem)){
      ##elem1=paste(elem,"*",sep="")
      ##tes1<-c(tes1,elem1)
    ##}else if(!sjmisc::is_empty(tryCatch({AD[AD$formula==elem,]$exactMass[1]},warning=function(cond){message("error in the database search info")}))){
     ## val= tryCatch({AD[AD$formula==elem,]$exactMass[1]},warning=function(cond){message("error in data base search")})
      ##tes1<-c(tes1,val)
    ##}else{
     ## x<-"Pass"
   ## }
 ## }
 ## return(paste(tes1,collapse = ""))
##}
########################################################################
########################################################################
comwithNu<-function(elem)
{
  tes1<-c()
  ################
  if(grepl("^[[:digit:]]+",elem) & !grepl("^[0-9]*[.])?[0-9]+",elem)){
    elem1<-stri_extract_first_regex(elem,"[0-9]+")
    elem2<-substring(elem,2)
    if(numbers_only(elem1)){
      nelem1=paste(elem1,"*",sep="")
      tes1<-c(tes1,nelem1)
    }else{
      tes1<-c(tes1,elem1)
    }
  ######################
    if(!sjmisc::is_empty(tryCatch({AD[AD$formula==elem2,]$exactMass[1]},warning=function(cond){message("error in the database search info")}))){
      ##print(elem2)
      val= tryCatch({AD[AD$formula==elem2,]$exactMass[1]},warning=function(cond){message("error in data base search")})
      tes1<-c(tes1,val)
      ##print(val)
    }else{
      tes1<-c(tes1,elem)
    }
  }## end of if
  ###############################
  return(paste(tes1,collapse = ""))
}
#########################################################################################

comNFMDA<-function(tes)
{
  tes1<-c()
  for(i in c(1:length(tes)))
  {

    if(tes[i]=="M")
    {
      tes1<-c(tes1,tes[i])
    }else if(tes[i]=="+"){
      tes1<-c(tes1,tes[i])
    }else if(tes[i]=="-"){
      tes1<-c(tes1,tes[i])
    }else if(tes[i]=="[M]"){
      tes1<-c(tes1,"M")
    }else{
      if(!sjmisc::is_empty(tryCatch({AD[AD$formula==tes[i],]$exactMass[1]},warning=function(cond){message("error in data base search")})))
      {
        val=tryCatch({AD[AD$formula==tes[i],]$exactMass[1]},warning=function(cond){message("error in data base search")})
        tes1<-c(tes1,val)
      }else if(grepl("^[[:digit:]]+", tes[i])){
        VL<-comwithNu(tes[i])
        tes1<-c(tes1,VL)

      }else{
        tes1<-c(tes1,tes[i])
      }

    }

  }
  tes2<-paste(tes1,collapse = "")
  if(tes2=="M+"|tes2=="M-"){
  return("M")
  }else{
    return(tes2)
  }
}

###################################################################################
FADINF<-function(addu1)
{
  faddu<-c()
cha<-stringr::word(addu1, 2, sep="]")
if(cha=="-")
{

  cha<-"1-"
  addu2<-tryCatch({qdapRegex::ex_between(addu1, "[", "]")[[1]]},warning=function(cond){message("error happened in precursortype extraction")})
  addu3<-tryCatch({AIN[AIN$V1==addu2 & AIN$V3==cha,]$V2},warning=function(cond){message("adduct match is empty")})
  if(!sjmisc::is_empty(addu3))
  {
    ##print(addu3)
    faddu<-c(faddu,addu3)

  }else{
    NAF<-unlist(strsplit(addu2, "(?=[+-])", perl = TRUE))
    NAF1<-comNFMDA(NAF)
    faddu<-c(faddu,NAF1)

  }

}else if(cha=="+"){
  cha<-"1+"
  addu2<-tryCatch({qdapRegex::ex_between(addu1, "[", "]")[[1]]},warning=function(cond){message("error happened in precursortype extraction")})
  addu3<-tryCatch({AIN[AIN$V1==addu2 & AIN$V3==cha,]$V2},warning=function(cond){message("adduct match is empty")})
  if(!sjmisc::is_empty(addu3))
  {
    ##print(addu1)
    faddu<-c(faddu,addu3)

  }else{
    NAF<-unlist(strsplit(addu2, "(?=[+-])", perl = TRUE))
    NAF1<-comNFMDA(NAF)
    faddu<-c(faddu,NAF1)
    ##print(NAF1)
  }

}else if(cha=="2-"){

  ##x<-"pass1"
  dIV<-"2"
  addu2<-tryCatch({qdapRegex::ex_between(addu1, "[", "]")[[1]]},warning=function(cond){message("error happened in precursortype extraction")})
  addu3<-tryCatch({AIN[AIN$V1==addu2 & AIN$V3==cha,]$V2},warning=function(cond){message("adduct match is empty")})
  if(!sjmisc::is_empty(addu3))
  {
    ##print(addu1)
    faddu<-c(faddu,addu3)

  }else{
    NAF<-unlist(strsplit(addu2, "(?=[+-])", perl = TRUE))
    NAF1<-comNFMDA(NAF)
    NAF2<-finADWC(NAF1,dIV)
    faddu<-c(faddu,NAF2)
    ##print(NAF2)
    ##print(NAF1)
  }

}else if(cha=="2+"){
  dIV<-"2"
  addu2<-tryCatch({qdapRegex::ex_between(addu1, "[", "]")[[1]]},warning=function(cond){message("error happened in precursortype extraction")})
  addu3<-tryCatch({AIN[AIN$V1==addu2 & AIN$V3==cha,]$V2},warning=function(cond){message("adduct match is empty")})
  if(!sjmisc::is_empty(addu3))
  {
    ##print(addu1)
    faddu<-c(faddu,addu3)

  }else{
    NAF<-unlist(strsplit(addu2, "(?=[+-])", perl = TRUE))
    NAF1<-comNFMDA(NAF)
    NAF2<-finADWC(NAF1,dIV)
    faddu<-c(faddu,NAF2)
    ##faddu<-c(faddu,NAF1)
    ##print(NAF1)
  }

}else if(cha=="3-"){
  ##x<-"pass1"
  dIV<-"3"
  addu2<-tryCatch({qdapRegex::ex_between(addu1, "[", "]")[[1]]},warning=function(cond){message("error happened in precursortype extraction")})
  addu3<-tryCatch({AIN[AIN$V1==addu2 & AIN$V3==cha,]$V2},warning=function(cond){message("adduct match is empty")})
  if(!sjmisc::is_empty(addu3))
  {
    ##print(addu1)
    faddu<-c(faddu,addu3)

  }else{
    NAF<-unlist(strsplit(addu2, "(?=[+-])", perl = TRUE))
    NAF1<-comNFMDA(NAF)
    NAF2<-finADWC(NAF1,dIV)
    faddu<-c(faddu,NAF2)
    ##faddu<-c(faddu,addu3)
    ##print(NAF1)
  }
}else if(cha=="3+"){
  dIV<-"3"
  addu2<-tryCatch({qdapRegex::ex_between(addu1, "[", "]")[[1]]},warning=function(cond){message("error happened in precursortype extraction")})
  addu3<-tryCatch({AIN[AIN$V1==addu2 & AIN$V3==cha,]$V2},warning=function(cond){message("adduct match is empty")})
  if(!sjmisc::is_empty(addu3))
  {
    ##print(addu1)
    faddu<-c(faddu,addu3)

  }else{
    NAF<-unlist(strsplit(addu2, "(?=[+-])", perl = TRUE))
    NAF1<-comNFMDA(NAF)
    NAF2<-finADWC(NAF1,dIV)
    faddu<-c(faddu,NAF2)

  }
}else{
  x<-"pass"

}
return(faddu)
}

#############################################################################################
#############################################################################################
#############################################################################################
MaKlist<-function(gFile)
{
  #############################################
  lines <- readLines(gFile)
  lst <-split(lines, cumsum(lines==""))
  lst1 <-lapply(lst, function(x) if (x[1] == "") x[-1] else x)
  LL<-sapply(lst, length)
  IR<-which(unname(LL) == 0)
  if(length(IR)>0){LL1<-LL[-(IR)]}else{LL1 <- LL}
  ###### MSP List#############################
  lst2 <-lst1[names(LL1)]
  lst2<-purrr::compact(lst2)
  #############################
  PMZL<-unname(sapply(lst2, function(x) grep("PRECURSORMZ",x)))
  PMZL1<-Filter(length,PMZL)
  ##### making addut information
  pattern <- "PRECURSORTYPE|ADDUCTIONNAME"
  PTYL<-unname(sapply(lst2, function(x) grep(pattern,x)))
  PTYL1<-Filter(length,PTYL)
  ##############################################
  ###### Adduct list ###########################
  ###########################################
  ###########################################
  faddu<-c()
  if(length(PTYL1)>0){
    for(i in c(1:length(PTYL1))){
      addu<-grep("PRECURSORTYPE:|ADDUCTIONNAME:",lst2[[i]])
      addu1<-stringr::str_trim(gsub('PRECURSORTYPE:|ADDUCTIONNAME:','',lst2[[i]][addu]))
      #####################
      cha<-stringr::word(addu1, 2, sep="]")
      if(cha=="-")
      {

        cha<-"1-"
        addu2<-tryCatch({qdapRegex::ex_between(addu1, "[", "]")[[1]]},warning=function(cond){message("error happened in precursortype extraction")})
        addu3<-tryCatch({AIN[AIN$V1==addu2 & AIN$V3==cha,]$V2},warning=function(cond){message("adduct match is empty")})
        if(!sjmisc::is_empty(addu3))
        {
          ##print(addu3)
          faddu<-c(faddu,addu3)

        }else{
          NAF<-unlist(strsplit(addu2, "(?=[+-])", perl = TRUE))
          NAF1<-comNFMDA(NAF)
          faddu<-c(faddu,NAF1)

        }

      }else if(cha=="+"){
        cha<-"1+"
        addu2<-tryCatch({qdapRegex::ex_between(addu1, "[", "]")[[1]]},warning=function(cond){message("error happened in precursortype extraction")})
        addu3<-tryCatch({AIN[AIN$V1==addu2 & AIN$V3==cha,]$V2},warning=function(cond){message("adduct match is empty")})
        if(!sjmisc::is_empty(addu3))
        {
          ##print(addu1)
          faddu<-c(faddu,addu3)

        }else{
          NAF<-unlist(strsplit(addu2, "(?=[+-])", perl = TRUE))
          NAF1<-comNFMDA(NAF)
          faddu<-c(faddu,NAF1)
          ##print(NAF1)
        }

      }else if(cha=="2-"){

        ##x<-"pass1"
        dIV<-"2"
        addu2<-tryCatch({qdapRegex::ex_between(addu1, "[", "]")[[1]]},warning=function(cond){message("error happened in precursortype extraction")})
        addu3<-tryCatch({AIN[AIN$V1==addu2 & AIN$V3==cha,]$V2},warning=function(cond){message("adduct match is empty")})
        if(!sjmisc::is_empty(addu3))
        {
          ##print(addu1)
          faddu<-c(faddu,addu3)

        }else{
          NAF<-unlist(strsplit(addu2, "(?=[+-])", perl = TRUE))
          NAF1<-comNFMDA(NAF)
          NAF2<-finADWC(NAF1,dIV)
          faddu<-c(faddu,NAF2)
          ##print(NAF2)
          ##print(NAF1)
        }

      }else if(cha=="2+"){
        dIV<-"2"
        addu2<-tryCatch({qdapRegex::ex_between(addu1, "[", "]")[[1]]},warning=function(cond){message("error happened in precursortype extraction")})
        addu3<-tryCatch({AIN[AIN$V1==addu2 & AIN$V3==cha,]$V2},warning=function(cond){message("adduct match is empty")})
        if(!sjmisc::is_empty(addu3))
        {
          ##print(addu1)
          faddu<-c(faddu,addu3)

        }else{
          NAF<-unlist(strsplit(addu2, "(?=[+-])", perl = TRUE))
          NAF1<-comNFMDA(NAF)
          NAF2<-finADWC(NAF1,dIV)
          faddu<-c(faddu,NAF2)
          ##faddu<-c(faddu,NAF1)
          ##print(NAF1)
        }

      }else if(cha=="3-"){
        ##x<-"pass1"
        dIV<-"3"
        addu2<-tryCatch({qdapRegex::ex_between(addu1, "[", "]")[[1]]},warning=function(cond){message("error happened in precursortype extraction")})
        addu3<-tryCatch({AIN[AIN$V1==addu2 & AIN$V3==cha,]$V2},warning=function(cond){message("adduct match is empty")})
        if(!sjmisc::is_empty(addu3))
        {
          ##print(addu1)
          faddu<-c(faddu,addu3)

        }else{
          NAF<-unlist(strsplit(addu2, "(?=[+-])", perl = TRUE))
          NAF1<-comNFMDA(NAF)
          NAF2<-finADWC(NAF1,dIV)
          faddu<-c(faddu,NAF2)
          ##faddu<-c(faddu,addu3)
          ##print(NAF1)
        }
      }else if(cha=="3+"){
        dIV<-"3"
        addu2<-tryCatch({qdapRegex::ex_between(addu1, "[", "]")[[1]]},warning=function(cond){message("error happened in precursortype extraction")})
        addu3<-tryCatch({AIN[AIN$V1==addu2 & AIN$V3==cha,]$V2},warning=function(cond){message("adduct match is empty")})
        if(!sjmisc::is_empty(addu3))
        {
          ##print(addu1)
          faddu<-c(faddu,addu3)

        }else{
          NAF<-unlist(strsplit(addu2, "(?=[+-])", perl = TRUE))
          NAF1<-comNFMDA(NAF)
          NAF2<-finADWC(NAF1,dIV)
          faddu<-c(faddu,NAF2)

        }
      }else{
        x<-"pass"

      }

    }
  }else(print("AIn is empty"))
  #### The list with added adduct value ########################################################
  lst3<-tryCatch({base::Map(c, lst2, faddu)},error=function(cond){message("faddu is empty")})
  ######### Precursor MZ information ###########################################################
  lst4<-c()
  for(i in 1:length(lst3)){
    LV<-tryCatch({lst3[[i]]},warning=function(cond){message("list is empty")})
    LaV<-tryCatch({BBmisc::getLast(LV)},warning=function(cond){message("last element is empty")})
    SR<-tryCatch({stringi::stri_startswith_fixed(LV, 'PRECURSORMZ:')},warning=function(cond){message("PRECURSORMZ is empty")})
    SR1<-tryCatch({which(SR)},warning=function(cond){message("index is empty")})
    SR2<-tryCatch({LV[[SR1]]},error=function(cond){message("index does not exist")})
    PV<-tryCatch({stringr::str_trim(stringr::str_replace(SR2, "PRECURSORMZ:", ""))},warning=function(cond){message("index does not exist")})
    S1r<-tryCatch({stringr::str_replace(LaV, "M", PV)},warning=function(cond){message("replacement does not exist")})
    S2r<-tryCatch({as.numeric(pander::evals(S1r)[[1]]$result)},error=function(cond){message("result is empty")})
    lst4<-c(lst4,S2r)
  }
  ####################################################################################################
  #### adding the adduct adjusted mass to make final list
  lst5<-tryCatch({base::Map(c, lst3, lst4)},error=function(cond){message("Either lst3 or lst4 is empty")})
  ######## Precursormz value #########################################################################
  fmass<-c()
  for(i in c(1:length(PMZL1))){
    LV1<-tryCatch({lst3[[i]]},warning=function(cond){message("list is empty")})
    S1R<-tryCatch({stringi::stri_startswith_fixed(LV1, 'PRECURSORMZ:') },warning=function(cond){message("PRECURSORMZ is empty")})
    S1R1<-tryCatch({which(S1R)},warning=function(cond){message("Index is empty")})
    S1R2<-tryCatch({LV1[S1R1]},error=function(cond){message("Element is empty")})
    PV1<-tryCatch({stringr::str_trim(stringr::str_replace(S1R2, "PRECURSORMZ:", ""))},warning=function(cond){message("PRECURSORMZ empty")})
    PV2<-tryCatch({as.numeric(PV1)},warning=function(cond){message("PRECURSORMZ does not exists")})
    fmass<-c(fmass,PV2)
  }
  #######################################################################################################
  lst7<-tryCatch({base::Map(c, lst5, fmass)},error=function(cond){message("Either lst5 or fmass is empty")})
  ######### Coping to new list
  lst8<-lst7
  ##### Adduct adjusted mass values #####################################################################
  AAMZV<-tryCatch({sapply(lst5, function(x) BBmisc::getLast(x))},error=function(cond){message("Either lst5 is empty")})
  AAMZV <- tryCatch({setNames(as.numeric(AAMZV), names(AAMZV))},error=function(cond){message("last numeric is empty")})
  ###############################
  ###### Retention time list#####
  ### Making the RT list#################################################################################
  RTL<-tryCatch({sapply(lst2, function(x) grep("RETENTIONTIME",x))},error=function(cond){message("retention time is empty")})
  IND<-tryCatch({unname(RTL)},error=function(cond){message("retention time is empty")})
  ############# Final Retention Test list ###############################################################
  FRTL <- c()
  if(length(IND)>0){
    for(i in c(1:length(IND))){
      FRTL <- c(FRTL,get(names(lst2)[i],lst2)[IND[i]])
    }}
  ########################################################################################################
  FRTL1<-tryCatch({as.numeric(sapply(strsplit(FRTL, ":"),`[`, 2))},error = function(cond){message("out is empty")})
  #######################################################################################################
  return(list(lst2, PMZL1, PTYL1,faddu,lst3,lst4,lst5,fmass,lst7,lst8,AAMZV,RTL,IND,FRTL,FRTL1))

}
##########################################################################################################
##########################################################################################################
NFFilter1<-function(InMEDA,InAdVA,InMSPL,InPMZ,InRTL)
{
  ################################
  FMa<-c()
  ################################
  if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & !startsWith(as.character(InMEDA[["InChI"]]),'not available') & !startsWith(as.character(InMEDA[["InChI"]]),'CAS:') & !startsWith(as.character(InMEDA[["InChI"]]),'InChI=')){
    if(tryCatch({webchem::is.inchikey(stringr::str_trim(as.character(InMEDA[["InChI"]])))},error=function(cond){message("inchikey..file must be empty")})){
      ########################################################################
      print("enter the Inchi key AREA..inchikey is not empty")
      ########################################################################
      IK<-as.character(InMEDA[["InChI"]])
      IK1<-tryCatch({get_cid(IK, from = "inchikey")},error=function(cond){message("Inchi name must be empty or rinchi not abe to fetch")})
      PCID<-tryCatch({IK1[[2]]},error=function(cond){message("Pubchem Id is empty")})
      ########################################
      PCID1<-tryCatch({webchem::pc_prop(as.numeric(PCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem Id is empty")})
      PCID2<-tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){message("Pubchem Id is empty")})
      ##########################################################################
      ##########################################################################
      if(!sjmisc::is_empty(PCID2))
      {
        FMa<-c(FMa,PCID2)
      }else{
	###########################################################
        if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),'InChI=')){
	  ##############################################
          IK1<-tryCatch({get_cid(IK, from = "inchi")},error=function(cond){message("Inchi name must be empty or rinchi not abe to fetch")})
          PCID<-tryCatch({IK1[[2]]},error=function(cond){message("Pubchem Id is empty")})
	  #######################################
	  #######################################
	  PCID1<-tryCatch({webchem::pc_prop(as.numeric(PCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem Id is empty")})
	  PCID2<-tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){message("Pubchem Id is empty")})
	  ###################################
          if(!sjmisc::is_empty(PCID2)){
            FMa<-c(FMa,PCID2)
          }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["InChI"]]))) & startsWith(as.character(InMEDA[["InChI"]]),'CAS:')){
	    ##############################################
            CV<-stringr::str_trim(as.character(InMEDA[["InChI"]]))
            CV1<-stringr::str_replace(CV,pattern='CAS:',replacement ="")
            CV2<-stringr::str_trim(as.character(CV1))
	    #############################################
            PCID<-tryCatch({get_cid(CV2, from = "xref/rn",match="first")},error=function(cond){message("Pubchem Id is empty")})
	    #########################################
	    ########################################
	    PCID1<-tryCatch({webchem::pc_prop(as.numeric(PCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem Id is empty")})
	    PCID2<-tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){message("Pubchem Id is empty")})
	    #############################################
            if(!sjmisc::is_empty(PCID2)){
              FMa<-c(FMa,PCID2)
            }else if(!sjmisc::is_empty(as.character(InMEDA[["SMILES"]])) & !startsWith(as.character(InMEDA[["SMILES"]]),'not available')){
	      ######################################
	      print("entering smiles area in InchiKey")
	       ##################################
              IK<-as.character(InMEDA[["SMILES"]])
	      ######################################
	      tes<-tryCatch({rcdk::parse.smiles(IK)},error=function(cond){message("smiles not abe to parse")})
	      tes1<-tryCatch({tes[[1]]},error=function(cond){message("smiles parse information is empty")})
	      tes2<-tryCatch({rcdk::get.exact.mass(tes1)},error=function(cond){message("smiles not abe to fetch")})
	      PCID2<-tryCatch({tes2},error=function(cond){message("smiles not abe to fetch")})
	      #######################################
              #######################################
              if(!sjmisc::is_empty(PCID2)){
                FMa<-c(FMa,PCID2)
              }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(InMEDA[["PubChem CID"]]),'not available')){
                FPUCID<-stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))
	        ####################################
	        PCID1<-tryCatch({webchem::pc_prop(as.numeric(FPUCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem CId is empty..did not get exact mass")})
		PCID2<-tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){message("Pubchem CID is empty..did not get exact mass")})
		######################################
                ######################################
                if(!sjmisc::is_empty(PCID2)){
                  FMa<-c(FMa,PCID2)
                }else if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){
		 ########################################################
		 EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem Id is empty")})
		 EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
		 #########################################################
		 #########################################################
                  if(!sjmisc::is_empty(EM1)){
                    FMa<-c(FMa,EM1)
                  }else{
                    FMa<-c(FMa,PCID2)
                  }
               ###########################
                }else{
                  ## formula not found and --smiles and pubchem failed to get PubchemID
                  FMa<-c(FMa,PCID2)
                }
              }else{
                ### PubchemId is not found ..entering the else loop
                if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){
		  #########################################
		  EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem Id is empty")})
	          EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
		  ##########################################
	          ############################################
                  if(!sjmisc::is_empty(EM1)){
                    ##print("formula area.....")
                    ##print(EM1)
                    FMa<-c(FMa,EM1)
                  }else{
                    ##print("formula area empty.....")
                    ##print(PCID2)
                    FMa<-c(FMa,PCID2)
                  }
                }
                ###############################
              }
              ####################  PCID2 ...end ###
            }else{
              ################################
              ## smiles end ... not found
              ################################
              if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(InMEDA[["PubChem CID"]]),'not available')){
		##########################################
                FPUCID<-stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))
	        ##########################################
	        ##########################################
	        PCID1<-tryCatch({webchem::pc_prop(as.numeric(PCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem Id is empty")})
		PCID2<-tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){message("Pubchem Id is empty")})
                ##############################
                if(!sjmisc::is_empty(PCID2)){

                  FMa<-c(FMa,PCID2)
                }else if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){
		 ##############################################
		  EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem Id is empty")})
		  EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
		  ###############################
		  ###############################
                  if(!sjmisc::is_empty(EM1)){

                    FMa<-c(FMa,EM1)
                  }else{

                    FMa<-c(FMa,PCID2)
                  }
            ###################################
                }else{

                  FMa<-c(FMa,PCID2)
                }
             #########################################
              }else{
                ###### PubchemID is not found... so getting exact mass from
                if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){
		  ########################################
		  EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem Id is empty")})
	          EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
		  ##########################################
	          ##########################################
                  if(!sjmisc::is_empty(EM1)){
                    ##print("formula area.....")
                    ##print(EM1)
                    FMa<-c(FMa,EM1)
                  }else{
                    ##print("formula area empty.....")
                    ##print(PCID2)
                    FMa<-c(FMa,PCID2)
                  }
                }
                #######################################
              }
              #############################
            } ##end if loop Pubchem CID
            ###########
          } ### end of else ...smiles end ... not found
          ########
        }else{
          ### Inchikey end start of smiles
          if(!sjmisc::is_empty(as.character(InMEDA[["SMILES"]])) & !startsWith(as.character(InMEDA[["SMILES"]]),'not available')){
	    print("enter the smiles area ..smiles")
	    ###########################################
            IK<-as.character(InMEDA[["SMILES"]])
	    ###########################################
	    tes<-tryCatch({rcdk::parse.smiles(IK)},error=function(cond){message("smiles not abe to parse")})
	    tes1<-tryCatch({tes[[1]]},error=function(cond){message("smiles parse information is empty")})
	    tes2<-tryCatch({rcdk::get.exact.mass(tes1)},error=function(cond){message("smiles not abe to fetch")})
	    PCID2<-tryCatch({tes2},error=function(cond){message("smiles not abe to fetch")})
	    ################################################
            ################################################
            if(!sjmisc::is_empty(PCID2)){
              FMa<-c(FMa,PCID2)
            }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(InMEDA[["PubChem CID"]]),'not available')){
              FPUCID<-stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))
	      ######################################
	      ######################################
	      PCID1<-tryCatch({webchem::pc_prop(as.numeric(FPUCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem Id is empty")})
	      PCID2<-tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){message("Pubchem Id is empty")})
              ##############################
              if(!sjmisc::is_empty(PCID2)){
                FMa<-c(FMa,PCID2)
              }else if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){
		###############################
		EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem Id is empty")})
		EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
		###############################
                if(!sjmisc::is_empty(EM1)){
                  FMa<-c(FMa,EM1)
                }else{
                  FMa<-c(FMa,PCID2)
                }
                ###########################
              }else{
                ## formula not found and --smiles and pubchem failed to get PubchemID
                FMa<-c(FMa,PCID2)
              }
            }else{
              ### PubchemId is not found ..entering the else loop
              if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){
		############################
		EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem Id is empty")})
                EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
    		###################################		
	        ###################################
                if(!sjmisc::is_empty(EM1)){
                  ##print("formula area.....")
                  ##print(EM1)
                  FMa<-c(FMa,EM1)
                }else{
                  ##print("formula area empty.....")
                  ##print(PCID2)
                  FMa<-c(FMa,PCID2)
                }
              }
              ####################################
            }
          ####################  PCID2 ...end ###
          }else{
            ## smiles end ... not found
            #######################################################
            if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(InMEDA[["PubChem CID"]]),'not available')){
	      ##########################################################
              FPUCID<-stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))
	      ############################
	      PCID1<-tryCatch({webchem::pc_prop(as.numeric(FPUCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem Id is empty")})
	      PCID2<-tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){message("Pubchem Id is empty")})
	      #####################################################
              #####################################################
              if(!sjmisc::is_empty(PCID2)){

                FMa<-c(FMa,PCID2)
              }else if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){

	        ###################################
		EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem Id is empty")})
	        EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
	        ##################################
                if(!sjmisc::is_empty(EM1)){

                  FMa<-c(FMa,EM1)
                }else{

                  FMa<-c(FMa,PCID2)
                }
                #########################
              }else{

                FMa<-c(FMa,PCID2)
              }
              ####################################
            }else{
              ###### PubchemID is not found... so getting exact mass from
              if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){
		#######################################
		EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem Id is empty")})
	    	EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
	        ##################################
                if(!sjmisc::is_empty(EM1)){
                  ##print("formula area.....")
                  ##print(EM1)
                  FMa<-c(FMa,EM1)
                }else{
                  ##print("formula area empty.....")
                  ##print(PCID2)
                  FMa<-c(FMa,PCID2)
                }
	       ###########################
              } ## formula
            } ## end of else loop
          } ## end of else
        }### InCHi.. end
      } ## end of main else
    } ### end of inchikey
  }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["InChI"]]))) & startsWith(as.character(InMEDA[["InChI"]]),'CAS:')){
    #####################################
    print("enter the cas function area.... cas is avilable")
    ####################################
    CV<-stringr::str_trim(as.character(InMEDA[["InChI"]]))
    CV1<-stringr::str_replace(CV,pattern='CAS:',replacement ="")
    CV2<-stringr::str_trim(as.character(CV1))
    ###################################
    PCID<-tryCatch({get_cid(CV2, from = "xref/rn",match="first")},error=function(cond){message("Pubchem Id is empty")})
    ###################################
    PCID1<-tryCatch({webchem::pc_prop(as.numeric(PCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem Id is empty")})
    PCID2<-tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){message("Pubchem Id is empty")})
    ##############################################
    ##############################################
    if(!sjmisc::is_empty(PCID2)){

      FMa<-c(FMa,PCID2)
    }else if(!sjmisc::is_empty(as.character(InMEDA[["SMILES"]])) & !startsWith(as.character(InMEDA[["SMILES"]]),'not available')){
      ##print("enter the smiles..cas")
      ##########################################
      IK<-as.character(InMEDA[["SMILES"]])
      tes<-tryCatch({rcdk::parse.smiles(IK)},error=function(cond){message("smiles not abe to parse")})
      tes1<-tryCatch({tes[[1]]},error=function(cond){message("smiles parse information is empty")})
      tes2<-tryCatch({rcdk::get.exact.mass(tes1)},error=function(cond){message("smiles not abe to fetch")})
      PCID2<-tryCatch({tes2},error=function(cond){message("smiles not abe to fetch")})
      ##########################################
      ##########################################
      if(!sjmisc::is_empty(PCID2)){

        FMa<-c(FMa,PCID2)
      }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(InMEDA[["PubChem CID"]]),'not available')){

        FPUCID<-stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))
        ############################################
        ############################################
        PCID1<-tryCatch({webchem::pc_prop(as.numeric(FPUCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem Id is empty")})
	PCID2<-tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){message("Pubchem Id is empty")})
        ###########################################
        if(!sjmisc::is_empty(PCID2)){

          FMa<-c(FMa,PCID2)
        }else if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){
	  #######################################
	  EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem Id is empty")})
	  EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
	  #########################################
	  #########################################
          if(!sjmisc::is_empty(EM1)){

            FMa<-c(FMa,EM1)
          }else{

            FMa<-c(FMa,PCID2)
          }

        }else{
          ### formula not found ..exact mass
          FMa<-c(FMa,PCID2)
        }
        #######################################
      }else{
        ##print("enter the else area ..2")
        #### PubchemID is not found
	###################################################
        if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){
          ##########################################
	  EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem Id is empty")})
          EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
	  ###########################################
          ##########################################
          if(!sjmisc::is_empty(EM1)){
            ##print("formula area.....")
            ##print(EM1)
            FMa<-c(FMa,EM1)
          }else{
            ##print("formula area empty.....")
            ##print(PCID2)
            FMa<-c(FMa,PCID2)
          }
        }
        ##########################
      }
      ############## Smiles is empty and Pubchem CID
    }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(InMEDA[["PubChem CID"]]),'not available')){

      FPUCID<-stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))
      ###############################
      PCID1<-tryCatch({webchem::pc_prop(as.numeric(FPUCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem Id is empty")})
      PCID2<-tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){message("Pubchem Id is empty")})
      ##################################
      ##################################
      if(!sjmisc::is_empty(PCID2)){

        FMa<-c(FMa,PCID2)
      }else if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){
        #####################################
	EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem Id is empty")})
        EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
	####################################
        ####################################
        if(!sjmisc::is_empty(EM1)){

          FMa<-c(FMa,EM1)
        }else{

          FMa<-c(FMa,PCID2)
        }
	##############################

      }else{
        print("entering Pubchem CID else loop---not got excat mass from pubchem CID")
        print(FPUCID)
        FMa<-c(FMa,PCID2)
      }
      #######################################
    }else{
      ##print("enter the else area ..2")
      #### PubchemID is not found
      if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){
	##############################
	EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem Id is empty")})
        EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
	########################################
        ###################################
        if(!sjmisc::is_empty(EM1)){
          ##print("formula area.....")
          ##print(EM1)
          FMa<-c(FMa,EM1)
        }else{
          ##print("formula area empty.....")
          ##print(PCID2)
          FMa<-c(FMa,PCID2)
        }
	################################
      }else{
        FMa<-c(FMa,PCID2)
      }
    }
  }else if(!sjmisc::is_empty(as.character(InMEDA[["SMILES"]])) & !startsWith(as.character(InMEDA[["SMILES"]]),'not available')){
	  print("enter smiles area --2")
    #######################################
    IK<-as.character(InMEDA[["SMILES"]])
    IK1<-tryCatch({get_cid(IK, from = "smiles")},error=function(cond){message("smiles not abe to fetch")})
    PCID<-tryCatch({IK1[[2]]},error=function(cond){message("Pubchem Id is empty")})
    ########################### changing the smiles part to get exact mass
    tes<-tryCatch({rcdk::parse.smiles(IK)},error=function(cond){message("smiles not abe to parse")})
    tes1<-tryCatch({tes[[1]]},error=function(cond){message("smiles parse information is empty")})
    tes2<-tryCatch({rcdk::get.exact.mass(tes1)},error=function(cond){message("smiles not abe to fetch")})
    PCID2<-tryCatch({tes2},error=function(cond){message("smiles not abe to fetch")})
    #######################################
    #######################################
    if(!sjmisc::is_empty(PCID2)){
      print("enter the if loop ..smiles")
      ##print("checking if entering this area")
      print(PCID2)
      FMa<-c(FMa,PCID2)
    }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(InMEDA[["PubChem CID"]]),'not available')){
	    ########################################
	    print("enter the pubchem CID ..smiles area ..meaning smiles are there and not able to get exact mass..pubchem CID is avilable")
      ##########################################
      FPUCID<-stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))
      ##########################################
      ###########################################
      PCID1<-tryCatch({webchem::pc_prop(as.numeric(FPUCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem Id is empty")})
      PCID2<-tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){message("Pubchem Id is empty")})
      #######################################
      print("entering the pubchem CID")
      print(FPUCID)
      ########################################
      if(!sjmisc::is_empty(PCID2)){
        FMa<-c(FMa,PCID2)
      }else if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){
	      print("enter the if pubchem CID ..else if..")
	##############################
	EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem Id is empty")})
	EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
	##################################
	#################################
        ##print(EM1)
        if(!sjmisc::is_empty(EM1)){
          FMa<-c(FMa,EM1)
        }else{
		print("enter the pubchem CID area ..else part in else if..that means formula is empty")
          FMa<-c(FMa,PCID2)
        }
       ##################################

      }else{
	      print("enter the pubchem CID area...else ...that means neither formula ...nothing is avilable for Pubchem CID")
        FMa<-c(FMa,PCID2)
      }
      ####################################

    }else{
	    #### Pubchem CID #########################
	    print("entering the else part---smiles ... will check formula...for exact mass")
	    ##########################################
	    if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){
		    #######################
		    EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem Id is empty")})
		    EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
		    ##################################
		    ##EM1<-OrgMassSpecR::MolecularWeight(formula = OrgMassSpecR::ListFormula(as.character(InMEDA[["Formula"]])))
		    if(!sjmisc::is_empty(EM1)){
			    FMa<-c(FMa,EM1)
			    }else{
				    FMa<-c(FMa,PCID2)
		    }


	    }else{
		    FMa<-c(FMa,PCID2)
	    }
	    ###########################################
    }
  ###########################################
  }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(InMEDA[["PubChem CID"]]),'not available')){
    ####################################
    print("entering the Pubchem CID area")
    ################################
    FPUCID<-stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))
    ##############################
    PCID1<-tryCatch({webchem::pc_prop(as.numeric(FPUCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem Id is empty")})
    PCID2<-tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){message("Pubchem Id is empty")})
    ###################################
    print("the exact mass of Pubchem CID")
    print(PCID2)
    ###################################
    if(!sjmisc::is_empty(PCID2)){
      print("entering the if loop")
      FMa<-c(FMa,PCID2)
    }else if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){
      ##############################
      EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem Id is empty")})
      EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
      ##############################
      print("entering the else loop of formula..Pubchem CID")
      ##############################
      if(!sjmisc::is_empty(EM1)){
        FMa<-c(FMa,EM1)
      }else{
        FMa<-c(FMa,EM1)
      }

    }else{
      print("entering the else of Pubchem CID")
      FMa<-c(FMa,EM1)
    }
    ################################

  }else{
     #################################	  
    if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){
      #################################
      EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem Id is empty")})
      EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
      #########################
      if(!sjmisc::is_empty(EM1)){
        FMa<-c(FMa,EM1)
      }else{
        FMa<-c(FMa,EM1)
      }
    }else{
      print("this is entering the formula else in final else loop")
      print(EM1)
      FMa<-c(FMa,EM1)
    }
  ##############################
  }

}

##print("enter the area before Ikfilter")
#########################################################################################################
#########################################################################################################
Ikfilter <- function(InKeyVal,InMSPL,InMEDA,InAdVA,InPMZ,InRTL){
  ####################
  print("entering the Inchikey area")
  ###################
  out<-c()
  #####################
  IINF<-tryCatch({webchem::cts_compinfo(InKeyVal)},error=function(cond){message("webchecm could not fetch the info")})
  #####################
  ##if(length(IINF)> 0 & !is.na(IINF)){
    ##!sjmisc::is_empty
  #############################
  if(length(IINF)> 0 & !sjmisc::is_empty(IINF)){	  
    ############################
    IK<-tryCatch({IINF[[1]][1]},error=function(cond){message("Inchiley value is empty")})
    PMZ<-tryCatch({IINF[[1]][4]},error=function(cond){message("PrecursorMZ value is empty")})
    FM<-tryCatch({IINF[[1]][5]},error=function(cond){message("Formula value is empty")})
    CID<-tryCatch({webchem::get_cid(IK, from = "inchikey")},error=function(cond){message("webchecm could not fetch the info")})
    CID1<-tryCatch({CID%>% dplyr::select(cid)},error=function(cond){message("CompoundID is empty; check previous step")})
    CID2<-as.character(CID1)
    CID3<-gsub("[[:punct:]]", "",CID2 )
    CID4<-unlist(strsplit(CID3, " "))
    CID5<-paste(CID4, collapse = ';')
    ###########################
    SM<-tryCatch({webchem::cir_query(IK,"smiles")},error=function(cond){message("webchecm could not fetch the info")})
    SM1<-tryCatch({SM[[1]][1]},error=function(cond){message("smiles Information fetch error")})
    ########################### Inchi value
    InchiV<-tryCatch({rinchi::get.inchi(SM1)},error=function(cond){message("webchecm could not fetch the info")})
    ############################
    AUIN<-as.character(InMEDA[["Adduct"]])
    ############################
    InKeyVal<-tryCatch({IK$inchikey},error=function(cond){message("InchiKey is empty")})
    ############################
    if(!sjmisc::is_empty(AUIN) & !sjmisc::is_empty(PMZ)){
      #if(length(AUIN) > 0 & !is.na(PMZ) ){
      ##########################
      AUIN1<-tryCatch({qdapRegex::ex_between(AUIN, "[", "]")[[1]]},error=function(cond){message("Adduct value is missing")})
      AUIN2<-tryCatch({FADINF(AUIN)},error=function(cond){message("adduct value matching is not found")})
      ##AUIN2<-tryCatch({InAdVA[InAdVA$V1==AUIN1,]$V8},warning=function(cond){message("Adduct value is missing")})
      AAMS<-tryCatch({stringr::str_replace(AUIN2, "M",as.character(PMZ$exactmass))},error=function(cond){message("Missing adduct replacement")})
      AAMS1<-tryCatch({as.numeric(pander::evals(AAMS)[[1]]$result)},error=function(cond){message("Error in adduct replacement step")})
      #########################
      PPm=AAMS1*(25/(1000000))
      #########################
      MPPmL=AAMS1-PPm
      MPPmU=AAMS1+PPm
      ###########################
      Tmass<-InPMZ[InPMZ >= MPPmL & InPMZ <= MPPmU]
      ITmass<-which(InPMZ %in% Tmass)
      ###########################
      ##VRT<-as.numeric(as.character(InMEDA[["RT (min)"]]))
      ##VRTL<-VRT-0.20
      ##VRTU<-VRT+0.20
      #############################
      ##TRTL<-InRTL[InRTL >=VRTL & InRTL <= VRTU]
      ##ITRTL<-which(InRTL %in% TRTL)
      ####################################

      ################################
      ##print("enter the line ...1195")
      ##if(length(ITRTL) >= 1){
      ###################################
        print("enter the line ...1197")
        if(length(ITmass) >= 1){
          print("enter the line ...1198")
          #####################
	  INLL<-ITmass
          ##INLL<-intersect(ITmass,ITRTL)
          #######################
          if(length(INLL) == 1){
	    ###################
            print("enter the line ...1199")
            ###################
            F1FPL<-InMSPL[INLL]
            #################
            SM1<-as.character(InMEDA[["SMILES"]])
            #####################
            F2FPL<-tryCatch({F1FPL},error=function(cond){message("List value is empty")})
            ######################
            Find<-tryCatch({length(F2FPL[[1]])},error=function(cond){message("List value is empty")})
            ########################
            FNA<-tryCatch({F1FPL[[1]]},error=function(cond){message("there is an error in list")})
            ######################## adding this new
	    ##print(FNA)
	    ###########################################
	    ###########################################
            PT <- c("PRECURSORTYPE:", "ADDUCTIONNAME:")
            NPT<-tryCatch({grep(paste(PT,collapse="|"), FNA, value=TRUE)},error=function(cond){message("List value is empty")})
            NPT1<-tryCatch({match(NPT,FNA)},error=function(cond){message("List value is empty")})
            PTV <-tryCatch({stringr::str_remove(FNA[NPT1],c("PRECURSORTYPE:","ADDUCTIONNAME:"))},error=function(cond){message("List value is empty")})
            PTV1<-tryCatch({PTV[1]},error=function(cond){message("List value is empty")})
            PTV2<-tryCatch({stringr::str_trim(PTV1)},error=function(cond){message("List value is empty")})
	    P1TV2<-tryCatch({stringr::str_trim(gsub("ADDUCTIONNAME:","",PTV2))},error=function(cond){message("List value is empty")})
            PTV3<-tryCatch({paste("PRECURSORTYPE:",P1TV2)},error=function(cond){message("List value is empty")})
            ###########################################################
	    print("enter my test...1")
	    print(P1TV2)
	    print(as.character(InMEDA[["Adduct"]]))
	    ###########################################################
            if(identical(P1TV2,as.character(InMEDA[["Adduct"]]))){
	      ########################
              print("enter the line ...95")
              ########################
              FNA1<-which(stringi::stri_detect_fixed(FNA,"NAME:"))
              FNA2<-InMEDA[["Name"]]
              FNA3<-as.character(FNA2)
              #######################
              FNAM<-paste("NAME:",FNA3,sep=" ")
              out<-c(out,FNAM)
              ########################
              FRA1<-which(stringi::stri_detect_fixed(FNA,"RETENTIONTIME:"))
              F1RA1<-FNA[FRA1]
              out<-c(out,F1RA1)
              ################################
              FMZ1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORMZ:"))
              F1MZ1<-FNA[FMZ1]
              out<-c(out,F1MZ1)
              ################################
              FPT1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORTYPE:"))
              F1PT1<-FNA[FPT1]
              out<-c(out,PTV3)
              #out<-c(out,F1PT1)
              #################################
              FIN1<-InMEDA[["Ionization mode"]]
              F1IN1<-as.character(FIN1)
              F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
              out<-c(out,F2IN1)
              ##################################################


	      ###################################################
	      IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(InKeyVal)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
              ##IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
              ##ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
	      #################################################
              if(!sjmisc::is_empty(IKCRV)){
		      ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
		      F1ONT<-paste("Ontology:",ONTV,sep=" ")
		      out<-c(out,F1ONT)
	      }else{
		      F1ONT<-paste("Ontology:","",sep=" ")
		      out<-c(out,F1ONT)
	      }

	      #############################################
              ##F1ONT<-paste("Ontology:",ONTV,sep=" ")
              ##out<-c(out,F1ONT)
              ############################################
              print("enter the line ...1878")
	      ############################################
              FINK<-paste("INCHIKEY:",tryCatch({IK$inchikey},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
              out<-c(out,FINK)
              FINCH<-paste("INCHI:",InchiV,sep=" ")
              out<-c(out,FINCH)
              FSIM<-paste("SMILES:",SM1,sep=" ")
              out<-c(out,FSIM)
              ##############################
              FFOR<-FM$formula
              FFOR1<-paste("FORMULA:",FFOR,sep=" ")
              out<-c(out,FFOR1)
              #############################
              FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
              FINS1<-FNA[FINS]
              out<-c(out,FINS1)
              ############################
              FAUT<-as.character(InMEDA[["Authors"]])
              FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
              out<-c(out,FAUT1)
              ##########################
              ##FLIC<-paste("LICENSE:",sep=" ")
	      FLIC<-paste("LICENSE:","CC BY",sep=" ")
              out<-c(out,FLIC)
              ###########################
              FCIE<-paste("COLLISIONENERGY:",as.character(InMEDA[["Collision energy"]]),sep=" ")
              out<-c(out,FCIE)
              #########################
              FINST<-as.character(InMEDA[["INSTRUMENT_TYPE"]])
              FINST1<-paste("INSTRUMENTTYPE:",FINST,sep=" ")
              out<-c(out,FINST1)
              ########################
              FINS<-as.character(InMEDA[["INSTRUMENT"]])
              FINS1<-paste("INSTRUMENT:",FINS,sep=" ")
              out<-c(out,FINS1)
              ####################
              ##FCOM<-paste("COMMENT:")
	      FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
              out<-c(out,FCOM)
              ##################
              FNPA<-which(stringi::stri_detect_fixed(FNA,"Num Peaks:"))
              F1NPA<-FNA[FNPA]
              out<-c(out,F1NPA)
              ###################
              Fpea<-FNA[(FNPA+1):Find]
              #########################
              if(is.na(Fpea))
              {
                Fpea1<-FNA[(FNPA+1)]


              }else{

                MV=AAMS1
                tes1<-unlist(strsplit(Fpea, "\t"))
                tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
                tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
                tes4<-which(tes2 > (3+MV))
                if(length(tes4)>1)
                {
                  tes5<-tes2[-tes4]
                  tes6<-tes3[-tes4]
                  tes7<-paste(tes5,tes6,sep="\t")
                  out<-c(out,tes7)
                }else{
                  out<-c(out,Fpea)

                }
              }
              ### adding this }
            }
          ###############################################################
          } else if(length(INLL) > 1){
            #####################################
            print("entering the 1569")
	    #####################################
            MONMS=InMSPL[INLL]
	    ##############################
	    #####FMZ1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORMZ:"))
	    ######PRECURSORMZ:
	    #####AAMS1
	    ##### INstead of VRT
	    ##### Commenting above three lines
	    ########################### commenting the original 
            ###TRA<-unname(rapply(MONMS, function(x) grep("RETENTIONTIME:",x, value=TRUE)))
            ###TRA1<-as.numeric(stringr::str_trim(stringr::str_replace(TRA, "RETENTIONTIME:", "")))           
            ###TRA2<-abs(VRT-TRA1)
	    #################### adding this new
	    TRA<-unname(rapply(MONMS, function(x) grep("PRECURSORMZ:",x, value=TRUE)))
	    TRA1<-as.numeric(stringr::str_trim(stringr::str_replace(TRA, "PRECURSORMZ:", "")))
	    TRA2<-abs(AAMS1-TRA1)
	    #########################
            TRA3<-which.min(TRA2)
            TRA4<-INLL[TRA3]
            TRA5<-InMSPL[TRA4]
            #######################
            F1FPL<-TRA5
            ######################
            SM1<-as.character(InMEDA[["SMILES"]])
            ######################
            #InMEDA[["SMILES"]]<-SM1
            #InMEDA[["PubChem CID"]]<-CID5
            #####################
            F2FPL<-F1FPL
            ######################
            Find<-tryCatch({length(F2FPL[[1]])},error=function(cond){message("List value is empty")})
            ########################
            FNA<-tryCatch({F1FPL[[1]]},error=function(cond){message("List value is empty")})
            ########################### adding this new
	    ##print(FNA)
	    ##############################################
	    ##############################################
            PT <- c("PRECURSORTYPE:", "ADDUCTIONNAME:")
	    NPT<-tryCatch({grep(paste(PT,collapse="|"), FNA, value=TRUE)},error=function(cond){message("List value is empty")})
	    NPT1<-tryCatch({match(NPT,FNA)},error=function(cond){message("List value is empty")})
	    PTV <-tryCatch({stringr::str_remove(FNA[NPT1],c("PRECURSORTYPE:","ADDUCTIONNAME:"))},error=function(cond){message("List value is empty")})
	    PTV1<-tryCatch({PTV[1]},error=function(cond){message("List value is empty")})
            PTV2<- tryCatch({stringr::str_trim(PTV1)},error=function(cond){message("List value is empty")})
	    P1TV2<-tryCatch({stringr::str_trim(gsub("ADDUCTIONNAME:","",PTV2))},error=function(cond){message("List value is empty")})
            PTV3 <- tryCatch({paste("PRECURSORTYPE:",P1TV2)},error=function(cond){message("List value is empty")})
            ################################################
	    print("enter my test ...2")
            print(P1TV2)
	    ##########################################
            print(as.character(InMEDA[["Adduct"]]))
	    ##########################################
            if(identical(P1TV2,as.character(InMEDA[["Adduct"]]))){
              ########################
              FNA1<-which(stringi::stri_detect_fixed(FNA,"NAME:"))
              FNA2<-InMEDA[["Name"]]
              FNA3<-as.character(FNA2)
              ######################
              FNAM<-paste("NAME:",FNA3,sep=" ")
              out<-c(out,FNAM)
              ######################
              FRA1<-which(stringi::stri_detect_fixed(FNA,"RETENTIONTIME:"))
              F1RA1<-FNA[FRA1]
              out<-c(out,F1RA1)
              ################################
              FMZ1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORMZ:"))
              F1MZ1<-FNA[FMZ1]
              out<-c(out,F1MZ1)
              ################################
              FPT1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORTYPE:"))
              F1PT1<-FNA[FPT1]
              out<-c(out,PTV3)
              #out<-c(out,F1PT1)
              #############################
              FIN1<-InMEDA[["Ionization mode"]]
              F1IN1<-as.character(FIN1)
              F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
              out<-c(out,F2IN1)
              #############################################


	      #############################################
	      IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(InKeyVal)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
              ##ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
	      ###########################################
	      if(!sjmisc::is_empty(IKCRV)){
		      ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
		      F1ONT<-paste("Ontology:",ONTV,sep=" ")
		      out<-c(out,F1ONT)
	      }else{
		      F1ONT<-paste("Ontology:","",sep=" ")
		      out<-c(out,F1ONT)
	      }
	      ###########################################
              ##F1ONT<-paste("Ontology:",ONTV,sep=" ")
              ##out<-c(out,F1ONT)
              ##########################################
              FINK<-paste("INCHIKEY:",tryCatch({IK$inchikey},error=function(cond){message("Inchikey value is empty")}),sep=" ")
              out<-c(out,FINK)
              FINCH<-paste("INCHI:",InchiV,sep=" ")
              out<-c(out,FINCH)
              FSIM<-paste("SMILES:",SM1,sep=" ")
              out<-c(out,FSIM)
              ##############################
              FFOR<-FM$formula
              FFOR1<-paste("FORMULA:",FFOR,sep=" ")
              out<-c(out,FFOR1)
              ###############################
              FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
              FINS1<-FNA[FINS]
              out<-c(out,FINS1)
              #############################
              FAUT<-as.character(InMEDA[["Authors"]])
              FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
              out<-c(out,FAUT1)
              #############################
              ##FLIC<-paste("LICENSE:",sep=" ")
	      FLIC<-paste("LICENSE:","CC BY",sep=" ")
              out<-c(out,FLIC)
              #############################
              FCIE<-paste("COLLISIONENERGY:",as.character(InMEDA[["Collision energy"]]),sep=" ")
              out<-c(out,FCIE)
              ############################
              FINST<-as.character(InMEDA[["INSTRUMENT_TYPE"]])
              FINST1<-paste("INSTRUMENTTYPE:",FINST,sep=" ")
              out<-c(out,FINST1)
              ##########################
              FINS<-as.character(InMEDA[["INSTRUMENT"]])
              FINS1<-paste("INSTRUMENT:",FINS,sep=" ")
              out<-c(out,FINS1)
              ########################
              ##FCOM<-paste("COMMENT:")
	      FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
              out<-c(out,FCOM)
              #######################
              FNPA<-which(stringi::stri_detect_fixed(FNA,"Num Peaks:"))
              F1NPA<-FNA[FNPA]
              out<-c(out,F1NPA)
              #######################
              Fpea<-FNA[(FNPA+1):Find]
              ########################
              if(is.na(Fpea))
              {
                Fpea1<-FNA[(FNPA+1)]


              }else{

                MV=AAMS1
                tes1<-unlist(strsplit(Fpea, "\t"))
                tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
                tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
                tes4<-which(tes2 > (3+MV))
                if(length(tes4)>1)
                {
                  tes5<-tes2[-tes4]
                  tes6<-tes3[-tes4]
                  tes7<-paste(tes5,tes6,sep="\t")
                  out<-c(out,tes7)
                }else{
                  out<-c(out,Fpea)

                }
              }
              ##############################
            } ### checking ...if the identical
            ################################
          } else{
            #print("entering the line 750")
            PASS<-RRV
          }
          ####### This is testing , if this works
          return(out)
          #################
        } ### this is mz closing brace
	####### this is my commented code and testing this works
      ##} ## this is RT closing braces
    }else{
      #####################################
      print("enter the line 317")
      ######################################
      NCIDV<-tryCatch({webchem::get_cid(InKeyVal, from = "inchikey")},error=function(cond){message("webchecm could not fetch the info")})
      NCIDV1<-tryCatch({as.numeric(NCIDV$cid)},error=function(cond){message("webchecm could not fetch the info")})
      FMWFS<-tryCatch({webchem::cs_compinfo(NCIDV1,c("Formula","MolecularWeight","MonoisotopicMass"), verbose = TRUE)},error=function(cond){message("webchecm could not fetch the info")})
      FMWFS1<-tryCatch({FMWFS$monoisotopicMass},error=function(cond){message("CID is empty")})
      ######################################
      #print(FMWFS1)
      ######################################
      if(!sjmisc::is_empty(AUIN) & !sjmisc::is_empty(FMWFS1)){
	############################################################
        print("enter the line 325")
        AUIN1<-tryCatch({qdapRegex::ex_between(AUIN, "[", "]")[[1]]},error=function(cond){message("Adduct value is missing")})
	AUIN2<-tryCatch({FADINF(AUIN)},error=function(cond){message("adduct value matching is not found")})
        ##AUIN2<-tryCatch({InAdVA[InAdVA$V1==AUIN1,]$V8},warning=function(cond){message("Adduct value is missing")})
        AAMS<-tryCatch({stringr::str_replace(AUIN2, "M",as.character(FMWFS1))},error=function(cond){message("Missing adduct replacement")})
        AAMS1<-tryCatch({as.numeric(pander::evals(AAMS)[[1]]$result)},error=function(cond){message("Error in adduct replacement step")})
        ###################################
        PPm=AAMS1*(25/(1000000))
        #########################
        MPPmL=AAMS1-PPm
        MPPmU=AAMS1+PPm
        ###########################
        Tmass<-InPMZ[InPMZ >= MPPmL & InPMZ <= MPPmU]
        ITmass<-which(InPMZ %in% Tmass)
        ###########################
        ##VRT<-as.numeric(as.character(InMEDA[["RT (min)"]]))
        ##VRTL<-VRT-0.20
        ##VRTU<-VRT+0.20
        #############################
        ##TRTL<-InRTL[InRTL >=VRTL & InRTL <= VRTU]
        ##ITRTL<-which(InRTL %in% TRTL)
        #############################
	##commenting this line and see if this works
	###############################
        ###if(length(ITRTL) >= 1){
	##################################
          print("enter the line ...1197")
          if(length(ITmass) >= 1){
            print("enter the line ...1198")
            ########################
	    INLL<-ITmass
            ##INLL<-intersect(ITmass,ITRTL)
            #######################
            if(length(INLL) == 1){
	      ########################
	      print("enter the line ...1199")
              #####################
              F1FPL<-InMSPL[INLL]
              ####################
              SM1<-as.character(InMEDA[["SMILES"]])
              #####################
              F2FPL<-tryCatch({F1FPL},error=function(cond){message("List value is empty")})
              ######################
              Find<-tryCatch({length(F2FPL[[1]])},error=function(cond){message("List value is empty")})
              ########################
              FNA<-tryCatch({F1FPL[[1]]},error=function(cond){message("there is an error in list")})
              ######################## adding this new ################
	     ## print(FNA)
	      ##########################################################
	      ##########################################################
              PT <- c("PRECURSORTYPE:", "ADDUCTIONNAME:")
              NPT<-tryCatch({grep(paste(PT,collapse="|"), FNA, value=TRUE)},error=function(cond){message("List value is empty")})
              NPT1<-tryCatch({match(NPT,FNA)},error=function(cond){message("List value is empty")})
              PTV <- tryCatch({stringr::str_remove(FNA[NPT1],c("PRECURSORTYPE:","ADDUCTIONNAME:"))},error=function(cond){message("List value is empty")})
              PTV1<- tryCatch({PTV[1]},error=function(cond){message("List value is empty")})
              PTV2<-tryCatch({stringr::str_trim(PTV1)},error=function(cond){message("List value is empty")})
	      P1TV2<-tryCatch({stringr::str_trim(gsub("ADDUCTIONNAME:","",PTV2))},error=function(cond){message("List value is empty")})
              PTV3<-tryCatch({paste("PRECURSORTYPE:",P1TV2)},error=function(cond){message("List value is empty")})
              #########################################################
	      print("enter my test...3")
              print(P1TV2)
              print(as.character(InMEDA[["Adduct"]]))
	      ##############################
              if(identical(P1TV2,as.character(InMEDA[["Adduct"]]))){
                ########################
                FNA1<-which(stringi::stri_detect_fixed(FNA,"NAME:"))
                FNA2<-InMEDA[["Name"]]
                FNA3<-as.character(FNA2)
                #######################
                FNAM<-paste("NAME:",FNA3,sep=" ")
                out<-c(out,FNAM)
                ########################
                FRA1<-which(stringi::stri_detect_fixed(FNA,"RETENTIONTIME:"))
                F1RA1<-FNA[FRA1]
                out<-c(out,F1RA1)
                ################################
                FMZ1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORMZ:"))
                F1MZ1<-FNA[FMZ1]
                out<-c(out,F1MZ1)
                ##########################
                FPT1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORTYPE:"))
                F1PT1<-FNA[FPT1]
                out<-c(out,PTV3)
                #out<-c(out,F1PT1)
                ###########################
                FIN1<-InMEDA[["Ionization mode"]]
                F1IN1<-as.character(FIN1)
                F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                out<-c(out,F2IN1)
                ##################################################
		##################################################
		IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(InKeyVal)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
               ##	ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
		########################################
		if(!sjmisc::is_empty(IKCRV)){
			ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
			F1ONT<-paste("Ontology:",ONTV,sep=" ")
		}else{
			F1ONT<-paste("Ontology:","",sep=" ")
			out<-c(out,F1ONT)
		}
		#################################################
                ##F1ONT<-paste("Ontology:",ONTV,sep=" ")
                ##out<-c(out,F1ONT)
                #################################################
                FINK<-paste("INCHIKEY:",tryCatch({IK$inchikey},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                out<-c(out,FINK)
                FINCH<-paste("INCHI:",InchiV,sep=" ")
                out<-c(out,FINCH)
                FSIM<-paste("SMILES:",SM1,sep=" ")
                out<-c(out,FSIM)
                #######################
                FFOR<-FM$formula
                FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                out<-c(out,FFOR1)
                #######################
                FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                FINS1<-FNA[FINS]
                out<-c(out,FINS1)
                #######################
                FAUT<-as.character(InMEDA[["Authors"]])
                FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
                out<-c(out,FAUT1)
                ##################
                ##FLIC<-paste("LICENSE:",sep=" ")
		FLIC<-paste("LICENSE:","CC BY",sep=" ")
                out<-c(out,FLIC)
                ##################
                FCIE<-paste("COLLISIONENERGY:",as.character(InMEDA[["Collision energy"]]),sep=" ")
                out<-c(out,FCIE)
                ##################
                FINST<-as.character(InMEDA[["INSTRUMENT_TYPE"]])
                FINST1<-paste("INSTRUMENTTYPE:",FINST,sep=" ")
                out<-c(out,FINST1)
                #####################
                FINS<-as.character(InMEDA[["INSTRUMENT"]])
                FINS1<-paste("INSTRUMENT:",FINS,sep=" ")
                out<-c(out,FINS1)
                ####################
                ##FCOM<-paste("COMMENT:")
		FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
                out<-c(out,FCOM)
                ##################
                FNPA<-which(stringi::stri_detect_fixed(FNA,"Num Peaks:"))
                F1NPA<-FNA[FNPA]
                out<-c(out,F1NPA)
                ###################
                Fpea<-FNA[(FNPA+1):Find]
                #########################
                if(is.na(Fpea))
                {
                  Fpea1<-FNA[(FNPA+1)]


                }else{

                  MV=AAMS1
                  tes1<-unlist(strsplit(Fpea, "\t"))
                  tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
                  tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
                  tes4<-which(tes2 > (3+MV))
                  if(length(tes4)>1)
                  {
                    tes5<-tes2[-tes4]
                    tes6<-tes3[-tes4]
                    tes7<-paste(tes5,tes6,sep="\t")
                    out<-c(out,tes7)
                  }else{
                    out<-c(out,Fpea)

                  }
                }

                ##################################
              }
            ############################################################################################
            } else if(length(INLL) > 1){
              ###################################
              print("entering the 1569")
	      ###################################
              MONMS=InMSPL[INLL]
	      ############# commenting the oronal code
              ##TRA<-unname(rapply(MONMS, function(x) grep("RETENTIONTIME:",x, value=TRUE)))
              ##TRA1<-as.numeric(stringr::str_trim(stringr::str_replace(TRA, "RETENTIONTIME:", "")))
              ##TRA2<-abs(VRT-TRA1)
              ##TRA3<-which.min(TRA2)
	      ################# adding this new
	      TRA<-unname(rapply(MONMS, function(x) grep("PRECURSORMZ:",x, value=TRUE)))
	      TRA1<-as.numeric(stringr::str_trim(stringr::str_replace(TRA, "PRECURSORMZ:", "")))
	      TRA2<-abs(AAMS1-TRA1)
	     ##############################
	      TRA3<-which.min(TRA2)
	     ############################## 
              TRA4<-INLL[TRA3]
              TRA5<-InMSPL[TRA4]
              #######################
              F1FPL<-TRA5
              ####################
              SM1<-as.character(InMEDA[["SMILES"]])
              #InMEDA[["SMILES"]]<-SM1
              #InMEDA[["PubChem CID"]]<-CID5
              #####################
              F2FPL<-F1FPL
              ######################
              Find<-tryCatch({length(F2FPL[[1]])},error=function(cond){message("List value is empty")})
              ########################
              FNA<-tryCatch({F1FPL[[1]]},error=function(cond){message("List value is empty")})
              ########################### adding this new
	      ##print(FNA)
	      ###########################################
	      ############################################
              PT <- c("PRECURSORTYPE:", "ADDUCTIONNAME:")
              NPT<-tryCatch({grep(paste(PT,collapse="|"), FNA, value=TRUE)},error=function(cond){message("List value is empty")})
              NPT1<-tryCatch({match(NPT,FNA)},error=function(cond){message("List value is empty")})
              PTV <-tryCatch({stringr::str_remove(FNA[NPT1],c("PRECURSORTYPE:","ADDUCTIONNAME:"))},error=function(cond){message("List value is empty")})
              PTV1<-tryCatch({PTV[1]},error=function(cond){message("List value is empty")})
              PTV2<-tryCatch({stringr::str_trim(PTV1)},error=function(cond){message("List value is empty")})
	      P1TV2<-tryCatch({stringr::str_trim(gsub("ADDUCTIONNAME:","",PTV2))},error=function(cond){message("List value is empty")})
              PTV3<-tryCatch({paste("PRECURSORTYPE:",PTV2)},error=function(cond){message("List value is empty")})
	      ####################################
	      print("enter my test...4")
              print(P1TV2)
              print(as.character(InMEDA[["Adduct"]]))
              ######################## adding this new
              if(identical(P1TV2,as.character(InMEDA[["Adduct"]]))){
                ##if(!sjmisc::is_empty(AUIN) & !sjmisc::is_empty(FMWFS1)){
                ########################
                FNA1<-which(stringi::stri_detect_fixed(FNA,"NAME:"))
                FNA2<-InMEDA[["Name"]]
                FNA3<-as.character(FNA2)
                ######################
                FNAM<-paste("NAME:",FNA3,sep=" ")
                out<-c(out,FNAM)
                ######################
                FRA1<-which(stringi::stri_detect_fixed(FNA,"RETENTIONTIME:"))
                F1RA1<-FNA[FRA1]
                out<-c(out,F1RA1)
                ################################
                FMZ1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORMZ:"))
                F1MZ1<-FNA[FMZ1]
                out<-c(out,F1MZ1)
                ################################
                FPT1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORTYPE:"))
                F1PT1<-FNA[FPT1]
                out<-c(out,PTV3)
                #out<-c(out,F1PT1)
                #############################
                FIN1<-InMEDA[["Ionization mode"]]
                F1IN1<-as.character(FIN1)
                F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                out<-c(out,F2IN1)
                ###################################
		IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(InKeyVal)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
                ##ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
		###################################
		if(!sjmisc::is_empty(IKCRV)){
			ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
			F1ONT<-paste("Ontology:",ONTV,sep=" ")
			out<-c(out,F1ONT)
		}else{
			F1ONT<-paste("Ontology:","",sep=" ")
			out<-c(out,F1ONT)
		}
		###################################
                ###IKCRV<-classyfireR::get_classification(InKeyVal)
                ###ONTV<-do.call(paste, c(as.list(IKCRV@classification$Classification), sep = ","))
		##############################################
                ##F1ONT<-paste("Ontology:",ONTV,sep=" ")
                ##out<-c(out,F1ONT)
                ###################################
                FINK<-paste("INCHIKEY:",tryCatch({IK$inchikey},error=function(cond){message("Inchikey value is empty")}),sep=" ")
                out<-c(out,FINK)
                FINCH<-paste("INCHI:",InchiV,sep=" ")
                out<-c(out,FINCH)
                FSIM<-paste("SMILES:",SM1,sep=" ")
                out<-c(out,FSIM)
                ##############################
                FFOR<-FM$formula
                FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                out<-c(out,FFOR1)
                ###############################
                FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                FINS1<-FNA[FINS]
                out<-c(out,FINS1)
                #############################
                FAUT<-as.character(InMEDA[["Authors"]])
                FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
                out<-c(out,FAUT1)
                #############################
                ##FLIC<-paste("LICENSE:",sep=" ")
		FLIC<-paste("LICENSE:","CC BY",sep=" ")
                out<-c(out,FLIC)
                #############################
                FCIE<-paste("COLLISIONENERGY:",as.character(InMEDA[["Collision energy"]]),sep=" ")
                out<-c(out,FCIE)
                ############################
                FINST<-as.character(InMEDA[["INSTRUMENT_TYPE"]])
                FINST1<-paste("INSTRUMENTTYPE:",FINST,sep=" ")
                out<-c(out,FINST1)
                ##########################
                FINS<-as.character(InMEDA[["INSTRUMENT"]])
                FINS1<-paste("INSTRUMENT:",FINS,sep=" ")
                out<-c(out,FINS1)
                ########################
                ##FCOM<-paste("COMMENT:")
		FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
                out<-c(out,FCOM)
                #######################
                FNPA<-which(stringi::stri_detect_fixed(FNA,"Num Peaks:"))
                F1NPA<-FNA[FNPA]
                out<-c(out,F1NPA)
                #######################
                Fpea<-FNA[(FNPA+1):Find]
                ########################
                if(is.na(Fpea))
                {
                  Fpea1<-FNA[(FNPA+1)]


                }else{

                  MV=AAMS1
                  tes1<-unlist(strsplit(Fpea, "\t"))
                  tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
                  tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
                  tes4<-which(tes2 > (3+MV))
                  if(length(tes4)>1)
                  {
                    tes5<-tes2[-tes4]
                    tes6<-tes3[-tes4]
                    tes7<-paste(tes5,tes6,sep="\t")
                    out<-c(out,tes7)
                  }else{
                    out<-c(out,Fpea)

                  }
                }
            ######################
              }
              #### need to add } to
           ################################
            } else{
              ##print("entering the line 750")
              PASS<-RRV
            }
            ####### This is testing , if this works
            return(out)
         #################
          } ### this is mz closing brace
	  #### commenting this line and see if this works
        ##} ## this is RT closing braces
      ######################
      }

    }

  }
}


#################################################################################################
#################################################################################################
NFFilter<-function(InMEDA,InAdVA,InMSPL,InPMZ,InRTL)
{
  ######################
  out<-c()
  ######################
  print("entering the line...189")
  ######################
  IV<-as.character(InMEDA[["InChI"]])
  EM<-as.character(InMEDA[["Exact mass"]])
  ######################
  EM1<-as.numeric(EM)
  ##################
  NEM1<-tryCatch({NFFilter1(InMEDA,InAdVA,InMSPL,InPMZ,InRTL)},error=function(cond){message("retention time is empty")})
  ##################
  ##print("enter my test area")
  ##print(NEM1)
  ######################
  if(!sjmisc::is_empty(NEM1)){
    #############################
    print("enter the line ..196")
    ##############################
    ADV<-as.character(InMEDA[["Adduct"]])
    ADV1<-tryCatch({qdapRegex::ex_between(ADV, "[", "]")[[1]]},error=function(cond){message("adduct value is missing")})
    ##ADV2<-tryCatch({InAdVA[InAdVA$V1==ADV1,]$V2},warning=function(cond){message("adduct value matching is not found")})
    ADV2<-tryCatch({FADINF(ADV)},error=function(cond){message("adduct value matching is not found")})
    ##ADV3<-tryCatch({stringr::str_replace(ADV2, "M",as.character(EM1))},warning=function(cond){message("adduct value replacement is not found")})
    ADV3<-tryCatch({stringr::str_replace(ADV2, "M",as.character(NEM1))},error=function(cond){message("adduct value replacement is not found")})
    ADV4<-tryCatch({as.numeric(pander::evals(ADV3)[[1]]$result)},error=function(cond){message("getting the result")})
    #####################
    #####################
    PPm=ADV4*(25/(1000000))
    #####################
    MPPmL=ADV4-PPm
    MPPmU=ADV4+PPm
    #####################
    Tmass<-InPMZ[InPMZ >= MPPmL & InPMZ <= MPPmU]
    ITmass<-which(InPMZ %in% Tmass)
    #####################
    #####################
    ##VRT<-as.numeric(as.character(InMEDA[["RT (min)"]]))
    ##VRTL<-VRT-0.20
    ##VRTU<-VRT+0.20
    ######################
    ##TRTL<-InRTL[InRTL >=VRTL & InRTL <= VRTU]
    ##ITRTL<-which(InRTL %in% TRTL)
    ######################
    ##INLL<-intersect(ITmass,ITRTL)
    ###################################
    INLL<-ITmass
    ########################################
    ##if(length(ITRTL) >= 1){
      print("entering the line 219")
      if(length(ITmass) >= 1){
        print("entering the line 221")
        if(length(INLL) == 1){
          ##################################
          print("entering the line 224")
          ##################################
          F1FPL<-InMSPL[INLL]
          F2FPL<-F1FPL
          ############################## added this new ###################
          #################################################################
          Find<-tryCatch({length(F2FPL[[1]])},error=function(cond){message("there is an error in list")})
          #######################################
          FNA<-tryCatch({F1FPL[[1]]},error=function(cond){message("there is an error in list")})
          ######################## added this new ###############
          #######################################################
          ###ADDUCTIONNAME:
          PT <- c("PRECURSORTYPE:", "ADDUCTIONNAME:")
          NPT<-tryCatch({grep(paste(PT,collapse="|"), FNA, value=TRUE)},error=function(cond){message("List value is empty")})
          NPT1<-tryCatch({match(NPT,FNA)},error=function(cond){message("List value is empty")})
          PTV <-tryCatch({stringr::str_remove(FNA[NPT1],c("PRECURSORTYPE:","ADDUCTIONNAME:"))},error=function(cond){message("List value is empty")})
          PTV1<-tryCatch({PTV[1]},error=function(cond){message("List value is empty")})
          PTV2<-tryCatch({stringr::str_trim(PTV1)},error=function(cond){message("List value is empty")})
          P1TV2<-tryCatch({stringr::str_trim(gsub("ADDUCTIONNAME:","",PTV2))},error=function(cond){message("List value is empty")})
          PTV3<-tryCatch({paste("PRECURSORTYPE:",PTV2)},error=function(cond){message("List value is empty")})
          ####################################################
          print("enter the test area...1")
          print(P1TV2)
          print(as.character(InMEDA[["Adduct"]]))
          ######################## adding this new#########
          if(identical(P1TV2,as.character(InMEDA[["Adduct"]]))){
            ##############################################
            print("entering the line 234")
            ##############################################
            FNA1<-which(stringi::stri_detect_fixed(FNA,"NAME:"))
            FNA2<-tryCatch({InMEDA[["Name"]]},warning=function(cond){message("name is empty")})
            FNA3<-as.character(FNA2)
            FNAM<-paste("NAME:",FNA3,sep=" ")
            out<-c(out,FNAM)
            ##########################
            FRA1<-which(stringi::stri_detect_fixed(FNA,"RETENTIONTIME:"))
            F1RA1<-FNA[FRA1]
            out<-c(out,F1RA1)
            ###########################
            FMZ1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORMZ:"))
            F1MZ1<-FNA[FMZ1]
            out<-c(out,F1MZ1)
            ############################
            FPT1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORTYPE:"))
            F1PT1<-FNA[FPT1]
            out<-c(out,PTV3)
            #out<-c(out,F1PT1)
            ###########################
            FIN1<-tryCatch({InMEDA[["Ionization mode"]]},error=function(cond){message("name is empty")})
            F1IN1<-as.character(FIN1)
            F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
            out<-c(out,F2IN1)
            ########################### need to add this new ## this is the ontology part
            ##print("check if problem in ontology part")
            FINIKOT<-MaKE.ONT.REC(InMEDA)
            out<-c(out,FINIKOT)
            ################################################
            ############ this is new #######################
            ################################################
            FSIM<-paste("SMILES:",as.character(InMEDA[["SMILES"]]),sep=" ")
            out<-c(out,FSIM)
            ###################################
            FFOR<-InMEDA[["Formula"]]
            FFOR1<-paste("FORMULA:",FFOR,sep=" ")
            out<-c(out,FFOR1)
            ################################
            FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
            FINS1<-FNA[FINS]
            out<-c(out,FINS1)
            ##############################
            FAUT<-as.character(InMEDA[["Authors"]])
            FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
            out<-c(out,FAUT1)
            #############################
            ##FLIC<-paste("LICENSE:",sep=" ")
            FLIC<-paste("LICENSE:","CC BY",sep=" ")
            out<-c(out,FLIC)
            ###########################
            FCIE<-paste("COLLISIONENERGY:",as.character(InMEDA[["Collision energy"]]),sep=" ")
            out<-c(out,FCIE)
            ##########################
            FINST<-as.character(InMEDA[["INSTRUMENT_TYPE"]])
            FINST1<-paste("INSTRUMENTTYPE:",FINST,sep=" ")
            out<-c(out,FINST1)
            #########################
            FINS<-as.character(InMEDA[["INSTRUMENT"]])
            FINS1<-paste("INSTRUMENT:",FINS,sep=" ")
            out<-c(out,FINS1)
            #########################
            ##FCOM<-paste("COMMENT:")
            FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
            out<-c(out,FCOM)
            ########################
            FNPA<-which(stringi::stri_detect_fixed(FNA,"Num Peaks:"))
            F1NPA<-FNA[FNPA]
            out<-c(out,F1NPA)
            #######################
            Fpea<-FNA[(FNPA+1):Find]
            ######################
            if(is.na(Fpea))
            {
              Fpea1<-FNA[(FNPA+1)]
              out<-c(out,Fpea1)

            }else{
              MV=ADV4
              tes1<-unlist(strsplit(Fpea, "\t"))
              tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
              tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
              tes4<-which(tes2 > (3+MV))
              if(length(tes4)>1)
              {
                tes5<-tes2[-tes4]
                tes6<-tes3[-tes4]
                tes7<-paste(tes5,tes6,sep="\t")
                out<-c(out,tes7)
              }else{
                out<-c(out,Fpea)
              }
            } ## end of else
         #####################################################
          } ### adduct match
        #########################################################
        }else if(length(INLL) > 1){
          ##############################
          print("enter the line ...511")
          ##############################
          MONMS=InMSPL[INLL]
	  ##################### Commenting the Original
          ##TRA<-unname(rapply(MONMS, function(x) grep("RETENTIONTIME:",x, value=TRUE)))
          ##TRA1<-as.numeric(stringr::str_trim(stringr::str_replace(TRA, "RETENTIONTIME:", "")))
          ##TRA2<-abs(VRT-TRA1)
	  ##################### Adding this new
	  TRA<-unname(rapply(MONMS, function(x) grep("PRECURSORMZ:",x, value=TRUE)))
	  TRA1<-as.numeric(stringr::str_trim(stringr::str_replace(TRA, "PRECURSORMZ:", "")))
	  TRA2<-abs(ADV4-TRA1)
	  #############################
          TRA3<-which.min(TRA2)
          TRA4<-INLL[TRA3]
          TRA5<-InMSPL[TRA4]
          ###########################
          F1FPL<-TRA5
          F2FPL<-F1FPL
          ##########################
          Find<-tryCatch({length(F2FPL[[1]])},error=function(cond){message("there is an error in list")})
          ########################
          FNA<-tryCatch({F1FPL[[1]]},error=function(cond){message("there is an error in list")})
          ####### adding this new ###################################################
          ###########################################################################
          PT <- c("PRECURSORTYPE:", "ADDUCTIONNAME:")
          NPT<-tryCatch({grep(paste(PT,collapse="|"), FNA, value=TRUE)},error=function(cond){message("List value is empty")})
          NPT1<-tryCatch({match(NPT,FNA)},error=function(cond){message("List value is empty")})
          PTV <-tryCatch({stringr::str_remove(FNA[NPT1],c("PRECURSORTYPE:","ADDUCTIONNAME:"))},error=function(cond){message("List value is empty")})
          PTV1<-tryCatch({PTV[1]},error=function(cond){message("List value is empty")})
          PTV2<-tryCatch({stringr::str_trim(PTV1)},error=function(cond){message("List value is empty")})
          P1TV2<-tryCatch({stringr::str_trim(gsub("ADDUCTIONNAME:","",PTV2))},error=function(cond){message("List value is empty")})
          PTV3<-tryCatch({paste("PRECURSORTYPE:",PTV2)},error=function(cond){message("List value is empty")})
          ##########################################################################
          print("enter the test area..2")
          print(P1TV2)
          print(as.character(InMEDA[["Adduct"]]))
          ##################################################
          if(identical(P1TV2,as.character(InMEDA[["Adduct"]]))){
            ################################################
            print("enter the line ...511")
            ################################################
            FNA1<-which(stringi::stri_detect_fixed(FNA,"NAME:"))
            FNA2<-InMEDA[["Name"]]
            FNA3<-as.character(FNA2)
            FNAM<-paste("NAME:",FNA3,sep=" ")
            out<-c(out,FNAM)
            ####################
            FRA1<-which(stringi::stri_detect_fixed(FNA,"RETENTIONTIME:"))
            F1RA1<-FNA[FRA1]
            out<-c(out,F1RA1)
            ##################
            FMZ1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORMZ:"))
            F1MZ1<-FNA[FMZ1]
            out<-c(out,F1MZ1)
            ###################
            FPT1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORTYPE:"))
            F1PT1<-FNA[FPT1]
            out<-c(out,PTV3)
            #out<-c(out,F1PT1)
            ############################################################
            FIN1<-InMEDA[["Ionization mode"]]
            F1IN1<-as.character(FIN1)
            F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
            out<-c(out,F2IN1)
            ##################### adding this new ####################
            FINIKOT<-MaKE.ONT.REC(InMEDA)
            out<-c(out,FINIKOT)
            ###################### this is end########################
            ############################################### this is new
            FSIM<-paste("SMILES:",as.character(InMEDA[["SMILES"]]),sep=" ")
            out<-c(out,FSIM)
            ##########################################################
            FFOR<-InMEDA[["Formula"]]
            FFOR1<-paste("FORMULA:",FFOR,sep=" ")
            out<-c(out,FFOR1)
            ##########################################################
            FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
            FINS1<-FNA[FINS]
            out<-c(out,FINS1)
            ###########################################
            FAUT<-as.character(InMEDA[["Authors"]])
            FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
            out<-c(out,FAUT1)
            ###########################################
            ##FLIC<-paste("LICENSE:",sep=" ")
            FLIC<-paste("LICENSE:","CC BY",sep=" ")
            out<-c(out,FLIC)
            ##########################################
            FCIE<-paste("COLLISIONENERGY:",as.character(InMEDA[["Collision energy"]]),sep=" ")
            out<-c(out,FCIE)
            #########################################
            FINST<-as.character(InMEDA[["INSTRUMENT_TYPE"]])
            FINST1<-paste("INSTRUMENTTYPE:",FINST,sep=" ")
            out<-c(out,FINST1)
            #########################################
            FINS<-as.character(InMEDA[["INSTRUMENT"]])
            FINS1<-paste("INSTRUMENT:",FINS,sep=" ")
            out<-c(out,FINS1)
            ########################################
            ##FCOM<-paste("COMMENT:")
            FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
            out<-c(out,FCOM)
            ######################################
            FNPA<-which(stringi::stri_detect_fixed(FNA,"Num Peaks:"))
            F1NPA<-FNA[FNPA]
            out<-c(out,F1NPA)
            ###################
            Fpea<-FNA[(FNPA+1):Find]
            ###################################
            if(is.na(Fpea))
            {
              Fpea1<-FNA[(FNPA+1)]
              out<-c(out,Fpea1)

            }else{
              MV=ADV4
              tes1<-unlist(strsplit(Fpea, "\t"))
              tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
              tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
              tes4<-which(tes2 > (3+MV))
              if(length(tes4)>1)
              {
                tes5<-tes2[-tes4]
                tes6<-tes3[-tes4]
                tes7<-paste(tes5,tes6,sep="\t")
                out<-c(out,tes7)
              }else{
                out<-c(out,Fpea)
              }
            } ### this is the else MV= ADV4...closing
        ############################## adding this new
          }
        ###############################
        } else {
          #print("entering the line 458")
          PASS<-RRV
        }  ## main else part
        ###testing if this works
        return(out)
        #####################
      } #ITmass
    #######will comment this line and check if this works
    ##} ##IITRTL
    #######################
  }else{
    ### that means not able to get mass from smiles,pubchem ID and
    ##########################################
    print("entering this line..line 767")
    ###########################################
    FM<-as.character(InMEDA[["Exact mass"]])
    ##########################################
    print("test area ...exact mass area is avilable")
    print(!sjmisc::is_empty(FM))
    ###########################################
    if(!sjmisc::is_empty(FM))
    {
      ###################################################################
      print("entering the line ... 772")
      ###################################################################
      ##TE<-FM
      ##EM1<-OrgMassSpecR::MolecularWeight(formula = OrgMassSpecR::ListFormula(FM))
      ##if(!sjmisc::is_empty( EM1)) {
      ####################################################################
      if(!sjmisc::is_empty(FM) & !sjmisc::is_empty(as.character(InMEDA[["Adduct"]]))) {
        #############################################
        print("entering the line ... 779")
        #############################################
        ADV<-as.character(InMEDA[["Adduct"]])
        ADV1<-tryCatch({qdapRegex::ex_between(ADV, "[", "]")[[1]]},error=function(cond){message("adduct value is missing")})
        ##ADV2<-tryCatch({InAdVA[InAdVA$V1==ADV1,]$V2},warning=function(cond){message("adduct value matching is not found")})
	ADV2<-tryCatch({FADINF(ADV)},error=function(cond){message("adduct value matching is not found")})
        ADV3<-tryCatch({stringr::str_replace(ADV2, "M",as.character(FM))},error=function(cond){message("adduct value replacement is not found")})
        ADV4<-tryCatch({as.numeric(pander::evals(ADV3)[[1]]$result)},error=function(cond){message("getting the result")})
        ######################
        PPm=ADV4*(25/(1000000))
        ######################
        MPPmL=ADV4-PPm
        MPPmU=ADV4+PPm
        #######################
        Tmass<-InPMZ[InPMZ >= MPPmL & InPMZ <= MPPmU]
        ITmass<-which(InPMZ %in% Tmass)
        #####################
        ##VRT<-as.numeric(as.character(InMEDA[["RT (min)"]]))
        ##VRTL<-VRT-0.20
        ##VRTU<-VRT+0.20
        ######################
        ##TRTL<-InRTL[InRTL >=VRTL & InRTL <= VRTU]
        ##ITRTL<-which(InRTL %in% TRTL)
        ######################
        ##INLL<-intersect(ITmass,ITRTL)
	############################
	INLL<-ITmass
        ########################################
        ##if(length(ITRTL) >= 1){
	#########################################
          print("entering the line 514")
          if(length(ITmass) >= 1){
            print("entering the line 517")
            if(length(INLL) == 1){
              ############################
              print("entering the line 519")
              #############################
              F1FPL<-InMSPL[INLL]
              F2FPL<-F1FPL
              #######################
              Find<-tryCatch({length(F2FPL[[1]])},error=function(cond){message("there is an error in list")})
              ########################
              FNA<-tryCatch({F1FPL[[1]]},error=function(cond){message("there is an error in list")})
              ######################## adding this new ###########
              ####################################################
              PT <- c("PRECURSORTYPE:", "ADDUCTIONNAME:")
              NPT<-tryCatch({grep(paste(PT,collapse="|"), FNA, value=TRUE)},error=function(cond){message("List value is empty")})
              NPT1<-tryCatch({match(NPT,FNA)},error=function(cond){message("List value is empty")})
              PTV <-tryCatch({stringr::str_remove(FNA[NPT1],c("PRECURSORTYPE:","ADDUCTIONNAME:"))},error=function(cond){message("List value is empty")})
              PTV1<-tryCatch({PTV[1]},error=function(cond){message("List value is empty")})
              PTV2<-tryCatch({stringr::str_trim(PTV1)},error=function(cond){message("List value is empty")})
              P1TV2<-tryCatch({stringr::str_trim(gsub("ADDUCTIONNAME:","",PTV2))},error=function(cond){message("List value is empty")})
              PTV3<-tryCatch({paste("PRECURSORTYPE:",PTV2)},error=function(cond){message("List value is empty")})
              #####################################################
              print("enter the test area..3")
	      ####################################################
              print(P1TV2)
              print(as.character(InMEDA[["Adduct"]]))
              ################################ adding this new ###
              if(identical(P1TV2,as.character(InMEDA[["Adduct"]]))){
                ##################################################
                print("entering the line 587")
                ##############################
                FNA1<-which(stringi::stri_detect_fixed(FNA,"NAME:"))
                FNA2<-InMEDA[["Name"]]
                FNA3<-as.character(FNA2)
                FNAM<-paste("NAME:",FNA3,sep=" ")
                out<-c(out,FNAM)
                ####################
                FRA1<-which(stringi::stri_detect_fixed(FNA,"RETENTIONTIME:"))
                F1RA1<-FNA[FRA1]
                out<-c(out,F1RA1)
                ##################
                FMZ1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORMZ:"))
                F1MZ1<-FNA[FMZ1]
                out<-c(out,F1MZ1)
                ###################
                FPT1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORTYPE:"))
                F1PT1<-FNA[FPT1]
                out<-c(out,PTV3)
                #out<-c(out,F1PT1)
                #################################
                FIN1<-InMEDA[["Ionization mode"]]
                F1IN1<-as.character(FIN1)
                F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                out<-c(out,F2IN1)
                ########################### adding this new ###
                FINIKOT<-MaKE.ONT.REC(InMEDA)
                out<-c(out,FINIKOT)
                ########################## this is the end #####
                ################################################
                FSIM<-paste("SMILES:",as.character(InMEDA[["SMILES"]]),sep=" ")
                out<-c(out,FSIM)
                ################################
                FFOR<-InMEDA[["Formula"]]
                FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                out<-c(out,FFOR1)
                ###############################
                FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                FINS1<-FNA[FINS]
                out<-c(out,FINS1)
                ##############################
                FAUT<-as.character(InMEDA[["Authors"]])
                FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
                out<-c(out,FAUT1)
                #############################
                ##FLIC<-paste("LICENSE:",sep=" ")
                FLIC<-paste("LICENSE:","CC BY",sep=" ")
                out<-c(out,FLIC)
                ###########################
                FCIE<-paste("COLLISIONENERGY:",as.character(InMEDA[["Collision energy"]]),sep=" ")
                out<-c(out,FCIE)
                ##########################
                FINST<-as.character(InMEDA[["INSTRUMENT_TYPE"]])
                FINST1<-paste("INSTRUMENTTYPE:",FINST,sep=" ")
                out<-c(out,FINST1)
                #########################
                FINS<-as.character(InMEDA[["INSTRUMENT"]])
                FINS1<-paste("INSTRUMENT:",FINS,sep=" ")
                out<-c(out,FINS1)
                #########################
                ##FCOM<-paste("COMMENT:")
                FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
                out<-c(out,FCOM)
                ########################
                FNPA<-which(stringi::stri_detect_fixed(FNA,"Num Peaks:"))
                F1NPA<-FNA[FNPA]
                out<-c(out,F1NPA)
                #######################
                Fpea<-FNA[(FNPA+1):Find]
                ######################
                if(is.na(Fpea))
                {
                  Fpea1<-FNA[(FNPA+1)]
                  out<-c(out,Fpea1)

                }else{
                  MV=ADV4
                  tes1<-unlist(strsplit(Fpea, "\t"))
                  tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
                  tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
                  tes4<-which(tes2 > (3+MV))
                  if(length(tes4)>1)
                  {
                    tes5<-tes2[-tes4]
                    tes6<-tes3[-tes4]
                    tes7<-paste(tes5,tes6,sep="\t")
                    out<-c(out,tes7)
                  }else{
                    out<-c(out,Fpea)
                  }
                } ## end of else
           #################################### adding this } new
              } ## adduct match
           #########################################################
            }else if(length(INLL) > 1){
              ##############################################
              print("enter the line ...1020")
              ##############################################
              MONMS=InMSPL[INLL]
	      ####################### Commenting the original
              ##TRA<-unname(rapply(MONMS, function(x) grep("RETENTIONTIME:",x, value=TRUE)))
              ##TRA1<-as.numeric(stringr::str_trim(stringr::str_replace(TRA, "RETENTIONTIME:", "")))
              ##TRA2<-abs(VRT-TRA1)
	      ###################### Adding this new
	      TRA<-unname(rapply(MONMS, function(x) grep("PRECURSORMZ:",x, value=TRUE)))
	      TRA1<-as.numeric(stringr::str_trim(stringr::str_replace(TRA, "PRECURSORMZ:", "")))
	      TRA2<-abs(ADV4-TRA1)
	      ##########################
              TRA3<-which.min(TRA2)
              TRA4<-INLL[TRA3]
              TRA5<-InMSPL[TRA4]
              ###########################
              F1FPL<-TRA5
              F2FPL<-F1FPL
              ##########################
              Find<-tryCatch({length(F2FPL[[1]])},error=function(cond){message("there is an error in list")})
              ##########################
              FNA<-tryCatch({F1FPL[[1]]},error=function(cond){message("there is an error in list")})
              ######################## adding this new
              ###########################################
              PT <- c("PRECURSORTYPE:", "ADDUCTIONNAME:")
              NPT<-tryCatch({grep(paste(PT,collapse="|"), FNA, value=TRUE)},error=function(cond){message("List value is empty")})
              NPT1<-tryCatch({match(NPT,FNA)},error=function(cond){message("List value is empty")})
              PTV <-tryCatch({stringr::str_remove(FNA[NPT1],c("PRECURSORTYPE:","ADDUCTIONNAME:"))},error=function(cond){message("List value is empty")})
              PTV1<-tryCatch({PTV[1]},error=function(cond){message("List value is empty")})
              PTV2<-tryCatch({stringr::str_trim(PTV1)},error=function(cond){message("List value is empty")})
              P1TV2<-tryCatch({stringr::str_trim(gsub("ADDUCTIONNAME:","",PTV2))},error=function(cond){message("List value is empty")})
              PTV3<-tryCatch({paste("PRECURSORTYPE:",P1TV2)},error=function(cond){message("List value is empty")})
              ########################################
              print("enter the test area..4")
	      ########################################
              print(P1TV2)
              print(as.character(InMEDA[["Adduct"]]))
              ################################## adding this new
              if(identical(P1TV2,as.character(InMEDA[["Adduct"]]))){
                ##################################################
                print("enter the line ...1020")
                #################################################
                FNA1<-which(stringi::stri_detect_fixed(FNA,"NAME:"))
                FNA2<-InMEDA[["Name"]]
                FNA3<-as.character(FNA2)
                FNAM<-paste("NAME:",FNA3,sep=" ")
                out<-c(out,FNAM)
                ####################
                FRA1<-which(stringi::stri_detect_fixed(FNA,"RETENTIONTIME:"))
                F1RA1<-FNA[FRA1]
                out<-c(out,F1RA1)
                ########################
                FMZ1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORMZ:"))
                F1MZ1<-FNA[FMZ1]
                out<-c(out,F1MZ1)
                ##########################
                FPT1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORTYPE:"))
                F1PT1<-FNA[FPT1]
                out<-c(out,PTV3)
                #out<-c(out,F1PT1)
                ###########################
                FIN1<-InMEDA[["Ionization mode"]]
                F1IN1<-as.character(FIN1)
                F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                out<-c(out,F2IN1)
                ########################### adding this new #####
                FINIKOT<-MaKE.ONT.REC(InMEDA)
                out<-c(out,FINIKOT)
                ###########################this is end ########
                ###############################################
                FSIM<-paste("SMILES:",as.character(InMEDA[["SMILES"]]),sep=" ")
                out<-c(out,FSIM)
                ############################################
                FFOR<-InMEDA[["Formula"]]
                FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                out<-c(out,FFOR1)
                ############################################
                FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                FINS1<-FNA[FINS]
                out<-c(out,FINS1)
                ###########################################
                FAUT<-as.character(InMEDA[["Authors"]])
                FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
                out<-c(out,FAUT1)
                ###########################################
                ##FLIC<-paste("LICENSE:",sep=" ")
                FLIC<-paste("LICENSE:","CC BY",sep=" ")
                out<-c(out,FLIC)
                ##########################################
                FCIE<-paste("COLLISIONENERGY:",as.character(InMEDA[["Collision energy"]]),sep=" ")
                out<-c(out,FCIE)
                #########################################
                FINST<-as.character(InMEDA[["INSTRUMENT_TYPE"]])
                FINST1<-paste("INSTRUMENTTYPE:",FINST,sep=" ")
                out<-c(out,FINST1)
                #########################################
                FINS<-as.character(InMEDA[["INSTRUMENT"]])
                FINS1<-paste("INSTRUMENT:",FINS,sep=" ")
                out<-c(out,FINS1)
                ######################
                ##FCOM<-paste("COMMENT:")
                FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
                out<-c(out,FCOM)
                ####################
                FNPA<-which(stringi::stri_detect_fixed(FNA,"Num Peaks:"))
                F1NPA<-FNA[FNPA]
                out<-c(out,F1NPA)
                ###################
                Fpea<-FNA[(FNPA+1):Find]
                ###################
                if(is.na(Fpea))
                {
                  Fpea1<-FNA[(FNPA+1)]
                  out<-c(out,Fpea1)

                }else{
                  MV=ADV4
                  tes1<-unlist(strsplit(Fpea, "\t"))
                  tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
                  tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
                  tes4<-which(tes2 > (3+MV))
                  if(length(tes4)>1)
                  {
                    tes5<-tes2[-tes4]
                    tes6<-tes3[-tes4]
                    tes7<-paste(tes5,tes6,sep="\t")
                    out<-c(out,tes7)
                  }else{
                    out<-c(out,Fpea)
                  }
                } ### this is the else MV= ADV4...closing
            ###########################
              }
            ##########################
            } else {
              #print("entering the line 764")
	      print("enter the final else part")
              PASS<-RRV
            }  ## main else part
            ###testing if this works
	    ### THis is original return area
	    #### Have to change this I guess
            ##return(out)
	    ####################
            #####################
          } #ITmass
	###############
        ##} #ITRTL
        ###################
      } ## adduct match and exacct mass present
    } else{
      ##print("enter the final else..2 part")	    
      PASS<-RRV
    }
  } ## end of else
  ############ adding this here
  return(out)
  ###########################
}###NFFilter


      
      
      
    
###########################################################################################

##print("enter the are before main loop")
###########################################################################################
for(i in 1:length(LmeCmu1))
###for(i in 1:2)
{
  ########################
  print("entering the main function")
  ########################
  
  Val=LmeCmu1[i]
  Val1=LmeCmu1[i+1]
  nVal=Val+1
  #if(!is.na(Val) & !is.na(Val1))
  if(!sjmisc::is_empty(Val) & !sjmisc::is_empty(Val1))
  {
    #####################
    print("entering the line ...1707")
    NRXF3<-RXF3[nVal:Val1,]
    ##NRXF3<-RXF3[Val:Val1,]
    FiNA<-SFileNam[i]
    BNFiNA<-basename(FiNA)
    BNFiNAEX<-tools::file_path_sans_ext(BNFiNA)
    BNFiNAEX1<-gsub("NRG[\\]Set01", "", BNFiNAEX)
    BNFiNAEX2<-gsub("NRG[\\]", "", BNFiNAEX1)
    BNFiNAEX3<-gsub("NRG[\\]NRG set 12[\\]", "", BNFiNAEX2)
    BNFiNAEX4<-gsub("NRG[\\]Set 11 NRG[\\]", "", BNFiNAEX3)
    BNFiNAEX5<-gsub("NRG[\\]Set 5 NRG[\\]", "", BNFiNAEX4)
    BNFiNAEX6<-gsub("NRG[\\]SET 6 NRG[\\]", "", BNFiNAEX5)
    BNFiNAEX7<-gsub("NRG[\\]SET 8 NRG[\\]", "", BNFiNAEX6)
    BNFiNAEX8<-gsub("NRG[\\]SET 9 NRG[\\]", "", BNFiNAEX7)
    BNFiNAEX9<-gsub("NRG[\\]Set01 NRG[\\]", "", BNFiNAEX8)
    BNFiNAEX10<-gsub("NRG[\\]Set 11 NRG[\\]", "", BNFiNAEX9)
    BNFiNAEX11<-gsub("NRG[\\]NRG set 12[\\]", "", BNFiNAEX10)
    BNFiNAEX12<-gsub("Set02[\\]", "", BNFiNAEX11)
    BNFiNAEX13<-gsub("Set 11", "", BNFiNAEX12)
    BNFiNAEX14<-gsub("NRG set 12[\\]", "", BNFiNAEX13)
    BNFiNAEX15<-gsub("Set 13", "", BNFiNAEX14)
    BNFiNAEX16<-gsub("Set02", "", BNFiNAEX15)
    BNFiNAEX17<-gsub("SET03", "", BNFiNAEX16)
    BNFiNAEX18<-gsub("Set 4", "", BNFiNAEX17)
    BNFiNAEX19<-gsub("Set 5", "", BNFiNAEX18)
    BNFiNAEX20<-gsub("SET 6", "", BNFiNAEX19)
    BNFiNAEX21<-gsub("SET 8", "", BNFiNAEX20)
    BNFiNAEX22<-gsub("SET 9", "", BNFiNAEX21)
    BNFiNAEX23<-stringr::str_trim(BNFiNAEX22)
    BNFiNAEX24<-gsub("NRG[\\]Set 10 NRG[\\]", "", BNFiNAEX23)
    BNFiNAEX25<-stringr::str_trim(BNFiNAEX24)
    ########### changing this############# 
    ##FINAMSP<-paste(BNFiNAEX23,"msp",sep=".")
    FINAMSP<-paste(BNFiNAEX25,"msp",sep=".")
    ####################################
    ##print(FINAMSP)
    ##print(Fi6)
    ####################################
    OUNA<-tools::file_path_sans_ext(FiNA)
    OUNA1<-paste(OUNA,"passed","msp",sep = ".")
    OUNA2<-paste(Fi9,OUNA1,sep="")
    ##########mkdir mz.25ppm.20RT
    OUNA3<-paste(Fi11,OUNA1,sep="")
    #########mkdir mz.40ppm.35RT
    OUNA4<-paste(Fi13,OUNA1,sep="")
    #########mkdir mz.50ppm.40RT
    OUNA5<-paste(Fi15,OUNA1,sep="")
    ##############################
    ##print(OUNA5)
    ###############FiNA1 ... INput ..msp file ###########################
    #####OUNA2 ... Output file
    ##FiNA1<-tryCatch({list.files(path =Fi6 , pattern = FINAMSP, recursive = TRUE, full.names = TRUE)},warning=function(cond){message("some mistake happened in file search files")})
    ######################################################################
    ###############FiNA1 ... INput ..msp file
    ##NFiNA<-tryCatch({list.files(path =N1Fi6, recursive = TRUE, full.names = TRUE)},warning=function(cond){message("some mistake happened in file search files")})
    #######################################################################
    NFiNA<-tryCatch({list.files(path =Fi6 , recursive = TRUE, full.names = TRUE)},warning=function(cond){message("some mistake happened in file search files")})
    NFiNA1<-tryCatch({match(FINAMSP,basename(NFiNA))},warning=function(cond){message("some mistake happened in file search files")})
    NFiNA2<-tryCatch({NFiNA[NFiNA1]},warning=function(cond){message("some mistake happened in file search files")})
    #####OUNA2 ... Output file ##########################################
    FiNA1<-NFiNA2
    #####################################################################
    print(FiNA1)
    print(FINAMSP)
    #####################################################################
    if((length(FiNA1) >= 1))
    {
      #####################################	    
      print("entering the line ...line 1723")
      #####################################
      CaIF<-tryCatch({MaKlist(FiNA1[1])},error=function(cond){message("some mistake happened in filename..file name must be empty")})
      lst2<-tryCatch({CaIF[[1]]},error=function(cond){message("some mistake happened in indexes")})
      fmass<-tryCatch({CaIF[[8]]},error=function(cond){message("some mistake happened in indexes")})
      FRTL1<-tryCatch({CaIF[[15]]},error=function(cond){message("some mistake happened in indexes")})
      ########################################
      if(file.exists(FiNA1[1]))
      {
	#######################################
        print("entering the line ...1732")
        ######################################
        if(length(NRXF3)>0){
          for (i in 1:length(NRXF3)){
            ################
            RRV<-NRXF3[i,]
            ################
            PIN<-ifelse(is.na(RRV[["InChI"]]), 1, 0)
            PSM<-ifelse(is.na(RRV[["SMILES"]]), 1, 0)
            PPC<-ifelse(is.na(RRV[["PubChem CID"]]), 1, 0)
            ####################
            if(PIN == 1 & PSM == 1 & PPC == 1){
              #########
              PASS<-RRV
              #########
            }
            ######################################################################################################################################
	    else{
              print("enter the line ...1749")
	      #####################################################
              if(!sjmisc::is_empty(stringr::str_trim(as.character(RRV[["InChI"]]))) & !startsWith(as.character(RRV[["InChI"]]),'not available') & !startsWith(as.character(RRV[["InChI"]]),'CAS:') & !startsWith(as.character(RRV[["InChI"]]),'InChI=')){
                print("enter the line ...1841")
		##################################################
                if(tryCatch({webchem::is.inchikey(stringr::str_trim(as.character(RRV[["InChI"]])))},warning=function(cond){message("not able to pass is.inchiKey")})){
                  IV<-stringr::str_trim(as.character(RRV[["InChI"]]))
	                #########################################
		         print("enter the line ...1867")
		              if(!sjmisc::is_empty(IV)){  
                    outn<-Ikfilter(IV,lst2,RRV,AIN,fmass,FRTL1)
                    len1<-length(outn)
                    ##############################################
                    if(len1 > 1)
                    {
                      print("enter the line ...845")
                      cat(sapply(outn, toString), file=OUNA3, sep="\n",append=TRUE)
                      cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
                      
                    }else{
                      print("enter the line ...1765")
                      FINKE<-NFFilter(RRV,AIN,lst2,fmass,FRTL1)
                      ### adding this code new 
                      len3<-length(FINKE)
                      if(len3 > 1){
                        print("enter the line ...856")
                        cat(sapply(FINKE, toString), file=OUNA3, sep="\n",append=TRUE)
                        cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
                      }else{
			      ### adding this new for error report
			      DN<-dirname(dirname(File1))
			      DN1<-paste(DN,"Error-Report",sep="/")
			      FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
			      write.table(unname(as.data.frame(RRV)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
		      }
                      
                    } # end of the else loop
                                      
                }else{
	        ###########################################	
                  print("enter the line ...1780")
                  FINKE<-NFFilter(RRV,AIN,lst2,fmass,FRTL1)
                  ### adding this code new 
                  len3<-length(FINKE)
                  if(len3 > 1){
                    print("enter the line ...873")
                    cat(sapply(FINKE, toString), file=OUNA3, sep="\n",append=TRUE)
                    cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
                  }else{
			  ### adding this new for error report
			  DN<-dirname(dirname(File1))
			  DN1<-paste(DN,"Error-Report",sep="/")
			  FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
			  write.table(unname(as.data.frame(RRV)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
		  }
		 ####################################### 
                }
		}
              ######################################################################################################################################
              }else if(!sjmisc::is_empty(as.character(RRV[["InChI"]])) & startsWith(as.character(RRV[["InChI"]]),'InChI=')){
                print("enter the line ...861")
                ##IV1<-stringr::str_trim(as.character(RRV[["InChI"]]))
	        #####################################################
	        IV1<-stringr::str_trim(as.character(RRV[["InChI"]]))
		#####################################################
		FSMV<-tryCatch({rinchi::parse.inchi(IV1)},error=function(cond){message("Inchi name must be empty or rinchi not abe to fetch")})
		FSMV1<-tryCatch({rcdk::get.smiles(FSMV[[1]])},error=function(cond){message("name is empty")})
		FSMV2<-tryCatch({rinchi::get.inchi.key(FSMV1)},error=function(cond){message("webchecm could not fetch the info")})
		#############################################################
		##FSMV<-tryCatch({res <- R.utils::withTimeout({chemspiderapi::post_convert(IV1,inputFormat = "InChI",outputFormat ="SMILES", apikey <- apikey)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
                ##FSMV<-tryCatch({chemspiderapi::post_convert(IV1,inputFormat = "InChI",outputFormat ="SMILES", apikey <- apikey)},warning=function(cond){message("webchecm could not fetch the info")})
                ##FSMV1<-tryCatch({unname(FSMV)},warning=function(cond){message("webchecm could not fetch the info")})
                ##FSMV2<-tryCatch({rinchi::get.inchi.key(FSMV1)},warning=function(cond){message("webchecm could not fetch the info")})
                ##########################################################
		            if(!sjmisc::is_empty(FSMV2)){
		##################################################		    
			    print("enter the line ...868")	    
                  outn1<-Ikfilter(FSMV2,lst2,RRV,AIN,fmass,FRTL1)
                  len2<-length(outn1)
                  ####################################################
                  if(len2 > 1)
                  {
                    print("enter the line ...878")
                    cat(sapply(outn1, toString), file=OUNA3, sep="\n",append=TRUE)
                    cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
                    
                  ##} # end of if loop
                }else{
		  print("enter the line ...880")	
                  FINKE1<-NFFilter(RRV,AIN,lst2,fmass,FRTL1)
                  ### adding this code new 
                  len3<-length(FINKE1)
                  if(len3 > 1){
                    print("enter the line ...888")
                    cat(sapply(FINKE1, toString), file=OUNA3, sep="\n",append=TRUE)
                    cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
                  }else{
			  ### adding this new
			  DN<-dirname(dirname(File1))
			  DN1<-paste(DN,"Error-Report",sep="/")
			  FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
			  write.table(unname(as.data.frame(RRV)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)

		  }
		  ####################################
                } # end of the else loop
	       
	      }else{
                  print("enter the line ...1780")
                  FINKE1<-NFFilter(RRV,AIN,lst2,fmass,FRTL1)
                  ### adding this code new 
                  len3<-length(FINKE1)
                  if(len3 > 1){
                    print("enter the line ...873")
                    cat(sapply(FINKE1, toString), file=OUNA3, sep="\n",append=TRUE)
                    cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
                  }else{
			  ### adding this new for error report
			  DN<-dirname(dirname(File1))
			  DN1<-paste(DN,"Error-Report",sep="/")
			  FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
			  write.table(unname(as.data.frame(RRV)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
		  }
		  ###################################

                }

              ###############################################################################################################################
              }else if(!sjmisc::is_empty(stringr::str_trim(as.character(RRV[["InChI"]]))) & startsWith(as.character(RRV[["InChI"]]),'CAS:')){
                print("enter the line ...884")
	        ###################################
                CV<-stringr::str_trim(as.character(RRV[["InChI"]]))
                CV1<-stringr::str_replace(CV,pattern='CAS:',replacement ="")
                CV2<-stringr::str_trim(as.character(CV1))
                ####################################
                ##FIV=tryCatch({aw_query(CV2, from = 'cas')[[1]]},warning=function(cond){message("webchecm could not fetch the info")})
		####################################
		FIV=tryCatch({webchem::aw_query(CV2, from = 'cas')},error=function(cond){message("some mistake happened in cas finding")})
		FIV1=tryCatch({FIV[[1]]},error=function(cond){message("some mistake happened in file search files")})
		FIV2=tryCatch({FIV1$inchikey},error=function(cond){message("some mistake happened in file search files")})
                ###################################
		##if(!sjmisc::is_empty(FIV2))
                ##if(!is.na(FIV))
                ##{
		##################################
                  ##FIV1<-FIV2
		#################################
		  if(!sjmisc::is_empty(FIV2)){
                  ##if(!is.na(FIV1)){
                    outn3<-Ikfilter(FIV2,lst2,RRV,AIN,fmass,FRTL1)
                    len4<-length(outn3)
                    if(len4 > 1)
                    {
                      print("enter the line ...912")
                      cat(sapply(outn3, toString), file=OUNA3, sep="\n",append=TRUE)
                      cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
                      
                    ##} # end of len4
                  }else{
                    FSMIL<-NFFilter(RRV,AIN,lst2,fmass,FRTL1)
                    len3<-length(FSMIL)
                    if(len3 > 1){
                      print("enter the line ...921")
                      cat(sapply(FSMIL, toString), file=OUNA3, sep="\n",append=TRUE)
                      cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
                    }else{
			    DN<-dirname(dirname(File1))
			    DN1<-paste(DN,"Error-Report",sep="/")
			    FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
			    write.table(unname(as.data.frame(RRV)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
		    }
		   ##########################
                  } # end of else loop
		  
		}else{
                  
                  FCAS<-NFFilter(RRV,AIN,lst2,fmass,FRTL1)
                  ### adding this code new 
                  len4<-length(FCAS)
                  if(len4 > 1){
                    print("enter the line ...935")
                    cat(sapply(FCAS, toString), file=OUNA3, sep="\n",append=TRUE)
                    cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
                  }else{
			  DN<-dirname(dirname(File1))
			  DN1<-paste(DN,"Error-Report",sep="/")
			  FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
			  write.table(unname(as.data.frame(RRV)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)

		  }
		  #############################
                }
              #########################################################################################################################
              }else if(!sjmisc::is_empty(as.character(RRV[["SMILES"]])) & !startsWith(as.character(RRV[["SMILES"]]),'not available')){
                print("enter the line ...915")
	        #####################################################
                F1SM<-stringr::str_trim(as.character(RRV[["SMILES"]]))
                F1SM1<-tryCatch({rinchi::get.inchi.key(F1SM)},error=function(cond){message("webchecm could not fetch the info")})
                ######################################################
		            if(!sjmisc::is_empty(F1SM1)){
                  outn4<-Ikfilter(F1SM1,lst2,RRV,AIN,fmass,FRTL1)
                  loutn5<-length(outn4)
                  ###################################################
                  if(loutn5 > 1)
                  {
                    print("enter the line ...954")
                    cat(sapply(outn4, toString), file=OUNA3, sep="\n",append=TRUE)
                    cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
                    
                  }else{
                    print("enter the line ...959")
                    FSMIL<-NFFilter(RRV,AIN,lst2,fmass,FRTL1)
                    len3<-length(FSMIL)
                    if(len3 > 1){
                      print("enter the line ...963")
                      cat(sapply(FSMIL, toString), file=OUNA3, sep="\n",append=TRUE)
                      cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
                    }else{
			    DN<-dirname(dirname(File1))
			    DN1<-paste(DN,"Error-Report",sep="/")
			    FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
			    write.table(unname(as.data.frame(RRV)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
		    }
		   ################### 
                  } ## end of else loop
	         
	         }else{
                  #print("enter the line ...970")
                  FSMIL<-NFFilter(RRV,AIN,lst2,fmass,FRTL1)
                  len3<-length(FSMIL)
                  if(len3 > 1){
                    print("enter the line ...974")
                    cat(sapply(FSMIL, toString), file=OUNA3, sep="\n",append=TRUE)
                    cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
                  }else{
			  DN<-dirname(dirname(File1))
			  DN1<-paste(DN,"Error-Report",sep="/")
			  FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
			  write.table(unname(as.data.frame(RRV)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
		  }
		 ########################### 
                }
                
              #############################################################################################################################
              }else if(!sjmisc::is_empty(stringr::str_trim(as.character(RRV[["PubChem CID"]]))) & !startsWith(as.character(RRV[["PubChem CID"]]),'not available')){
                print("enter the line ...939")
                FPUCID<-stringr::str_trim(as.character(RRV[["PubChem CID"]]))
		####################################
		print("entering the PubchemID test area")
		print(FPUCID)
		######################################################################
		####################################################################
                FINSM<-tryCatch({webchem::cs_convert(as.numeric(FPUCID),from = "csid", to = "smiles")},error=function(cond){message("webchecm could not fetch the info")})
		FIINK<-tryCatch({rinchi::get.inchi.key(FINSM)},error=function(cond){message("webchecm could not fetch the info")})
                ##################################################################
		  ##################################################################
		              if(!sjmisc::is_empty(FIINK)){
                    outn6<-Ikfilter(FIINK,lst2,RRV,AIN,fmass,FRTL1)
                    len7<-length(outn6)
                    if(len7 > 1)
                    {
                      print("enter the line ...994")
                      cat(sapply(outn6, toString), file=OUNA3, sep="\n",append=TRUE)
                      cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
                      
                    }else{
                      print("enter the line ...999")
                      FCIDR<-NFFilter(RRV,AIN,lst2,fmass,FRTL1)
                      len3<-length(FCIDR)
                      if(len3 > 1){
                        print("enter the line ...1003")
                        cat(sapply(FCIDR, toString), file=OUNA3, sep="\n",append=TRUE)
                        cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
                      }else{
			      DN<-dirname(dirname(File1))
			      DN1<-paste(DN,"Error-Report",sep="/")
			      FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
			      write.table(unname(as.data.frame(RRV)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
		      }
		      ###############################
                    } # end of else loop
		   }else{
                    #print("enter the line ...1011")
                    FSMIL<-NFFilter(RRV,AIN,lst2,fmass,FRTL1)
                    len3<-length(FSMIL)
                    if(len3 > 1){
                      print("enter the line ...1014")
                      cat(sapply(FSMIL, toString), file=OUNA3, sep="\n",append=TRUE)
                      cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
                    }else{
			    DN<-dirname(dirname(File1))
			    DN1<-paste(DN,"Error-Report",sep="/")
			    FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
			    write.table(unname(as.data.frame(RRV)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
		    }
		##################################################### 
                  }# end of else loop
                ##############################################################################
                ##}else{                                                                    ##
                  #print("enter the line ...1037")					    ##
                 ## FSMIL<-NFFilter(RRV,AIN,lst2,fmass,FRTL1)				    ##
                 ## len3<-length(FSMIL)							    ##
                  ##if(len3 > 1){							    ##
                   ## print("enter the line ...1014")					    ##
                   ## cat(sapply(FSMIL, toString), file=OUNA3, sep="\n",append=TRUE)	    ##
                    ##cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)           ##
                  ##} ### end of len3							    ##
                ##}									    ##
	       ##################################	                                    ##
              ##}									    ##
              ################################################################################
	      }else{
                PASS1 <-RRV
              }
              
	       }
        } ### big else loop
      } ## for loop NRXF3
    } ## End ...NRXF3
 #### adding this temporary### else part
  }
	 
###################################
} ## is_empty(Val and Val1)
}
