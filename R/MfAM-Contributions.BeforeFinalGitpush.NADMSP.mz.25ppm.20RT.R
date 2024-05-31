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
library(magrittr)
library("metaMS")
library("rcellminer")
library("R.utils")
library(RCurl)
library(Rdisop)
#########################################################
##BiocManager::install("RMassBank")
###BiocManager::install("metaMS")
##BiocManager::install("rcellminerData")
##BiocManager::install("rcellminer")
##BiocManager::install("zlibbioc")
##BiocManager::install("ChemmineOB")
##BiocManager::install("Rdisop")
##library(Rdisop)
#########################################################
##devtools::install_github("r-lib/xml2")
##install.packages("RcppEigen",dependencies = TRUE)
##install.packages('Cairo',dependencies = TRUE)
##install.packages("sjmisc")
########## Step1: Reading the API key ####################################
Sys.setenv(CHEMSPIDER_KEY = "qeTvHEEPlYAefWqUpv4dJGG8w1UuxV5G")
rr_auth("qeTvHEEPlYAefWqUpv4dJGG8w1UuxV5G")
apikey="qeTvHEEPlYAefWqUpv4dJGG8w1UuxV5G"
##########################################################################
args <- commandArgs(TRUE)
File1<-args[1]
##########################################################################
##mz.25ppm.20RT.R
##########################################################################
##File1<-"/mnt/ifs/data/IPB/Projects/2017_005_MS-databases/mFam contributions/maximilian_frey@uni-hohenheim.de/meta data/230181218_Frey_Compound_Spreadsheet_For_MSMS_v17_GB.xlsx"
##########################################################################
##print(File1)
##########################################################################
Fi <- unlist(strsplit(File1, "/"))
##Fi1 <-c(paste(Fi[-length(Fi)], collapse = "/"), last(Fi))
Fi1 <-c(paste(Fi[-length(Fi)], collapse = "/"), data.table::last(Fi))
Fi2<-Fi1[1]
Fi3<-unlist(strsplit(Fi2, "/"))
##Fi4<-c(paste(Fi3[-length(Fi3)], collapse = "/"), last(Fi3))
Fi4<-c(paste(Fi3[-length(Fi3)], collapse = "/"), data.table::last(Fi3))
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
##NFi6<-c(paste(NFi5[-length(NFi5)], collapse = "/"), last(NFi5))
NFi6<-c(paste(NFi5[-length(NFi5)], collapse = "/"), data.table::last(NFi5))
N1Fi6<-NFi6[1]
################################################################
##print(Fi6)
################################################################
####### Step2: Reading the All the required files
################################################################
AIN1<-read.table("/mnt/ifs/data/IPB/Projects/2017_005_MS-databases/mFam contributions/scripts/ADI.4.csv",sep=",",header=F,quote="",stringsAsFactors = FALSE)
AIN<-AIN1
####################### adding the new database table###########
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
#########################################################

#########################################################
######### Adding this new################################
Lmeda<-tryCatch({base::rle(as.character(RXF3[["File"]]))},error=function(cond){message("file structure is wrong")})
NLmeda<-tryCatch({Lmeda$lengths},error=function(cond){message("repeat length calculate mistakes")})
LmeCmu<-tryCatch({abs(cumsum(NLmeda))},error=function(cond){message("cumulative calculate mistakes happening")})
LmeCmu1<-tryCatch({c(0,LmeCmu)},error=function(cond){message("cumulative calculate mistakes happening")})
SFileNam<-tryCatch({Lmeda$values},error=function(cond){message("rle is not able to fetch the information properly")})
#########################################################
#########################################################
##### REST API Functions
##########################################################
##########################################################
PuInKtoSM<-function(getINK)
{
  url<- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"
  out<-tryCatch({jsonlite::fromJSON(paste0(url,getINK, "/JSON"))} ,error = function(x) {return(NA)})
  prop.names  <-tryCatch({out$PC_Compounds$props[[1]][[1]]},error = function(x) {return(NA)})
  prop.values <- tryCatch({out$PC_Compounds$props[[1]][[2]]},error = function(x) {return(NA)})
  sm <-tryCatch({grep("smiles", prop.names[,"label"], ignore.case = TRUE)},error = function(x) {return(NA)})
  csmiles<-c()
  if(length(sm) >= 1) {
    can <- tryCatch({grep("canonical", prop.names[,"name"], ignore.case = TRUE)},error= function(x) {return(NA)})
    can1<-tryCatch({prop.values[sm[1],"sval"]},warning= function(x) {return(NA)})
    csmiles<-c(csmiles,can1)
  }else{
    csmiles<-c(csmiles,NA)
  }
  return(csmiles)
  
}
############################################################
############################################################
PuInKtoSM1<-function(getINK)
{
  ###### This return Isomeric smiles	
  url<- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"
  out<-tryCatch({jsonlite::fromJSON(paste0(url,getINK, "/JSON"))}, error = function(x) {return(NA)})
  prop.names  <-tryCatch({out$PC_Compounds$props[[1]][[1]]}, error = function(x) {return(NA)})
  prop.values <- tryCatch({out$PC_Compounds$props[[1]][[2]]}, error = function(x) {return(NA)})
  sm <-tryCatch({grep("smiles", prop.names[,"label"], ignore.case = TRUE)}, error = function(x) {return(NA)})
  csmiles<-c()
  if(length(sm) == 2) {
    can <- tryCatch({grep("canonical", prop.names[,"name"], ignore.case = TRUE)}, error = function(x) {return(NA)})
    can1<-tryCatch({prop.values[sm[2],"sval"]},warning= function(x) {return(NA)})
    csmiles<-c(csmiles,can1)
  }else{
    csmiles<-c(csmiles,NA)
  }
  return(csmiles)
  
}
##########################################################################
###########################################################################
PuCIDtoEM<-function(getCID)
{
  ######################
  #######################
  URL="https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/"
  URL1=paste0(URL,getCID, "/JSON/?response_type=display")
############################
  ########################
  data=tryCatch({jsonlite::fromJSON(URL1)} ,error = function(x) {return(NA)})
  ########################
  data1<-tibble::enframe(unlist(data))
  data2<-as.data.frame(data1)
  ##########################
  inVa<-which(data2$name %in% "Record.Section.Section.Section.Information.Value.StringWithMarkup.String")
  ##########################
  TEST<-sapply(data2, "[", inVa)
  TEST1<-grep("InChI=",TEST)
  TEST2<-TEST1[1]
  InchI<-TEST[TEST2]
  inchikey<-TEST[TEST2+1]
  #########################
  url<- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"
  #########################
  out<-tryCatch({jsonlite::fromJSON(paste0(url,inchikey, "/JSON"))} ,error = function(x) {return(NA)})
  #######################
  EMV<-out$PC_Compounds$props[[1]][22,]$value$sval
  #######################
  EMV1<-c()
  ######################
  if(!sjmisc::is_empty(EMV))
  {
    EMV1<-c(EMV1,EMV)
  }else{
    EMV1<-c(EMV1,0)
  }
  #####################
  return(EMV1)
  ####################
}
###########################################################################
###########################################################################
PuCAStoOI<-function(getCAS)
{
  ###CAS: 328-50-7
  getCAS1<-stringr::str_replace(getCAS,pattern='CAS:',replacement ="")
  getCAS2<-stringr::str_trim(getCAS1)
  url<- "http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
  out<-tryCatch({jsonlite::fromJSON(paste0(url,getCAS2, "/property/InChIKey"))}, warning = function(x) {return(NA)})
  out1<-tryCatch({jsonlite::fromJSON(paste0(url,getCAS2, "/property/CanonicalSMILES"))}, warning = function(x) {return(NA)})
  out2<-tryCatch({jsonlite::fromJSON(paste0(url,getCAS2, "/property/InChI"))}, warning = function(x) {return(NA)})
  ###########################
  OIN<-tryCatch({out2$PropertyTable$Properties$InChI},error=function(cond){message("List value is empty")})
  OIK<-tryCatch({out$PropertyTable$Properties$InChIKey},error=function(cond){message("List value is empty")})
  OSM<-tryCatch({out1$PropertyTable$Properties$CanonicalSMILES},error=function(cond){message("List value is empty")})
  OCID<-tryCatch({out1$PropertyTable$Properties$CID},error=function(cond){message("List value is empty")})
  ############################
  ##returning Inchi,InchiKey,Smiles,CompoundID...in order
  ###############################
  return(c(OIN,OIK,OSM,OCID))
  ###############################
}
#########################################################################
#########################################################################
ConvPCIDtoOCN<-function(getPCID)
{
  url<- "https://www.metabolomicsworkbench.org/rest/compound/pubchem_cid/"
  out<-tryCatch({jsonlite::fromJSON(paste0(url,getPCID, "/all"))}, warning = function(x) {return(NA)})
  #########################
  OIK<-tryCatch({out$inchi_key},error=function(cond){message("Inchikey value is empty")})
  OSM<-tryCatch({out$smiles},error=function(cond){message("smiles value is empty")})
  OCID<-tryCatch({out$pubchem_cid},error=function(cond){message("Pubchem CID value is empty")})
  OEM<-tryCatch({out$exactmass},error=function(cond){message("Exact mass value is empty")})
  OFOR<-tryCatch({out$formula},error=function(cond){message("FORMULA value is empty")})
  ################################
  return(c(OIK,OSM,OCID,OEM,OFOR))
  ################################

}
##########################################################################
###########################################################################
ConvINKtoOID<-function(getINK)
{
  url<- "https://www.metabolomicsworkbench.org/rest/compound/inchi_key/"
  out<-tryCatch({jsonlite::fromJSON(paste0(url,getINK, "/all"))}, warning = function(x) {return(NA)})
  #########################
  OIK<-tryCatch({out$inchi_key},error=function(cond){message("Inchikey value is empty")})
  OSM<-tryCatch({out$smiles},error=function(cond){message("smiles value is empty")})
  OCID<-tryCatch({out$pubchem_cid},error=function(cond){message("Pubchem CID value is empty")})
  OEM<-tryCatch({out$exactmass},error=function(cond){message("Exact mass value is empty")})
  OFOR<-tryCatch({out$formula},error=function(cond){message("FORMULA value is empty")})
  ###########################
  return(c(OIK,OSM,OCID,OEM,OFOR))
  ###########################
  
}
##########################################################################
##########################################################################
ClassSmilesToOntolgy<-function(getSMILE)
{
  #####################
  url="https://gnps-structure.ucsd.edu/classyfire?smiles="
  url1=paste0(url,getSMILE)
  out=jsonlite::fromJSON(url1)
  res=do.call(paste, c(as.list(tryCatch({rev(out$ancestors)},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
  ###########
  return(res)
  ############
}
##########################################################################
##########################################################################
FuFtoRe<-function(InMEDA)
{
  if(!sjmisc::is_empty(as.character(InMEDA[["Name"]])) & !startsWith(as.character(InMEDA[["Name"]]),'not available')){
    GCID<-tryCatch({webchem::get_cid(InMEDA[["Name"]])},error=function(cond){message("name is empty")})
    GCID1<-tryCatch({GCID[[2]][1]},error=function(cond){message("name is empty")})
    PCID1<-tryCatch({webchem::pc_prop(as.numeric(GCID1), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem CId is empty..did not get exact mass")})
    PCID2<-tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){message("Pubchem CID is empty..did not get exact mass")})
    ############################
    FMa<-c()
    ############################
    if(!sjmisc::is_empty(PCID2)){
      FMa<-c(FMa,PCID2)
    }else{
      EM<-stringr::str_trim(as.character(InMEDA[["Exact mass"]]))
      if(!sjmisc::is_empty(EM)){
        FMa<-c(FMa,EM)
      }else{
        FMa<-c(FMa,0)
      }
    }
  }else{
    EM<-stringr::str_trim(as.character(InMEDA[["Exact mass"]]))
    ############################
    if(!sjmisc::is_empty(EM)){
      FMa<-c(FMa,EM)
    }else{
      FMa<-c(FMa,0)
    }
    ###########################     
  }
  #################
  return(FMa)
  #################
}
###########################################################################
###########################################################################
MaKE.ONT.REC<-function(InMEDA)
{
  out<-c()
  ################
  print("entering the ontology area")
  ################
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
      IKCRV<-tryCatch({classyfireR::get_classification(IK)},error=function(cond){message("Classifier could not fetch the information")})
      #######################################
      ##IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(IK)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
      ########################################
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
	######################################      
        IK1<-paste("INCHIKEY:",IK,sep=" ")
        tes2<-paste("Ontology:","",sep=" ")
        FINCH<-paste("INCHI:",IN,sep=" ")
        ##################################
        out<-c(out,tes2)
        out<-c(out,IK1)
        out<-c(out,FINCH)
	##################################
        
      }
    ####################################################################
    }else if(!sjmisc::is_empty(SM)){
      ##################################################################
      tes<-tryCatch({classyfireR::submit_query(label = 'query_test', input = SM, type = 'STRUCTURE')},error=function(cond){message("Classyfire is empty")})
      ##tes1<-do.call(paste, c(as.list(tryCatch({tes@classification$Classification},error=function(cond){message("Classifier could not fecth the information")})), sep = ","))
      tes1<-ClassSmilesToOntolgy(SM)
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
	############################################      
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
      ######################
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
      print("check... if this is the error area")
      ##print(IN)
      #######################################
      FINCH<-paste("INCHI:",IN,sep=" ")
      IK<-as.character(InMEDA[["InChI"]])
      IK1<-paste("INCHIKEY:",IK,sep=" ")
      ############################
      if(!sjmisc::is_empty(IK)){
        ########################
        print("enter the line 8")
        #############################
        IKCRV<-tryCatch({classyfireR::get_classification(IK)},error=function(cond){message("Classifier could not fetch the information")})
	#################################
        ##IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(IK)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
	#####################################
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
	  #################################	
          IK1<-paste("INCHIKEY:",IK,sep=" ")
          tes2<-paste("Ontology:","",sep=" ")
          FINCH<-paste("INCHI:",IN,sep=" ")
          ################
          out<-c(out,tes2)
          out<-c(out,IK1)
          out<-c(out,FINCH)
          #####################
        }## end of else
      }else if(!sjmisc::is_empty(SM)){
	#############################################      
        tes<-tryCatch({classyfireR::submit_query(label = 'query_test', input = SM, type = 'STRUCTURE')},error=function(cond){message("Classyfire is empty")})
        ##tes1<-do.call(paste, c(as.list(tryCatch({tes@classification$Classification},error=function(cond){message("Classifier could not fecth the information")})), sep = ","))
        tes1<-ClassSmilesToOntolgy(SM)
        ################################
        IK1<-paste("INCHIKEY:",IK,sep=" ")
        tes2<-paste("Ontology:",tes1,sep=" ")
        FINCH<-paste("INCHI:",IN,sep=" ")
        ################################
        out<-c(out,tes2)
        out<-c(out,IK1)
        out<-c(out,FINCH)
        #################################
        
      }else{
        F1ONT<-paste("Ontology:","",sep=" ")
        ##FINCH<-paste("INCHI:",IN,sep=" ")
        ##############################################
        if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){
	 ####################################################	
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
          IKCRV<-tryCatch({classyfireR::get_classification(IK)},error=function(cond){message("Classifier could not fetch the information")})
	  ################################
          ##IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(IK)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
	  #####################################
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
            ###############################
            if(!sjmisc::is_empty(F1SM1)){
	      ################################	    
              IK<-F1SM1
	      IKCRV<-tryCatch({classyfireR::get_classification(IK)},error=function(cond){message("Classifier could not fetch the information")})
	      ############################################
              ##IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(IK)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
	      ####################################
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
		##############################
              }##IKCRV
            }else{
              ## not able to get inchikey from smiles too check pubchemID
              if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(InMEDA[["PubChem CID"]]),'not available')){
                ###################
                FPUCID<-stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))
                FPUCID1<-as.numeric(FPUCID)
                FINSM<-tryCatch({webchem::pc_prop(FPUCID1)},error=function(cond){message("webchecm could not fetch the info")})
                FIINK<-tryCatch({FINSM$InChIKey},error=function(cond){message("webchecm could not fetch the info")})
                #############################
                IK<-FIINK
                #############################
                if(!sjmisc::is_empty(IK)){
                  ##########################
		  IKCRV<-tryCatch({classyfireR::get_classification(IK)},error=function(cond){message("Classifier could not fetch the information")})	
		#####################################
                  ##IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(IK)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
		######################################
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
                  ######################
                }## end of else
                
                
              }else{
                F1ONT<-paste("Ontology:","",sep=" ")
                ##FINCH<-paste("INCHI:",SM,sep=" ")
                ##############################################
                if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){
		  #####################################	
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
	      ##################################	    
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
	  ######################################	
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
    ###################################
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
      IKCRV<-tryCatch({classyfireR::get_classification(IK)},error=function(cond){message("Classifier could not fetch the information")})	    #################################  
      ##IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(IK)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
    ######################################
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
	#############################################      
        IK<-as.character(InMEDA[["SMILES"]])
        F1SM<-stringr::str_trim(as.character(InMEDA[["SMILES"]]))
        F1SM1<-tryCatch({rinchi::get.inchi.key(F1SM)},error=function(cond){message("webchecm could not fetch the info")})
        ##############################
        if(!sjmisc::is_empty(F1SM1)){
	  ############################	
          IK<-F1SM1
	  IKCRV<-tryCatch({classyfireR::get_classification(IK)},error=function(cond){message("Classifier could not fetch the information")})
	  ###################################
          ##IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(IK)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
	  ############################################
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
	    ############################################	  
            FPUCID<-stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))
            FPUCID1<-as.numeric(FPUCID)
            FINSM<-tryCatch({webchem::pc_prop(FPUCID1)},error=function(cond){message("webchecm could not fetch the info")})
            FIINK<-tryCatch({FINSM$InChIKey},error=function(cond){message("webchecm could not fetch the info")})
            #############################
            IK<-FIINK
            #############################
            if(!sjmisc::is_empty(IK)){
              ##########################
	      IKCRV<-tryCatch({classyfireR::get_classification(IK)},error=function(cond){message("Classifier could not fetch the information")})
              ################################	    
              ##IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(IK)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
	      #######################################
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
		###############################      
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
              ######################
            }## end of else
            
            
          }else{
	    ##############################################	  
            F1ONT<-paste("Ontology:","",sep=" ")
            ##FINCH<-paste("INCHI:",SM,sep=" ")
            ##############################################
            if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){
	      #################################################	    
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
	  ##########################################	
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
    ###################################
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
        ##################################
	IKCRV<-tryCatch({classyfireR::get_classification(IK)},error=function(cond){message("Classifier could not fetch the information")}) 
        ####################################     
        ##IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(IK)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
        ####################################
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
	  #########################
        }else{
	  #################################	
          IK1<-paste("INCHIKEY:",IK,sep=" ")
          tes2<-paste("Ontology:","",sep=" ")
          FINCH<-paste("INCHI:",IN,sep=" ")
          out<-c(out,tes2)
          out<-c(out,IK1)
          out<-c(out,FINCH)
	  #################################
        }###else
      }else if(!sjmisc::is_empty(gSMI)){
        ##################################################
        print("enter the line 184")
        ##################################################
        tes<-tryCatch({classyfireR::submit_query(label = 'query_test', input = gSMI, type = 'STRUCTURE')},error=function(cond){message("adduct value is missing")})
        ##tes1<-do.call(paste, c(as.list(tryCatch({tes@classification$Classification},error=function(cond){message("Classifier could not fecth the information")})), sep = ","))
        tes1<-ClassSmilesToOntolgy(SM)
        ###############################################
        tes2<-paste("Ontology:",tes1,sep=" ")
        FINCH<-paste("INCHI:",IN,sep=" ")
        #######################
        out<-c(out,tes2)
        out<-c(out,IK1)
        out<-c(out,FINCH)
        #################
      }else{
	###############################################################      
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
	  ##################################	
          FINCH<-paste("INCHI:","",sep=" ")
          IK<-paste("INCHIKEY:",as.character(InMEDA[["InChI"]]),sep=" ")
          #################
          out<-c(out,F1ONT)
          out<-c(out,IK)
          out<-c(out,FINCH)
          ####################
        }##else..if
      ###############################################################  
      }###else...else if ..if
    }else{
      ######check the smiles is not empty
      ###PCSM<-SMV
      SMV<-PCSM
      if(!sjmisc::is_empty(SMV)){
	#############################      
        IK<-tryCatch({rinchi::get.inchi.key(SMV)},error=function(cond){message("rinchi could not fetch inchikey missing")})
        SMV1<-tryCatch({rinchi::get.inchi(SMV)},error=function(cond){message("rinchi could not fetch inchi missing")})
        IK1<-paste("INCHIKEY:",IK,sep=" ")
        FINCH<-paste("INCHI:",SMV1,sep=" ")
        #########################
        if(!sjmisc::is_empty(IK)){
	  #############################
	  IKCRV<-tryCatch({classyfireR::get_classification(IK)},error=function(cond){message("Classifier could not fetch the information")})	        ####################################
          ##IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(IK)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
	  ##################################################
          if(!sjmisc::is_empty(IKCRV)){
            ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
            ###########################
            IK1<-paste("INCHIKEY:",IK,sep=" ")
            tes2<-paste("Ontology:",ONTV,sep=" ")
            ###########################
            out<-c(out,tes2)
            out<-c(out,IK1)
            out<-c(out,FINCH)
	    ##################
            
          }else{
	    ########################	  
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
	##########################################	
          tes<-tryCatch({classyfireR::submit_query(label = 'query_test', input = SMV, type = 'STRUCTURE')},error=function(cond){message("Classfire not able to fetch empty")})
          tes1<-ClassSmilesToOntolgy(SM)
          ##tes1<-do.call(paste, c(as.list(tryCatch({tes@classification$Classification},error=function(cond){message("Classifier could not fecth the information")})), sep = ","))
          tes2<-paste("Ontology:",tes1,sep=" ")
          teIK2<-SMV1
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
            ####################
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
#########################################################################################
#########################################################################################
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
      ###################################
    }###end of else ..else if ..if 
  }### end of the for loop
  return(paste(tes1,collapse = ""))
}

#########################################################################################
#########################################################################################
comwithNu<-function(elem)
{
  tes1<-c()
  ##########################
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
  #########################
  }## end of if
  ###############################
  return(paste(tes1,collapse = ""))
  ##############################
}

#########################################################################################
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
  ######## Precursormz value #########################################################################
  fmass<-c()
  for(i in c(1:length(PMZL1))){
    LV1<-tryCatch({lst2[[i]]},warning=function(cond){message("list is empty")})
    S1R<-tryCatch({stringi::stri_startswith_fixed(LV1, 'PRECURSORMZ:') },warning=function(cond){message("PRECURSORMZ is empty")})
    S1R1<-tryCatch({which(S1R)},warning=function(cond){message("Index is empty")})
    S1R2<-tryCatch({LV1[S1R1]},warning=function(cond){message("Element is empty")})
    PV1<-tryCatch({stringr::str_trim(stringr::str_replace(S1R2, "PRECURSORMZ:", ""))},warning=function(cond){message("PRECURSORMZ empty")})
    PV2<-tryCatch({as.numeric(PV1)},warning=function(cond){message("PRECURSORMZ does not exists")})
    fmass<-c(fmass,PV2)
  }
  #######################################################################################################
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
  ### return(list(lst2, PMZL1, PTYL1,faddu,lst3,lst4,lst5,fmass,lst7,lst8,AAMZV,RTL,IND,FRTL,FRTL1))
  #####################################
  return(list(lst2, fmass,FRTL1))
  
}

##############################################################################################
##############################################################################################

##########################################################################################################
#########################################################################################################
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
      IK1<-tryCatch({webchem::get_cid(IK, from = "inchikey")},error=function(cond){message("Inchi name must be empty or rinchi not abe to fetch")})
      PCID<-tryCatch({IK1[[2]][1]},error=function(cond){message("Pubchem Id is empty")})
      ########################################################################
      PCID1<-tryCatch({webchem::pc_prop(as.numeric(PCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem Id is empty")})
      PCID2<-tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){message("Pubchem Id is empty")})
      ##########################################################################
      ##########################################################################
      if(!sjmisc::is_empty(PCID2))
      {
        #####################################      
        print("enter the if loop...inkikey")
        ####################################
        FMa<-c(FMa,PCID2)
      }else{
        ###########################################################      
        print("enter the else part")      
        ###########################################################
        if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),'InChI=')){
          ##############################################
          print("enter the Inchi part ...in else loop...Inchikey")
          ##############################################
          IK1<-tryCatch({webchem::get_cid(IK, from = "inchi")},error=function(cond){message("Inchi name must be empty or rinchi not abe to fetch")})
          PCID<-tryCatch({IK1[[2]][1]},error=function(cond){message("Pubchem Id is empty")})
          #######################################
          #######################################
          PCID1<-tryCatch({webchem::pc_prop(as.numeric(PCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem Id is empty")})
          PCID2<-tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){message("Pubchem Id is empty")})
          ###################################
          if(!sjmisc::is_empty(PCID2)){
            ###########################################	  
            print("enter the if loop ...inchi area")
            ###########################################	  
            FMa<-c(FMa,PCID2)
          }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["InChI"]]))) & startsWith(as.character(InMEDA[["InChI"]]),'CAS:')){
            ##############################################
            print("enter the Inchi part ...in else loop...CAS area")
            ###############################################
            CV<-stringr::str_trim(as.character(InMEDA[["InChI"]]))
            CV1<-stringr::str_replace(CV,pattern='CAS:',replacement ="")
            CV2<-stringr::str_trim(as.character(CV1))
            #############################################
            PCID<-tryCatch({webchem::get_cid(CV2, from = "xref/rn",match="first")},error=function(cond){message("Pubchem Id is empty")})
            ##########################################
            ##########################################
            PCID1<-tryCatch({webchem::pc_prop(as.numeric(PCID[[2]][1]), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES","InChI","InChIKey"))},error=function(cond){message("Pubchem Id is empty")})
            ##PCID1<-tryCatch({webchem::pc_prop(as.numeric(PCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem Id is empty")})
            PCID2<-tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){message("Pubchem Id is empty")})
            #############################################
            #############################################
            if(!sjmisc::is_empty(PCID2)){
              ########################################	    
              print("enter the if loop ...cas area")
              ########################################	    
              FMa<-c(FMa,PCID2)
            }else if(!sjmisc::is_empty(as.character(InMEDA[["SMILES"]])) & !startsWith(as.character(InMEDA[["SMILES"]]),'not available')){
              ######################################
              print("entering smiles area in InchiKey")
              ######################################
              IK<-as.character(InMEDA[["SMILES"]])
              #######################################
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
                #####################################
                PCID1<-tryCatch({webchem::pc_prop(as.numeric(FPUCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem CId is empty..did not get exact mass")})
                PCID2<-tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){message("Pubchem CID is empty..did not get exact mass")})
                ######################################
                ######################################
                if(!sjmisc::is_empty(PCID2)){
                  ############################################	
                  print("entering the if loop..Pubchem CID")
                  ###########################################
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
                    
                    #############################
                    PMA<-tryCatch({ConvPCIDtoOCN(FPUCID)},error=function(cond){message("Pubchem value is empty")})
                    FMa<-c(FMa,ifelse(!sjmisc::is_empty(tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")})),tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")}),0))
                    ########################################
                    ########################################
                  }
               ######################################
                }else{
                  ## formula not found and --smiles and pubchem failed to get PubchemID
                  ##FMa<-c(FMa,PCID2)
                  ############################
                  PMA<-FuFtoRe(InMEDA)
                  FMa<-c(FMa,PMA)
                  ############################
                }
	      ##################################	
              }else{
                ### PubchemId is not found ..entering the else loop
                if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){
                  ############################################
                  EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem Id is empty")})
                  EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
                  ############################################
                  ############################################
                  if(!sjmisc::is_empty(EM1)){
                    
			  FMa<-c(FMa,EM1)
                  }else{
                    ##print("formula area empty.....")
                    ##print(PCID2)
                    ##FMa<-c(FMa,PCID2)
                    #############################
                    PMA<-FuFtoRe(InMEDA)
                    FMa<-c(FMa,PMA)
                    #############################
                  }
                }
             ###############################
              }
            ####################  PCID2 ...end ###
            }else{
              #####################################
              ## smiles end ... not found
              #####################################
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
                    
                    ##################################
                    ### adding this new code
                    ########################
                    PMA<-tryCatch({ConvPCIDtoOCN(FPUCID)},error=function(cond){message("Pubchem value is empty")})
                    FMa<-c(FMa,ifelse(!sjmisc::is_empty(tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")})),tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")}),0))
                    #####################################
                    ####################################
                  }
                  ###################################
                }else{
                  
                  ##FMa<-c(FMa,PCID2)
                  ###############################
                  PMA<-FuFtoRe(InMEDA)
                  FMa<-c(FMa,PMA)
                  #################################
                }
                #########################################
              }else{
                ###### PubchemID is not found... so getting exact mass from
                if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){
                  #########################################
                  EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem Id is empty")})
                  EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
                  ##########################################
                  ##########################################
                  if(!sjmisc::is_empty(EM1)){
                    ##print("formula area.....")
                    ##print(EM1)
                    FMa<-c(FMa,EM1)
                  }else{
                    ##FMa<-c(FMa,PCID2)
                    #############################
                    PMA<-FuFtoRe(InMEDA)
                    FMa<-c(FMa,PMA)
                    ##############################
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
            ###########################################
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
              #####################################
              ####################################
              if(!sjmisc::is_empty(PCID2)){
                FMa<-c(FMa,PCID2)
              }else if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){
                ##################################
                EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem Id is empty")})
                EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
                ##################################
                if(!sjmisc::is_empty(EM1)){
                  FMa<-c(FMa,EM1)
                }else{
                  ######################
                  ##PMA<-PuCIDtoEM(as.numeric(FPUCID))
                  ##FMa<-c(FMa,PMA)
                  ######################
                  PMA<-tryCatch({ConvPCIDtoOCN(FPUCID)},error=function(cond){message("Pubchem value is empty")})
                  FMa<-c(FMa,ifelse(!sjmisc::is_empty(tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")})),tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")}),0))
                  #############################
                  #############################
                }
                ################################
              }else{
                ## formula not found and --smiles and pubchem failed to get PubchemID
                ##FMa<-c(FMa,PCID2)
                ##############################
                PMA<-FuFtoRe(InMEDA)
                FMa<-c(FMa,PMA)
                ###############################
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
                  ##FMa<-c(FMa,PCID2)
                  #########################
                  PMA<-FuFtoRe(InMEDA)
                  FMa<-c(FMa,PMA)
                  #############################
                }
              }
              ####################################
            }
            ####################  PCID2 ...end ###
          }else{
            ## smiles end ... not found
            ############################################################
            if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(InMEDA[["PubChem CID"]]),'not available')){
              ##########################################################
              FPUCID<-stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))
              ##########################################################
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
                  ### Pubchem CID is there ...adding this new
                  ##FMa<-c(FMa,PCID2)
                  ####################
                  ##PMA<-PuCIDtoEM(as.numeric(FPUCID))
                  ##FMa<-c(FMa,PMA) 
                  ### adding the new code here
                  ############################
                  PMA<-tryCatch({ConvPCIDtoOCN(FPUCID)},error=function(cond){message("Pubchem value is empty")})
                  FMa<-c(FMa,ifelse(!sjmisc::is_empty(tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")})),tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")}),0))

                  #################################
                  ##################################
                }
                #######################################
              }else{
                ### Trying to get the exact mass from name provided
                ###############################
                PMA<-FuFtoRe(InMEDA)
                FMa<-c(FMa,PMA)
                ##################################
                
              }
              ####################################
            }else{
              ###### PubchemID is not found... so getting exact mass from
              if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){
                #######################################
                EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem Id is empty")})
                EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
                ######################################
                if(!sjmisc::is_empty(EM1)){
                  ##print("formula area.....")
                  ##print(EM1)
                  FMa<-c(FMa,EM1)
                }else{
                  ##print("formula area empty.....")
                  ##print(PCID2)
                  ###FMa<-c(FMa,PCID2)
                  #####################
                  PMA<-FuFtoRe(InMEDA)
                  FMa<-c(FMa,PMA)
                  #######################
                }
                ###########################
              } ## formula
            } ## end of else loop
          } ## end of else
        }### InCHi.. end
      } ## end of main else
    } ### end of inchikey
  }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["InChI"]]))) & startsWith(as.character(InMEDA[["InChI"]]),'CAS:')){
    ############################################
    print("enter the cas function area.... cas is avilable")
    ############################################
    CV<-stringr::str_trim(as.character(InMEDA[["InChI"]]))
    CV1<-stringr::str_replace(CV,pattern='CAS:',replacement ="")
    CV2<-stringr::str_trim(as.character(CV1))
    ############################################
    PCID<-tryCatch({webchem::get_cid(CV2, from = "xref/rn",match="first")},error=function(cond){message("Pubchem Id is empty")})
    #############################################
    PCID1<-tryCatch({webchem::pc_prop(as.numeric(PCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem Id is empty")})
    PCID2<-tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){message("Pubchem Id is empty")})
    ##############################################
    ##############################################
    if(!sjmisc::is_empty(PCID2)){
      
      FMa<-c(FMa,PCID2)
    }else if(!sjmisc::is_empty(as.character(InMEDA[["SMILES"]])) & !startsWith(as.character(InMEDA[["SMILES"]]),'not available')){
      ##########################################
      ##print("enter the smiles..in cas area")
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
        ###############################################################
        FPUCID<-stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))
        ############################################
        ############################################
        PCID1<-tryCatch({webchem::pc_prop(as.numeric(FPUCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem Id is empty")})
        PCID2<-tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){message("Pubchem Id is empty")})
        ###########################################
        if(!sjmisc::is_empty(PCID2)){
          
          FMa<-c(FMa,PCID2)
          
        }else if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){
          ########################################
          ##########################################
          EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem Id is empty")})
          EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
          #########################################
          #########################################
          if(!sjmisc::is_empty(EM1)){
            
            FMa<-c(FMa,EM1)
          }else{
            ### Pubchem CID is found
            ##FMa<-c(FMa,PCID2)
            ##############################
            ##PMA<-PuCIDtoEM(as.numeric(FPUCID))
            ##FMa<-c(FMa,PMA)
            ### adding the new code here 
            ###################################
            PMA<-tryCatch({ConvPCIDtoOCN(FPUCID)},error=function(cond){message("Pubchem value is empty")})
            FMa<-c(FMa,ifelse(!sjmisc::is_empty(tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")})),tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")}),0))
            ###################################
            ##################################
          }
          
        }else{
          ### formula not found ..exact mass
          ##FMa<-c(FMa,PCID2)
          ###################################
          PMA<-FuFtoRe(InMEDA)
          FMa<-c(FMa,PMA)
          #################################
        }
        #######################################
      }else{
        ##print("enter the else area ..2")
        #### PubchemID is not found
        ###################################################
        if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){
          ################################################
          EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem Id is empty")})
          EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
          ###########################################
          ###########################################
          if(!sjmisc::is_empty(EM1)){
            ##print("formula area.....")
            ##print(EM1)
            ###################
            FMa<-c(FMa,EM1)
            ####################
          }else{
            ##print("formula area empty.....")
            ##print(PCID2)
            ##FMa<-c(FMa,PCID2)
            #####################
            PMA<-FuFtoRe(InMEDA)
            FMa<-c(FMa,PMA)
            ######################
          }
        }
        #######################################
      }
      ############## Smiles is empty and Pubchem CID
    }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(InMEDA[["PubChem CID"]]),'not available')){
      ########################################################
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
          
          ##FMa<-c(FMa,PCID2)
          ##################################
          ##PMA<-PuCIDtoEM(as.numeric(FPUCID))
          ##FMa<-c(FMa,PMA)
          ###################################
          PMA<-tryCatch({ConvPCIDtoOCN(FPUCID)},error=function(cond){message("Pubchem value is empty")})
          FMa<-c(FMa,ifelse(!sjmisc::is_empty(tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")})),tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")}),0))
          #################################
          ##################################
        }
        ####################################
        ###################################
        
      }else{
        print("entering Pubchem CID else loop---not got excat mass from pubchem CID")
        ##print(FPUCID)
        ##FMa<-c(FMa,PCID2)
        #######################
        PMA<-FuFtoRe(InMEDA)
        FMa<-c(FMa,PMA)
        ######################
        
      }
      #######################################
    }else{
      ##print("enter the else area ..pubcehm CID not found")
      #### PubchemID is not found
      if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){
        #########################################
        EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem Id is empty")})
        EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
        ########################################
        #######################################
        if(!sjmisc::is_empty(EM1)){
          
          FMa<-c(FMa,EM1)
        }else{
          ##print("formula area empty.....")
          ##print(PCID2)
          ##FMa<-c(FMa,PCID2)
          ###########################
          PMA<-FuFtoRe(InMEDA)
          FMa<-c(FMa,PMA)
          ###########################
          
        }
        ################################
      }else{
        ### print formula not found 
        ##FMa<-c(FMa,PCID2)
        ###############################
        PMA<-FuFtoRe(InMEDA)
        FMa<-c(FMa,PMA)
        ###############################
        
      }
    }###end of else loop 
  }else if(!sjmisc::is_empty(as.character(InMEDA[["SMILES"]])) & !startsWith(as.character(InMEDA[["SMILES"]]),'not available')){
    ########################################
    print("enter smiles area --2")
    #######################################
    IK<-as.character(InMEDA[["SMILES"]])
    IK1<-tryCatch({webchem::get_cid(IK, from = "smiles")},error=function(cond){message("smiles not abe to fetch")})
    PCID<-tryCatch({IK1[[2]][1]},error=function(cond){message("Pubchem Id is empty")})
    tes<-tryCatch({rcdk::parse.smiles(IK)},error=function(cond){message("smiles not abe to parse")})
    tes2<-tryCatch({rcdk::get.exact.mass(tes1)},error=function(cond){message("smiles not abe to fetch")})
    PCID2<-tryCatch({tes2},error=function(cond){message("smiles not abe to fetch")})
    #######################################
    #######################################
    if(!sjmisc::is_empty(PCID2)){
      print("enter the if loop ..smiles")
      
      FMa<-c(FMa,PCID2)
      
    }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(InMEDA[["PubChem CID"]]),'not available')){
      ############################################
      print("enter the pubchem CID ..smiles area ..meaning smiles are there and not able to get exact mass..pubchem CID is avilable")
      #############################################
      FPUCID<-stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))
      ###########################################
      ###########################################
      PCID1<-tryCatch({webchem::pc_prop(as.numeric(FPUCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem Id is empty")})
      PCID2<-tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){message("Pubchem Id is empty")})
      ##########################################
      ########################################
      if(!sjmisc::is_empty(PCID2)){
        FMa<-c(FMa,PCID2)
      }else if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){
        print("enter the if pubchem CID ..else if..")
        ##################################
        EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem Id is empty")})
        EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
        ##################################
        ##################################
        if(!sjmisc::is_empty(EM1)){
          FMa<-c(FMa,EM1)
        }else{
          print("enter the pubchem CID area ..else part in else if..that means formula is empty")
          
          ###################################
          PMA<-tryCatch({ConvPCIDtoOCN(FPUCID)},error=function(cond){message("Pubchem value is empty")})
          FMa<-c(FMa,ifelse(!sjmisc::is_empty(tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")})),tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")}),0))
          ###################################
          ##################################
        }
        ##################################
        ##################################
      }else{
        ### formula not aviable ...else for ..else if else
        print("enter the pubchem CID area...else ...that means neither formula ...nothing is avilable for Pubchem CID")
        ##FMa<-c(FMa,PCID2)
        #######################################
        ##PMA<-PuCIDtoEM(as.numeric(FPUCID))
        ##FMa<-c(FMa,PMA)
        ####################################
        PMA<-tryCatch({ConvPCIDtoOCN(FPUCID)},error=function(cond){message("Pubchem value is empty")})
        FMa<-c(FMa,ifelse(!sjmisc::is_empty(tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")})),tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")}),0))
        ######################################
        ######################################
      }
      ####################################
      
    }else{
      #### Pubchem CID #########################
      print("entering the else part---smiles ... will check formula...for exact mass")
      ##########################################
      if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){
        #################################
        EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem Id is empty")})
        EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
        ##################################
        
        if(!sjmisc::is_empty(EM1)){
          FMa<-c(FMa,EM1)
        }else{
          
          ##########################
          PMA<-FuFtoRe(InMEDA)
          FMa<-c(FMa,PMA)
          ##########################
        }
        
        
      }else{
        ### print formula not found 
        ##FMa<-c(FMa,PCID2)
        #######################
        PMA<-FuFtoRe(InMEDA)
        FMa<-c(FMa,PMA)
        #######################
      }
      ###########################################
    }
    ###########################################
  }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(InMEDA[["PubChem CID"]]),'not available')){
    ####################################
    print("entering the Pubchem CID area")
    #######################################
    FPUCID<-stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))
    ####################################
    PCID1<-tryCatch({webchem::pc_prop(as.numeric(FPUCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem Id is empty")})
    PCID2<-tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){message("Pubchem Id is empty")})
    ###################################
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
       
        ###################################
        ### changing this function too
        ##PMA<-PuCIDtoEM(as.numeric(FPUCID))
        #####################################
        PMA<-tryCatch({ConvPCIDtoOCN(FPUCID)},error=function(cond){message("Pubchem value is empty")})
        FMa<-c(FMa,ifelse(!sjmisc::is_empty(tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")})),tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")}),0))
        #############################
        ######################################
      }
      
    }else{
      print("entering the else of Pubchem CID")
      print("entering this new")
      ##############################
      ###########################
      PMA<-FuFtoRe(InMEDA)
      FMa<-c(FMa,PMA)
      ############################
      ############################
    }
    ################################
    
  }else{
    ###################################
    print("pubchem CID is  not available ...so checking formula")
    #################################	  
    if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){
      #################################
      EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem Id is empty")})
      EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
      #########################
      if(!sjmisc::is_empty(EM1)){
        print("this is entering the formula if loop...")      
        FMa<-c(FMa,EM1)
      }else{
        print("adding this new code here")
        if(!sjmisc::is_empty(as.character(InMEDA[["Name"]])) & !startsWith(as.character(InMEDA[["Name"]]),'not available')){
          GCID<-tryCatch({webchem::get_cid(InMEDA[["Name"]])},error=function(cond){message("name is empty")})
          GCID1<-tryCatch({GCID[[2]][1]},error=function(cond){message("name is empty")})
          PCID1<-tryCatch({webchem::pc_prop(as.numeric(GCID1), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem CId is empty..did not get exact mass")})
          PCID2<-tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){message("Pubchem CID is empty..did not get exact mass")})
          ########################
          FMa<-c(FMa,PCID2)
          ########################
        }else{
          ### adding this new part here 
          ###########################
          ###########################
          PMA<-FuFtoRe(InMEDA)
          FMa<-c(FMa,PMA)
          ##############################
          ##############################
        }### this is the end of else
       
      }
      print("this is entering the formula else...that means formula not avilable... in final else loop")
      
      #######################################
      #######################################
      PMA<-FuFtoRe(InMEDA)
      FMa<-c(FMa,PMA)
      ##########################################
      ##########################################
      
      
    }### this is the end of else...formula not avilable
    ##############################
  } ### end of the else
  ################################
  ################################
  return(FMa)
  ################################
  ################################
  
}
##########################################################################################################
##########################################################################################################

####################################################################################
####################################################################################
##print("enter the area before Ikfilter")
####################################################################################
gETSmiles<-function(InMEDA)
{
  ######################
  out<-c()
  #######################
  if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI="))
  {
    IN<-as.character(InMEDA[["InChI"]])
    mol <-tryCatch({rinchi::parse.inchi(IN)},error=function(cond){message("name is empty")})
    SM<-tryCatch({rcdk::get.smiles(mol[[1]])},error=function(cond){message("name is empty")})
    ########################
    ########################
    ### adding the new conditions here 
    if(!sjmisc::is_empty(SM))
    {
      out<-c(out,SM)
    }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["SMILES"]])))){
      SM<-as.character(InMEDA[["SMILES"]])
      out<-c(out,SM)
    }else{
      out<-c(out,"NA")
    }
  ################################
  }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["InChI"]]))) & startsWith(as.character(InMEDA[["InChI"]]),'CAS:')){
    ##########################	  
    CV<-stringr::str_trim(as.character(InMEDA[["InChI"]]))
    CV1<-stringr::str_replace(CV,pattern='CAS:',replacement ="")
    CV2<-stringr::str_trim(as.character(CV1))
    ##############################
    PCID<-tryCatch({webchem::get_cid(CV2, from = "xref/rn",match="first")},error=function(cond){message("Pubchem Id is empty")})
    PCID1<-tryCatch({webchem::pc_prop(as.numeric(PCID[[2]][1]), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES","InChI","InChIKey"))},error=function(cond){message("Pubchem Id is empty")})
    #############################
    gSMI<-tryCatch({PCID1$CanonicalSMILES},error=function(cond){message("smiles is not found")})
    IK<-tryCatch({PCID1$InChIKey},error=function(cond){message("some mistake happened in file search files")})
    IN<-tryCatch({PCID1$InChI},error=function(cond){message("some mistake happened in file search files")})
    ################
    ##out<-c(out,gSMI)
    ################
    if(!sjmisc::is_empty(gSMI))
    {
      out<-c(out,gSMI)
    }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["SMILES"]])))){
      ## Adding this new
      SM<-tryCatch({PuCAStoOI(CV2)},error=function(cond){message("CAS value is empty")})
      out<-c(out,ifelse(!sjmisc::is_empty(tryCatch({ConvPCIDtoOCN(PC)[2]},error=function(cond){message("Pubchem value is empty")})),tryCatch({ConvPCIDtoOCN(PC)[2]},error=function(cond){message("Pubchem value is empty")}),as.character(InMEDA[["SMILES"]])))
      ###########################
      ##out<-c(out,tryCatch({SM[3]},error=function(cond){message("CAS value is empty")}))
      ##### This is the old code
      ##SM<-as.character(InMEDA[["SMILES"]])
      ##out<-c(out,SM)
      ######################
    }else{
      ##out<-c(out,"NA")
      SM<-as.character(InMEDA[["SMILES"]])
      out<-c(out,SM)
    }
 ###################################################
  }else if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]]))){
    if(tryCatch({webchem::is.inchikey(as.character(InMEDA[["InChI"]]))},error=function(cond){message("inchikey validation failed")})){
      tes<-tryCatch({webchem::get_cid(stringr::str_trim(as.character(InMEDA[["InChI"]])), from = "inchikey")},error=function(cond){message("webchem not able to get cid from Inchikey")})
      tes1<-tryCatch({tes$cid},error=function(cond){message("Inchikey to CID did not convert")})
      tes2<-tryCatch({webchem::pc_prop(as.numeric(tes1[1]), properties = c("MolecularFormula", "MolecularWeight","CanonicalSMILES","InChI","InChIKey"))},error=function(cond){message("Inchikey to CID did not convert so did not get properties")})
      IN<-tryCatch({tes2$InChI},error=function(cond){message("Inchi to Inchikey failed because of CID not converting")})
      SM<-tryCatch({tes2$CanonicalSMILES},error=function(cond){message("Inchi to Inchikey failed because of CID not converting")})
      out<-c(out,SM)
    }else{
      if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["SMILES"]])))){
        
        ######This is old code
        SM<-as.character(InMEDA[["SMILES"]])
        out<-c(out,SM)
        ####################
      }else{
        print("entering the else part")
        ######################
        out<-c(out,"NA")
        ####################
      }##inside else part
    }###end of else inchikey checking
  }else if(!sjmisc::is_empty(as.character(InMEDA[["PubChem CID"]]))){
 ################################################################
    tes<-tryCatch({webchem::pc_prop(PCID1, properties = c("MolecularFormula", "MolecularWeight","CanonicalSMILES","InChI","InChIKey"))},error=function(cond){message("Did not get properties from Pubchem CID")})
    gSMI<-tryCatch({tes$CanonicalSMILES},error=function(cond){message("smiles is not found")})
    IK<-tryCatch({tes$InChIKey},error=function(cond){message("some mistake happened in file search files")})
    IN<-tryCatch({tes$InChI},error=function(cond){message("some mistake happened in file search files")})
    ##out<-c(out,gSMI)
    if(!sjmisc::is_empty(gSMI))
    {
      out<-c(out,gSMI)
    }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["SMILES"]])))){
      SM<-as.character(InMEDA[["SMILES"]])
      out<-c(out,SM)
    }else{
      ### This is old code I am adding new code
      ##out<-c(out,"NA")
      ####################################
      PC<-as.character(InMEDA[["PubChem CID"]])
      out<-c(out,ifelse(!sjmisc::is_empty(tryCatch({ConvPCIDtoOCN(PC)[2]},error=function(cond){message("Pubchem value is empty")})),tryCatch({ConvPCIDtoOCN(PC)[2]},error=function(cond){message("Pubchem value is empty")}),"NA"))
      ########################
    }
    
    
  }else if(!sjmisc::is_empty(as.character(InMEDA[["SMILES"]]))){
    SM<-as.character(InMEDA[["SMILES"]])
    out<-c(out,SM)
  }else{
    out<-c(out,"NA")
  }
  ############
  return(out)
  ############
}



#########################################################################################################
#########################################################################################################
Ikfilter <- function(InKeyVal,InMSPL,InMEDA,InAdVA,InPMZ,InRTL){
  ####################
  ##print("checking if problem is in Inchikey area")
  ###################
  out<-c()
  ###################################
  IINF<-tryCatch({webchem::cts_compinfo(InKeyVal)},error=function(cond){message("webchecm could not fetch the info")})
  ####################################
  ##if(length(IINF)> 0 & !is.na(IINF)){
  #####################################
  if(length(IINF)> 0 & !sjmisc::is_empty(IINF)){	  
    ##!sjmisc::is_empty	  
    ###########################
    IK<-tryCatch({IINF[[1]][1]},error=function(cond){message("Inchikey value is empty")})
    IN<-tryCatch({IINF[[1]][2]},error=function(cond){message("Inchi value is empty")})
    PMZ<-tryCatch({IINF[[1]][4]},error=function(cond){message("PrecursorMZ value is empty")})
    FM<-tryCatch({IINF[[1]][5]},error=function(cond){message("Formula value is empty")})
    ############################################
    ##SM<-tryCatch({PuInKtoSM(IK)},error=function(cond){message("inchikey to smiles conversion failed")})
    ##########################################
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
    ##InKeyVal<-IK
    InKeyVal<-tryCatch({IK$inchikey},error=function(cond){message("Inchikey is empty")})
    ############################
    ##if(!sjmisc::is_empty(AUIN) & !sjmisc::is_empty(PMZ)){
    if(!sjmisc::is_empty(PMZ)){
      ################################
      print("enter the function ...PMZ(Exact mass) and AUIN(Adduct) not empty ...Inchikey function")
      #if(length(AUIN) > 0 & !is.na(PMZ) ){
      ##########################
      AUIN1<-tryCatch({qdapRegex::ex_between(AUIN, "[", "]")[[1]]},error=function(cond){message("Adduct value is missing")})
      AUIN2<-tryCatch({FADINF(AUIN)},error=function(cond){message("adduct value matching is not found")})      
      ##AUIN2<-tryCatch({InAdVA[InAdVA$V1==AUIN1,]$V8},warning=function(cond){message("Adduct value is missing")})
      AAMS<-tryCatch({stringr::str_replace(AUIN2, "M",as.character(PMZ$exactmass))},error=function(cond){message("Missing adduct replacement")})
      AAMS1<-tryCatch({as.numeric(pander::evals(AAMS)[[1]]$result)},error=function(cond){message("Error in adduct replacement step")})
      #########################
      PPm=AAMS1*(mz_Tol/(1000000))
      ###PPm=AAMS1*(25/(1000000))
      #########################
      MPPmL=AAMS1-PPm
      MPPmU=AAMS1+PPm
      ###########################
      Tmass<-InPMZ[InPMZ >= MPPmL & InPMZ <= MPPmU]
      ITmass<-which(InPMZ %in% Tmass)
      ##ITmass<-match(Tmass,InPMZ)
      ###########################
      VRT<-as.numeric(as.character(InMEDA[["RT (min)"]]))
      VRTL<-VRT-RT_Tol
      VRTU<-VRT+RT_Tol
      #############################
      TRTL<-InRTL[InRTL >=VRTL & InRTL <= VRTU]
      ITRTL<-which(InRTL %in% TRTL)
      ###############################
      ##ITRTL<-match(TRTL,InRTL)
      ##################################
      ##print("enter temp test...function 1")
      ##print(PMZ$exactmass)
      ##print(MPPmL)
      ##print(MPPmU)
      ##print(VRT)
      ##print(VRTL)
      ##print(VRTU)
      ################################
      print("enter the line ...1195")
      if(length(ITRTL) >= 1){
        print("enter the line ...1197")
        if(length(ITmass) >= 1){
          print("enter the line ...1198")
	  ##############################
	  ##print(Tmass)
	  ##print(ITmass)
	  ##print(lst2[ITmass])
	  ##print(intersect(ITmass,ITRTL))
	  ##print(ITRTL)
          #############################
          INLL<-intersect(ITmass,ITRTL)
          #####################
          if(length(INLL) == 1){
	    ###################
            print("enter the line ...1199")
            ###################
            F1FPL<-InMSPL[INLL]
            ###################
            ##SM1<-as.character(InMEDA[["SMILES"]])
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
#             PT <- c("PRECURSORTYPE:", "ADDUCTIONNAME:")
#             NPT<-tryCatch({grep(paste(PT,collapse="|"), FNA, value=TRUE)},error=function(cond){message("List value is empty")})
#             NPT1<-tryCatch({match(NPT,FNA)},error=function(cond){message("List value is empty")})
#             PTV <-tryCatch({stringr::str_remove(FNA[NPT1],c("PRECURSORTYPE:","ADDUCTIONNAME:"))},error=function(cond){message("List value is empty")})
#             PTV1<-tryCatch({PTV[1]},error=function(cond){message("List value is empty")})
#             PTV2<-tryCatch({stringr::str_trim(PTV1)},error=function(cond){message("List value is empty")})
# 	    P1TV2<-tryCatch({stringr::str_trim(gsub("ADDUCTIONNAME:","",PTV2))},error=function(cond){message("List value is empty")})            
# 	    PTV3<-tryCatch({paste("PRECURSORTYPE:",P1TV2)},error=function(cond){message("List value is empty")})
#             ###########################################################
	    ####################################################################
	    ##print("enter my test...1")
	    ##print(FNA)
	    ##print(P1TV2)
	    ##print(as.character(InMEDA[["Adduct"]]))
	    ##print(AAMS1)
	    ##print(MPPmL)
	    ##print(MPPmU)
	    ##print(VRT)
	    ##print(VRTL)
	    ##print(VRTU)
      ####################################################################      
	    ################################################################# commenting this
	     ###   if(!sjmisc::is_empty(P1TV2) || !sjmisc::is_empty(as.character(InMEDA[["Adduct"]]))){
	    ####################################################################
	    ####################################################################         
            ##if(identical(P1TV2,as.character(InMEDA[["Adduct"]]))){
	      ##############################
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
              ###out<-c(out,PTV3)
              out<-c(out,F1PT1)
              #################################
              FIN1<-InMEDA[["Ionization mode"]]
              F1IN1<-as.character(FIN1)
              F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
              out<-c(out,F2IN1)
              ##################################################
	    
	      IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
	      ###################################################
	      ###IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(InKeyVal)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
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
              ############################################
              ##print("enter the line ...1878")
	      ############################################
              FINK<-paste("INCHIKEY:",tryCatch({IK$inchikey},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
	      ##FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
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
	      FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:","",sep=""))
              out<-c(out,FINS2)
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
              } ### end of else Fea
	   ###################################################
            ### adding this }
	   ####################################################
    #####################################################        
	      ##  }
    ######################################################
    #####################################################
            ### commenting all this 
#             else{
# 		    ###############################
# 		    print("enter the part 1..... P1TV2...Adduct")
# 		    ##########################
# 		    ##print(as.character(InMEDA[["Adduct"]]))
# 		    ##print(INLL)
# 		    ##print(FNA)
# 		    ##print(as.character(InMEDA[["Adduct"]]) == "[M]+")
# 		    ##"[M]+"
# 		    #########################################
# 		    #########################################
# 		    if(as.character(InMEDA[["Adduct"]]) == "[M]+" || as.character(InMEDA[["Adduct"]]) == "[M]-"){
# 		    ######################################	    
# 		    ##if(as.character(InMEDA[["Adduct"]]) == "[M]+"){
# 		    ###########################
# 		    print("enter the part 1..if loop")
# 		    ###############################
# 		    FNA1<-which(stringi::stri_detect_fixed(FNA,"NAME:"))
# 		    FNA2<-InMEDA[["Name"]]
# 		    FNA3<-as.character(FNA2)
# 		    ###########################
# 		    FNAM<-paste("NAME:",FNA3,sep=" ")
# 		    out<-c(out,FNAM)
# 		    ########################
# 		    FRA1<-which(stringi::stri_detect_fixed(FNA,"RETENTIONTIME:"))
# 		    F1RA1<-FNA[FRA1]
# 		    out<-c(out,F1RA1)
# 		    ########################
# 		    FMZ1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORMZ:"))
# 		    F1MZ1<-FNA[FMZ1]
# 		    out<-c(out,F1MZ1)
# 		    ########################
# 		    FPT1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORTYPE:"))
# 		    F1PT1<-FNA[FPT1]
# 		    out<-c(out,PTV3)
# 		    #################################
# 		    FIN1<-InMEDA[["Ionization mode"]]
# 		    F1IN1<-as.character(FIN1)
# 		    F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
# 		    out<-c(out,F2IN1)
# 		    ##################################################
# 		    IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
# 		    ###################################################
# 		    ##IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(InKeyVal)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
# 		    #################################################
# 		    if(!sjmisc::is_empty(IKCRV)){
#   			ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
#   			F1ONT<-paste("Ontology:",ONTV,sep=" ")
#   			out<-c(out,F1ONT)
# 		    }else{
#   			F1ONT<-paste("Ontology:","",sep=" ")
#   			out<-c(out,F1ONT)
# 		    }
# 
# 		    #############################################
# 		    #############################################
# 		    FINK<-paste("INCHIKEY:",tryCatch({IK$inchikey},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
# 		    out<-c(out,FINK)
# 		    FINCH<-paste("INCHI:",InchiV,sep=" ")
# 		    out<-c(out,FINCH)
# 		    FSIM<-paste("SMILES:",SM1,sep=" ")
# 		    out<-c(out,FSIM)
# 		    ############################
# 		    FFOR<-FM$formula
# 		    FFOR1<-paste("FORMULA:",FFOR,sep=" ")
# 		    out<-c(out,FFOR1)
# 		    #############################
# 		    FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
# 		    FINS1<-FNA[FINS]
# 		    FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:","",sep=""))
# 		    out<-c(out,FINS2)
# 		    ########################
# 		    FAUT<-as.character(InMEDA[["Authors"]])
# 		    FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
# 		    out<-c(out,FAUT1)
# 		    #######################
# 		    FLIC<-paste("LICENSE:","CC BY",sep=" ")
# 		    out<-c(out,FLIC)
# 		    #####################
# 		    FCIE<-paste("COLLISIONENERGY:",as.character(InMEDA[["Collision energy"]]),sep=" ")
# 		    out<-c(out,FCIE)
# 		    #########################
# 		    FINST<-as.character(InMEDA[["INSTRUMENT_TYPE"]])
# 		    FINST1<-paste("INSTRUMENTTYPE:",FINST,sep=" ")
# 		    out<-c(out,FINST1)
# 		    ########################
# 		    FINS<-as.character(InMEDA[["INSTRUMENT"]])
# 		    FINS1<-paste("INSTRUMENT:",FINS,sep=" ")
# 		    out<-c(out,FINS1)
# 		    ####################
# 		    FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
# 		    out<-c(out,FCOM)
# 		    ##################
# 		    FNPA<-which(stringi::stri_detect_fixed(FNA,"Num Peaks:"))
# 		    F1NPA<-FNA[FNPA]
# 		    out<-c(out,F1NPA)
# 		    ###################
# 		    Fpea<-FNA[(FNPA+1):Find]
# 		    ###########################
# 		    if(is.na(Fpea))
# 		    {
#   			Fpea1<-FNA[(FNPA+1)]
# 		    }else{
# 			    MV=AAMS1
# 			    tes1<-unlist(strsplit(Fpea, "\t"))
# 			    tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
# 			    tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
# 			    tes4<-which(tes2 > (3+MV))
# 			    if(length(tes4)>1)
# 			    {
# 			     	tes5<-tes2[-tes4]
# 			        tes6<-tes3[-tes4]
# 			        tes7<-paste(tes5,tes6,sep="\t")
# 				out<-c(out,tes7)
# 			    }else{
# 				out<-c(out,Fpea)
# 			    }
# 			   }    
# 			    ###end of Fea        
# 		   ########################################################
# 		     }  ### end of "[M]+"
# 	########################################################
#   #########################################################            
#        #### end of else ### THis is need to commented 
#             }
  ###############################################################          
  ###############################################################  
          } else if(length(INLL) > 1){
            #####################################		  
            print("entering the 1569")
	    print("you are choosing RT time ..Inchikey")
	    #####################################
            MONMS=InMSPL[INLL]
            TRA<-unname(rapply(MONMS, function(x) grep("RETENTIONTIME:",x, value=TRUE)))
            TRA1<-as.numeric(stringr::str_trim(stringr::str_replace(TRA, "RETENTIONTIME:", "")))
            TRA2<-abs(VRT-TRA1)
            TRA3<-which.min(TRA2)
            TRA4<-INLL[TRA3]
            TRA5<-InMSPL[TRA4]
            #######################
            F1FPL<-TRA5
            ######################
            ##SM1<-as.character(InMEDA[["SMILES"]])
            ######################
            #InMEDA[["SMILES"]]<-SM1
            #InMEDA[["PubChem CID"]]<-CID5
            #####################
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
#             PT <- c("PRECURSORTYPE:", "ADDUCTIONNAME:")
# 	    NPT<-tryCatch({grep(paste(PT,collapse="|"), FNA, value=TRUE)},error=function(cond){message("List value is empty")})
# 	    NPT1<-tryCatch({match(NPT,FNA)},error=function(cond){message("List value is empty")})
# 	    PTV <-tryCatch({stringr::str_remove(FNA[NPT1],c("PRECURSORTYPE:","ADDUCTIONNAME:"))},error=function(cond){message("List value is empty")})
# 	    PTV1<-tryCatch({PTV[1]},error=function(cond){message("List value is empty")})
#             PTV2<- tryCatch({stringr::str_trim(PTV1)},error=function(cond){message("List value is empty")})
# 	    P1TV2<-tryCatch({stringr::str_trim(gsub("ADDUCTIONNAME:","",PTV2))},error=function(cond){message("List value is empty")})
#             PTV3 <- tryCatch({paste("PRECURSORTYPE:",P1TV2)},error=function(cond){message("List value is empty")})
#             ################################################
# 	    print("enter my test ...2")
#             print(P1TV2)
# 	    ##########################################
#             print(as.character(InMEDA[["Adduct"]]))
	    ####################################################
      ######################################################      
        #####    if(!sjmisc::is_empty(P1TV2) || !sjmisc::is_empty(as.character(InMEDA[["Adduct"]]))){
	    ####################################################
      #####################################################        
            ##if(identical(P1TV2,as.character(InMEDA[["Adduct"]]))){
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
              ####out<-c(out,PTV3)
              out<-c(out,F1PT1)
              #############################
              FIN1<-InMEDA[["Ionization mode"]]
              F1IN1<-as.character(FIN1)
              F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
              out<-c(out,F2IN1)
              #############################################
	    
	      IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
	      #############################################
	      ##IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(InKeyVal)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
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
	      ##FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}),sep=" ")
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
	      FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:","",sep=""))
              out<-c(out,FINS2)
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
	 ##############################################################
        ##############################################################
            ######}
        #############################################################
        ##############################################################    
#             else{
# 		    ##############################
# 		    print("enter the part 2")
# 		    ##############################
# 		    ##print(as.character(InMEDA[["Adduct"]]))
# 		    ##print(INLL)
#                     ##print(FNA)
# 		    ##print(as.character(InMEDA[["Adduct"]]) == "[M]+")
# 		    #####################################
# 		    #####################################
# 		    if(as.character(InMEDA[["Adduct"]]) == "[M]+" || as.character(InMEDA[["Adduct"]]) == "[M]-"){
# 		    #########################################	    
# 		    ##if(as.character(InMEDA[["Adduct"]]) == "[M]+"){
# 		    ###############################
# 		    print("enter the part 2...if loop")
# 		    ##################################
# 		    FNA1<-which(stringi::stri_detect_fixed(FNA,"NAME:"))
# 		    FNA2<-InMEDA[["Name"]]
# 		    FNA3<-as.character(FNA2)
# 		    #######################
# 		    FNAM<-paste("NAME:",FNA3,sep=" ")
# 		    out<-c(out,FNAM)
# 		    #########################
# 		    FRA1<-which(stringi::stri_detect_fixed(FNA,"RETENTIONTIME:"))
# 		    F1RA1<-FNA[FRA1]
# 		    out<-c(out,F1RA1)
# 		    ################################
# 		    FMZ1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORMZ:"))
# 		    F1MZ1<-FNA[FMZ1]
# 		    out<-c(out,F1MZ1)
# 		    ################################
# 		    FPT1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORTYPE:"))
# 		    F1PT1<-FNA[FPT1]
# 		    out<-c(out,PTV3)
# 		    #################################
# 		    FIN1<-InMEDA[["Ionization mode"]]
# 		    F1IN1<-as.character(FIN1)
# 		    F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
# 		    out<-c(out,F2IN1)
# 		    ###################################################
# 		    IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
# 		    ###################################################
# 		    ##IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(InKeyVal)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
# 		    #################################################
# 		    if(!sjmisc::is_empty(IKCRV)){
# 			    ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
# 			    F1ONT<-paste("Ontology:",ONTV,sep=" ")
# 			    out<-c(out,F1ONT)
# 	            }else{
# 			    F1ONT<-paste("Ontology:","",sep=" ")
# 			    out<-c(out,F1ONT)
# 		    }
# 		    ################################
# 		    FINK<-paste("INCHIKEY:",tryCatch({IK$inchikey},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
# 		    out<-c(out,FINK)
# 		    FINCH<-paste("INCHI:",InchiV,sep=" ")
# 		    out<-c(out,FINCH)
# 		    FSIM<-paste("SMILES:",SM1,sep=" ")
# 		    out<-c(out,FSIM)
# 		    ##############################
# 		    FFOR<-FM$formula
# 		    FFOR1<-paste("FORMULA:",FFOR,sep=" ")
# 		    out<-c(out,FFOR1)
# 		    #############################
# 		    FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
# 		    FINS1<-FNA[FINS]
# 		    FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:","",sep=""))
# 		    out<-c(out,FINS2)
# 		    ############################
# 		    FAUT<-as.character(InMEDA[["Authors"]])
# 		    FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
# 		    out<-c(out,FAUT1)
# 		    ##########################
# 		    FLIC<-paste("LICENSE:","CC BY",sep=" ")
# 		    out<-c(out,FLIC)
# 		    ###########################
# 		    FCIE<-paste("COLLISIONENERGY:",as.character(InMEDA[["Collision energy"]]),sep=" ")
# 		    out<-c(out,FCIE)
# 		    ############################
# 		    FINST<-as.character(InMEDA[["INSTRUMENT_TYPE"]])
# 		    FINST1<-paste("INSTRUMENTTYPE:",FINST,sep=" ")
# 		    out<-c(out,FINST1)
# 		    ########################
# 		    FINS<-as.character(InMEDA[["INSTRUMENT"]])
# 		    FINS1<-paste("INSTRUMENT:",FINS,sep=" ")
# 		    out<-c(out,FINS1)
# 		    ########################
# 		    FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
# 		    out<-c(out,FCOM)
# 		    #########################
# 		    FNPA<-which(stringi::stri_detect_fixed(FNA,"Num Peaks:"))
# 		    F1NPA<-FNA[FNPA]
# 		    out<-c(out,F1NPA)
# 		    ###################
# 		    Fpea<-FNA[(FNPA+1):Find]
# 		    #########################
# 		    if(is.na(Fpea))
# 		     {
# 			Fpea1<-FNA[(FNPA+1)]
# 		     }else{
# 		        MV=AAMS1
# 		        tes1<-unlist(strsplit(Fpea, "\t"))
# 			tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
# 			tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
# 			tes4<-which(tes2 > (3+MV))
# 			if(length(tes4)>1)
# 			{
# 			    tes5<-tes2[-tes4]
# 		            tes6<-tes3[-tes4]
# 			    tes7<-paste(tes5,tes6,sep="\t")
# 			    out<-c(out,tes7)
# 
# 			}else{
# 			    out<-c(out,Fpea)
# 		        }
# 		      }
# 		    #############################################
# 	             ##}
#            ######################################################
# 	    }
# 	   ##########################################################
#     ############################################################          
# 	    } 
    ##########################################################
  #############################################################          
          } else{
            #print("entering the line 750")
            PASS<-RRV
          }
          ####### This is testing , if this works
          return(out)
          #################
        } ### this is mz closing brace
      } ## this is RT closing braces
    }else{
      ########################################
      ########################################	    
      print("enter the line 317......either PMZ or Adduct(AUIN) empty")
      print("entering the else part ---inchikey")
      ########################################
      NCIDV<-tryCatch({webchem::get_cid(InKeyVal, from = "inchikey")},error=function(cond){message("webchecm could not fetch the info")})
      NCIDV1<-tryCatch({as.numeric(NCIDV[[2]][1])},error=function(cond){message("webchecm could not fetch the info")})
      ########################################
      ##NCIDV1<-tryCatch({as.numeric(NCIDV$cid)},error=function(cond){message("webchecm could not fetch the info")})
      ##FMWFS<-tryCatch({webchem::cs_compinfo(NCIDV1,c("Formula","MolecularWeight","MonoisotopicMass"), verbose = TRUE)},error=function(cond){message("webchecm could not fetch the info")})
      ##FMWFS1<-tryCatch({FMWFS$monoisotopicMass},error=function(cond){message("CID is empty")})
      #######################################
      FMWFS<-tryCatch({webchem::pc_prop(NCIDV1, properties = c("MolecularFormula", "ExactMass","CanonicalSMILES","InChI","InChIKey"))},error=function(cond){message("webchem conversion from CID is not sucessful")})
      FMWFS1<-tryCatch({FMWFS$ExactMass},error=function(cond){message("CID is empty")})
      SM<-tryCatch({FMWFS$CanonicalSMILES},error=function(cond){message("CID is empty")})
      IN<-tryCatch({FMWFS$InChI},error=function(cond){message("CID is empty")})
      IK<-tryCatch({FMWFS$InChIKey},error=function(cond){message("CID is empty")})
      ######################################
      #print(FMWFS1)
      ######################################
      if(!sjmisc::is_empty(FMWFS1)){
      ##if(!sjmisc::is_empty(AUIN) & !sjmisc::is_empty(FMWFS1)){
	############################################################
        print("enter the line 325")
        ############################################################
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
	##ITmass<-match(Tmass,InPMZ)
        ###########################
        VRT<-as.numeric(as.character(InMEDA[["RT (min)"]]))
        VRTL<-VRT-0.20 
        VRTU<-VRT+0.20
	##################################
	print("enter temp test ---function 1 ...else part")
	##print(FMWFS1)
	##print(AAMS1)
	##print(MPPmL)
	##print(MPPmU)
	##print(VRT)
	##print(VRTL)
	##print(VRTU)
        #############################
        TRTL<-InRTL[InRTL >=VRTL & InRTL <= VRTU]
        ITRTL<-which(InRTL %in% TRTL)
	#############################
	##ITRTL<-match(TRTL,InRTL)
        #############################
	print("enter the line ...1196")
        if(length(ITRTL) >= 1){
          print("enter the line ...1197")
          if(length(ITmass) >= 1){
            print("enter the line ...1198")
            ###########################
            INLL<-intersect(ITmass,ITRTL)
            ##########################
            if(length(INLL) == 1){
	      ########################
	      print("enter the line ...1199")
              ######################
              F1FPL<-InMSPL[INLL]
              #####################
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
#               PT <- c("PRECURSORTYPE:", "ADDUCTIONNAME:")
#               NPT<-tryCatch({grep(paste(PT,collapse="|"), FNA, value=TRUE)},error=function(cond){message("List value is empty")})
#               NPT1<-tryCatch({match(NPT,FNA)},error=function(cond){message("List value is empty")})
#               PTV <- tryCatch({stringr::str_remove(FNA[NPT1],c("PRECURSORTYPE:","ADDUCTIONNAME:"))},error=function(cond){message("List value is empty")})
#               PTV1<- tryCatch({PTV[1]},error=function(cond){message("List value is empty")})
#               PTV2<-tryCatch({stringr::str_trim(PTV1)},error=function(cond){message("List value is empty")})
# 	      P1TV2<-tryCatch({stringr::str_trim(gsub("ADDUCTIONNAME:","",PTV2))},error=function(cond){message("List value is empty")})
#               PTV3<-tryCatch({paste("PRECURSORTYPE:",P1TV2)},error=function(cond){message("List value is empty")})
#               #########################################################
# 	      print("enter my test...3")
#               print(P1TV2)
#               print(as.character(InMEDA[["Adduct"]]))
	      #############################################
        #############################################      
            ###  if(!sjmisc::is_empty(P1TV2) || !sjmisc::is_empty(as.character(InMEDA[["Adduct"]]))){
	      ###############################################
        ###############################################      
              ##if(identical(P1TV2,as.character(InMEDA[["Adduct"]]))){
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
                ##out<-c(out,PTV3)
                out<-c(out,F1PT1)
                ###########################
                FIN1<-InMEDA[["Ionization mode"]]
                F1IN1<-as.character(FIN1)
                F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                out<-c(out,F2IN1)
                ##################################################
		IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
		##################################################
		##IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(InKeyVal)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
		########################################
		if(!sjmisc::is_empty(IKCRV)){
			ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
			F1ONT<-paste("Ontology:",ONTV,sep=" ")
			out<-c(out,F1ONT)
		}else{
			F1ONT<-paste("Ontology:","",sep=" ")
			out<-c(out,F1ONT)
		}
		#################################################
                #################################################
                FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                out<-c(out,FINK)
                ##FINCH<-paste("INCHI:",InchiV,sep=" ")
		FINCH<-paste("INCHI:",IN,sep=" ")
                out<-c(out,FINCH)
                FSIM<-paste("SMILES:",SM,sep=" ")
                out<-c(out,FSIM)
                #######################
                FFOR<-FM$formula
                FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                out<-c(out,FFOR1)
                #######################
                FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                FINS1<-FNA[FINS]
		FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:","",sep=""))
                out<-c(out,FINS2)
                #######################
                FAUT<-as.character(InMEDA[["Authors"]])
                FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
                out<-c(out,FAUT1)
                ##########################
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
            #############################################################################    
            #############################################################################
              ###}
            #############################################################################
            ##############################################################################  
              # else{
		      #######################################
# 		      print("enter the part 3")
# 		      #########################################
# 		      ##print(as.character(InMEDA[["Adduct"]]))
# 		      ##print(INLL)
#                       ##print(FNA)
# 		      ##print(as.character(InMEDA[["Adduct"]]) == "[M]+")
# 		      ######################################
# 		      ##if(as.character(InMEDA[["Adduct"]]) == "[M]+"){
# 		      #########################################
# 		      #########################################
# 		      if(as.character(InMEDA[["Adduct"]]) == "[M]+" || as.character(InMEDA[["Adduct"]]) == "[M]-"){
# 		      #######################################
# 		      print("enter the part 3..if loop")
# 		      #######################################
# 		      FNA1<-which(stringi::stri_detect_fixed(FNA,"NAME:"))
# 		      FNA2<-InMEDA[["Name"]]
# 		      FNA3<-as.character(FNA2)
# 		      #######################
# 		      FNAM<-paste("NAME:",FNA3,sep=" ")
# 		      out<-c(out,FNAM)
# 		      #######################
# 		      FRA1<-which(stringi::stri_detect_fixed(FNA,"RETENTIONTIME:"))
# 		      F1RA1<-FNA[FRA1]
# 		      out<-c(out,F1RA1)
# 		      ################################
# 		      FMZ1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORMZ:"))
# 		      F1MZ1<-FNA[FMZ1]
# 		      out<-c(out,F1MZ1)
# 		      ################################
# 		      FPT1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORTYPE:"))
# 		      F1PT1<-FNA[FPT1]
# 		      out<-c(out,PTV3)
# 		      #################################
# 		      FIN1<-InMEDA[["Ionization mode"]]
# 		      F1IN1<-as.character(FIN1)
# 		      F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
# 		      out<-c(out,F2IN1)
# 		      #################################################
# 		      IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
# 		      #################################################
# 		      ##IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(InKeyVal)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
# 		      #################################################
# 		      if(!sjmisc::is_empty(IKCRV)){
# 		        ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
# 		        F1ONT<-paste("Ontology:",ONTV,sep=" ")
# 		        out<-c(out,F1ONT)
# 		      }else{
# 		        F1ONT<-paste("Ontology:","",sep=" ")
# 		        out<-c(out,F1ONT)
# 		      }
# 		      
# 		      #############################################
# 		      ############################################
# 		      ##print("enter the line ...1878")
# 		      ############################################
# 		      FINK<-paste("INCHIKEY:",tryCatch({IK$inchikey},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
# 		      ##FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
# 		      out<-c(out,FINK)
# 		      FINCH<-paste("INCHI:",InchiV,sep=" ")
# 		      out<-c(out,FINCH)
# 		      FSIM<-paste("SMILES:",SM1,sep=" ")
# 		      out<-c(out,FSIM)
# 		      ##############################
# 		      FFOR<-FM$formula
# 		      FFOR1<-paste("FORMULA:",FFOR,sep=" ")
# 		      out<-c(out,FFOR1)
# 		      #############################
# 		      FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
# 		      FINS1<-FNA[FINS]
# 		      FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:","",sep=""))
# 		      out<-c(out,FINS2)
# 		      ############################
# 		      FAUT<-as.character(InMEDA[["Authors"]])
# 		      FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
# 		      out<-c(out,FAUT1)
# 		      ##########################
# 		      ##FLIC<-paste("LICENSE:",sep=" ")
# 		      FLIC<-paste("LICENSE:","CC BY",sep=" ")
# 		      out<-c(out,FLIC)
# 		      ###########################
# 		      FCIE<-paste("COLLISIONENERGY:",as.character(InMEDA[["Collision energy"]]),sep=" ")
# 		      out<-c(out,FCIE)
# 		      #########################
# 		      FINST<-as.character(InMEDA[["INSTRUMENT_TYPE"]])
# 		      FINST1<-paste("INSTRUMENTTYPE:",FINST,sep=" ")
# 		      out<-c(out,FINST1)
# 		      ########################
# 		      FINS<-as.character(InMEDA[["INSTRUMENT"]])
# 		      FINS1<-paste("INSTRUMENT:",FINS,sep=" ")
# 		      out<-c(out,FINS1)
# 		      ####################
# 		      FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
# 		      out<-c(out,FCOM)
# 		      ##################
# 		      FNPA<-which(stringi::stri_detect_fixed(FNA,"Num Peaks:"))
# 		      F1NPA<-FNA[FNPA]
# 		      out<-c(out,F1NPA)
# 		      ###################
# 		      Fpea<-FNA[(FNPA+1):Find]
# 		      #########################
# 		      if(is.na(Fpea))
# 		      {
# 		        Fpea1<-FNA[(FNPA+1)]
# 		        
# 		        
# 		      }else{
# 		        
# 		        MV=AAMS1
# 		        tes1<-unlist(strsplit(Fpea, "\t"))
# 		        tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
# 		        tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
# 		        tes4<-which(tes2 > (3+MV))
# 		        if(length(tes4)>1)
# 		        {
# 		          tes5<-tes2[-tes4]
# 		          tes6<-tes3[-tes4]
# 		          tes7<-paste(tes5,tes6,sep="\t")
# 		          out<-c(out,tes7)
# 		        }else{
# 		          out<-c(out,Fpea)
# 		          
# 		        }
# 		      }
# 	    ###########################################
# 		      }
# 	    ###########################################
#     ###########################################
#               }
    #############################################          
    ############################################################################################  
            } else if(length(INLL) > 1){
              ###################################		    
              print("entering the 1569")
	      ###################################
              MONMS=InMSPL[INLL]
              TRA<-unname(rapply(MONMS, function(x) grep("RETENTIONTIME:",x, value=TRUE)))
              TRA1<-as.numeric(stringr::str_trim(stringr::str_replace(TRA, "RETENTIONTIME:", "")))
              TRA2<-abs(VRT-TRA1)
              TRA3<-which.min(TRA2)
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
#               PT <- c("PRECURSORTYPE:", "ADDUCTIONNAME:")
#               NPT<-tryCatch({grep(paste(PT,collapse="|"), FNA, value=TRUE)},error=function(cond){message("List value is empty")})
#               NPT1<-tryCatch({match(NPT,FNA)},error=function(cond){message("List value is empty")})
#               PTV <-tryCatch({stringr::str_remove(FNA[NPT1],c("PRECURSORTYPE:","ADDUCTIONNAME:"))},error=function(cond){message("List value is empty")})
#               PTV1<-tryCatch({PTV[1]},error=function(cond){message("List value is empty")})
#               PTV2<-tryCatch({stringr::str_trim(PTV1)},error=function(cond){message("List value is empty")})
# 	      P1TV2<-tryCatch({stringr::str_trim(gsub("ADDUCTIONNAME:","",PTV2))},error=function(cond){message("List value is empty")})
#               PTV3<-tryCatch({paste("PRECURSORTYPE:",P1TV2)},error=function(cond){message("List value is empty")})
# 	      ##PTV3<-tryCatch({paste("PRECURSORTYPE:",PTV2)},error=function(cond){message("List value is empty")})
# 	      #########################################
# 	      print("enter my test...4")
#               print(P1TV2)
#               print(as.character(InMEDA[["Adduct"]]))
              ######################## adding this new
              ###########################################
            ###  if(!sjmisc::is_empty(P1TV2) || !sjmisc::is_empty(as.character(InMEDA[["Adduct"]]))){
	      ##############################################
                ####################################
              ##if(identical(P1TV2,as.character(InMEDA[["Adduct"]]))){
                ##if(!sjmisc::is_empty(AUIN) & !sjmisc::is_empty(FMWFS1)){
                ###########################################
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
                ##out<-c(out,PTV3)
                out<-c(out,F1PT1)
                #############################
                FIN1<-InMEDA[["Ionization mode"]]
                F1IN1<-as.character(FIN1)
                F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                out<-c(out,F2IN1)
                ###################################
		IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
		###################################
		##IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(InKeyVal)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
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
                FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}),sep=" ")
                out<-c(out,FINK)
		FINCH<-paste("INCHI:",IN,sep=" ")
                ##FINCH<-paste("INCHI:",InchiV,sep=" ")
                out<-c(out,FINCH)
                FSIM<-paste("SMILES:",SM,sep=" ")
                out<-c(out,FSIM)
                ##############################
                FFOR<-tryCatch({FM$formula},error=function(cond){message("Formula is empty")})
                FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                out<-c(out,FFOR1)
                ###############################
                FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                FINS1<-FNA[FINS]
		FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:","",sep=""))
                out<-c(out,FINS2)
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
         ##########################################
            #####  }
          ###########################################
          ############################################    
#               else{
#           ########################################
#           print("enter the part 4")
# 	  ##print(as.character(InMEDA[["Adduct"]]))
# 	  ##print(INLL)
#           ##print(FNA)
# 	  ##print(as.character(InMEDA[["Adduct"]]) == "[M]+")
#           ########################
# 	  ##if(as.character(InMEDA[["Adduct"]]) == "[M]+"){
# 	  ####################################
# 	  ####################################
# 	  if(as.character(InMEDA[["Adduct"]]) == "[M]+" || as.character(InMEDA[["Adduct"]]) == "[M]-"){
#           ########################
# 	  print("enter the part 4..if loop")
#           ##########################
#           FNA1<-which(stringi::stri_detect_fixed(FNA,"NAME:"))
#           FNA2<-InMEDA[["Name"]]
#           FNA3<-as.character(FNA2)
#           #######################
#           FNAM<-paste("NAME:",FNA3,sep=" ")
#           out<-c(out,FNAM)
#           ########################
#           FRA1<-which(stringi::stri_detect_fixed(FNA,"RETENTIONTIME:"))
#           F1RA1<-FNA[FRA1]
#           out<-c(out,F1RA1)
#           ################################
#           FMZ1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORMZ:"))
#           F1MZ1<-FNA[FMZ1]
#           out<-c(out,F1MZ1)
#           ################################
#           FPT1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORTYPE:"))
#           F1PT1<-FNA[FPT1]
#           out<-c(out,PTV3)
#           #out<-c(out,F1PT1)
#           #################################
#           FIN1<-InMEDA[["Ionization mode"]]
#           F1IN1<-as.character(FIN1)
#           F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
#           out<-c(out,F2IN1)
#           ##################################################
#           
#           IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
#           ###################################################
#           ##IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(InKeyVal)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
#           #################################################
#           if(!sjmisc::is_empty(IKCRV)){
#             ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
#             F1ONT<-paste("Ontology:",ONTV,sep=" ")
#             out<-c(out,F1ONT)
#           }else{
#             F1ONT<-paste("Ontology:","",sep=" ")
#             out<-c(out,F1ONT)
#           }
#           
#           #############################################
#           ############################################
#           ##print("enter the line ...1878")
#           ############################################
#           FINK<-paste("INCHIKEY:",tryCatch({IK$inchikey},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
#           ##FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
#           out<-c(out,FINK)
#           FINCH<-paste("INCHI:",InchiV,sep=" ")
#           out<-c(out,FINCH)
#           FSIM<-paste("SMILES:",SM1,sep=" ")
#           out<-c(out,FSIM)
#           ##############################
#           FFOR<-FM$formula
#           FFOR1<-paste("FORMULA:",FFOR,sep=" ")
#           out<-c(out,FFOR1)
#           #############################
#           FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
#           FINS1<-FNA[FINS]
# 	  FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:","",sep=""))
#           out<-c(out,FINS2)
#           ############################
#           FAUT<-as.character(InMEDA[["Authors"]])
#           FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
#           out<-c(out,FAUT1)
#           ##########################
#           ##FLIC<-paste("LICENSE:",sep=" ")
#           FLIC<-paste("LICENSE:","CC BY",sep=" ")
#           out<-c(out,FLIC)
#           ###########################
#           FCIE<-paste("COLLISIONENERGY:",as.character(InMEDA[["Collision energy"]]),sep=" ")
#           out<-c(out,FCIE)
#           #########################
#           FINST<-as.character(InMEDA[["INSTRUMENT_TYPE"]])
#           FINST1<-paste("INSTRUMENTTYPE:",FINST,sep=" ")
#           out<-c(out,FINST1)
#           ########################
#           FINS<-as.character(InMEDA[["INSTRUMENT"]])
#           FINS1<-paste("INSTRUMENT:",FINS,sep=" ")
#           out<-c(out,FINS1)
#           ####################
#           FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
#           out<-c(out,FCOM)
#           ##################
#           FNPA<-which(stringi::stri_detect_fixed(FNA,"Num Peaks:"))
#           F1NPA<-FNA[FNPA]
#           out<-c(out,F1NPA)
#           ###################
#           Fpea<-FNA[(FNPA+1):Find]
#           #########################
#           if(is.na(Fpea))
#           {
#             Fpea1<-FNA[(FNPA+1)]
#             
#             
#           }else{
#             
#             MV=AAMS1
#             tes1<-unlist(strsplit(Fpea, "\t"))
#             tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
#             tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
#             tes4<-which(tes2 > (3+MV))
#             if(length(tes4)>1)
#             {
#               tes5<-tes2[-tes4]
#               tes6<-tes3[-tes4]
#               tes7<-paste(tes5,tes6,sep="\t")
#               out<-c(out,tes7)
#             }else{
#               out<-c(out,Fpea)
#               
#             }
#           }
#         ############################
# 	  }             
# 	#################################
#   #################################
# 	      }
  #### need to add } to #####################  
################################  
            } else{
              #print("entering the line 750")
              PASS<-RRV
            }
        ####### This is testing , if this works
            return(out)
        #################
          } ### this is mz closing brace
        } ## this is RT closing braces
      #############################################
      #############################################  
         }
      ##########################################
      ##########################################
      else{
        ### THis is the new function I am adding from REST API conversion ..Inchikey to exact value 
        #######################################
        #######################################
        FMWFS1<-ifelse(!sjmisc::is_empty(tryCatch({ConvINKtoOID(InKeyVal)[4]},error=function(cond){message("Inchikey conversion is empty")})),tryCatch({ConvINKtoOID(InKeyVal)[4]},error=function(cond){message("iNCHIKEY conversion is empty")}),"NA")
        ##if(!sjmisc::is_empty(AUIN) & !sjmisc::is_empty(FMWFS1)){
        ##########################################################
        ###########################################################  
          if(!sjmisc::is_empty(FMWFS1)){
          ############################################################
          print("enter the line 325")
          ############################################################
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
          ##ITmass<-match(Tmass,InPMZ)
          ###########################
          VRT<-as.numeric(as.character(InMEDA[["RT (min)"]]))
          VRTL<-VRT-0.20 
          VRTU<-VRT+0.20
          ##################################
          print("enter temp test ---function 1 ...else part")
          ##print(FMWFS1)
          ##print(AAMS1)
          ##print(MPPmL)
          ##print(MPPmU)
          ##print(VRT)
          ##print(VRTL)
          ##print(VRTU)
          #############################
          TRTL<-InRTL[InRTL >=VRTL & InRTL <= VRTU]
          ITRTL<-which(InRTL %in% TRTL)
          #############################
          ##ITRTL<-match(TRTL,InRTL)
          #############################
          print("enter the line ...1196")
          if(length(ITRTL) >= 1){
            print("enter the line ...1197")
            if(length(ITmass) >= 1){
              print("enter the line ...1198")
              ###########################
              INLL<-intersect(ITmass,ITRTL)
              ##########################
              if(length(INLL) == 1){
                ########################
                print("enter the line ...1199")
                ######################
                F1FPL<-InMSPL[INLL]
                #####################
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
                # PT <- c("PRECURSORTYPE:", "ADDUCTIONNAME:")
                # NPT<-tryCatch({grep(paste(PT,collapse="|"), FNA, value=TRUE)},error=function(cond){message("List value is empty")})
                # NPT1<-tryCatch({match(NPT,FNA)},error=function(cond){message("List value is empty")})
                # PTV <- tryCatch({stringr::str_remove(FNA[NPT1],c("PRECURSORTYPE:","ADDUCTIONNAME:"))},error=function(cond){message("List value is empty")})
                # PTV1<- tryCatch({PTV[1]},error=function(cond){message("List value is empty")})
                # PTV2<-tryCatch({stringr::str_trim(PTV1)},error=function(cond){message("List value is empty")})
                # P1TV2<-tryCatch({stringr::str_trim(gsub("ADDUCTIONNAME:","",PTV2))},error=function(cond){message("List value is empty")})
                # PTV3<-tryCatch({paste("PRECURSORTYPE:",P1TV2)},error=function(cond){message("List value is empty")})
                # #########################################################
                # print("enter my test...3")
                # print(P1TV2)
                # print(as.character(InMEDA[["Adduct"]]))
                #################################################
                ##if(!sjmisc::is_empty(P1TV2) || !sjmisc::is_empty(as.character(InMEDA[["Adduct"]]))){
              ###############################################
              ################################################  
                  ##if(identical(P1TV2,as.character(InMEDA[["Adduct"]]))){
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
                  ##out<-c(out,PTV3)
                  out<-c(out,F1PT1)
                  ###########################
                  FIN1<-InMEDA[["Ionization mode"]]
                  F1IN1<-as.character(FIN1)
                  F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                  out<-c(out,F2IN1)
                  ##################################################
                  IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
                  ##################################################
                  ##IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(InKeyVal)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
                  ########################################
                  if(!sjmisc::is_empty(IKCRV)){
                    ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
                    F1ONT<-paste("Ontology:",ONTV,sep=" ")
                    out<-c(out,F1ONT)
                  }else{
                    F1ONT<-paste("Ontology:","",sep=" ")
                    out<-c(out,F1ONT)
                  }
                  #################################################
                  #################################################
                  FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                  out<-c(out,FINK)
                  ##FINCH<-paste("INCHI:",InchiV,sep=" ")
                  FINCH<-paste("INCHI:",IN,sep=" ")
                  out<-c(out,FINCH)
                  FSIM<-paste("SMILES:",SM,sep=" ")
                  out<-c(out,FSIM)
                  #######################
                  FFOR<-FM$formula
                  FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                  out<-c(out,FFOR1)
                  #######################
                  FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                  FINS1<-FNA[FINS]
                  FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:","",sep=""))
                  out<-c(out,FINS2)
                  #######################
                  FAUT<-as.character(InMEDA[["Authors"]])
                  FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
                  out<-c(out,FAUT1)
                  ##########################
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
            #############################################################################    
            #############################################################################
               ### }
              #########################################################################
              #########################################################################  
            #     else{
            #       #######################################
            #       print("enter the part 3")
            #       #########################################
            #       ##print(as.character(InMEDA[["Adduct"]]))
            #       ##print(INLL)
            #       ##print(FNA)
            #       ##print(as.character(InMEDA[["Adduct"]]) == "[M]+")
            #       ######################################
            #       ##if(as.character(InMEDA[["Adduct"]]) == "[M]+"){
            #       #########################################
            #       #########################################
            #       if(as.character(InMEDA[["Adduct"]]) == "[M]+" || as.character(InMEDA[["Adduct"]]) == "[M]-"){
            #         #######################################
            #         print("enter the part 3..if loop")
            #         #######################################
            #         FNA1<-which(stringi::stri_detect_fixed(FNA,"NAME:"))
            #         FNA2<-InMEDA[["Name"]]
            #         FNA3<-as.character(FNA2)
            #         #######################
            #         FNAM<-paste("NAME:",FNA3,sep=" ")
            #         out<-c(out,FNAM)
            #         #######################
            #         FRA1<-which(stringi::stri_detect_fixed(FNA,"RETENTIONTIME:"))
            #         F1RA1<-FNA[FRA1]
            #         out<-c(out,F1RA1)
            #         ################################
            #         FMZ1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORMZ:"))
            #         F1MZ1<-FNA[FMZ1]
            #         out<-c(out,F1MZ1)
            #         ################################
            #         FPT1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORTYPE:"))
            #         F1PT1<-FNA[FPT1]
            #         out<-c(out,PTV3)
            #         #################################
            #         FIN1<-InMEDA[["Ionization mode"]]
            #         F1IN1<-as.character(FIN1)
            #         F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
            #         out<-c(out,F2IN1)
            #         #################################################
            #         IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
            #         #################################################
            #         ##IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(InKeyVal)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
            #         #################################################
            #         if(!sjmisc::is_empty(IKCRV)){
            #           ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
            #           F1ONT<-paste("Ontology:",ONTV,sep=" ")
            #           out<-c(out,F1ONT)
            #         }else{
            #           F1ONT<-paste("Ontology:","",sep=" ")
            #           out<-c(out,F1ONT)
            #         }
            #         
            #         #############################################
            #         ############################################
            #         ##print("enter the line ...1878")
            #         ############################################
            #         FINK<-paste("INCHIKEY:",tryCatch({IK$inchikey},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
            #         ##FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
            #         out<-c(out,FINK)
            #         FINCH<-paste("INCHI:",InchiV,sep=" ")
            #         out<-c(out,FINCH)
            #         FSIM<-paste("SMILES:",SM1,sep=" ")
            #         out<-c(out,FSIM)
            #         ##############################
            #         FFOR<-FM$formula
            #         FFOR1<-paste("FORMULA:",FFOR,sep=" ")
            #         out<-c(out,FFOR1)
            #         #############################
            #         FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
            #         FINS1<-FNA[FINS]
            #         FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:","",sep=""))
            #         out<-c(out,FINS2)
            #         ############################
            #         FAUT<-as.character(InMEDA[["Authors"]])
            #         FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
            #         out<-c(out,FAUT1)
            #         ##########################
            #         ##FLIC<-paste("LICENSE:",sep=" ")
            #         FLIC<-paste("LICENSE:","CC BY",sep=" ")
            #         out<-c(out,FLIC)
            #         ###########################
            #         FCIE<-paste("COLLISIONENERGY:",as.character(InMEDA[["Collision energy"]]),sep=" ")
            #         out<-c(out,FCIE)
            #         #########################
            #         FINST<-as.character(InMEDA[["INSTRUMENT_TYPE"]])
            #         FINST1<-paste("INSTRUMENTTYPE:",FINST,sep=" ")
            #         out<-c(out,FINST1)
            #         ########################
            #         FINS<-as.character(InMEDA[["INSTRUMENT"]])
            #         FINS1<-paste("INSTRUMENT:",FINS,sep=" ")
            #         out<-c(out,FINS1)
            #         ####################
            #         FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
            #         out<-c(out,FCOM)
            #         ##################
            #         FNPA<-which(stringi::stri_detect_fixed(FNA,"Num Peaks:"))
            #         F1NPA<-FNA[FNPA]
            #         out<-c(out,F1NPA)
            #         ###################
            #         Fpea<-FNA[(FNPA+1):Find]
            #         #########################
            #         if(is.na(Fpea))
            #         {
            #           Fpea1<-FNA[(FNPA+1)]
            #           
            #           
            #         }else{
            #           
            #           MV=AAMS1
            #           tes1<-unlist(strsplit(Fpea, "\t"))
            #           tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
            #           tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
            #           tes4<-which(tes2 > (3+MV))
            #           if(length(tes4)>1)
            #           {
            #             tes5<-tes2[-tes4]
            #             tes6<-tes3[-tes4]
            #             tes7<-paste(tes5,tes6,sep="\t")
            #             out<-c(out,tes7)
            #           }else{
            #             out<-c(out,Fpea)
            #             
            #           }
            #         }
            #         ###########################################
            #       }
            # ###########################################
            # ###########################################
            #     }
            ####################################################      
            ############################################################################################  
              } else if(length(INLL) > 1){
                ###################################		    
                print("entering the 1569")
                ###################################
                MONMS=InMSPL[INLL]
                TRA<-unname(rapply(MONMS, function(x) grep("RETENTIONTIME:",x, value=TRUE)))
                TRA1<-as.numeric(stringr::str_trim(stringr::str_replace(TRA, "RETENTIONTIME:", "")))
                TRA2<-abs(VRT-TRA1)
                TRA3<-which.min(TRA2)
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
                # PT <- c("PRECURSORTYPE:", "ADDUCTIONNAME:")
                # NPT<-tryCatch({grep(paste(PT,collapse="|"), FNA, value=TRUE)},error=function(cond){message("List value is empty")})
                # NPT1<-tryCatch({match(NPT,FNA)},error=function(cond){message("List value is empty")})
                # PTV <-tryCatch({stringr::str_remove(FNA[NPT1],c("PRECURSORTYPE:","ADDUCTIONNAME:"))},error=function(cond){message("List value is empty")})
                # PTV1<-tryCatch({PTV[1]},error=function(cond){message("List value is empty")})
                # PTV2<-tryCatch({stringr::str_trim(PTV1)},error=function(cond){message("List value is empty")})
                # P1TV2<-tryCatch({stringr::str_trim(gsub("ADDUCTIONNAME:","",PTV2))},error=function(cond){message("List value is empty")})
                # PTV3<-tryCatch({paste("PRECURSORTYPE:",P1TV2)},error=function(cond){message("List value is empty")})
                # ##PTV3<-tryCatch({paste("PRECURSORTYPE:",PTV2)},error=function(cond){message("List value is empty")})
                # ####################################
                # print("enter my test...4")
                # print(P1TV2)
                # print(as.character(InMEDA[["Adduct"]]))
                ######################## adding this new
                ################################################
              ###  if(!sjmisc::is_empty(P1TV2) || !sjmisc::is_empty(as.character(InMEDA[["Adduct"]]))){
                  ##############################################
                  ###############################################
                  ##if(identical(P1TV2,as.character(InMEDA[["Adduct"]]))){
                  ##if(!sjmisc::is_empty(AUIN) & !sjmisc::is_empty(FMWFS1)){
                  ###########################################
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
                  ###out<-c(out,PTV3)
                  out<-c(out,F1PT1)
                  #############################
                  FIN1<-InMEDA[["Ionization mode"]]
                  F1IN1<-as.character(FIN1)
                  F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                  out<-c(out,F2IN1)
                  ###################################
                  IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
                  ###################################
                  ##IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(InKeyVal)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
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
                  FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}),sep=" ")
                  out<-c(out,FINK)
                  FINCH<-paste("INCHI:",IN,sep=" ")
                  ##FINCH<-paste("INCHI:",InchiV,sep=" ")
                  out<-c(out,FINCH)
                  FSIM<-paste("SMILES:",SM,sep=" ")
                  out<-c(out,FSIM)
                  ##############################
                  FFOR<-tryCatch({FM$formula},error=function(cond){message("Formula is empty")})
                  FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                  out<-c(out,FFOR1)
                  ###############################
                  FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                  FINS1<-FNA[FINS]
                  FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:","",sep=""))
                  out<-c(out,FINS2)
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
            ##################################################      
            ##########################################
               ### }
                
      #           else{
      #     ########################################
      #     ###############################################        
      #             print("enter the part 4")
      #             ##print(as.character(InMEDA[["Adduct"]]))
      #             ##print(INLL)
      #             ##print(FNA)
      #             ##print(as.character(InMEDA[["Adduct"]]) == "[M]+")
      #             ########################
      #             ##if(as.character(InMEDA[["Adduct"]]) == "[M]+"){
      #             ####################################
      #             ####################################
      #             if(as.character(InMEDA[["Adduct"]]) == "[M]+" || as.character(InMEDA[["Adduct"]]) == "[M]-"){
      #               ########################
      #               print("enter the part 4..if loop")
      #               ##########################
      #               FNA1<-which(stringi::stri_detect_fixed(FNA,"NAME:"))
      #               FNA2<-InMEDA[["Name"]]
      #               FNA3<-as.character(FNA2)
      #               #######################
      #               FNAM<-paste("NAME:",FNA3,sep=" ")
      #               out<-c(out,FNAM)
      #               ########################
      #               FRA1<-which(stringi::stri_detect_fixed(FNA,"RETENTIONTIME:"))
      #               F1RA1<-FNA[FRA1]
      #               out<-c(out,F1RA1)
      #               ################################
      #               FMZ1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORMZ:"))
      #               F1MZ1<-FNA[FMZ1]
      #               out<-c(out,F1MZ1)
      #               ################################
      #               FPT1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORTYPE:"))
      #               F1PT1<-FNA[FPT1]
      #               out<-c(out,PTV3)
      #               #out<-c(out,F1PT1)
      #               #################################
      #               FIN1<-InMEDA[["Ionization mode"]]
      #               F1IN1<-as.character(FIN1)
      #               F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
      #               out<-c(out,F2IN1)
      #               ##################################################
      #               
      #               IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
      #               ###################################################
      #               ##IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(InKeyVal)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
      #               #################################################
      #               if(!sjmisc::is_empty(IKCRV)){
      #                 ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
      #                 F1ONT<-paste("Ontology:",ONTV,sep=" ")
      #                 out<-c(out,F1ONT)
      #               }else{
      #                 F1ONT<-paste("Ontology:","",sep=" ")
      #                 out<-c(out,F1ONT)
      #               }
      #               
      #               #############################################
      #               ############################################
      #               ##print("enter the line ...1878")
      #               ############################################
      #               FINK<-paste("INCHIKEY:",tryCatch({IK$inchikey},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
      #               ##FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
      #               out<-c(out,FINK)
      #               FINCH<-paste("INCHI:",InchiV,sep=" ")
      #               out<-c(out,FINCH)
      #               FSIM<-paste("SMILES:",SM1,sep=" ")
      #               out<-c(out,FSIM)
      #               ##############################
      #               FFOR<-FM$formula
      #               FFOR1<-paste("FORMULA:",FFOR,sep=" ")
      #               out<-c(out,FFOR1)
      #               #############################
      #               FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
      #               FINS1<-FNA[FINS]
      #               FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:","",sep=""))
      #               out<-c(out,FINS2)
      #               ############################
      #               FAUT<-as.character(InMEDA[["Authors"]])
      #               FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
      #               out<-c(out,FAUT1)
      #               ##########################
      #               ##FLIC<-paste("LICENSE:",sep=" ")
      #               FLIC<-paste("LICENSE:","CC BY",sep=" ")
      #               out<-c(out,FLIC)
      #               ###########################
      #               FCIE<-paste("COLLISIONENERGY:",as.character(InMEDA[["Collision energy"]]),sep=" ")
      #               out<-c(out,FCIE)
      #               #########################
      #               FINST<-as.character(InMEDA[["INSTRUMENT_TYPE"]])
      #               FINST1<-paste("INSTRUMENTTYPE:",FINST,sep=" ")
      #               out<-c(out,FINST1)
      #               ########################
      #               FINS<-as.character(InMEDA[["INSTRUMENT"]])
      #               FINS1<-paste("INSTRUMENT:",FINS,sep=" ")
      #               out<-c(out,FINS1)
      #               ####################
      #               FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
      #               out<-c(out,FCOM)
      #               ##################
      #               FNPA<-which(stringi::stri_detect_fixed(FNA,"Num Peaks:"))
      #               F1NPA<-FNA[FNPA]
      #               out<-c(out,F1NPA)
      #               ###################
      #               Fpea<-FNA[(FNPA+1):Find]
      #               #########################
      #               if(is.na(Fpea))
      #               {
      #                 Fpea1<-FNA[(FNPA+1)]
      #                 
      #                 
      #               }else{
      #                 
      #                 MV=AAMS1
      #                 tes1<-unlist(strsplit(Fpea, "\t"))
      #                 tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
      #                 tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
      #                 tes4<-which(tes2 > (3+MV))
      #                 if(length(tes4)>1)
      #                 {
      #                   tes5<-tes2[-tes4]
      #                   tes6<-tes3[-tes4]
      #                   tes7<-paste(tes5,tes6,sep="\t")
      #                   out<-c(out,tes7)
      #                 }else{
      #                   out<-c(out,Fpea)
      #                   
      #                 }
      #               }
      #               ############################
      #             }             
      # ##################################
      #   ###################################
      #           }
          #### need to add } to   
        ################################  
              } else{
                #print("entering the line 750")
                PASS<-RRV
              }
              ####### This is testing , if this works
              return(out)
              #################
            } ### this is mz closing brace
          } ## this is RT closing braces
    #######################################
    ########################################      
        }else{
          PASS1 <-RRV
        } 
  ########################################
  ########################################
      }
  ####################################
  ##########################################    
      ##If I have to add further function 
      ######################################
    }
    
  } 
}
###################################################################################################################################################

##print("enter the area for NFfilter")
##############################################################################################################################################
##############################################################################################################################################
NFFilter<-function(InMEDA,InAdVA,InMSPL,InPMZ,InRTL)
{
  ###########################
  out<-c()
  ############################
  print("entering the line...189")
  ###########################
  IV<-as.character(InMEDA[["InChI"]])
  EM<-as.character(InMEDA[["Exact mass"]])
  #############################
  EM1<-as.numeric(EM)
  #############################
  NEM1<-tryCatch({NFFilter1(InMEDA,InAdVA,InMSPL,InPMZ,InRTL)},error=function(cond){message("retention time is empty")})
  #############################
  print("enter my test area")
  print(NEM1)
  ###############################
  if(!sjmisc::is_empty(NEM1)){
    #############################
    print("enter the line ..196")
    ####################################
    ####################################
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
    #########################
    ##ITmass<-match(Tmass,InPMZ)
    #####################
    #####################
    VRT<-as.numeric(as.character(InMEDA[["RT (min)"]]))
    VRTL<-VRT-0.20
    VRTU<-VRT+0.20
    ######################
    TRTL<-InRTL[InRTL >=VRTL & InRTL <= VRTU]
    ITRTL<-which(InRTL %in% TRTL)
    ##ITRTL<-match(TRTL,InRTL)
    #########################
    ##print("entering ...second function")
    ##print(NEM1)
    ##print(ADV4)
    ##print(MPPmL)
    ##print(MPPmU)
    ##print(Tmass)
    ##print(ITmass)
    ##print(VRT)
    ##print(VRTL)
    ##print(VRTU)
    ##print(TRTL)
    ##print(ITRTL)
    ######################################
    INLL<-intersect(ITmass,ITRTL)
    ########################################
    print("entering the line 218")
    if(length(ITRTL) >= 1){
      print("entering the line 219")
      if(length(ITmass) >= 1){
        print("entering the line 220")
        if(length(INLL) == 1){
          ##################################
          print("entering the line 221")
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
	  #######################################################
#           PT <- c("PRECURSORTYPE:", "ADDUCTIONNAME:")
#           NPT<-tryCatch({grep(paste(PT,collapse="|"), FNA, value=TRUE)},error=function(cond){message("List value is empty")})
#           NPT1<-tryCatch({match(NPT,FNA)},error=function(cond){message("List value is empty")})
#           PTV <-tryCatch({stringr::str_remove(FNA[NPT1],c("PRECURSORTYPE:","ADDUCTIONNAME:"))},error=function(cond){message("List value is empty")})
#           PTV1<-tryCatch({PTV[1]},error=function(cond){message("List value is empty")})
#           PTV2<-tryCatch({stringr::str_trim(PTV1)},error=function(cond){message("List value is empty")})
#           P1TV2<-tryCatch({stringr::str_trim(gsub("ADDUCTIONNAME:","",PTV2))},error=function(cond){message("List value is empty")})
# 	  PTV3<-tryCatch({paste("PRECURSORTYPE:",P1TV2)},error=function(cond){message("List value is empty")})
#           ##PTV3<-tryCatch({paste("PRECURSORTYPE:",PTV2)},error=function(cond){message("List value is empty")})
#           ####################################################
#           print("enter the test area...1")
# 	  print(PTV1)
#           print(P1TV2)
#           print(as.character(InMEDA[["Adduct"]]))
          ######################## adding this new############
          #####################################################
          ###if(!sjmisc::is_empty(P1TV2) || !sjmisc::is_empty(as.character(InMEDA[["Adduct"]]))){
          #####################################################
          ##if(identical(P1TV2,as.character(InMEDA[["Adduct"]]))){
            ##############################################
            print("entering the line 234")
            ##############################################
            FNA1<-which(stringi::stri_detect_fixed(FNA,"NAME:"))
            FNA2<-tryCatch({InMEDA[["Name"]]},error=function(cond){message("name is empty")})
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
            ###out<-c(out,PTV3)
            out<-c(out,F1PT1)
            ###########################
            FIN1<-tryCatch({InMEDA[["Ionization mode"]]},error=function(cond){message("name is empty")})
            F1IN1<-as.character(FIN1)
            F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
            out<-c(out,F2IN1)
            ########################### need to add this new ## this is the ontology part
            ##print("check if problem in ontology part")
	    ########### Adding this new ####################
	    ################################################
            FINIKOT<-MaKE.ONT.REC(InMEDA)
	    ##############################
	    ##############################
	    ##print("entering the Ontology test area...checking what it prints...1")
	    ##print(FINIKOT)
	    #####################
	    if(!sjmisc::is_empty(FINIKOT)){
	    #####################
            out<-c(out,FINIKOT)
	    ####################
	    }else{
		    ###########################
		    InKeyVal1<-which(stringi::stri_detect_fixed(FNA,"INCHIKEY:"))
		    InKeyVal2<-FNA[InKeyVal1]
		    InKeyVal3<-gsub("INCHIKEY:","",InKeyVal2)
		    InKeyVal4<-str_trim(gsub("INCHIKEY:","",InKeyVal3))
		    InKeyVal<-InKeyVal4
		    ################################
		    IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},error=function(cond){message("Classifier could not fecth the information")})
		    ################################
		    if(!sjmisc::is_empty(IKCRV)){
  			ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
  			F1ONT<-paste("Ontology:",ONTV,sep=" ")
  			out<-c(out,F1ONT)
		   }else{
  			F1ONT<-paste("Ontology:","",sep=" ")
  			out<-c(out,F1ONT)
		    }
		    ###############################
		    FINK<-paste("INCHIKEY:",tryCatch({InKeyVal},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
		    out<-c(out,FINK)
		    ############################
		    INV<-which(stringi::stri_detect_fixed(FNA,"INCHI:"))
		    INV1<-FNA[INV]
		    ###########################
		    INV2<-gsub("INCHI:","",INV1)
		    INV3<-str_trim(INV2)
		    FINCH<-paste("INCHI:",INV3,sep=" ")
		    out<-c(out,FINCH)
		    ######################################

	    }
            ################################################
            ############ this is new #######################
            ################################################
	    FSIM<-paste("SMILES:",gETSmiles(InMEDA),sep=" ")
            out<-c(out,FSIM)
            ###################################
            FFOR<-InMEDA[["Formula"]]
            FFOR1<-paste("FORMULA:",FFOR,sep=" ")
            out<-c(out,FFOR1)
	    ################################
            ################################
            FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
            FINS1<-FNA[FINS]
	    FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:","",sep=""))
            out<-c(out,FINS2)
            ###############
	    ##print("checking what it is printing")
	    ##print(FINS1)
	    ##print(FINS2)
	    ########################
            ##print("checking the intensity is working")
            ##print(FNA)
            ##print(FINS1)
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
	#########################################################
        ################### Need to comment before thsi ###################
         #### }
        ###################################################
    ###########################################################      
#           else{
# 	     #################################	  
# 	     print("enter the else part...1")
#             ####################################
# 	    ####################################
# 	    if(as.character(InMEDA[["Adduct"]]) == "[M]+" || as.character(InMEDA[["Adduct"]]) == "[M]-"){
# 	    #####################################
#             ####################################
#             FNA1<-which(stringi::stri_detect_fixed(FNA,"NAME:"))
#             FNA2<-InMEDA[["Name"]]
#             FNA3<-as.character(FNA2)
#             #######################
#             FNAM<-paste("NAME:",FNA3,sep=" ")
#             out<-c(out,FNAM)
#             ########################
#             FRA1<-which(stringi::stri_detect_fixed(FNA,"RETENTIONTIME:"))
#             F1RA1<-FNA[FRA1]
#             out<-c(out,F1RA1)
#             ################################
#             FMZ1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORMZ:"))
#             F1MZ1<-FNA[FMZ1]
#             out<-c(out,F1MZ1)
#             ################################
#             FPT1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORTYPE:"))
#             F1PT1<-FNA[FPT1]
#             out<-c(out,PTV3)
#             #out<-c(out,F1PT1)
#             #################################
#             FIN1<-InMEDA[["Ionization mode"]]
#             F1IN1<-as.character(FIN1)
#             F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
#             out<-c(out,F2IN1)
#             ##################################################
#             ##################################################
# 	    ###################################
# 	    if(!sjmisc::is_empty(FINIKOT)){
# 		    #############################
# 		    out<-c(out,FINIKOT)
# 		    #############################
# 	    }else{
# 		    ################################
# 		    InKeyVal1<-which(stringi::stri_detect_fixed(FNA,"INCHIKEY:"))
# 		    InKeyVal2<-FNA[InKeyVal1]
# 		    InKeyVal3<-gsub("INCHIKEY:","",InKeyVal2)
# 		    InKeyVal4<-str_trim(gsub("INCHIKEY:","",InKeyVal3))
# 		    InKeyVal<-InKeyVal4
# 		    ################################
# 		    IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
# 		    #############################
# 		    if(!sjmisc::is_empty(IKCRV)){
# 			    ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
# 			    F1ONT<-paste("Ontology:",ONTV,sep=" ")
# 			    out<-c(out,F1ONT)
# 	            }else{
# 			    F1ONT<-paste("Ontology:","",sep=" ")
# 			    out<-c(out,F1ONT)
# 			    }
# 		    ##########################
# 		    FINK<-paste("INCHIKEY:",tryCatch({InKeyVal},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
# 		    out<-c(out,FINK)
# 		    ################################
# 		    INV<-which(stringi::stri_detect_fixed(FNA,"INCHI:"))
# 		    INV1<-FNA[INV]
# 		    ###################
# 		    INV2<-gsub("INCHI:","",INV1)
# 		    INV3<-str_trim(INV2)
# 		    FINCH<-paste("INCHI:",INV3,sep=" ")
# 		    out<-c(out,FINCH)
# 		    #############################
# 
# 
# 	    }
# 	    ###########################
# 	    ###########################
# 	    FSIM<-paste("SMILES:",gETSmiles(InMEDA),sep=" ")
# 	    out<-c(out,FSIM)
# 	    #######################
# 	    FFOR<-InMEDA[["Formula"]]
# 	    FFOR1<-paste("FORMULA:",FFOR,sep=" ")
# 	    out<-c(out,FFOR1)
# 	    #############################
#             #############################
#             FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
#             FINS1<-FNA[FINS]
# 	    FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:","",sep=""))
#             out<-c(out,FINS2)
#             ############################
#             FAUT<-as.character(InMEDA[["Authors"]])
#             FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
#             out<-c(out,FAUT1)
#             ##########################
#             ##FLIC<-paste("LICENSE:",sep=" ")
# 	    #############################
#             FLIC<-paste("LICENSE:","CC BY",sep=" ")
#             out<-c(out,FLIC)
#             ###########################
#             FCIE<-paste("COLLISIONENERGY:",as.character(InMEDA[["Collision energy"]]),sep=" ")
#             out<-c(out,FCIE)
#             #########################
#             FINST<-as.character(InMEDA[["INSTRUMENT_TYPE"]])
#             FINST1<-paste("INSTRUMENTTYPE:",FINST,sep=" ")
#             out<-c(out,FINST1)
#             ########################
#             FINS<-as.character(InMEDA[["INSTRUMENT"]])
#             FINS1<-paste("INSTRUMENT:",FINS,sep=" ")
#             out<-c(out,FINS1)
#             ####################
#             FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
#             out<-c(out,FCOM)
#             ##################
#             FNPA<-which(stringi::stri_detect_fixed(FNA,"Num Peaks:"))
#             F1NPA<-FNA[FNPA]
#             out<-c(out,F1NPA)
#             ###################
#             Fpea<-FNA[(FNPA+1):Find]
#             #########################
#             if(is.na(Fpea))
#             {
#               Fpea1<-FNA[(FNPA+1)]
#               
#               
#             }else{
#               
#               MV=AAMS1
#               tes1<-unlist(strsplit(Fpea, "\t"))
#               tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
#               tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
#               tes4<-which(tes2 > (3+MV))
#               if(length(tes4)>1)
#               {
#                 tes5<-tes2[-tes4]
#                 tes6<-tes3[-tes4]
#                 tes7<-paste(tes5,tes6,sep="\t")
#                 out<-c(out,tes7)
#               }else{
#                 out<-c(out,Fpea)
#                 
#               }
#             }
#           ###################################
# 	    }
#   ################################################
#   #################################################          
#           }
  #########################################################        
	#########################################################  
        #########################################################
        }else if(length(INLL) > 1){
          ##############################
          print("enter the line ...511")
          ##############################
          MONMS=InMSPL[INLL]
          TRA<-unname(rapply(MONMS, function(x) grep("RETENTIONTIME:",x, value=TRUE)))
          TRA1<-as.numeric(stringr::str_trim(stringr::str_replace(TRA, "RETENTIONTIME:", "")))
          TRA2<-abs(VRT-TRA1)
          TRA3<-which.min(TRA2)
          TRA4<-INLL[TRA3]
          TRA5<-InMSPL[TRA4]
          ###########################
          F1FPL<-TRA5
          F2FPL<-F1FPL
          ##########################
          Find<-tryCatch({length(F2FPL[[1]])},error=function(cond){message("there is an error in list")})
          ###########################
          FNA<-tryCatch({F1FPL[[1]]},error=function(cond){message("there is an error in list")})
          ####### adding this new ###################################################
          ###########################################################################
#           PT <- c("PRECURSORTYPE:", "ADDUCTIONNAME:")
#           NPT<-tryCatch({grep(paste(PT,collapse="|"), FNA, value=TRUE)},error=function(cond){message("List value is empty")})
#           NPT1<-tryCatch({match(NPT,FNA)},error=function(cond){message("List value is empty")})
#           PTV <-tryCatch({stringr::str_remove(FNA[NPT1],c("PRECURSORTYPE:","ADDUCTIONNAME:"))},error=function(cond){message("List value is empty")})
#           PTV1<-tryCatch({PTV[1]},error=function(cond){message("List value is empty")})
#           PTV2<-tryCatch({stringr::str_trim(PTV1)},error=function(cond){message("List value is empty")})
#           P1TV2<-tryCatch({stringr::str_trim(gsub("ADDUCTIONNAME:","",PTV2))},error=function(cond){message("List value is empty")})
# 	  PTV3<-tryCatch({paste("PRECURSORTYPE:",P1TV2)},error=function(cond){message("List value is empty")})
#           ##PTV3<-tryCatch({paste("PRECURSORTYPE:",PTV2)},error=function(cond){message("List value is empty")})
#           ##########################################################################
#           print("enter the test area..2")
#           print(P1TV2)
#           print(as.character(InMEDA[["Adduct"]]))
          ##################################################
          #####################################################
          ###if(!sjmisc::is_empty(P1TV2) || !sjmisc::is_empty(as.character(InMEDA[["Adduct"]]))){
          ##################################################		  
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
            ####out<-c(out,PTV3)
            out<-c(out,F1PT1)
            ############################################################
            FIN1<-InMEDA[["Ionization mode"]]
            F1IN1<-as.character(FIN1)
            F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
            out<-c(out,F2IN1)
	    ##########################################################
            ##################### adding this new ####################
            FINIKOT<-MaKE.ONT.REC(InMEDA)
            ##out<-c(out,FINIKOT)
            ###################### this is end########################
            ############################################### this is new
            ##FSIM<-paste("SMILES:",as.character(InMEDA[["SMILES"]]),sep=" ")
            #####################################################
            #####################################################
            if(!sjmisc::is_empty(FINIKOT)){
              #####################
              out<-c(out,FINIKOT)
              ####################
            }else{
              
              ###############################################
              InKeyVal1<-which(stringi::stri_detect_fixed(FNA,"INCHIKEY:"))
              InKeyVal2<-FNA[InKeyVal1]
              InKeyVal3<-gsub("INCHIKEY:","",InKeyVal2)
              InKeyVal4<-str_trim(gsub("INCHIKEY:","",InKeyVal3))
              InKeyVal<-InKeyVal4
              ################################
              ###############################
              IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
              #################################################
              if(!sjmisc::is_empty(IKCRV)){
                ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
                F1ONT<-paste("Ontology:",ONTV,sep=" ")
                out<-c(out,F1ONT)
              }else{
                F1ONT<-paste("Ontology:","",sep=" ")
                out<-c(out,F1ONT)
              }
              ##################################
              FINK<-paste("INCHIKEY:",tryCatch({InKeyVal},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
              out<-c(out,FINK)
              ################################
              INV<-which(stringi::stri_detect_fixed(FNA,"INCHI:"))
              INV1<-FNA[INV]
              ###################
              INV2<-gsub("INCHI:","",INV1)
              INV3<-str_trim(INV2)
              FINCH<-paste("INCHI:",INV3,sep=" ")
              out<-c(out,FINCH)
              ###################
            }
            #########################################################
	    ##########################################################
	    FSIM<-paste("SMILES:",gETSmiles(InMEDA),sep=" ")
            out<-c(out,FSIM)
            ##########################################################
            FFOR<-InMEDA[["Formula"]]
            FFOR1<-paste("FORMULA:",FFOR,sep=" ")
            out<-c(out,FFOR1)
	    ##########################################################
            ##########################################################
            FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
            FINS1<-FNA[FINS]
	    FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:","",sep=""))
            out<-c(out,FINS2)
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
        ############################## adding this new### commenting below this
	 #############################################
        ###  }
    #############################################
    ###################################################      
#           else{
# 	   ####################################
# 	    print("enter the else part...2")
#             ###########################################
# 	    ###########################################
# 	    if(as.character(InMEDA[["Adduct"]]) == "[M]+" || as.character(InMEDA[["Adduct"]]) == "[M]-"){
#             ############################################		    
#             ############################################
#             FNA1<-which(stringi::stri_detect_fixed(FNA,"NAME:"))
#             FNA2<-InMEDA[["Name"]]
#             FNA3<-as.character(FNA2)
#             #######################
#             FNAM<-paste("NAME:",FNA3,sep=" ")
#             out<-c(out,FNAM)
#             ########################
#             FRA1<-which(stringi::stri_detect_fixed(FNA,"RETENTIONTIME:"))
#             F1RA1<-FNA[FRA1]
#             out<-c(out,F1RA1)
#             ################################
#             FMZ1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORMZ:"))
#             F1MZ1<-FNA[FMZ1]
#             out<-c(out,F1MZ1)
#             ################################
#             FPT1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORTYPE:"))
#             F1PT1<-FNA[FPT1]
#             out<-c(out,PTV3)
#             #out<-c(out,F1PT1)
#             #################################
#             FIN1<-InMEDA[["Ionization mode"]]
#             F1IN1<-as.character(FIN1)
#             F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
#             out<-c(out,F2IN1)
# 	    ###########################################################
#             ############### Adding this new####################################
# 	    FINIKOT<-MaKE.ONT.REC(InMEDA)
# 	    ##out<-c(out,FINIKOT)
# 	    ###################
#             ######################################################
#             if(!sjmisc::is_empty(FINIKOT)){
#               #####################
#               out<-c(out,FINIKOT)
#               ####################
#             }else{
#               
#               ###############################################
#               InKeyVal1<-which(stringi::stri_detect_fixed(FNA,"INCHIKEY:"))
#               InKeyVal2<-FNA[InKeyVal1]
#               InKeyVal3<-gsub("INCHIKEY:","",InKeyVal2)
#               InKeyVal4<-str_trim(gsub("INCHIKEY:","",InKeyVal3))
#               InKeyVal<-InKeyVal4
#               ################################
#               ###############################
#               IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
#               #################################################
#               if(!sjmisc::is_empty(IKCRV)){
#                 ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
#                 F1ONT<-paste("Ontology:",ONTV,sep=" ")
#                 out<-c(out,F1ONT)
#               }else{
#                 F1ONT<-paste("Ontology:","",sep=" ")
#                 out<-c(out,F1ONT)
#               }
#               ##################################
#               FINK<-paste("INCHIKEY:",tryCatch({InKeyVal},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
#               out<-c(out,FINK)
#               ################################
#               INV<-which(stringi::stri_detect_fixed(FNA,"INCHI:"))
#               INV1<-FNA[INV]
#               ###################
#               INV2<-gsub("INCHI:","",INV1)
#               INV3<-str_trim(INV2)
#               FINCH<-paste("INCHI:",INV3,sep=" ")
#               out<-c(out,FINCH)
#               ###################
#             }
#       ##############################
#             #########################
# 	    FSIM<-paste("SMILES:",gETSmiles(InMEDA),sep=" ")
# 	    out<-c(out,FSIM)
# 	    ####################################
# 	    FFOR<-InMEDA[["Formula"]]
# 	    FFOR1<-paste("FORMULA:",FFOR,sep=" ")
# 	    out<-c(out,FFOR1)
#             
#             ####################################################
#             ###################################################
#             FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
#             FINS1<-FNA[FINS]
# 	    FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:","",sep=""))
#             out<-c(out,FINS2)
#             ############################
#             FAUT<-as.character(InMEDA[["Authors"]])
#             FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
#             out<-c(out,FAUT1)
#             ##########################
#             ##FLIC<-paste("LICENSE:",sep=" ")
#             FLIC<-paste("LICENSE:","CC BY",sep=" ")
#             out<-c(out,FLIC)
#             ###########################
#             FCIE<-paste("COLLISIONENERGY:",as.character(InMEDA[["Collision energy"]]),sep=" ")
#             out<-c(out,FCIE)
#             #########################
#             FINST<-as.character(InMEDA[["INSTRUMENT_TYPE"]])
#             FINST1<-paste("INSTRUMENTTYPE:",FINST,sep=" ")
#             out<-c(out,FINST1)
#             ########################
#             FINS<-as.character(InMEDA[["INSTRUMENT"]])
#             FINS1<-paste("INSTRUMENT:",FINS,sep=" ")
#             out<-c(out,FINS1)
#             ####################
#             FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
#             out<-c(out,FCOM)
#             ##################
#             FNPA<-which(stringi::stri_detect_fixed(FNA,"Num Peaks:"))
#             F1NPA<-FNA[FNPA]
#             out<-c(out,F1NPA)
#             ###################
#             Fpea<-FNA[(FNPA+1):Find]
#             #########################
#             if(is.na(Fpea))
#             {
#               Fpea1<-FNA[(FNPA+1)]
#               
#               
#             }else{
#               
#               MV=AAMS1
#               tes1<-unlist(strsplit(Fpea, "\t"))
#               tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
#               tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
#               tes4<-which(tes2 > (3+MV))
#               if(length(tes4)>1)
#               {
#                 tes5<-tes2[-tes4]
#                 tes6<-tes3[-tes4]
#                 tes7<-paste(tes5,tes6,sep="\t")
#                 out<-c(out,tes7)
#               }else{
#                 out<-c(out,Fpea)
#                 
#               }
#             }
#          ############################
# 	     }
# #########################################
# ############### Comenting below this 
# 	  }
 #########################################
  #########################################          
        } else {
          #print("entering the line 458")
          PASS<-RRV
        }  ## main else part
        ###testing if this works
        return(out)
        #####################
      } #ITmass
    } ##IITRTL
  }else{
    ##############################################	  
    ##############################################	  
    ### that means not able to get mass from smiles,pubchem ID and
    ##########################################
    print("entering this line..line 771")
    ###########################################
    FM<-as.character(InMEDA[["Exact mass"]])
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
	##ITmass<-match(Tmass,InPMZ)
        #####################
        VRT<-as.numeric(as.character(InMEDA[["RT (min)"]]))
        VRTL<-VRT-0.20
        VRTU<-VRT+0.20
	#####################
	print("enter the function 2 ..else part")
	##print(FM)
	##print(ADV4)
	##print(MPPmL)
	##print(MPPmU)
	##print(VRT)
	##print(VRTL)
	##print(VRTU)
        ######################
        TRTL<-InRTL[InRTL >=VRTL & InRTL <= VRTU]
        ITRTL<-which(InRTL %in% TRTL)
	#############################
	##ITRTL<-match(TRTL,InRTL)
        ######################
        INLL<-intersect(ITmass,ITRTL)
        ########################################
	print("entering the line 513")
        if(length(ITRTL) >= 1){
          print("entering the line 515")
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
#               PT <- c("PRECURSORTYPE:", "ADDUCTIONNAME:")
#               NPT<-tryCatch({grep(paste(PT,collapse="|"), FNA, value=TRUE)},error=function(cond){message("List value is empty")})
#               NPT1<-tryCatch({match(NPT,FNA)},error=function(cond){message("List value is empty")})
#               PTV <-tryCatch({stringr::str_remove(FNA[NPT1],c("PRECURSORTYPE:","ADDUCTIONNAME:"))},error=function(cond){message("List value is empty")})
#               PTV1<-tryCatch({PTV[1]},error=function(cond){message("List value is empty")})
#               PTV2<-tryCatch({stringr::str_trim(PTV1)},error=function(cond){message("List value is empty")})
#               P1TV2<-tryCatch({stringr::str_trim(gsub("ADDUCTIONNAME:","",PTV2))},error=function(cond){message("List value is empty")})
# 	      PTV3<-tryCatch({paste("PRECURSORTYPE:",P1TV2)},error=function(cond){message("List value is empty")})
#               ##PTV3<-tryCatch({paste("PRECURSORTYPE:",PTV2)},error=function(cond){message("List value is empty")})
#               #####################################################
#               print("enter the test area..3")
# 	      ####################################################
#               print(P1TV2)
#               print(as.character(InMEDA[["Adduct"]]))
              ################################ adding this new ###
              ##if(identical(P1TV2,as.character(InMEDA[["Adduct"]]))){
	      ###############################################
        ####################################################      
           ###   if(!sjmisc::is_empty(P1TV2) || !sjmisc::is_empty(as.character(InMEDA[["Adduct"]]))){
            ##################################################
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
                ###out<-c(out,PTV3)
                out<-c(out,F1PT1)
                #################################
                FIN1<-InMEDA[["Ionization mode"]]
                F1IN1<-as.character(FIN1)
                F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                out<-c(out,F2IN1)
		###############################################
		###############################################
		InKeyVal1<-which(stringi::stri_detect_fixed(FNA,"INCHIKEY:"))
	        InKeyVal2<-FNA[InKeyVal1]
		InKeyVal3<-gsub("INCHIKEY:","",InKeyVal2)
		InKeyVal4<-str_trim(gsub("INCHIKEY:","",InKeyVal3))
		InKeyVal<-InKeyVal4
		################################
		###############################
		IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
		#################################################
		if(!sjmisc::is_empty(IKCRV)){
		      ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
		      F1ONT<-paste("Ontology:",ONTV,sep=" ")
		      out<-c(out,F1ONT)
	        }else{
		      F1ONT<-paste("Ontology:","",sep=" ")
		      out<-c(out,F1ONT)
	        }
               ##################################
		FINK<-paste("INCHIKEY:",tryCatch({InKeyVal},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
		out<-c(out,FINK)
                ################################
		INV<-which(stringi::stri_detect_fixed(FNA,"INCHI:"))
                INV1<-FNA[INV]
		###################
		INV2<-gsub("INCHI:","",INV1)
		INV3<-str_trim(INV2)
	        FINCH<-paste("INCHI:",INV3,sep=" ")
                out<-c(out,FINCH)
                ###################
	        ################################
                ########################### adding this new ###
                ##FINIKOT<-MaKE.ONT.REC(InMEDA)
                ##out<-c(out,FINIKOT)
                ########################## this is the end #####
                ################################################
		FSIM<-paste("SMILES:",gETSmiles(InMEDA),sep=" ")
                out<-c(out,FSIM)
                ################################
                FFOR<-InMEDA[["Formula"]]
                FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                out<-c(out,FFOR1)
                ###############################
		##############################
                FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                FINS1<-FNA[FINS]
		FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:","",sep=""))
                out<-c(out,FINS2)
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
	    #################################### Before this comment this 
             ####### }
          ######################################
          ############## Commenting this ##############    
#               else{
# 		  #########################################    
# 		      print("enter the else part...3")
#                 ##########################################
# 		###########################################      
# 		if(as.character(InMEDA[["Adduct"]]) == "[M]+" || as.character(InMEDA[["Adduct"]]) == "[M]-"){
#                 ###########################################
# 		###########################################
# 			print("enter the else part...3..if loop")
# 		###########################
#                 FNA1<-which(stringi::stri_detect_fixed(FNA,"NAME:"))
#                 FNA2<-InMEDA[["Name"]]
#                 FNA3<-as.character(FNA2)
#                 #######################
#                 FNAM<-paste("NAME:",FNA3,sep=" ")
#                 out<-c(out,FNAM)
#                 ########################
#                 FRA1<-which(stringi::stri_detect_fixed(FNA,"RETENTIONTIME:"))
#                 F1RA1<-FNA[FRA1]
#                 out<-c(out,F1RA1)
#                 ################################
#                 FMZ1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORMZ:"))
#                 F1MZ1<-FNA[FMZ1]
#                 out<-c(out,F1MZ1)
#                 ################################
#                 FPT1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORTYPE:"))
#                 F1PT1<-FNA[FPT1]
#                 out<-c(out,PTV3)
#                 #out<-c(out,F1PT1)
#                 #################################
#                 FIN1<-InMEDA[["Ionization mode"]]
#                 F1IN1<-as.character(FIN1)
#                 F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
#                 out<-c(out,F2IN1)
#                 #################Adding this new#################################
#                 ##################################################
# 		##FINIKOT<-MaKE.ONT.REC(InMEDA)
# 		##out<-c(out,FINIKOT)
# 		#################################
# 		InKeyVal1<-which(stringi::stri_detect_fixed(FNA,"INCHIKEY:"))
# 		InKeyVal2<-FNA[InKeyVal1]
# 		InKeyVal3<-gsub("INCHIKEY:","",InKeyVal2)
# 		InKeyVal4<-str_trim(gsub("INCHIKEY:","",InKeyVal3))
# 		InKeyVal<-InKeyVal4
# 		########################
# 		IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
# 		######################
# 		if(!sjmisc::is_empty(IKCRV)){
# 			ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
# 			F1ONT<-paste("Ontology:",ONTV,sep=" ")
# 			out<-c(out,F1ONT)
# 		}else{
# 			F1ONT<-paste("Ontology:","",sep=" ")
# 			out<-c(out,F1ONT)
# 		}
# 		######################
# 		FINK<-paste("INCHIKEY:",tryCatch({InKeyVal},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
# 		out<-c(out,FINK)
# 		#####################
# 		INV<-which(stringi::stri_detect_fixed(FNA,"INCHI:"))
# 		INV1<-FNA[INV]
# 		###################
# 		INV2<-gsub("INCHI:","",INV1)
# 		INV3<-str_trim(INV2)
# 		FINCH<-paste("INCHI:",INV3,sep=" ")
# 		out<-c(out,FINCH)
# 		##########################
# 		###########################
# 		FSIM<-paste("SMILES:",gETSmiles(InMEDA),sep=" ")
# 		out<-c(out,FSIM)
# 		#######################
# 		FFOR<-InMEDA[["Formula"]]
# 		FFOR1<-paste("FORMULA:",FFOR,sep=" ")
# 		out<-c(out,FFOR1) 
#                 ###################################################		
#                 ###################################################
# 		#############################
#                 FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
#                 FINS1<-FNA[FINS]
# 		FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:","",sep=""))
#                 out<-c(out,FINS2)
#                 ############################
#                 FAUT<-as.character(InMEDA[["Authors"]])
#                 FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
#                 out<-c(out,FAUT1)
#                 ##########################
#                 ##FLIC<-paste("LICENSE:",sep=" ")
#                 FLIC<-paste("LICENSE:","CC BY",sep=" ")
#                 out<-c(out,FLIC)
#                 ###########################
#                 FCIE<-paste("COLLISIONENERGY:",as.character(InMEDA[["Collision energy"]]),sep=" ")
#                 out<-c(out,FCIE)
#                 #########################
#                 FINST<-as.character(InMEDA[["INSTRUMENT_TYPE"]])
#                 FINST1<-paste("INSTRUMENTTYPE:",FINST,sep=" ")
#                 out<-c(out,FINST1)
#                 ########################
#                 FINS<-as.character(InMEDA[["INSTRUMENT"]])
#                 FINS1<-paste("INSTRUMENT:",FINS,sep=" ")
#                 out<-c(out,FINS1)
#                 ####################
#                 FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
#                 out<-c(out,FCOM)
#                 ##################
#                 FNPA<-which(stringi::stri_detect_fixed(FNA,"Num Peaks:"))
#                 F1NPA<-FNA[FNPA]
#                 out<-c(out,F1NPA)
#                 ###################
#                 Fpea<-FNA[(FNPA+1):Find]
#                 #########################
#                 if(is.na(Fpea))
#                 {
#                   Fpea1<-FNA[(FNPA+1)]
#                   
#                   
#                 }else{
#                   
#                   MV=AAMS1
#                   tes1<-unlist(strsplit(Fpea, "\t"))
#                   tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
#                   tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
#                   tes4<-which(tes2 > (3+MV))
#                   if(length(tes4)>1)
#                   {
#                     tes5<-tes2[-tes4]
#                     tes6<-tes3[-tes4]
#                     tes7<-paste(tes5,tes6,sep="\t")
#                     out<-c(out,tes7)
#                   }else{
#                     out<-c(out,Fpea)
#                     
#                   }
#                 }
#              ############################
# 		 }
#     ########################### Commenthing this below
#     ############################################
# 	      }
    #########################################################
   ###########################################################              
            }else if(length(INLL) > 1){
              ##############################################
              print("enter the line ...1020")
              ##############################################
              MONMS=InMSPL[INLL]
              TRA<-unname(rapply(MONMS, function(x) grep("RETENTIONTIME:",x, value=TRUE)))
              TRA1<-as.numeric(stringr::str_trim(stringr::str_replace(TRA, "RETENTIONTIME:", "")))
              TRA2<-abs(VRT-TRA1)
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
#               PT <- c("PRECURSORTYPE:", "ADDUCTIONNAME:")
#               NPT<-tryCatch({grep(paste(PT,collapse="|"), FNA, value=TRUE)},error=function(cond){message("List value is empty")})
#               NPT1<-tryCatch({match(NPT,FNA)},error=function(cond){message("List value is empty")})
#               PTV <-tryCatch({stringr::str_remove(FNA[NPT1],c("PRECURSORTYPE:","ADDUCTIONNAME:"))},error=function(cond){message("List value is empty")})
#               PTV1<-tryCatch({PTV[1]},error=function(cond){message("List value is empty")})
#               PTV2<-tryCatch({stringr::str_trim(PTV1)},error=function(cond){message("List value is empty")})
#               P1TV2<-tryCatch({stringr::str_trim(gsub("ADDUCTIONNAME:","",PTV2))},error=function(cond){message("List value is empty")})
#               PTV3<-tryCatch({paste("PRECURSORTYPE:",P1TV2)},error=function(cond){message("List value is empty")})
# 	      ##PTV3<-tryCatch({paste("PRECURSORTYPE:",P1TV2)},error=function(cond){message("List value is empty")})
#               ########################################
#               print("enter the test area..4")
# 	      ########################################
#               print(P1TV2)
#               print(as.character(InMEDA[["Adduct"]]))
        ################################## adding this new##########
        ###############################################################      
	     ### if(as.character(InMEDA[["Adduct"]]) == "[M]+" || as.character(InMEDA[["Adduct"]]) == "[M]-"){
	      ###################################################	      
              ##if(identical(P1TV2,as.character(InMEDA[["Adduct"]]))){
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
                ###out<-c(out,PTV3)
                out<-c(out,F1PT1)
                ###########################
                FIN1<-InMEDA[["Ionization mode"]]
                F1IN1<-as.character(FIN1)
                F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                out<-c(out,F2IN1)
		#################################################
                ########################### adding this new #####
                ##FINIKOT<-MaKE.ONT.REC(InMEDA)
                ##out<-c(out,FINIKOT)
		######################################
		######################################
		InKeyVal1<-which(stringi::stri_detect_fixed(FNA,"INCHIKEY:"))
		InKeyVal2<-FNA[InKeyVal1]
		InKeyVal3<-gsub("INCHIKEY:","",InKeyVal2)
		InKeyVal4<-str_trim(gsub("INCHIKEY:","",InKeyVal3))
		InKeyVal<-InKeyVal4
		##############################
		IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
		##################################
		if(!sjmisc::is_empty(IKCRV)){
  			ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
  			F1ONT<-paste("Ontology:",ONTV,sep=" ")
  			out<-c(out,F1ONT)
		}else{
  			F1ONT<-paste("Ontology:","",sep=" ")
  			out<-c(out,F1ONT)
		}
		###############################
		INV<-which(stringi::stri_detect_fixed(FNA,"INCHI:"))
		INV1<-FNA[INV]
		######################
		INV2<-gsub("INCHI:","",INV1)
		INV3<-str_trim(INV2)
		FINCH<-paste("INCHI:",INV3,sep=" ")
		out<-c(out,FINCH)
		######################################
		#############################
		FSIM<-paste("SMILES:",gETSmiles(InMEDA),sep=" ")
                out<-c(out,FSIM)
                ############################################
                FFOR<-InMEDA[["Formula"]]
                FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                out<-c(out,FFOR1)
                ############################################
		############################################
                FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                FINS1<-FNA[FINS]
		FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:","",sep=""))
                out<-c(out,FINS2)
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
      ##########################################
	    ############################ This is commenting
	    ###  }
    ##############################################
    ###########################################          
  ##            else{
		 #####################################     
# 		      print("enter the else part...4")
#                 #####################################
# 		#####################################
# 		if(as.character(InMEDA[["Adduct"]]) == "[M]+" || as.character(InMEDA[["Adduct"]]) == "[M]-"){
#                 #####################################
# 		#####################################
#                 FNA1<-which(stringi::stri_detect_fixed(FNA,"NAME:"))
#                 FNA2<-InMEDA[["Name"]]
#                 FNA3<-as.character(FNA2)
#                 #######################
#                 FNAM<-paste("NAME:",FNA3,sep=" ")
#                 out<-c(out,FNAM)
#                 ########################
#                 FRA1<-which(stringi::stri_detect_fixed(FNA,"RETENTIONTIME:"))
#                 F1RA1<-FNA[FRA1]
#                 out<-c(out,F1RA1)
#                 ################################
#                 FMZ1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORMZ:"))
#                 F1MZ1<-FNA[FMZ1]
#                 out<-c(out,F1MZ1)
#                 ################################
#                 FPT1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORTYPE:"))
#                 F1PT1<-FNA[FPT1]
#                 out<-c(out,PTV3)
#                 #out<-c(out,F1PT1)
#                 #################################
#                 FIN1<-InMEDA[["Ionization mode"]]
#                 F1IN1<-as.character(FIN1)
#                 F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
#                 out<-c(out,F2IN1)
# 		########################################
#                 ################ Adding this new #######
#                 ##FINIKOT<-MaKE.ONT.REC(InMEDA)
# 		##out<-c(out,FINIKOT)
# 		##########################################
# 		##########################################
# 		InKeyVal1<-which(stringi::stri_detect_fixed(FNA,"INCHIKEY:"))
# 		InKeyVal2<-FNA[InKeyVal1]
# 		InKeyVal3<-gsub("INCHIKEY:","",InKeyVal2)
# 		InKeyVal4<-str_trim(gsub("INCHIKEY:","",InKeyVal3))
# 		InKeyVal<-InKeyVal4
# 		###########################################
# 		IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
# 		#########################
# 		if(!sjmisc::is_empty(IKCRV)){
#   			ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
#   			F1ONT<-paste("Ontology:",ONTV,sep=" ")
#   			out<-c(out,F1ONT)
# 		}else{
#   			F1ONT<-paste("Ontology:","",sep=" ")
#   			out<-c(out,F1ONT)
# 		}
# 		####################
# 		FINK<-paste("INCHIKEY:",tryCatch({InKeyVal},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
# 		out<-c(out,FINK)
# 		#####################
# 		INV<-which(stringi::stri_detect_fixed(FNA,"INCHI:"))
# 		INV1<-FNA[INV]
# 		#######################
# 		INV2<-gsub("INCHI:","",INV1)
# 		INV3<-str_trim(INV2)
# 		FINCH<-paste("INCHI:",INV3,sep=" ")
# 		out<-c(out,FINCH)
# 		########################
# 		#####################
# 		FSIM<-paste("SMILES:",gETSmiles(InMEDA),sep=" ")
# 		out<-c(out,FSIM)
# 		###############
# 		FFOR<-InMEDA[["Formula"]]
# 		FFOR1<-paste("FORMULA:",FFOR,sep=" ")
# 		out<-c(out,FFOR1)
#                 ###################################################		
#                 ###################################################
#                 ##################################
# 		##################################
#                 FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
#                 FINS1<-FNA[FINS]
# 		FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:","",sep=""))
#                 out<-c(out,FINS2)
#                 ############################
#                 FAUT<-as.character(InMEDA[["Authors"]])
#                 FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
#                 out<-c(out,FAUT1)
#                 ##########################
#                 ##FLIC<-paste("LICENSE:",sep=" ")
#                 FLIC<-paste("LICENSE:","CC BY",sep=" ")
#                 out<-c(out,FLIC)
#                 ###########################
#                 FCIE<-paste("COLLISIONENERGY:",as.character(InMEDA[["Collision energy"]]),sep=" ")
#                 out<-c(out,FCIE)
#                 #########################
#                 FINST<-as.character(InMEDA[["INSTRUMENT_TYPE"]])
#                 FINST1<-paste("INSTRUMENTTYPE:",FINST,sep=" ")
#                 out<-c(out,FINST1)
#                 ########################
#                 FINS<-as.character(InMEDA[["INSTRUMENT"]])
#                 FINS1<-paste("INSTRUMENT:",FINS,sep=" ")
#                 out<-c(out,FINS1)
#                 ####################
#                 FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
#                 out<-c(out,FCOM)
#                 ##################
#                 FNPA<-which(stringi::stri_detect_fixed(FNA,"Num Peaks:"))
#                 F1NPA<-FNA[FNPA]
#                 out<-c(out,F1NPA)
#                 ###################
#                 Fpea<-FNA[(FNPA+1):Find]
#                 #########################
#                 if(is.na(Fpea))
#                 {
#                   Fpea1<-FNA[(FNPA+1)]
#                   
#                   
#                 }else{
#                   
#                   MV=AAMS1
#                   tes1<-unlist(strsplit(Fpea, "\t"))
#                   tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
#                   tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
#                   tes4<-which(tes2 > (3+MV))
#                   if(length(tes4)>1)
#                   {
#                     tes5<-tes2[-tes4]
#                     tes6<-tes3[-tes4]
#                     tes7<-paste(tes5,tes6,sep="\t")
#                     out<-c(out,tes7)
#                   }else{
#                     out<-c(out,Fpea)
#                     
#                   }
#                 }
#              ############################
# 		 }
# #######################################
# #################### Need to coment this ####
# 	      }
    ##########################
  #################################              
            # } else {
            #   #print("entering the line 764")
            #   PASS<-InMEDA
            # }  ## main else part
            ###testing if this works
            return(out)
         #####################
          } #ITmass
        } #ITRTL
    ###############################
    ###########################    
      } ## adduct match and exacct mass present
    ###############################
    #######################################  
    } else{
      
      ########################################
      #### this entered the FM does not exists 
      FM2<-FuFtoRe(InMEDA)
      ##PASS<-RRV
      ###############################################################
      ###############################################################
      if(!sjmisc::is_empty(FM2))
      {
        ###################################################################
        print("entering the line ... 772")
        ###################################################################
        ##TE<-FM
        ##EM1<-OrgMassSpecR::MolecularWeight(formula = OrgMassSpecR::ListFormula(FM))
        ##if(!sjmisc::is_empty( EM1)) {
        ####################################################################
        if(!sjmisc::is_empty(FM2)) {
        ###if(!sjmisc::is_empty(FM2) & !sjmisc::is_empty(as.character(InMEDA[["Adduct"]]))) {
          ########################################################
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
          ##ITmass<-match(Tmass,InPMZ)
          #####################
          VRT<-as.numeric(as.character(InMEDA[["RT (min)"]]))
          VRTL<-VRT-0.20
          VRTU<-VRT+0.20
          #####################
          print("enter the function 2 ..else part")
          ##print(FM)
          ##print(ADV4)
          ##print(MPPmL)
          ##print(MPPmU)
          ##print(VRT)
          ##print(VRTL)
          ##print(VRTU)
          ######################
          TRTL<-InRTL[InRTL >=VRTL & InRTL <= VRTU]
          ITRTL<-which(InRTL %in% TRTL)
          #############################
          ##ITRTL<-match(TRTL,InRTL)
          ######################
          INLL<-intersect(ITmass,ITRTL)
          ########################################
          print("entering the line 513")
          if(length(ITRTL) >= 1){
            print("entering the line 515")
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
                # PT <- c("PRECURSORTYPE:", "ADDUCTIONNAME:")
                # NPT<-tryCatch({grep(paste(PT,collapse="|"), FNA, value=TRUE)},error=function(cond){message("List value is empty")})
                # NPT1<-tryCatch({match(NPT,FNA)},error=function(cond){message("List value is empty")})
                # PTV <-tryCatch({stringr::str_remove(FNA[NPT1],c("PRECURSORTYPE:","ADDUCTIONNAME:"))},error=function(cond){message("List value is empty")})
                # PTV1<-tryCatch({PTV[1]},error=function(cond){message("List value is empty")})
                # PTV2<-tryCatch({stringr::str_trim(PTV1)},error=function(cond){message("List value is empty")})
                # P1TV2<-tryCatch({stringr::str_trim(gsub("ADDUCTIONNAME:","",PTV2))},error=function(cond){message("List value is empty")})
                # PTV3<-tryCatch({paste("PRECURSORTYPE:",P1TV2)},error=function(cond){message("List value is empty")})
                # ##PTV3<-tryCatch({paste("PRECURSORTYPE:",PTV2)},error=function(cond){message("List value is empty")})
                # #####################################################
                # print("enter the test area..3")
                # ####################################################
                # print(P1TV2)
                # print(as.character(InMEDA[["Adduct"]]))
                # ################################ adding this new ###
                ##if(identical(P1TV2,as.character(InMEDA[["Adduct"]]))){
                ###############################################
                ###if(!sjmisc::is_empty(P1TV2) || !sjmisc::is_empty(as.character(InMEDA[["Adduct"]]))){
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
                  ####out<-c(out,PTV3)
                  out<-c(out,F1PT1)
                  #################################
                  FIN1<-InMEDA[["Ionization mode"]]
                  F1IN1<-as.character(FIN1)
                  F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                  out<-c(out,F2IN1)
                  ###############################################
                  ###############################################
                  InKeyVal1<-which(stringi::stri_detect_fixed(FNA,"INCHIKEY:"))
                  InKeyVal2<-FNA[InKeyVal1]
                  InKeyVal3<-gsub("INCHIKEY:","",InKeyVal2)
                  InKeyVal4<-str_trim(gsub("INCHIKEY:","",InKeyVal3))
                  InKeyVal<-InKeyVal4
                  ################################
                  ###############################
                  IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
                  #################################################
                  if(!sjmisc::is_empty(IKCRV)){
                    ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
                    F1ONT<-paste("Ontology:",ONTV,sep=" ")
                    out<-c(out,F1ONT)
                  }else{
                    F1ONT<-paste("Ontology:","",sep=" ")
                    out<-c(out,F1ONT)
                  }
                  ##################################
                  FINK<-paste("INCHIKEY:",tryCatch({InKeyVal},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                  out<-c(out,FINK)
                  ################################
                  INV<-which(stringi::stri_detect_fixed(FNA,"INCHI:"))
                  INV1<-FNA[INV]
                  ###################
                  INV2<-gsub("INCHI:","",INV1)
                  INV3<-str_trim(INV2)
                  FINCH<-paste("INCHI:",INV3,sep=" ")
                  out<-c(out,FINCH)
                  ###################
                  ################################
                  ########################### adding this new ###
                  ##FINIKOT<-MaKE.ONT.REC(InMEDA)
                  ##out<-c(out,FINIKOT)
                  ########################## this is the end #####
                  ################################################
                  FSIM<-paste("SMILES:",gETSmiles(InMEDA),sep=" ")
                  out<-c(out,FSIM)
                  ################################
                  FFOR<-InMEDA[["Formula"]]
                  FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                  out<-c(out,FFOR1)
                  ###############################
                  ##############################
                  FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                  FINS1<-FNA[FINS]
                  FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:","",sep=""))
                  out<-c(out,FINS2)
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
              ##########################################
               ### }
              ##########################################
              #########################################  
            #     else{
            #       #########################################    
            #       print("enter the else part...3")
            #       ##########################################
            #       ###########################################      
            #       if(as.character(InMEDA[["Adduct"]]) == "[M]+" || as.character(InMEDA[["Adduct"]]) == "[M]-"){
            #         ###########################################
            #         ###########################################
            #         print("enter the else part...3..if loop")
            #         ###########################
            #         FNA1<-which(stringi::stri_detect_fixed(FNA,"NAME:"))
            #         FNA2<-InMEDA[["Name"]]
            #         FNA3<-as.character(FNA2)
            #         #######################
            #         FNAM<-paste("NAME:",FNA3,sep=" ")
            #         out<-c(out,FNAM)
            #         ########################
            #         FRA1<-which(stringi::stri_detect_fixed(FNA,"RETENTIONTIME:"))
            #         F1RA1<-FNA[FRA1]
            #         out<-c(out,F1RA1)
            #         ################################
            #         FMZ1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORMZ:"))
            #         F1MZ1<-FNA[FMZ1]
            #         out<-c(out,F1MZ1)
            #         ################################
            #         FPT1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORTYPE:"))
            #         F1PT1<-FNA[FPT1]
            #         out<-c(out,PTV3)
            #         #out<-c(out,F1PT1)
            #         #################################
            #         FIN1<-InMEDA[["Ionization mode"]]
            #         F1IN1<-as.character(FIN1)
            #         F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
            #         out<-c(out,F2IN1)
            #         #################Adding this new#################################
            #         ##################################################
            #         ##FINIKOT<-MaKE.ONT.REC(InMEDA)
            #         ##out<-c(out,FINIKOT)
            #         #################################
            #         InKeyVal1<-which(stringi::stri_detect_fixed(FNA,"INCHIKEY:"))
            #         InKeyVal2<-FNA[InKeyVal1]
            #         InKeyVal3<-gsub("INCHIKEY:","",InKeyVal2)
            #         InKeyVal4<-str_trim(gsub("INCHIKEY:","",InKeyVal3))
            #         InKeyVal<-InKeyVal4
            #         ########################
            #         IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
            #         ######################
            #         if(!sjmisc::is_empty(IKCRV)){
            #           ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
            #           F1ONT<-paste("Ontology:",ONTV,sep=" ")
            #           out<-c(out,F1ONT)
            #         }else{
            #           F1ONT<-paste("Ontology:","",sep=" ")
            #           out<-c(out,F1ONT)
            #         }
            #         ######################
            #         FINK<-paste("INCHIKEY:",tryCatch({InKeyVal},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
            #         out<-c(out,FINK)
            #         #####################
            #         INV<-which(stringi::stri_detect_fixed(FNA,"INCHI:"))
            #         INV1<-FNA[INV]
            #         ###################
            #         INV2<-gsub("INCHI:","",INV1)
            #         INV3<-str_trim(INV2)
            #         FINCH<-paste("INCHI:",INV3,sep=" ")
            #         out<-c(out,FINCH)
            #         ##########################
            #         ###########################
            #         FSIM<-paste("SMILES:",gETSmiles(InMEDA),sep=" ")
            #         out<-c(out,FSIM)
            #         #######################
            #         FFOR<-InMEDA[["Formula"]]
            #         FFOR1<-paste("FORMULA:",FFOR,sep=" ")
            #         out<-c(out,FFOR1) 
            #         ###################################################		
            #         ###################################################
            #         #############################
            #         FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
            #         FINS1<-FNA[FINS]
            #         FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:","",sep=""))
            #         out<-c(out,FINS2)
            #         ############################
            #         FAUT<-as.character(InMEDA[["Authors"]])
            #         FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
            #         out<-c(out,FAUT1)
            #         ##########################
            #         ##FLIC<-paste("LICENSE:",sep=" ")
            #         FLIC<-paste("LICENSE:","CC BY",sep=" ")
            #         out<-c(out,FLIC)
            #         ###########################
            #         FCIE<-paste("COLLISIONENERGY:",as.character(InMEDA[["Collision energy"]]),sep=" ")
            #         out<-c(out,FCIE)
            #         #########################
            #         FINST<-as.character(InMEDA[["INSTRUMENT_TYPE"]])
            #         FINST1<-paste("INSTRUMENTTYPE:",FINST,sep=" ")
            #         out<-c(out,FINST1)
            #         ########################
            #         FINS<-as.character(InMEDA[["INSTRUMENT"]])
            #         FINS1<-paste("INSTRUMENT:",FINS,sep=" ")
            #         out<-c(out,FINS1)
            #         ####################
            #         FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
            #         out<-c(out,FCOM)
            #         ##################
            #         FNPA<-which(stringi::stri_detect_fixed(FNA,"Num Peaks:"))
            #         F1NPA<-FNA[FNPA]
            #         out<-c(out,F1NPA)
            #         ###################
            #         Fpea<-FNA[(FNPA+1):Find]
            #         #########################
            #         if(is.na(Fpea))
            #         {
            #           Fpea1<-FNA[(FNPA+1)]
            #           
            #           
            #         }else{
            #           
            #           MV=AAMS1
            #           tes1<-unlist(strsplit(Fpea, "\t"))
            #           tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
            #           tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
            #           tes4<-which(tes2 > (3+MV))
            #           if(length(tes4)>1)
            #           {
            #             tes5<-tes2[-tes4]
            #             tes6<-tes3[-tes4]
            #             tes7<-paste(tes5,tes6,sep="\t")
            #             out<-c(out,tes7)
            #           }else{
            #             out<-c(out,Fpea)
            #             
            #           }
            #         }
            #         ############################
            #       }
            # #############################################
            # ######################### This is add####
            #     }
            #########################################################      
            #########################################################
              }else if(length(INLL) > 1){
                ##############################################
                print("enter the line ...1020")
                ##############################################
                MONMS=InMSPL[INLL]
                TRA<-unname(rapply(MONMS, function(x) grep("RETENTIONTIME:",x, value=TRUE)))
                TRA1<-as.numeric(stringr::str_trim(stringr::str_replace(TRA, "RETENTIONTIME:", "")))
                TRA2<-abs(VRT-TRA1)
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
                # PT <- c("PRECURSORTYPE:", "ADDUCTIONNAME:")
                # NPT<-tryCatch({grep(paste(PT,collapse="|"), FNA, value=TRUE)},error=function(cond){message("List value is empty")})
                # NPT1<-tryCatch({match(NPT,FNA)},error=function(cond){message("List value is empty")})
                # PTV <-tryCatch({stringr::str_remove(FNA[NPT1],c("PRECURSORTYPE:","ADDUCTIONNAME:"))},error=function(cond){message("List value is empty")})
                # PTV1<-tryCatch({PTV[1]},error=function(cond){message("List value is empty")})
                # PTV2<-tryCatch({stringr::str_trim(PTV1)},error=function(cond){message("List value is empty")})
                # P1TV2<-tryCatch({stringr::str_trim(gsub("ADDUCTIONNAME:","",PTV2))},error=function(cond){message("List value is empty")})
                # PTV3<-tryCatch({paste("PRECURSORTYPE:",P1TV2)},error=function(cond){message("List value is empty")})
                # ##PTV3<-tryCatch({paste("PRECURSORTYPE:",P1TV2)},error=function(cond){message("List value is empty")})
                # ########################################
                # print("enter the test area..4")
                # ########################################
                # print(P1TV2)
                # print(as.character(InMEDA[["Adduct"]]))
                ################################## adding this new
               #### if(as.character(InMEDA[["Adduct"]]) == "[M]+" || as.character(InMEDA[["Adduct"]]) == "[M]-"){
                  ###################################################	      
                  ##if(identical(P1TV2,as.character(InMEDA[["Adduct"]]))){
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
                  ####out<-c(out,PTV3)
                  out<-c(out,F1PT1)
                  ###########################
                  FIN1<-InMEDA[["Ionization mode"]]
                  F1IN1<-as.character(FIN1)
                  F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                  out<-c(out,F2IN1)
                  #################################################
                  ########################### adding this new #####
                  ##FINIKOT<-MaKE.ONT.REC(InMEDA)
                  ##out<-c(out,FINIKOT)
                  ######################################
                  ######################################
                  InKeyVal1<-which(stringi::stri_detect_fixed(FNA,"INCHIKEY:"))
                  InKeyVal2<-FNA[InKeyVal1]
                  InKeyVal3<-gsub("INCHIKEY:","",InKeyVal2)
                  InKeyVal4<-str_trim(gsub("INCHIKEY:","",InKeyVal3))
                  InKeyVal<-InKeyVal4
                  ##############################
                  IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
                  ##################################
                  if(!sjmisc::is_empty(IKCRV)){
                    ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
                    F1ONT<-paste("Ontology:",ONTV,sep=" ")
                    out<-c(out,F1ONT)
                  }else{
                    F1ONT<-paste("Ontology:","",sep=" ")
                    out<-c(out,F1ONT)
                  }
                  ###############################
                  INV<-which(stringi::stri_detect_fixed(FNA,"INCHI:"))
                  INV1<-FNA[INV]
                  ######################
                  INV2<-gsub("INCHI:","",INV1)
                  INV3<-str_trim(INV2)
                  FINCH<-paste("INCHI:",INV3,sep=" ")
                  out<-c(out,FINCH)
                  ######################################
                  #############################
                  FSIM<-paste("SMILES:",gETSmiles(InMEDA),sep=" ")
                  out<-c(out,FSIM)
                  ############################################
                  FFOR<-InMEDA[["Formula"]]
                  FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                  out<-c(out,FFOR1)
                  ############################################
                  ############################################
                  FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                  FINS1<-FNA[FINS]
                  FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:","",sep=""))
                  out<-c(out,FINS2)
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
            ##########################################
              ##########################################
                ####}
              ############################################
              ############################################  
          #       else{
          #         #####################################     
          #         print("enter the else part...4")
          #         #####################################
          #         #####################################
          #         if(as.character(InMEDA[["Adduct"]]) == "[M]+" || as.character(InMEDA[["Adduct"]]) == "[M]-"){
          #           #####################################
          #           #####################################
          #           FNA1<-which(stringi::stri_detect_fixed(FNA,"NAME:"))
          #           FNA2<-InMEDA[["Name"]]
          #           FNA3<-as.character(FNA2)
          #           #######################
          #           FNAM<-paste("NAME:",FNA3,sep=" ")
          #           out<-c(out,FNAM)
          #           ########################
          #           FRA1<-which(stringi::stri_detect_fixed(FNA,"RETENTIONTIME:"))
          #           F1RA1<-FNA[FRA1]
          #           out<-c(out,F1RA1)
          #           ################################
          #           FMZ1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORMZ:"))
          #           F1MZ1<-FNA[FMZ1]
          #           out<-c(out,F1MZ1)
          #           ################################
          #           FPT1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORTYPE:"))
          #           F1PT1<-FNA[FPT1]
          #           out<-c(out,PTV3)
          #           #out<-c(out,F1PT1)
          #           #################################
          #           FIN1<-InMEDA[["Ionization mode"]]
          #           F1IN1<-as.character(FIN1)
          #           F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
          #           out<-c(out,F2IN1)
          #           ########################################
          #           ################ Adding this new #######
          #           ##FINIKOT<-MaKE.ONT.REC(InMEDA)
          #           ##out<-c(out,FINIKOT)
          #           ##########################################
          #           ##########################################
          #           InKeyVal1<-which(stringi::stri_detect_fixed(FNA,"INCHIKEY:"))
          #           InKeyVal2<-FNA[InKeyVal1]
          #           InKeyVal3<-gsub("INCHIKEY:","",InKeyVal2)
          #           InKeyVal4<-str_trim(gsub("INCHIKEY:","",InKeyVal3))
          #           InKeyVal<-InKeyVal4
          #           ###########################################
          #           IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
          #           #########################
          #           if(!sjmisc::is_empty(IKCRV)){
          #             ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
          #             F1ONT<-paste("Ontology:",ONTV,sep=" ")
          #             out<-c(out,F1ONT)
          #           }else{
          #             F1ONT<-paste("Ontology:","",sep=" ")
          #             out<-c(out,F1ONT)
          #           }
          #           ####################
          #           FINK<-paste("INCHIKEY:",tryCatch({InKeyVal},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
          #           out<-c(out,FINK)
          #           #####################
          #           INV<-which(stringi::stri_detect_fixed(FNA,"INCHI:"))
          #           INV1<-FNA[INV]
          #           #######################
          #           INV2<-gsub("INCHI:","",INV1)
          #           INV3<-str_trim(INV2)
          #           FINCH<-paste("INCHI:",INV3,sep=" ")
          #           out<-c(out,FINCH)
          #           ########################
          #           #####################
          #           FSIM<-paste("SMILES:",gETSmiles(InMEDA),sep=" ")
          #           out<-c(out,FSIM)
          #           ###############
          #           FFOR<-InMEDA[["Formula"]]
          #           FFOR1<-paste("FORMULA:",FFOR,sep=" ")
          #           out<-c(out,FFOR1)
          #           ###################################################		
          #           ###################################################
          #           ##IKCRV<-tryCatch({res <- R.utils::withTimeout({classyfireR::get_classification(InKeyVal)}, timeout=1.08, onTimeout="warning")}, warning=function(ex) {message("Classifier not able to fetch information")})
          #           ##out<-c(out,FFOR1)
          #           ##################################
          #           ##################################
          #           FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
          #           FINS1<-FNA[FINS]
          #           FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:","",sep=""))
          #           out<-c(out,FINS2)
          #           ############################
          #           FAUT<-as.character(InMEDA[["Authors"]])
          #           FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
          #           out<-c(out,FAUT1)
          #           ##########################
          #           ##FLIC<-paste("LICENSE:",sep=" ")
          #           FLIC<-paste("LICENSE:","CC BY",sep=" ")
          #           out<-c(out,FLIC)
          #           ###########################
          #           FCIE<-paste("COLLISIONENERGY:",as.character(InMEDA[["Collision energy"]]),sep=" ")
          #           out<-c(out,FCIE)
          #           #########################
          #           FINST<-as.character(InMEDA[["INSTRUMENT_TYPE"]])
          #           FINST1<-paste("INSTRUMENTTYPE:",FINST,sep=" ")
          #           out<-c(out,FINST1)
          #           ########################
          #           FINS<-as.character(InMEDA[["INSTRUMENT"]])
          #           FINS1<-paste("INSTRUMENT:",FINS,sep=" ")
          #           out<-c(out,FINS1)
          #           ####################
          #           FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
          #           out<-c(out,FCOM)
          #           ##################
          #           FNPA<-which(stringi::stri_detect_fixed(FNA,"Num Peaks:"))
          #           F1NPA<-FNA[FNPA]
          #           out<-c(out,F1NPA)
          #           ###################
          #           Fpea<-FNA[(FNPA+1):Find]
          #           #########################
          #           if(is.na(Fpea))
          #           {
          #             Fpea1<-FNA[(FNPA+1)]
          #             
          #             
          #           }else{
          #             
          #             MV=AAMS1
          #             tes1<-unlist(strsplit(Fpea, "\t"))
          #             tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
          #             tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
          #             tes4<-which(tes2 > (3+MV))
          #             if(length(tes4)>1)
          #             {
          #               tes5<-tes2[-tes4]
          #               tes6<-tes3[-tes4]
          #               tes7<-paste(tes5,tes6,sep="\t")
          #               out<-c(out,tes7)
          #             }else{
          #               out<-c(out,Fpea)
          #               
          #             }
          #           }
          #           ############################
          #         }
          # ########################## THis is e##
          # #######################################
          #       }
            #################################      
            ##################################
              } else {
                #print("entering the line 764")
                PASS<-RRV
              }  ## main else part
              ###testing if this works
              return(out)
              #####################
            } #ITmass
          } #ITRTL
          ###################
        } ## adduct match and exacct mass present
      }## FM2 is not empty 
      ###############################################################
      ###############################################################
    }###this is end of else ...where it is checking from NAME
  } ## end of else

}###NFFilter

}

###########################################################################################
###########################################################################################
##print("enter the are before main loop")
##print(length(LmeCmu1))
##print(LmeCmu1)
##print(length(LmeCmu1))

##print(dim(RXF3))
###########################################################################################
###########################################################################################
for(i in 1:length(LmeCmu1))
##for(i in 1:1)
{
  ########################
  print("entering the main function")
  ########################
  
  Val=LmeCmu1[i]
  Val1=LmeCmu1[i+1]
  nVal=Val+1
  #if(!is.na(Val) & !is.na(Val1))
  ########################
  print(Val)
  print(Val1)
  print(nVal)
  ##############################
  if(!sjmisc::is_empty(Val) & !sjmisc::is_empty(Val1))
  {
    ############################################
    print("entering the line ...1707")
    #############################################
    NRXF3<-RXF3[nVal:Val1,]
    ##print(NRXF3)
    ##############################################
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
    #######################################################################
    NFiNA<-tryCatch({list.files(path =Fi6 , recursive = TRUE, full.names = TRUE)},warning=function(cond){message("some mistake happened in file search files")})
    NFiNA1<-tryCatch({match(FINAMSP,basename(NFiNA))},warning=function(cond){message("some mistake happened in file search files")})
    NFiNA2<-tryCatch({NFiNA[NFiNA1]},warning=function(cond){message("some mistake happened in file search files")})
    #####OUNA2 ... Output file ##########################################
    FiNA1<-NFiNA2
    #####################################################################
    ######################################################################
    ##print(FiNA1)
    #####################################################################
    if((length(FiNA1) >= 1))
    {
      #####################################	    
      print("entering the line ...line 1723")
      #####################################
      CaIF<-tryCatch({MaKlist(FiNA1[1])},error=function(cond){message("some mistake happened in filename..file name must be empty")})
      lst2<-tryCatch({CaIF[[1]]},error=function(cond){message("some mistake happened in indexes")})
      fmass<-tryCatch({CaIF[[2]]},error=function(cond){message("some mistake happened in indexes")})
      FRTL1<-tryCatch({CaIF[[3]]},error=function(cond){message("some mistake happened in indexes")})
      ########################################
      if(file.exists(FiNA1[1]))
      {
	#####################################
        print("entering the line ...1732")
        ######################################
        if(length(NRXF3)>0){
          for (i in 1:length(NRXF3)){
            ###################################
	  
            RRV<-NRXF3[i,]
	    ####################################
	    ##print("my test.... adduct")
	    ##print(RRV)
            ####################################
	    ##PIN<-ifelse(sjmisc::is_empty(RRV[["InChI"]]), 1, 0)
	    ##PSM<-ifelse(sjmisc::is_empty(RRV[["SMILES"]]), 1, 0)
	    ##PPC<-ifelse(sjmisc::is_empty(RRV[["PubChem CID"]]), 1, 0)
	    ############################################
            PIN<-ifelse(is.na(RRV[["InChI"]]), 1, 0)
            PSM<-ifelse(is.na(RRV[["SMILES"]]), 1, 0)
            PPC<-ifelse(is.na(RRV[["PubChem CID"]]), 1, 0)
	    ############################################
	    ##print(PIN)
	    ##print(PSM)
	    ##print(PPC)
            ###################################
            if(PIN == 1 & PSM == 1 & PPC == 1){
              ##################################
	      print("enter the if loop as Inchi and smiles and pubchemid all are empty")
	      ##print(NFFilter(RRV,AIN,lst2,fmass,FRTL1))
	      ##print(RRV)
	      ##print(length(LmeCmu1))
	      ##print(NFFilter1(RRV,AIN,lst2,fmass,FRTL1))
	      ##print(length(lst2))
	      ##print(length(fmass))
	      ##print(length(FRTL1))
	      ##print(NFFilter(RRV,AIN,lst2,fmass,FRTL1))
	      #################### adding this new
	      if(sjmisc::is_empty(RRV)){
		te<-"PASS"
	      }else{
		      #####################################
		      #####################################
		      RRV1<-NRXF3[i,]
		      #########################################
		      #########################################
		      FINKE<-NFFilter(RRV1,AIN,lst2,fmass,FRTL1)
		      if(!sjmisc::is_empty(FINKE) & length(FINKE) > 1)
		      {
			      print("enter the line ...1680")
			      cat(sapply(FINKE, toString), file=OUNA3, sep="\n",append=TRUE)
			      cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
		      }else{
			      print("enter the line ...1691")
			      print("enter the else part")
			      #################################################
			      DN<-dirname(dirname(File1))
			      FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
			      write.table(unname(as.data.frame(RRV1)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
			      ###################################################

		      }###end of if ..else
		    #############################################################################
		    #############################################################################
	      }
	    #####################################################################################
	    #####################################################################################
            }else{
	     #######################################################
	      print("entering the else part")
	      RRV2<-NRXF3[i,]
	      #####################
	      ##print(FiNA1)
	      ##print(i)
	      ##print(RRV2[["Name"]])
	      ######################################################
              ##print(RRV2)
              ##print(RRV2[["Adduct"]])	      
	      ######################################################
              print("enter the line ...1749")
	      #####################################################
              if(!sjmisc::is_empty(stringr::str_trim(as.character(RRV2[["InChI"]]))) & !startsWith(as.character(RRV2[["InChI"]]),'not available') & !startsWith(as.character(RRV2[["InChI"]]),'CAS:') & !startsWith(as.character(RRV2[["InChI"]]),'InChI=')){
		##################################################
                print("enter the line ...1871")
	        print(stringr::str_trim(as.character(RRV2[["InChI"]])))
		print(webchem::is.inchikey(stringr::str_trim(as.character(RRV2[["InChI"]]))))
		##################################################
                if(tryCatch({webchem::is.inchikey(stringr::str_trim(as.character(RRV2[["InChI"]])))},error=function(cond){message("not able to pass is.inchiKey")})){
		  ###############################################
                  IV<-stringr::str_trim(as.character(RRV2[["InChI"]]))
		  ###########################################
		   ##print("test what is adduct is reading")
		   print("entering the inchikey pass area")
		   ##print(RRV2[["Adduct"]])
	           ###############################################
		         print("enter the line ...1873")
              ###############################################
              ################################################    
		              if(!sjmisc::is_empty(IV)){  
                    outn<-Ikfilter(IV,lst2,RRV2,AIN,fmass,FRTL1)
                    len1<-length(outn)
                    ##############################################
                    if(len1 > 1)
                    {
                      print("enter the line ...1878")
                      cat(sapply(outn, toString), file=OUNA3, sep="\n",append=TRUE)
                      cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
                      
                    }else{
                      print("enter the line ...1879")
                      FINKE<-NFFilter(RRV2,AIN,lst2,fmass,FRTL1)
                      ### adding this code new 
                      len3<-length(FINKE)
                      if(len3 > 1){
                        print("enter the line ...1880")
                        cat(sapply(FINKE, toString), file=OUNA3, sep="\n",append=TRUE)
                        cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
                      }else{
			      print("enter the line ...1891")
			      print("enter the else part...Error report")
			      ################################
			      DN<-dirname(dirname(File1))
			      DN1<-paste(DN,"Error-Report",sep="/")
			      FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
			      write.table(unname(as.data.frame(RRV2)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
			      #################################
		      }
        ################################################### 
                    } # end of the else loop
          #############################################                                
                }else{
	        ###################################################################	
                  print("enter the line ...1780")
                  FINKE<-NFFilter(RRV2,AIN,lst2,fmass,FRTL1)
                  ### adding this code new 
                  len3<-length(FINKE)
		  #################################################################
                  if(len3 > 1){
                    print("enter the line ...1783")
                    cat(sapply(FINKE, toString), file=OUNA3, sep="\n",append=TRUE)
                    cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
                  }else{
			  print("enter the line ...1785")
			  print("enter the else part")
			  ############################################################
			  DN<-dirname(dirname(File1))
			  DN1<-paste(DN,"Error-Report",sep="/")
			  FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
			  write.table(unname(as.data.frame(RRV2)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
			  ###########################################################
		  }
		############################################################## 
                } ## end of else
    ################################################################              
    ################################################################              
		}else{
			print("enter the else loop ...it did not pass is.inchikey")
			############################################################
			FINKE<-NFFilter(RRV2,AIN,lst2,fmass,FRTL1)
                  	### adding this code new
                  	len3<-length(FINKE)
                  	#################################################################
                  	if(len3 > 1){
                    		print("enter the line ...1783")
                    		cat(sapply(FINKE, toString), file=OUNA3, sep="\n",append=TRUE)
                    		cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
                  	}else{
                          	print("enter the line ...1785")
                          	print("enter the else part")
                          	############################################################
                          	DN<-dirname(dirname(File1))
                          	DN1<-paste(DN,"Error-Report",sep="/")
                          	FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
                          	write.table(unname(as.data.frame(RRV2)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
                          ###########################################################
                  	    }

	################################################################
		}### end of else ...where it did not pass Inchikey
	#####################################################################################################	
              ######################################################################################################################################
              }else if(!sjmisc::is_empty(as.character(RRV2[["InChI"]])) & startsWith(as.character(RRV2[["InChI"]]),'InChI=')){
		#################################################      
                print("enter the line ...861")
	        print("entering the inchi pass area")
                ##IV1<-stringr::str_trim(as.character(RRV[["InChI"]]))
	        #####################################################
	        IV1<-stringr::str_trim(as.character(RRV2[["InChI"]]))
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
                  outn1<-Ikfilter(FSMV2,lst2,RRV2,AIN,fmass,FRTL1)
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
                  FINKE1<-NFFilter(RRV2,AIN,lst2,fmass,FRTL1)
                  ### adding this code new 
                  len3<-length(FINKE1)
                  if(len3 > 1){
                    print("enter the line ...888")
                    cat(sapply(FINKE1, toString), file=OUNA3, sep="\n",append=TRUE)
                    cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
                  }else{
			  print("enter the line ...889")
			  print("enter the else part")
			  ###################################
			  DN<-dirname(dirname(File1))
			  DN1<-paste(DN,"Error-Report",sep="/")
			  FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
			  write.table(unname(as.data.frame(RRV2)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
			  ########################################

		  }
		##################################################
                } # end of the else loop
	       
	      }else{
		  ################################################    
                  print("enter the line ...870...inchi to inchikey conversion failed")
                  FINKE1<-NFFilter(RRV2,AIN,lst2,fmass,FRTL1)
                  ### adding this code new 
                  len3<-length(FINKE1)
                  if(len3 > 1){
                    print("enter the line ...873")
                    cat(sapply(FINKE1, toString), file=OUNA3, sep="\n",append=TRUE)
                    cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
                  }else{
			  print("enter the line ...875")
			  print("enter the else part")
			  #####################################
			  DN<-dirname(dirname(File1))
			  DN1<-paste(DN,"Error-Report",sep="/")
			  FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
			  write.table(unname(as.data.frame(RRV2)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
			  #######################################
		  }
		  ################################################

                }

              ###############################################################################################################################
              }else if(!sjmisc::is_empty(stringr::str_trim(as.character(RRV2[["InChI"]]))) & startsWith(as.character(RRV2[["InChI"]]),'CAS:')){
		####################################
                print("enter the line ...884")
	        print("enter the cas pass area")
	        ###################################
                CV<-stringr::str_trim(as.character(RRV2[["InChI"]]))
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
                    outn3<-Ikfilter(FIV2,lst2,RRV2,AIN,fmass,FRTL1)
                    len4<-length(outn3)
                    if(len4 > 1)
                    {
                      print("enter the line ...912")
                      cat(sapply(outn3, toString), file=OUNA3, sep="\n",append=TRUE)
                      cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
                      
                    ##} # end of len4
                  }else{
                    FSMIL<-NFFilter(RRV2,AIN,lst2,fmass,FRTL1)
                    len3<-length(FSMIL)
                    if(len3 > 1){
                      print("enter the line ...921")
                      cat(sapply(FSMIL, toString), file=OUNA3, sep="\n",append=TRUE)
                      cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
                    }else{
			    print("enter the line ...923")
			    print("enter the else part")
			    #####################################
			    DN<-dirname(dirname(File1))
			    DN1<-paste(DN,"Error-Report",sep="/")
			    FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
			    write.table(unname(as.data.frame(RRV2)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
			    ########################################
		    }### end of inner else
		   ##########################
                  } # end of else loop
		  
		}else{
		              ###################################
		              print("CAS to inchikey conversion failed")
		  IV<-ifelse(!sjmisc::is_empty(tryCatch({PuCAStoOI(CV2)[2]},error=function(cond){message("CAS is empty")})),tryCatch({PuCAStoOI(CV2)[2]},error=function(cond){message("CAS value is empty")}),"NA")
                  ######################################################
		  ###################################################################
		  if(!sjmisc::is_empty(IV)){  
		    outn<-Ikfilter(IV,lst2,RRV2,AIN,fmass,FRTL1)
		    len1<-length(outn)
		    ##############################################
		    if(len1 > 1)
		    {
		      print("enter the line ...1878")
		      cat(sapply(outn, toString), file=OUNA3, sep="\n",append=TRUE)
		      cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
		      
		    }else{
		      print("enter the line ...1879")
		      FINKE<-NFFilter(RRV2,AIN,lst2,fmass,FRTL1)
		      ### adding this code new 
		      len3<-length(FINKE)
		      if(len3 > 1){
		        print("enter the line ...1880")
		        cat(sapply(FINKE, toString), file=OUNA3, sep="\n",append=TRUE)
		        cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
		      }else{
		        print("enter the line ...1891")
		        print("enter the else part...Error report")
		        ################################
		        DN<-dirname(dirname(File1))
		        DN1<-paste(DN,"Error-Report",sep="/")
		        FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
		        write.table(unname(as.data.frame(RRV2)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
		        #################################
		      }
		      ################################################### 
		    } # end of the else loop
		    #############################################                                
		  }else{
		    ###################################################################	
		    print("enter the line ...1780")
		    FINKE<-NFFilter(RRV2,AIN,lst2,fmass,FRTL1)
		    ### adding this code new 
		    len3<-length(FINKE)
		    #################################################################
		    if(len3 > 1){
		      print("enter the line ...1783")
		      cat(sapply(FINKE, toString), file=OUNA3, sep="\n",append=TRUE)
		      cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
		    }else{
		      print("enter the line ...1785")
		      print("enter the else part")
		      ############################################################
		      DN<-dirname(dirname(File1))
		      DN1<-paste(DN,"Error-Report",sep="/")
		      FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
		      write.table(unname(as.data.frame(RRV2)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
		      ###########################################################
		    }
		    ############################################################## 
		  } ## end of else
		  #######################################################
		  ######################################################
#                   FCAS<-NFFilter(RRV2,AIN,lst2,fmass,FRTL1)
#                   ### adding this code new 
#                   len4<-length(FCAS)
#                   if(len4 > 1){
#                     print("enter the line ...935")
#                     cat(sapply(FCAS, toString), file=OUNA3, sep="\n",append=TRUE)
#                     cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
#                   }else{
# 			  print("enter the line ...937")
# 			  print("enter the else part")
# 			  #######################################
# 			  DN<-dirname(dirname(File1))
# 			  DN1<-paste(DN,"Error-Report",sep="/")
# 			  FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
# 			  write.table(unname(as.data.frame(RRV2)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
# 			  ##########################################
# 
# 		  }####end of else
		 #############################################
                }### end of main else 
              #########################################################################################################################
              }else if(!sjmisc::is_empty(as.character(RRV2[["SMILES"]])) & !startsWith(as.character(RRV2[["SMILES"]]),'not available')){
                print("enter the line ...915")
	        print("enter the smiles pass area")
	        #####################################################
                F1SM<-stringr::str_trim(as.character(RRV2[["SMILES"]]))
                F1SM1<-tryCatch({rinchi::get.inchi.key(F1SM)},error=function(cond){message("webchecm could not fetch the info")})
                ######################################################
		            if(!sjmisc::is_empty(F1SM1)){
                  outn4<-Ikfilter(F1SM1,lst2,RRV2,AIN,fmass,FRTL1)
                  loutn5<-length(outn4)
                  ###################################################
                  if(loutn5 > 1)
                  {
                    print("enter the line ...954")
                    cat(sapply(outn4, toString), file=OUNA3, sep="\n",append=TRUE)
                    cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
                    
                  }else{
                    print("enter the line ...959")
                    FSMIL<-NFFilter(RRV2,AIN,lst2,fmass,FRTL1)
                    len3<-length(FSMIL)
                    if(len3 > 1){
                      print("enter the line ...963")
                      cat(sapply(FSMIL, toString), file=OUNA3, sep="\n",append=TRUE)
                      cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
                    }else{
			    print("enter the line ...967")
			    print("enter the else part")
			    ##################################
			    DN<-dirname(dirname(File1))
			    DN1<-paste(DN,"Error-Report",sep="/")
			    FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
			    write.table(unname(as.data.frame(RRV2)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
			    ##########################################
		    }
	       ################################################################## 
                  } ## end of else loop
	         
	         }else{
	         ###############################################		 
                  #print("enter the line ...970")
	                 print("smiles to inchikey conversion is empty so checking other code")
                  FSMIL<-NFFilter(RRV2,AIN,lst2,fmass,FRTL1)
                  len3<-length(FSMIL)
                  if(len3 > 1){
                    print("enter the line ...971")
                    cat(sapply(FSMIL, toString), file=OUNA3, sep="\n",append=TRUE)
                    cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
                  }else{
			  print("enter the line ...974")
			  print("enter the else part")
			  ######################################
			  DN<-dirname(dirname(File1))
			  DN1<-paste(DN,"Error-Report",sep="/")
			  FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
			  write.table(unname(as.data.frame(RRV2)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
			  #########################################
		  }
	      ###################################################### 
                }### end of else ..smiles to inchikey conversion failed
                
              #############################################################################################################################
              }else if(!sjmisc::is_empty(stringr::str_trim(as.character(RRV2[["PubChem CID"]]))) & !startsWith(as.character(RRV2[["PubChem CID"]]),'not available')){
                print("enter the line ...939")
	        print("entering the pubchem pass area")
	        #############################################
                FPUCID<-stringr::str_trim(as.character(RRV2[["PubChem CID"]]))
		FPUCID1<-as.numeric(FPUCID)
		##############################################
		print("entering the PubchemID test area")
		print(FPUCID1)
		######################################################################
		######################################################################
		FINSM<-tryCatch({webchem::pc_prop(FPUCID1)},error=function(cond){message("webchecm could not fetch the info")})
		FIINK<-tryCatch({FINSM$InChIKey},error=function(cond){message("webchecm could not fetch the info")})
                #####################################################################
		#####################################################################
		              if(!sjmisc::is_empty(FIINK)){
                    outn6<-Ikfilter(FIINK,lst2,RRV2,AIN,fmass,FRTL1)
                    len7<-length(outn6)
                    if(len7 > 1)
                    {
                      print("enter the line ...994")
                      cat(sapply(outn6, toString), file=OUNA3, sep="\n",append=TRUE)
                      cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
                      
                    }else{
                      print("enter the line ...999")
                      FCIDR<-NFFilter(RRV2,AIN,lst2,fmass,FRTL1)
                      len3<-length(FCIDR)
                      if(len3 > 1){
                        print("enter the line ...1003")
                        cat(sapply(FCIDR, toString), file=OUNA3, sep="\n",append=TRUE)
                        cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
                      }else{
			      print("enter the line ...1005")
			      print("enter the else part")
			      ######################################
			      DN<-dirname(dirname(File1))
			      DN1<-paste(DN,"Error-Report",sep="/")
			      FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
			      write.table(unname(as.data.frame(RRV2)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
			      ################################################
		      }
		  #######################################################################
                    } # end of else loop
		   }else{
		    #####################################################################	   
                    #print("enter the line ...1011")
		                print("entering the else part...FIINK...pubchem CID to inchikey failed")
		          ############################################
		     #############entering new code here ##############
		     IV<-ifelse(!sjmisc::is_empty(tryCatch({ConvPCIDtoOCN(FPUCID1)[1]},error=function(cond){message("Pubchem is empty")})),tryCatch({ConvPCIDtoOCN(FPUCID1)[1]},error=function(cond){message("CAS value is empty")}),"NA")
		     if(!sjmisc::is_empty(IV)){  
		       outn<-Ikfilter(IV,lst2,RRV2,AIN,fmass,FRTL1)
		       len1<-length(outn)
		       ##############################################
		       if(len1 > 1)
		       {
		         print("enter the line ...1878")
		         cat(sapply(outn, toString), file=OUNA3, sep="\n",append=TRUE)
		         cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
		         
		       }else{
		         print("enter the line ...1879")
		         FINKE<-NFFilter(RRV2,AIN,lst2,fmass,FRTL1)
		         ### adding this code new 
		         len3<-length(FINKE)
		         if(len3 > 1){
		           print("enter the line ...1880")
		           cat(sapply(FINKE, toString), file=OUNA3, sep="\n",append=TRUE)
		           cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
		         }else{
		           print("enter the line ...1891")
		           print("enter the else part...Error report")
		           ################################
		           DN<-dirname(dirname(File1))
		           DN1<-paste(DN,"Error-Report",sep="/")
		           FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
		           write.table(unname(as.data.frame(RRV2)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
		           #################################
		         }
		         ################################################### 
		       } # end of the else loop
		       #############################################                                
		     }else{
		       ###################################################################	
		       print("enter the line ...1780")
		       FINKE<-NFFilter(RRV2,AIN,lst2,fmass,FRTL1)
		       ### adding this code new 
		       len3<-length(FINKE)
		       #################################################################
		       if(len3 > 1){
		         print("enter the line ...1783")
		         cat(sapply(FINKE, toString), file=OUNA3, sep="\n",append=TRUE)
		         cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
		       }else{
		         print("enter the line ...1785")
		         print("enter the else part")
		         ############################################################
		         DN<-dirname(dirname(File1))
		         DN1<-paste(DN,"Error-Report",sep="/")
		         FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
		         write.table(unname(as.data.frame(RRV2)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
		         ###########################################################
		       }
		       ############################################################## 
		     } ## end of else
		     ###################################################
		     ##################################################
#                     FSMIL<-NFFilter(RRV2,AIN,lst2,fmass,FRTL1)
#                     len3<-length(FSMIL)
#                     if(len3 > 1){
#                       print("enter the line ...1014")
#                       cat(sapply(FSMIL, toString), file=OUNA3, sep="\n",append=TRUE)
#                       cat(sapply("", toString), file=OUNA3, sep="\n",append=TRUE)
#                     }else{
# 			    print("enter the line ...1016")
# 			    print("enter the else part")
# 			    ###################################
# 			    DN<-dirname(dirname(File1))
# 			    DN1<-paste(DN,"Error-Report",sep="/")
# 			    FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
# 			    write.table(unname(as.data.frame(RRV2)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
# 			    ############################################
# 		    }
		############################################################################## 
                  }# end of else loop
                ##############################################################################
                ################################################################################
	      }else{
		print("entering the final else part")
                PASS1 <-RRV
              }
           #####################################################################################   
	   }
        } ### big else loop
      } ## for loop NRXF3
    } ## End ...NRXF3

  } ## FiNA is closing one
 ###################################
} ## is_empty(Val and Val1)
}
