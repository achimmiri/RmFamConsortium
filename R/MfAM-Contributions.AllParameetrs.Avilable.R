##install.packages("tidyverse")
##install.packages("readxl")
##install.packages("openxlsx")
##library("openxlsx")
##install.packages("devtools")
##devtools::install_github("ropensci/webchem")
##install.packages("webchem",dependencies = TRUE)
##install.packages("remotes")
##remotes::install_github("phively/wranglR")
##install.packages("janitor")
##install.packages("matlib")
##devtools::install_github("r-cas/ryacas")
##devtools::install_github("tidyverse/glue")
##remotes::install_github("jasenfinch/mzAnnotation")
##remotes::install_github("xia-lab/MetaboAnalystR")
##BiocManager::install("RMassBank")
###BiocManager::install("metaMS")
##BiocManager::install("rcellminerData")
##BiocManager::install("rcellminer")
##BiocManager::install("zlibbioc")
##BiocManager::install("ChemmineOB")
##BiocManager::install("Rdisop")
##install.packages("GMCM")
###################################################################
###################################################################
library(devtools)
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
library(RCurl)
library(Rdisop)
library(RChemMass)
library(mzAnnotation)
library("readxl")
library("metaMS")
library("rcellminer")
library("R.utils")
#########################################################
#########################################################
#######library(GMCM)
#########################################################
#########################################################


########## Step1: Reading the API key ####################################
##########################################################################
Sys.setenv(CHEMSPIDER_KEY = apikey)
rr_auth(apikey)
###########################################################################
###########################################################################
args <- commandArgs(TRUE)
File1<-args[1]
mz_Tol=as.numeric(args[2])
RT_Tol=(as.numeric(args[3])/100)
DEFAULT_MZ_TOLERANCE=as.numeric(args[4])
##########################################################################
###########################################################################
Fi <- unlist(strsplit(File1, "/"))
Fi1 <-c(paste(Fi[-length(Fi)], collapse = "/"), data.table::last(Fi))
Fi2<-Fi1[1]
Fi3<-unlist(strsplit(Fi2, "/"))
Fi4<-c(paste(Fi3[-length(Fi3)], collapse = "/"), data.table::last(Fi3))
Fi5<-Fi4[1]
ifelse(!dir.exists(file.path(Fi5,"exported as raw msp")), dir.create(file.path(Fi5,"exported as raw msp")), FALSE)
Fi6<-paste(Fi5,"raw data","exported as raw msp",sep="/")
### Input################
Fi7<-paste(Fi6,"/",sep="")
ifelse(!dir.exists(file.path(Fi5,"converted to msp")), dir.create(file.path(Fi5,"exported as raw msp")), FALSE)
Fi8<-paste(Fi5,"converted to msp",sep="/")
### Output##################
Fi9<-paste(Fi8,"/",sep="")
############################
Fi10<-paste("mz",".",mz_Tol,"ppm",".",RT_Tol,"RT", sep="")
Fi11<-paste0(Fi9,Fi10,"/")
############################
NFi5<-unlist(strsplit(Fi5, "/"))
NFi6<-c(paste(NFi5[-length(NFi5)], collapse = "/"), data.table::last(NFi5))
N1Fi6<-NFi6[1]
################################################################
################################################################
####### Step2: Reading the All the required files
################################################################
##AIN1<-read.table("/mnt/ifs/data/IPB/Projects/2017_005_MS-databases/mFam contributions/scripts/ADI.4.csv",sep=",",header=F,quote="",stringsAsFactors = FALSE)
##AIN2<-read.table("/mnt/ifs/data/IPB/Projects/2017_005_MS-databases/mFam contributions/scripts/ADI.4.csv",sep=",",header=F,stringsAsFactors = FALSE)
#####################################################################
AIN1<-read.csv2("/mnt/ifs/data/IPB/Projects/2017_005_MS-databases/mFam contributions/scripts/ADI.4.csv",sep=",",header=F,quote="",stringsAsFactors = FALSE)
AIN<-AIN1
####################################################################
AIN$V1<-gsub(pattern = "\"", replacement = "", x = AIN$V1)
AIN$V2<-gsub(pattern = "\"", replacement = "", x = AIN$V2)
AIN$V3<-gsub(pattern = "\"", replacement = "", x = AIN$V3)
AIN$V4<-gsub(pattern = "\"", replacement = "", x = AIN$V4)
AIN$V5<-gsub(pattern = "\"", replacement = "", x = AIN$V5)
AIN$V6<-gsub(pattern = "\"", replacement = "", x = AIN$V6)
AIN$V7<-gsub(pattern = "\"", replacement = "", x = AIN$V7)
AIN$V8<-gsub(pattern = "\"", replacement = "", x = AIN$V8)
####################### adding the new database table####################
AD<-read.table("/mnt/ifs/data/IPB/Projects/2017_005_MS-databases/mFam contributions/scripts/Database_Dec2017.txt",sep="\t",header=TRUE,fill = TRUE,stringsAsFactors = FALSE)
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
##########################################################

#########################################################
######### Adding this new################################
Lmeda<-tryCatch({base::rle(as.character(RXF3[["File"]]))},error=function(cond){message("file structure is wrong")})
NLmeda<-tryCatch({Lmeda$lengths},error=function(cond){message("repeat length calculate mistakes")})
###THis is the original one
LmeCmu<-tryCatch({abs(cumsum(NLmeda))},error=function(cond){message("cumulative calculate mistakes happening")})
#######LmeCmu<-tryCatch({abs(cummean(NLmeda))},error=function(cond){message("cumulative calculate mistakes happening")})
#######LmeCmu<-tryCatch({abs(dplyr::cummean(NLmeda))},error=function(cond){message("cumulative calculate mistakes happening")})
#######LmeCmu<-tryCatch({abs(GMCM:::cummean(NLmeda))},error=function(cond){message("cumulative calculate mistakes happening")})
#######print("adding this temporary")
#######print(LmeCmu)
LmeCmu1<-tryCatch({c(0,LmeCmu)},error=function(cond){message("cumulative calculate mistakes happening")})
SFileNam<-tryCatch({Lmeda$values},error=function(cond){message("rle is not able to fetch the information properly")})
#########################################################
#########################################################
##### REST API Functions
##########################################################
##########################################################
PuInKtoSM<-function(getINK)
{
  ###### This Functions return canonical smiles
  url<- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"
  out<-tryCatch({jsonlite::fromJSON(paste0(url,getINK, "/JSON"))} ,error = function(x) {return(NA)})
  prop.names  <-tryCatch({out$PC_Compounds$props[[1]][[1]]},error = function(x) {return(NA)})
  prop.values <- tryCatch({out$PC_Compounds$props[[1]][[2]]},error = function(x) {return(NA)})
  sm <-tryCatch({grep("smiles", prop.names[,"label"], ignore.case = TRUE)},error = function(x) {return(NA)})
  csmiles<-c()
  if(length(sm) >= 1 & !sjmisc::is_empty(sm)) {
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
  ###### This Functions return Isomeric smiles	
  url<- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"
  out<-tryCatch({jsonlite::fromJSON(paste0(url,getINK, "/JSON"))}, error = function(x) {return(NA)})
  prop.names  <-tryCatch({out$PC_Compounds$props[[1]][[1]]}, error = function(x) {return(NA)})
  prop.values <- tryCatch({out$PC_Compounds$props[[1]][[2]]}, error = function(x) {return(NA)})
  sm <-tryCatch({grep("smiles", prop.names[,"label"], ignore.case = TRUE)}, error = function(x) {return(NA)})
  csmiles<-c()
  if(length(sm) == 2 & !sjmisc::is_empty(sm)) {
    can <- tryCatch({grep("canonical", prop.names[,"name"], ignore.case = TRUE)}, error = function(x) {return(NA)})
    can1<-tryCatch({prop.values[sm[2],"sval"]},error= function(x) {return(NA)})
    csmiles<-c(csmiles,can1)
  }else{
    csmiles<-c(csmiles,NA)
  }
  return(csmiles)
  
}
#############################################################################
#############################################################################
ConvSMItoOID2<-function(getSMI)
{
url<- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/SMILES/"
out<-tryCatch({jsonlite::fromJSON(paste0(url,getSMI, "/JSON"))} ,error = function(x) {return(NA)})
prop.names  <-tryCatch({out$PC_Compounds$props[[1]][[1]]},error = function(x) {return(NA)})
prop.values <- tryCatch({out$PC_Compounds$props[[1]][[2]]},error = function(x) {return(NA)})
IUN<-tryCatch({prop.values[10,"sval"]},error = function(x) {return(NA)})
MF<-tryCatch({prop.values[17,"sval"]},error = function(x) {return(NA)})
MI<-tryCatch({prop.values[22,"sval"]},error = function(x) {return(0)})
MW<-tryCatch({prop.values[18 ,"sval"]},error = function(x) {return(0)})
IK<-tryCatch({prop.values[14 ,"sval"]},error = function(x) {return(NA)})
IN<-tryCatch({prop.values[13 ,"sval"]},error = function(x) {return(NA)})
return(c(IUN,MF,MI,MW,IK,IN))
}
##############################################################################
##############################################################################
ConvPCIDtoOCN<-function(getPCID)
{

  url<- "https://www.metabolomicsworkbench.org/rest/compound/pubchem_cid/"
  out<-tryCatch({jsonlite::fromJSON(paste0(url,getPCID, "/all"))}, error = function(x) {return(NA)})
  ##############################
  ##############################
  OIK<-tryCatch({out$inchi_key},error=function(cond){return(NA)})
  OSM<-tryCatch({out$smiles},error=function(cond){return(NA)})
  OCID<-tryCatch({out$pubchem_cid},error=function(cond){return(NA)})
  OEM<-tryCatch({out$exactmass},error=function(cond){return(0)})
  OFOR<-tryCatch({out$formula},error=function(cond){return(NA)})
  return(c(OIK[1],OSM[1],OCID[1],OEM[1],OFOR[1]))
  ###############################
  ###############################

}
##############################################################################
##############################################################################
ConvCIDtoOID1<-function(getCID)
{
  #################################
  url<-"http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
  out<-tryCatch({jsonlite::fromJSON(paste0(url,getCID,"/property/CanonicalSMILES,MonoisotopicMass,InChI,InChIKey,MolecularFormula"))},error = function(x) {return(NA)})
  ############################
  OIK<-tryCatch({out$PropertyTable$Properties$InChIKey},error=function(cond){return(NA)})
  OIN<-tryCatch({out$PropertyTable$Properties$InChI},error=function(cond){return(NA)})
  OSM<-tryCatch({out$PropertyTable$Properties$CanonicalSMILES},error=function(cond){return(NA)})
  OCID<-tryCatch({out$PropertyTable$Properties$CID},error=function(cond){return(NA)})
  EXM<-tryCatch({out$PropertyTable$Properties$MonoisotopicMass},error=function(cond){return(0)})
  OMF<-tryCatch({out$PropertyTable$Properties$MolecularFormula},error=function(cond){return(NA)})
  #############################
  return(c(tryCatch({OIN[1]},error = function(x) {return(NA)}),tryCatch({OIK[1]},error = function(x) {return(NA)}),tryCatch({OSM[1]},error = function(x) {return(NA)}),tryCatch({OCID[1]},error = function(x) {return(NA)}),tryCatch({EXM[1]},error = function(x) {return(0)}),tryCatch({OMF[1]},error = function(x) {return(NA)})))
  #############################
}

##############################################################################
##############################################################################
PuNAMEtoOI<-function(getNAME)
{
  ###################################	
  url<- "http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
  out<-tryCatch({jsonlite::fromJSON(paste0(url,getNAME, "/property/InChIKey"))},error = function(x) {return(NA)})
  out1<-tryCatch({jsonlite::fromJSON(paste0(url,getNAME, "/property/CanonicalSMILES"))}, error = function(x) {return(NA)})
  out2<-tryCatch({jsonlite::fromJSON(paste0(url,getNAME, "/property/InChI"))}, error = function(x) {return("NA")})
  ############### InChI
  OIN<-tryCatch({out2$PropertyTable$Properties$InChI},error=function(cond){return(NA)})
  ############### InchIkey
  OIK<-tryCatch({out$PropertyTable$Properties$InChIKey},error=function(cond){return(NA)})
  ############## SMILES
  OSM<-tryCatch({out1$PropertyTable$Properties$CanonicalSMILES},error=function(cond){return(NA)})
  ############# Compound CID
  OCID<-tryCatch({out1$PropertyTable$Properties$CID},error=function(cond){return(NA)})
  ########### Exact mass
  EXM<-tryCatch({ConvPCIDtoOCN(out1$PropertyTable$Properties$CID)[4]},error=function(cond){return(NA)})
  ###############################
  ###############################
  return(c(tryCatch({OIN[1]},error = function(x) {return(NA)}),tryCatch({OIK[1]},error = function(x) {return(NA)}),tryCatch({OSM[1]},error = function(x) {return(NA)}),tryCatch({OCID[1]},error = function(x) {return(NA)}),tryCatch({EXM[1]},error = function(x) {return(0)})))
  ###############################
  ###############################
}
#############################################################################
#############################################################################
getCactus <- function(identifier,representation){
  identifier <- gsub('#', '%23', identifier)
  ret <- tryCatch(httr::GET(paste("https://cactus.nci.nih.gov/chemical/structure/",
                                  URLencode(identifier), "/", representation, sep = "")),
                  error = function(e) NA)
  if (all(is.na(ret)))
    return(NA)
  if (ret["status_code"] == 404)
    return(NA)
  ret <- tryCatch({httr::content(ret)},error = function(x) {return(NA)})
  return(tryCatch({unlist(strsplit(ret, "\n"))},error = function(x) {return(NA)}))

}

############################################################################
############################################################################
PuInKtoIN<-function(getINK)
{
  url<- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"
  out<-tryCatch({jsonlite::fromJSON(paste0(url,getINK, "/JSON"))} ,error = function(x) {return(NA)})
  prop.names  <-tryCatch({out$PC_Compounds$props[[1]][[1]]},error = function(x) {return(NA)})
  prop.values <- tryCatch({out$PC_Compounds$props[[1]][[2]]},error = function(x) {return(NA)})
  InchiVal=tryCatch({prop.values[13,"sval"]},error = function(x) {return(NA)})
  return(tryCatch({InchiVal[1]},error = function(x) {return(NA)}))
}
###############################################################################
################################################################################
ConvINKtoOID<-function(getINK)
{
  ### ####This function return INCHIKEY to other identifiers like Inchi,Inchikey,Smiles,CompoundID
  url<- "https://www.metabolomicsworkbench.org/rest/compound/inchi_key/"
  out<-tryCatch({jsonlite::fromJSON(paste0(url,getINK, "/all"))}, error = function(x) {return(NA)})
  ###########################
  OIK<-tryCatch({out$inchi_key},error=function(cond){return(NA)})
  OSM<-tryCatch({out$smiles},error=function(cond){return(NA)})
  OCID<-tryCatch({out$pubchem_cid},error=function(cond){return(NA)})
  OEM<-tryCatch({out$exactmass},error=function(cond){return(NA)})
  OFOR<-tryCatch({out$formula},error=function(cond){return(NA)})
  ###########################
  return(c(tryCatch({OIK[1]},error = function(x) {return(NA)}),tryCatch({OSM[1]},error = function(x) {return(NA)}),tryCatch({OCID[1]},error = function(x) {return(NA)}),tryCatch({OEM[1]},error = function(x) {return(0)}),tryCatch({OFOR[1]},error = function(x) {return(NA)})))

  ###########################
}

#############################################################################
#############################################################################
ConvINKtoOID1<-function(getINK)
{
  ################################
  url<-"http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"
  out<-tryCatch({jsonlite::fromJSON(paste0(url,getINK, "/property/CanonicalSMILES,MonoisotopicMass,InChI,InChIKey"))}, error = function(x) {return(NA)})
  #################################
  OIK<-tryCatch({out$PropertyTable$Properties$InChIKey},error=function(cond){return(NA)})
  OIN<-tryCatch({out$PropertyTable$Properties$InChI},error=function(cond){return(NA)})
  OSM<-tryCatch({out$PropertyTable$Properties$CanonicalSMILES},error=function(cond){return(NA)})
  OCID<-tryCatch({out$PropertyTable$Properties$CID},error=function(cond){return(NA)})
  EXM<-tryCatch({out$PropertyTable$Properties$MonoisotopicMass},error=function(cond){return(0)})
  #################################
  return(c(tryCatch({OIN[1]},error = function(x) {return(NA)}),tryCatch({OIK[1]},error = function(x) {return(NA)}),tryCatch({OSM[1]},error = function(x) {return(NA)}),tryCatch({OCID[1]},error = function(x) {return(NA)}),tryCatch({EXM[1]},error = function(x) {return(0)})))
  #################################
}
#############################################################################
#############################################################################
PuSmilesToEM<-function(getSMILES)
{
  For<-tryCatch({RChemMass::MolFormFromSmiles.rcdk(getSMILES)},error=function(cond){return(NA)})
  EM<-tryCatch({Rdisop::getMolecule(For)},error=function(cond){return(NA)})
  EM1<-tryCatch({EM$exactmass},error=function(cond){return(0)})
  return(EM1)
}
###############################################################################
###############################################################################
MolFormFromSmiles.rcdk <- function(smiles) {

  mol <- tryCatch({rcdk::parse.smiles(smiles)[[1]]},error=function(cond){return(NA)})
  tryCatch({rcdk::convert.implicit.to.explicit(mol)},error=function(cond){return(NA)})
  charge1 <- tryCatch({rcdk::get.total.charge(mol)},error=function(cond){return(NA)})
  formula <- tryCatch({rcdk::get.mol2formula(mol, charge=charge1)},error=function(cond){return(NA)})
  return(tryCatch({formula@string},error=function(cond){return(NA)}))
}
##############################################################################
##############################################################################
ClassSmilesToOntolgy<-function(getSMILE)
{
  ### ####This function return smiles to ontology
  #####################
  url="https://gnps-structure.ucsd.edu/classyfire?smiles="
  url1=paste0(url,getSMILE)
  out=tryCatch({jsonlite::fromJSON(url1)},error = function(x) {return(NA)})
  res=do.call(paste, c(as.list(tryCatch({rev(out$ancestors)},error=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
  ####################
  return(res)
  #####################
}
#############################################################################
#############################################################################
PuSMtoCID<-function(getSMILES)
{
  url<-"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/"
  out<-tryCatch({jsonlite::fromJSON(paste0(url,getSMILES, "/cids"))} ,error = function(x) {return(NA)})
  return(tryCatch({out[[1]]$CID},error = function(x) {return(NA)}))

}

#############################################################################
#############################################################################
PuSMItoINKandInchi<-function(getSMI)
{
  url<- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/SMILES/"
  out<-tryCatch({jsonlite::fromJSON(paste0(url,getSMI, "/JSON"))} ,error = function(x) {return(NA)})
  prop.names  <-tryCatch({out$PC_Compounds$props[[1]][[1]]},error = function(x) {return(NA)})
  prop.values <- tryCatch({out$PC_Compounds$props[[1]][[2]]},error = function(x) {return(NA)})
  InchiVal=tryCatch({prop.values[13,"sval"]},error = function(x) {return(NA)})
  InchikeyVal=tryCatch({prop.values[14,"sval"]},error = function(x) {return(NA)})
  return(c(InchiVal,InchikeyVal))
}
##############################################################################
##############################################################################
ConvSMItoOID1<-function(getSMI)
{
  url<-"http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/"
  out<-tryCatch({jsonlite::fromJSON(paste0(url,getSMI,"/property/CanonicalSMILES,MonoisotopicMass,InChI,InChIKey"))},error = function(x) {return(NA)})
  ########################
  OIK<-tryCatch({out$PropertyTable$Properties$InChIKey},error=function(cond){return(NA)})
  OIN<-tryCatch({out$PropertyTable$Properties$InChI},error=function(cond){return(NA)})
  OSM<-tryCatch({out$PropertyTable$Properties$CanonicalSMILES},error=function(cond){return(NA)})
  OCID<-tryCatch({out$PropertyTable$Properties$CID},error=function(cond){return(NA)})
  EXM<-tryCatch({out$PropertyTable$Properties$MonoisotopicMass},error=function(cond){return(0)})
  ########################
  return(c(tryCatch({OIN[1]},error = function(x) {return(NA)}),tryCatch({OIK[1]},error = function(x) {return(NA)}),tryCatch({OSM[1]},error = function(x) {return(NA)}),tryCatch({OCID[1]},error = function(x) {return(NA)}),tryCatch({EXM[1]},error = function(x) {return(0)})))
  ########################
}
##############################################################################
##############################################################################
ConvHMDBtoOCN<-function(getHMDB)
{
  url<- "https://www.metabolomicsworkbench.org/rest/compound/hmdb_id/"
  out<-tryCatch({jsonlite::fromJSON(paste0(url,getHMDB, "/all"))}, error = function(x) {return(NA)})
  #########################
  #########################
  OIK<-tryCatch({out$inchi_key},error=function(cond){return(NA)})
  OSM<-tryCatch({out$smiles},error=function(cond){return(NA)})
  OCID<-tryCatch({out$pubchem_cid},error=function(cond){return(NA)})
  OEM<-tryCatch({out$exactmass},error=function(cond){return(NA)})
  OFOR<-tryCatch({out$formula},error=function(cond){return(NA)})
  return(c(OIK[1],OSM[1],OCID[1],OEM[1],OFOR[1]))
  ############################
  ############################
}
##############################################################################
##############################################################################
PuCIDtoEM<-function(getCID)
{
  #######################
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
  inVa<-tryCatch({which(data2$name %in% "Record.Section.Section.Section.Information.Value.StringWithMarkup.String")},error = function(x) {return(NA)})
  ##########################
  TEST<-tryCatch({sapply(data2, "[", inVa)},error = function(x) {return(NA)})
  TEST1<-tryCatch({grep("InChI=",TEST)},error = function(x) {return(NA)})
  TEST2<-tryCatch({TEST1[1]},error = function(x) {return(NA)})
  InchI<-tryCatch({TEST[TEST2]},error = function(x) {return(NA)})
  inchikey<-tryCatch({TEST[TEST2+1]},error = function(x) {return(NA)})
  #########################
  url<- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"
  #########################
  out<-tryCatch({jsonlite::fromJSON(paste0(url,inchikey, "/JSON"))},error = function(x) {return(NA)})
  #########################
  EMV<-tryCatch({out$PC_Compounds$props[[1]][22,]$value$sval},error = function(x) {return(NA)})
  ########################
  EMV1<-c()
  ########################
  if(!sjmisc::is_empty(EMV))
  {
    EMV1<-c(EMV1,EMV)
  }else{
    EMV1<-c(EMV1,0)
  }
  ########################
  return(EMV1)
  ########################
  ########################
}
################################################################################
##############################################################################
PuCAStoOI<-function(getCAS)
{
  ###CAS: 328-50-7
  getCAS1<-stringr::str_replace(getCAS,pattern='CAS:',replacement ="")
  getCAS2<-stringr::str_trim(getCAS1)
  ################################
  url<- "http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
  out<-tryCatch({jsonlite::fromJSON(paste0(url,getCAS2, "/property/InChIKey"))}, error = function(x) {return(NA)})
  out1<-tryCatch({jsonlite::fromJSON(paste0(url,getCAS2, "/property/CanonicalSMILES"))}, error = function(x) {return(NA)})
  out2<-tryCatch({jsonlite::fromJSON(paste0(url,getCAS2, "/property/InChI"))}, error = function(x) {return(NA)})
  #################################https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/
  OIN<-tryCatch({out2$PropertyTable$Properties$InChI},error=function(cond){return(NA)})
  OIK<-tryCatch({out$PropertyTable$Properties$InChIKey},error=function(cond){return(NA)})
  OSM<-tryCatch({out1$PropertyTable$Properties$CanonicalSMILES},error=function(cond){return(NA)})
  OCID<-tryCatch({out1$PropertyTable$Properties$CID},error=function(cond){return(NA)})
  ################################
  ################################
  return(c(tryCatch({OIN[1]},error=function(cond){return(NA)}),tryCatch({OSM[1]},error=function(cond){return(NA)}),tryCatch({OCID[1]},error=function(cond){return(NA)})))
}
###############################################################################
###############################################################################
getOntoSM<-function(SM)
{
tes1<-tryCatch({ClassSmilesToOntolgy(SM)},error=function(cond){message("smiles to ontology is failing")})
ts<-tes1
ts1<-unlist(strsplit(ts, ";"))
ts2<-ts1[order(grepl("^O",ts1),ts1,decreasing =T)]
ts3=ts2[grepl("^O",ts2)]
t1s3<-ts3[order(grepl("^Organic compounds",ts3),ts3,decreasing =T)]
ts4=ts2[!grepl("^O",ts2)]
ts5=sort(ts4[!grepl("^[0-9]+",ts4)])
ts6=ts4[grepl("^[0-9]+",ts4)]
ts7=c(ts5,ts6)
ts8=c(t1s3,ts7)
ts9=paste(ts8, collapse = ";")
tes2<-paste("Ontology:",ts9,sep=" ")
return(tes2)
}
##############################################################################
##############################################################################
###print("this function remove the mz if it has no intensity value")
PeakVali<-function(Fpea)
{

  PV1=c()
  for(i in 1:length(Fpea))
  {
    
    PV=Fpea[i]
   
    tes1<-unlist(strsplit(PV, "\t|\t\t|;|,"))
   
   
    if(length(tes1)==2)
    {
     
     
      PV1=c(PV1,PV)
    }
    
  }


  return(PV1)
}

###################################################################################
###################################################################################

###print("This is the peak centroiding function that I wrote")
cenPeaks <- function(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE )
{
################################
mzV1<-c(0,mzV)
###########################
###########################
av<-abs(mzV[2:length(mzV)] - mzV[1:(length(mzV)-1)]) <= MZ_TOLERANCE
###########################
###########################
Ntes2=mzV[which(!av)]
Ntes3=intenV[which(!av)]
###########################
CVL=split(which(av), cumsum(c(1, diff(which(av)) != 1)))
###########################
###########################
NMZV1=c()
INDV1=c()
############################
TV=c()
for(i in 1:length(CVL))
{
  LV=CVL[i]
  LV1=purrr::flatten_dbl(LV)
  LV2=tail(LV1,n=1)
  LV3=LV2+1
  NLV=unname(unlist(c(LV,LV3), recursive = FALSE))
  TV=append(TV,NLV)
  #############################
  #############################
  weightsV=intenV[NLV]
  dataV=mzV[NLV]
  weight_sum = sum(weightsV)
  #######################################
  #######################################

  NDV=sum(dataV)/length(dataV)
  NIV=sum(weightsV)/length(weightsV)

  NMZV1=append(NMZV1,NDV)
  INDV1=append(INDV1,NIV)


}## end of for loop

NCVLIN=seq(1,length(mzV))
DIFFIN=setdiff(NCVLIN,TV)

for(i in 1:length(DIFFIN))
{
  Val=DIFFIN[i]
  weightsV=intenV[Val]
  dataV=mzV[Val]

  ##print(dataV)
  NMZV1=append(NMZV1,dataV)
  INDV1=append(INDV1,weightsV)

}


N1tes2=NMZV1
N1tes3=INDV1


N2tes2=N1tes2[order(N1tes2)]
N2tes3=N1tes3[order(N1tes2)]

N2tes4=N2tes2[!is.na(N2tes2)]
N2tes5=N2tes3[!is.na(N2tes2)]

return(list(N2tes4,as.numeric(N2tes5)))

}

##############################################################################
##############################################################################
Centroid<-function(Np,tes2,tes3,tes4,NTES,NTES4,NTES5,DEFAULT_MZ_TOLERANCE,Fpea2)
{

  out<-c()

  if((Np >= 60) & (length(NTES) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))
  {

   	NTES1<-c(NTES4,tes4)


        if((Np >= 60) & (length(NTES1) > 1))
  	{


		tes5<-tes2[-NTES1]
  		tes6<-tes3[-NTES1]

  		mzV<-tes5
  		intenV<-tes6



  		if(length(mzV) > 1)
  		{


    			nMZV=tryCatch({cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]},error=function(cond){return(0)})
    			nINV=tryCatch({cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]},error=function(cond){return(0)})
    			tes7<-paste(nMZV,nINV,sep="\t")
    			F1NPA<-paste0("Num Peaks: ",length(tes7))

    			out<-c(out,F1NPA)
    			out<-c(out,tes7)

  		}else{

    			tes5<-tes2
    			tes6<-tes3

    			nMZV=tryCatch({cenPeaks(tes2,tes3, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]},error=function(cond){return(0)})
    			nINV=tryCatch({cenPeaks(tes2,tes3, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]},error=function(cond){return(0)})

    			if(length(nMZV) > 1 & length(nINV) > 1)
    			{

      				tes7<-paste(nMZV,nINV,sep="\t")

      				F1NPA<-paste0("Num Peaks: ",length(tes7))
      				out<-c(out,F1NPA)
      				out<-c(out,tes7)

    			}else{

    				tes7<-paste(tes5,tes6,sep="\t")
    				F1NPA<-paste0("Num Peaks: ",length(tes7))
    				out<-c(out,F1NPA)
    				out<-c(out,tes7)
    			}## end of if else

  			}###end of else



  		}else{




			if((length(tes2[-NTES]) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))
  			{
    				tes5<-tes2[-NTES]
    				tes6<-tes3[-NTES]
    				tes7<-paste(tes5,tes6,sep="\t")
    				F1NPA<-paste0("Num Peaks: ",length(tes7))
    				out<-c(out,F1NPA)
    				out<-c(out,tes7)

  			}else{

    				tes7<-paste(tes2,tes3,sep="\t")
				F1NPA<-paste0("Num Peaks: ",length(tes7))
    				out<-c(out,F1NPA)
    				out<-c(out,tes7)

  			}## end of if else loop

  			}### end of the else





		  }else{

		  if((length(NTES) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))
                  {
			  if(length(tes2[-NTES]) > 1)
    			  {

      				tes5<-tes2[-NTES]
      				tes6<-tes3[-NTES]
      				tes7<-paste(tes5,tes6,sep="\t")
      				### Adding this new here
      				F1NPA<-paste0("Num Peaks: ",length(tes7))
      				out<-c(out,F1NPA)
      				#############################
      				out<-c(out,tes7)
      				##############################
    			 }else{

      				tes5<-tes2
      				tes6<-tes3
      				tes7<-paste(tes5,tes6,sep="\t")
      				F1NPA<-paste0("Num Peaks: ",length(tes7))
      				out<-c(out,F1NPA)
      				out<-c(out,tes7)

    			}### if,else loop
                    ############################
		    ############################
                  }else{

			  ##if((Np >= 80))
			  if((Np >= 60))
			  {
				  mzV<-tes2
				  intenV<-tes3

				  nMZV=cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]
				  nINV=cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]

				  tes7<-paste(nMZV,nINV,sep="\t")
				  F1NPA<-paste0("Num Peaks: ",length(tes7))

				  out<-c(out,F1NPA)
				  out<-c(out,tes7)

			  }else{

		             	########################
		    		F1NPA<-paste0("Num Peaks: ",tryCatch({length(Fpea2)},error=function(cond){message("Fpea is empty")}))
	            		out<-c(out,F1NPA)
	            		out<-c(out,Fpea2)
		    		########################


                     	}## end of Np >= 80 and else
                 	############################

                  }### end of else


		  }


  return(list(out[1],out[2:length(out)]))
}

###############################################################################
###############################################################################
### 2. Converting the Identifiers in the Metadata Sheet
###############################################################################
###############################################################################
FuFtoRe<-function(InMEDA)
{
  ####This function returns exact mass given the Name of the compound
  if(!sjmisc::is_empty(as.character(InMEDA[["Name"]])) & !startsWith(as.character(InMEDA[["Name"]]),'not available')){
    getNAME<-stringr::str_trim(as.character(InMEDA[["Name"]]))
    GCID<-tryCatch({webchem::get_cid(InMEDA[["Name"]])},error=function(cond){return(NA)})
    GCID1<-tryCatch({GCID[[2]][1]},error=function(cond){message("Name value is empty")})
    #####G1CID1<-stringr::str_trim(GCID1)
    G1CID1<-stringr::str_trim(gsub("[[:punct:]]", "",GCID1))
    PCID1<-tryCatch({webchem::pc_prop(as.numeric(G1CID1), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){return(NA)})
    PCID2<-tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){return(NA)})
    ############################
    FMa<-c()
    ############################
    if(!sjmisc::is_empty(PCID2)){
      FMa<-c(FMa,PCID2)
    }else{
      EMV<-tryCatch({PuNAMEtoOI(getNAME)},error = function(x) {return(NA)})
      if(!sjmisc::is_empty(EMV)){
        FMa<-c(FMa,EMV[5])
      }else{
        if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){
          EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem CId is empty")})
          EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
          if(!sjmisc::is_empty(EM1)){
            FMa<-c(FMa,EM1)
          }else{
            EM<-stringr::str_trim(as.character(InMEDA[["Exact mass"]]))
            if(!sjmisc::is_empty(EM)){
              FMa<-c(FMa,EM)
            }else{
              FMa<-c(FMa,0)
            }
          }###end of else
        }##end of ifloop
      }
    }###end of else PCID2

}else{
  if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){
    EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem CId value is empty")})
    EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem CId is value1 empty")})
    if(!sjmisc::is_empty(EM1)){
      FMa<-c(FMa,EM1)
    }else{
      EM<-stringr::str_trim(as.character(InMEDA[["Exact mass"]]))
      if(!sjmisc::is_empty(EM)){
        FMa<-c(FMa,EM)
      }else{
        FMa<-c(FMa,0)
      }
    }###end of else
  }##end of ifloop
}
  return(FMa)
}
#########################################################################################
#########################################################################################
CONcidtoEM<-function(InMEDA)
{
  #################################
  FMa<-c()
  ################################
if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(InMEDA[["PubChem CID"]]),'not available')){
  ########################################
  ##FPUCID<-stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))
  FPUCID<-stringr::str_trim(gsub("[[:punct:]]", "",as.character(InMEDA[["PubChem CID"]])))
  ########################################
  PCID1<-tryCatch({webchem::pc_prop(as.numeric(FPUCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem CId is empty..did not get exact mass")})
  PCID2<-tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){message("Pubchem CID is empty..did not get exact mass")})
  ######################################
  ######################################
  if(!sjmisc::is_empty(PCID2)){
    ############################################	
    print("entering the if loop..Pubchem CID")
    ###########################################
    FMa<-c(FMa,PCID2)
  }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(InMEDA[["PubChem CID"]]),'not available') & tryCatch({ConvCIDtoOID1(INMEDA[["PubChem CID"]])[5]},error = function(x) {return(NA)}) != 0){
    PCID2<-tryCatch({ConvCIDtoOID1(INMEDA[["PubChem CID"]])[5]},error = function(x) {return(NA)})
    FMa<-c(FMa,PCID2)
    ######################################
  }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(InMEDA[["PubChem CID"]]),'not available') & tryCatch({ConvPCIDtoOCN(INMEDA[["PubChem CID"]])[4]},error=function(cond){return(NA)}) != 0){
    PCID2<-tryCatch({ConvPCIDtoOCN(INMEDA[["PubChem CID"]])[4]},error=function(cond){return(NA)})
    FMa<-c(FMa,PCID2)
  }else{
    ### PubchemId is not found ..entering the else loop
    if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["Formula"]]))) & !startsWith(stringr::str_trim(as.character(InMEDA[["Formula"]])),'not available')){
      ############################################
      ############################################	    
      EM<-tryCatch({Rdisop::getMolecule(stringr::str_trim(as.character(InMEDA[["Formula"]])))},error=function(cond){message("Pubchem Id is empty")})
      EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
      ############################################
      ############################################
      if(!sjmisc::is_empty(EM1)){
        FMa<-c(FMa,EM1)
      }else{
        #############################
        PMA<-tryCatch({FuFtoRe(InMEDA)},error=function(cond){message("NAME to exact mass...FuFtoRe")})
        FMa<-c(FMa,PMA)
        #############################
      }###end of else
    }else{
      PMA<-tryCatch({FuFtoRe(InMEDA)},error=function(cond){message("NAME to exact mass...FuFtoRe")})
      FMa<-c(FMa,PMA)
    }
  }##end of else
}else{
  if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["Formula"]]))) & !startsWith(stringr::str_trim(as.character(InMEDA[["Formula"]])),'not available')){
    ############################################
    ############################################	  
    EM<-tryCatch({Rdisop::getMolecule(stringr::str_trim(as.character(InMEDA[["Formula"]])))},error=function(cond){message("Pubchem CId is empty")})
    EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
    ############################################
    ############################################
    if(!sjmisc::is_empty(EM1)){
      FMa<-c(FMa,EM1)
    }else{
      #############################
      PMA<-tryCatch({FuFtoRe(InMEDA)},error=function(cond){message("FUFtoRe value is empty...1")})
      FMa<-c(FMa,PMA)
      #############################
    }###end of else
  }else{
    PMA<-tryCatch({FuFtoRe(InMEDA)},error=function(cond){message("FUFtoRe value is empty...2")})
    FMa<-c(FMa,PMA)
  }
}###end of else
  #################################
  return(FMa)
  ################################
}
##############################################################################
##############################################################################
CONSMItoEM<-function(InMEDA)
{
  ###############################
  FMa<-c()
  ###############################
if(!sjmisc::is_empty(as.character(InMEDA[["SMILES"]])) & !startsWith(as.character(InMEDA[["SMILES"]]),'not available')){
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
  }else if(!sjmisc::is_empty(as.character(InMEDA[["SMILES"]])) & !startsWith(as.character(InMEDA[["SMILES"]]),'not available') & !sjmisc::is_empty(tryCatch({ConvSMItoOID1(InMEDA[["SMILES"]])[5]}, error = function(x) {return(0)}))){
     print("entering the smiles pass in CONSMItoEM")
    PCID2<-ifelse(!sjmisc::is_empty(tryCatch({ConvSMItoOID1(InMEDA[["SMILES"]])[5]}, error = function(x) {return(0)})),tryCatch({ConvSMItoOID1(InMEDA[["SMILES"]])[5]}, error = function(x) {return(0)}),tryCatch({mzAnnotation::smileToAccurateMass(as.character(InMEDA[["SMILES"]]))}, error = function(x) {return(0)}))
    ##PCID2<-tryCatch({ConvSMItoOID1(InMEDA[["SMILES"]])[5]}, error = function(x) {return(0)})
    FMa<-c(FMa,PCID2)
  }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(InMEDA[["PubChem CID"]]),'not available')){
    PCID2<-tryCatch({CONcidtoEM(InMEDA)},error=function(cond){message("Pubchem CID to exact mass failed")})
    FMa<-c(FMa,PCID2)
  }else{
    PMA<-tryCatch({FuFtoRe(InMEDA)},error=function(cond){message("NAME to exact mass failed")})
    FMa<-c(FMa,PMA)
  }
}else{
  ### Smiles not found
  PCID2<-tryCatch({CONcidtoEM(InMEDA)},error=function(cond){message("pubchem CID to exact mass failed")})
  FMa<-c(FMa,PCID2)
}
  ############################
  return(FMa)
  ############################
}
#############################################################################
#############################################################################
getOntFromSMI<-function(IK,InMEDA)
{
  out<-c()
  if(!sjmisc::is_empty(IK)){
  ##################################
  IKCRV<-tryCatch({classyfireR::get_classification(IK)},error=function(cond){message("Classifier could not fetch the information")})
  ####################################
  if(!sjmisc::is_empty(IKCRV)){
    ####################
    ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
    ##################################
    IK1<-paste("INCHIKEY:",IK,sep=" ")
    tes2<-paste("Ontology:",ONTV,sep=" ")
    FINCH<-paste("INCHI:",tryCatch({PuInKtoIN(IK)},error = function(x) {return(NA)}),sep=" ")
    #################################
    out<-c(out,tes2)
    out<-c(out,IK1)
    out<-c(out,FINCH)
    #########################
  }else{
    SM<-tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)})
    IN<-paste("INCHI:",tryCatch({rinchi::get.inchi(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)}),sep=" ")
    IK<-paste("INCHIKEY:",tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)}) ,sep=" ")
    tes1<-tryCatch({ClassSmilesToOntolgy(SM)},error=function(cond){message("Classyfire smiles to Ontology empty")})
    #####################################
    tes2<-getOntoSM(SM)
    ######################################
    out<-c(out,tes2)
    out<-c(out,IK)
    out<-c(out,IN)
  }
  }
  return(out)
}
#########################################################################
#########################################################################
###  This is a function that will select Inchi or INCHIKEY value
#########################################################################
SelectInchiorInchiKey<-function(InMEDA,IK,InchiV)
{
  out<-c()
if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & !startsWith(as.character(InMEDA[["InChI"]]),'not available') & !startsWith(as.character(InMEDA[["InChI"]]),'CAS:') & !startsWith(as.character(InMEDA[["InChI"]]),'InChI=')){
  if(tryCatch({webchem::is.inchikey(stringr::str_trim(as.character(InMEDA[["InChI"]])))},error=function(cond){message("inchikey..file must be empty")})){
    FINK<-ifelse(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])),paste("INCHIKEY:",as.character(InMEDA[["InChI"]])),paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" "))
    FINCH<-paste("INCHI:",InchiV,sep=" ")
    out<-c(out,FINK)
    out<-c(out,FINCH)
  }else{
    FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
    FINCH<-ifelse(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])),paste("INCHI:",as.character(InMEDA[["InChI"]])),paste("INCHI:",InchiV,sep=" "))
    out<-c(out,FINK)
    out<-c(out,FINCH)
  }
}else{
  ####print("enter the else value")
  FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
  out<-c(out,FINK)
  FINCH<-paste("INCHI:",InchiV,sep=" ")
  out<-c(out,FINCH)

}
  return(out)
}
########################################################################
########################################################################
CONcastoEM<-function(InMEDA)
{
  ########################
  FMa<-c()
  ######################
if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["InChI"]]))) & startsWith(as.character(InMEDA[["InChI"]]),'CAS:')){
  ##############################################
  print("enter the Inchi/CAS part ...in loop...CAS area")
  ###############################################
  CV<-stringr::str_trim(as.character(InMEDA[["InChI"]]))
  CV1<-stringr::str_replace(CV,pattern='CAS:',replacement ="")
  CV2<-stringr::str_trim(as.character(CV1))
  SM<-tryCatch({webchem::cir_query(CV2, "smiles")[[2]]},error=function(cond){message("CAS not abe to fetch smiles")})
  EMS<-tryCatch({mzAnnotation::smileToAccurateMass(SM)}, error = function(x) {return(NA)})
  #############################################
  PCID<-tryCatch({webchem::get_cid(CV2, from = "xref/rn",match="first")},error=function(cond){message("Pubchem Id is empty")})
  ##########################################
  ##########################################
  PCID1<-tryCatch({webchem::pc_prop(as.numeric(PCID[[2]][1]), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES","InChI","InChIKey"))},error=function(cond){message("Pubchem CID is empty")})
  PCID2<-tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){message("Pubchem CID is empty")})
  #############################################
  print("smiles to exact mass value")
  prit(EMS)
  #############################################
  if(!sjmisc::is_empty(PCID2)){
    ###########################################	    
    print("enter the if loop ...cas area")
    ########################################	    
    FMa<-c(FMa,PCID2)
  }else if(!sjmisc::is_empty(EMS)){
	  #######################################################
	  print("entering else if smile to exact value is there")
	  #######################################################
	  FMa<-c(FMa,EMS)
  }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["InChI"]]))) & startsWith(as.character(InMEDA[["InChI"]]),'CAS:') & !sjmisc::is_empty(tryCatch({PuCAStoOI(CV2)[5]}, error = function(x) {return(0)}))){
    ############################################################	  
    print("Enter the else if loop INCHI/CAS value is avilable")
    ##############################################################  
    PCID2<-tryCatch({PuCAStoOI(InMEDA[["InChI"]])[5]}, error = function(x) {return(0)})
    FMa<-c(FMa,PCID2)
  }else if(!sjmisc::is_empty(as.character(InMEDA[["SMILES"]])) & !startsWith(as.character(InMEDA[["SMILES"]]),'not available')){
    ######################################
    print("entering smiles avilable area in InchiKey/CAS pass area")
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
    }else if(!sjmisc::is_empty(as.character(InMEDA[["SMILES"]])) & !startsWith(as.character(InMEDA[["SMILES"]]),'not available')){
	########################################################################    
        print("entering the else if loop ... Inchi/CAS part...smiles avilable")
        ########################################################################
      PCID2<-ifelse(!sjmisc::is_empty(tryCatch({ConvSMItoOID1(InMEDA[["SMILES"]])[5]}, error = function(x) {return(0)})),tryCatch({ConvSMItoOID1(InMEDA[["SMILES"]])[5]}, error = function(x) {return(0)}),tryCatch({PuSmilesToEM(InMEDA[["SMILES"]])}, error = function(x) {return(0)}))
      FMa<-c(FMa,PCID2)
    }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(InMEDA[["PubChem CID"]]),'not available')){
      PCID2<-tryCatch({CCONcidtoEM(InMEDA)},error=function(cond){message("CID not abe to fetch")})
      FMa<-c(FMa,PCID2)
    }else{
	    ############################################################
	    print("enter the else part ...CAS to exact mass ..smiles")
	    ############################################################
            PMA<-tryCatch({FuFtoRe(InMEDA)},error=function(cond){message("smiles not abe to fetch")})
            FMa<-c(FMa,PMA)
    }
  }else{
    ###########################################################	  
    print("smiles is not avilable..CAS to exact mass..") 
    ###########################################################
    PCID2<-tryCatch({CONcidtoEM(InMEDA)},error=function(cond){message("CID is not able to convert exact mass")})
    FMa<-c(FMa,PMA)
  }
}else{
  #######################################################	
  print("CAS not found.... CAS to exact mass...")
  #######################################################
  PCID2<-tryCatch({CONSMItoEM(InMEDA)},error=function(cond){message("convert smiles to exact mass failing")})
  FMa<-c(FMa,PCID2)
}
  ###################
  return(FMa)
  ###################
}
##############################################################################
##############################################################################
CONinctoEM<-function(InMEDA)
{
  #########################
  FMa<-c()
  ##########################
if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),'InChI=')){
  ##############################################
  print("enter the Inchi part ...in else loop...Inchikey")
  IK<-as.character(InMEDA[["InChI"]])
  ##############################################
  IK1<-tryCatch({webchem::get_cid(IK, from = "inchi")},error=function(cond){message("Inchi name must be empty or rinchi not abe to fetch")})
  PCID<-tryCatch({IK1[[2]][1]},error=function(cond){return(NA)})
  #######################################
  #######################################
  #########PCID1<-tryCatch({webchem::pc_prop(as.numeric(PCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem Id is empty")})
  #######################################
  #######################################
  PCID1<-tryCatch({webchem::pc_prop(as.numeric(stringr::str_trim(gsub("[[:punct:]]", "",PCID))), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){return(NA)})
  PCID2<-tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){return(NA)})
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
    PCID1<-tryCatch({webchem::pc_prop(as.numeric(PCID[[2]][1]), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES","InChI","InChIKey"))},error=function(cond){message("Pubchem CID is empty")})
    #############################################
    ##############################################
    FPUCID1<-tryCatch({as.numeric(PCID[[2]][1])},error=function(cond){message("Pubchem CID is empty")})
    PCID2<-ifelse(!sjmisc::is_empty(as.numeric(PCID1$ExactMass)),as.numeric(PCID1$ExactMass),ifelse(!sjmisc::is_empty(tryCatch({as.numeric(ConvCIDtoOID1(FPUCID1)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")})),tryCatch({as.numeric(ConvCIDtoOID1(FPUCID1)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}),tryCatch({PuCIDtoEM(FPUCID1)},error=function(cond){message("Pubchem CID is empty..did not get exact mass1")})))
    #############################################
    #############################################
    if(!sjmisc::is_empty(PCID2) & (PCID2 !=0)){
      ########################################
      print("enter the if loop ...cas area")
      ########################################
      FMa<-c(FMa,PCID2)
    }else if(!sjmisc::is_empty(as.character(InMEDA[["SMILES"]])) & !startsWith(as.character(InMEDA[["SMILES"]]),'not available')){
      ######################################
      PCID2<-tryCatch({CONSMItoEM(InMEDA)},error=function(cond){message("convert smiles to exact mass is empty")})
      FMa<-c(FMa,PCID2)
     }else{
      #####################################
       PCID2<-tryCatch({CONcidtoEM(InMEDA)},error=function(cond){message("convert cid to exact mass is empty")})
       FMa<-c(FMa,PCID2)
    ##################################
    } ##end if loop Pubchem CID
  #####################################
  }else{
    PMA<-tryCatch({FuFtoRe(InMEDA)},error=function(cond){message("FuFtoRe exact mass is empty..3")})
    FMa<-c(FMa,PMA)
  }
}else{
  PCID2<-tryCatch({CONSMItoEM(InMEDA)},error=function(cond){message("convert smiles to exact mass is empty")})
  FMa<-c(FMa,PCID2)
}
  ############################
  return(FMa)
  #############################
}
##############################################################################
##############################################################################
CONinktoEM<-function(InMEDA)
{
  #########################
  FMa<-c()
  ##########################
if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & !startsWith(as.character(InMEDA[["InChI"]]),'not available') & !startsWith(as.character(InMEDA[["InChI"]]),'CAS:') & !startsWith(as.character(InMEDA[["InChI"]]),'InChI=')){
    #######################################################################	
  if(tryCatch({webchem::is.inchikey(stringr::str_trim(as.character(InMEDA[["InChI"]])))},error=function(cond){message("inchikey..file must be empty")})){
    ########################################################################
    print("enter the Inchi key AREA..inchikey value is not empty")
    ########################################################################
    IK<-as.character(InMEDA[["InChI"]])
    ######################
    ######################
    IK1<-tryCatch({webchem::get_cid(IK, from = "inchikey")},error=function(cond){message("Inchi name must be empty or rinchi not abe to fetch")})
    PCID<-tryCatch({IK1[[2]][1]},error=function(cond){message("Pubchem Id is empty")})
    #########################################################################
    PCID1<-tryCatch({webchem::pc_prop(as.numeric(stringr::str_trim(gsub("[[:punct:]]", "",PCID))), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem Id is empty")})
    #########################################################################
    ##########################################################################
    PCID2<-tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){message("Pubchem Id is empty")})
    ##########################################################################
    ##########################################################################
    if(!sjmisc::is_empty(PCID2))
    {
      #####################################
      print("enter the if loop...inkikey")
      #####################################
      FMa<-c(FMa,PCID2)
    }else{

      if(tryCatch({webchem::is.inchikey(stringr::str_trim(as.character(InMEDA[["InChI"]])))},error=function(cond){message("inchikey..file must be empty")}) & tryCatch({ConvINKtoOID(InMEDA[["InChI"]])[4]}, error = function(x) {return(0)}) != 0){
        PCID2<-tryCatch({ConvINKtoOID(InMEDA[["InChI"]])[4]}, error = function(x) {return(0)})
        FMa<-c(FMa,PCID2)
        }else if(tryCatch({webchem::is.inchikey(stringr::str_trim(as.character(InMEDA[["InChI"]])))},error=function(cond){message("inchikey..file must be empty")}) & tryCatch({ConvINKtoOID1(InMEDA[["InChI"]])[5]}, error = function(x) {return(0)})){
        PCID2<-tryCatch({ConvINKtoOID1(InMEDA[["InChI"]])[5]}, error = function(x) {return(0)})
        FMa<-c(FMa,PCID2)
        }else{
          PCID2<-tryCatch({CONcastoEM(InMEDA)},error=function(cond){message("CAS to exact mass is empty...1")})
          FMa<-c(FMa,PCID2)
        }



  } ### end of inchikey
  }else{
    PCID2<-tryCatch({CONcastoEM(InMEDA)},error=function(cond){message("CAS to exact mass is empty.....2")})
    FMa<-c(FMa,PCID2)
    ## No inchikey found

  }## end of else
}else{
  ### Inchikey and inchi is not found
  PCID2<-tryCatch({CONSMItoEM(InMEDA)},error=function(cond){message("SMILES to exact mass is empty")})
  FMa<-c(FMa,PCID2)
}
  ###########################
  return(FMa)
  ###########################
}
###############################################################################
###############################################################################
FuFtoRe1<-function(InMEDA)
{
  #############################
  FMa<-c()
  ############################
  if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & !startsWith(as.character(InMEDA[["InChI"]]),'not available') & !startsWith(as.character(InMEDA[["InChI"]]),'CAS:') & !startsWith(as.character(InMEDA[["InChI"]]),'InChI=')){
    if(tryCatch({webchem::is.inchikey(stringr::str_trim(as.character(InMEDA[["InChI"]])))},error=function(cond){message("inchikey..file must be empty")})){
      PCID2<-tryCatch({CONinktoEM(InMEDA)},error=function(cond){message("Inchikey to exact mass is empty")})
      FMa<-c(FMa,PCID2)
    }else{
      PCID2<-tryCatch({CONinctoEM(InMEDA)},error=function(cond){message("Inchi to exact mass is empty")})
      FMa<-c(FMa,PCID2)
    }

  }else{
    PCID2<- tryCatch({CONcastoEM(InMEDA)},error=function(cond){message("CAS to exact mass is empty")})
    
    if(!sjmisc::is_empty(PCID2)){

            FMa<-c(FMa,PCID2)

      }else{


	      if(!sjmisc::is_empty(as.character(InMEDA[["SMILES"]])) & !startsWith(as.character(InMEDA[["SMILES"]]),'not available')){
		         
		         PCID2<-ifelse(!sjmisc::is_empty(tryCatch({CONSMItoEM(InMEDA)},error=function(cond){message("SMILES to exact mass is empty")})),tryCatch({CONSMItoEM(InMEDA)},error=function(cond){message("SMILES to exact mass is empty")}),ifelse(!sjmisc::is_empty(tryCatch({mzAnnotation::smileToAccurateMass(as.character(InMEDA[["SMILES"]]))}, error = function(x) {return(0)})), tryCatch({mzAnnotation::smileToAccurateMass(as.character(InMEDA[["SMILES"]]))}, error = function(x) {return(0)}), tryCatch({PuSmilesToEM(InMEDA[["SMILES"]])},error = function(x) {return(0)})))
			 
			 FMa<-c(FMa,PCID2)

               }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(InMEDA[["PubChem CID"]]),'not available')){
                      PCID2<-tryCatch({CONcastoEM(InMEDA)},error=function(cond){message("CAS to exact mass is empty.....3")})
	              FMa<-c(FMa,PCID2)
                }else{
			PCID2<-tryCatch({FuFtoRe(InMEDA)},error=function(cond){message("NAME to exact mass is empty")})
                     }

           }###end of the else  
  }
 #############################
  return(FMa)
 ############################
}
########################################################################
########################################################################
gETSmiles<-function(InMEDA)
{
  ######################
  out<-c()
  #######################
  if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["SMILES"]]))))
  {
	  SM <- stringr::str_trim(as.character(InMEDA[["SMILES"]]))
	  out<-c(out,SM)

  }else if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){
    IN <- stringr::str_trim(as.character(InMEDA[["InChI"]]))
    mol <-tryCatch({rinchi::parse.inchi(IN)},error=function(cond){message("parese inchi conversion issue")})
    SM <- tryCatch({rcdk::get.smiles(mol[[1]])},error=function(cond){message("get smiles conversion is a problem")})
    ##print(SM)
    ########################
    ########################
    ### adding the new conditions here 
    if(!sjmisc::is_empty(SM))
    {
      out<-c(out,SM)
    }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["SMILES"]])))){
      SM <- stringr::str_trim(as.character(InMEDA[["SMILES"]]))
      out<-c(out,SM)
    }else{
      SM1<-ifelse(!sjmisc::is_empty(tryCatch({webchem::cs_convert(IN, from = "inchi", to = "smiles")},error=function(cond){message("cs_convert IN to smiles is failing")})),tryCatch({webchem::cs_convert(IN, from = "inchi", to = "smiles")},error=function(cond){message("cs_convert IN to smiles is failing in gETSmiles")}),tryCatch({getCactus(IN,"smiles")},error=function(cond){message("get Cactus function from IN to smile is failing")}))	    
      
      out<-c(out,SM1)
    }
    ################################
    ################################
  }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["InChI"]]))) & startsWith(as.character(InMEDA[["InChI"]]),'CAS:')){
    ###############################	  
    CV<-stringr::str_trim(as.character(InMEDA[["InChI"]]))
    CV1<-stringr::str_replace(CV,pattern='CAS:',replacement ="")
    CV2<-stringr::str_trim(as.character(CV1))
    ##############################
    PCID<-tryCatch({webchem::get_cid(CV2, from = "xref/rn",match="first")},error=function(cond){return(NA)})
    PCID1<-tryCatch({webchem::pc_prop(as.numeric(PCID[[2]][1]), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES","InChI","InChIKey"))},error=function(cond){return(NA)})
    ###############################
    gSMI<-tryCatch({PCID1$CanonicalSMILES},error=function(cond){return(NA)})
    IK<-tryCatch({PCID1$InChIKey},error=function(cond){return(NA)})
    IN<-tryCatch({PCID1$InChI},error=function(cond){return(NA)})
    SM<-tryCatch({getCactus(CV2, "smiles")},error=function(cond){message("CAS to smiles conversion is a problem")})
    ################################
    ################################
    if(!sjmisc::is_empty(gSMI))
    {
      out<-c(out,gSMI)
    }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["SMILES"]])))){
      ## Adding this new
      SM<-tryCatch({PuCAStoOI(CV2)},error=function(cond){message("CAS value is empty")})
      out<-c(out,ifelse(!sjmisc::is_empty(tryCatch({ConvPCIDtoOCN(PC)[2]},error=function(cond){message("Pubchem value is empty")})),tryCatch({ConvPCIDtoOCN(PC)[2]},error=function(cond){message("Pubchem value is empty")}),as.character(InMEDA[["SMILES"]])))
      ###############################
    }else{
      ##out<-c(out,"NA")
      SM<-ifelse(!sjmisc::is_empty(SM),SM,tryCatch({webchem::cir_query(CV2,"smiles")[[2]]},error=function(cond){message("webchem cir query value empty")}))	    
      out<-c(out,SM)
    }
    #################################
  }else if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]]))){
      if(tryCatch({webchem::is.inchikey(as.character(stringr::str_trim(InMEDA[["InChI"]])))},error=function(cond){message("inchikey validation failed")})){	  
      ############################################
      IK<-stringr::str_trim(as.character(InMEDA[["InChI"]]))	    
      tes<-tryCatch({webchem::get_cid(stringr::str_trim(as.character(InMEDA[["InChI"]])), from = "inchikey")},error=function(cond){message("webchem not able to get cid from Inchikey")})
      tes1<-tryCatch({tes$cid},error=function(cond){message("Inchikey to CID did not convert")})
      tes2<-tryCatch({webchem::pc_prop(as.numeric(tes1[1]), properties = c("MolecularFormula", "MolecularWeight","CanonicalSMILES","InChI","InChIKey"))},error=function(cond){message("Inchikey to CID did not convert so did not get properties")})
      IN<-tryCatch({tes2$InChI},error=function(cond){message("Inchi to Inchikey failed because of CID not converting")})
      SM<-tryCatch({tes2$CanonicalSMILES},error=function(cond){message("Inchi to Inchikey failed because of CID not converting")})
      SM1=ifelse(!sjmisc::is_empty(SM),SM,ifelse(!sjmisc::is_empty(tryCatch({PuInKtoSM(InMEDA[["InChI"]])},error=function(cond){message("Inchikey to smile conversion")})),tryCatch({PuInKtoSM(InMEDA[["InChI"]])},error=function(cond){message("Inchikey to smile conversion")}),ifelse(!sjmisc::is_empty(tryCatch({getCactus(InMEDA[["InChI"]], "smiles")},error=function(cond){message("Inchikey to smile conversion")})),tryCatch({getCactus(InMEDA[["InChI"]], "smiles")},error=function(cond){message("Inchikey to smile conversion")}),tryCatch({PuNAMEtoOI(InMEDA[["Name"]])[3]},error=function(cond){message("name to smile conversion")}))))
      out<-c(out,SM1)
    
     }else{

      if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["SMILES"]])))){
        
        SM<-stringr::str_trim(as.character(InMEDA[["SMILES"]]))
        out<-c(out,SM)
        ####################
      }else{
        #############################################
	print("entering the else part..getSMILES..function ..inchikey and smiles failed")
        ######################
        PCID1<-stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))
	SM<-ifelse(!sjmisc::is_empty(PCID1),tryCatch({webchem::cs_convert(as.numeric(PCID1),from="csid",to="smiles")},error=function(cond){message("pubchem CID to smile conversion")}),tryCatch({ConvCIDtoOID1(as.numeric(PCID1))[3]},error=function(cond){message("pubchem CID to smile conversion")}))
        out<-c(out,SM)
        ####################
      }##inside else part
    }###end of else inchikey checking
  }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]])))){
    ####################################################
    PCID1<-stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))
    ####################################################
    tes<-tryCatch({webchem::pc_prop(PCID1, properties = c("MolecularFormula", "MolecularWeight","CanonicalSMILES","InChI","InChIKey"))},error=function(cond){message("Did not get properties from Pubchem CID")})
    gSMI<-tryCatch({tes$CanonicalSMILES},error=function(cond){message("smiles is not found")})
    IK<-tryCatch({tes$InChIKey},error=function(cond){message("some mistake happened in file search files")})
    IN<-tryCatch({tes$InChI},error=function(cond){message("some mistake happened in file search files")})
    ##out<-c(out,gSMI)
    if(!sjmisc::is_empty(gSMI))
    {
      out<-c(out,gSMI)
    }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["SMILES"]])))){
      SM<-stringr::str_trim(as.character(InMEDA[["SMILES"]]))
      out<-c(out,SM)
    }else{
      #########################################	    
      ### This is old code I am adding new code
      ##out<-c(out,"NA")
      ##########################################
      PC<-as.numeric(stringr::str_trim(as.character(InMEDA[["PubChem CID"]])))
      out<-c(out,ifelse(!sjmisc::is_empty(tryCatch({ConvPCIDtoOCN(PC)[2]},error=function(cond){message("Pubchem value is empty")})),tryCatch({ConvPCIDtoOCN(PC)[2]},error=function(cond){message("Pubchem value is empty")}),ifelse(sjmisc::is_empty(tryCatch({ConvCIDtoOID1(PC)[3]},error=function(cond){message("Pubchem CID is empty")})),!sjmisc::is_empty(tryCatch({ConvCIDtoOID1(PC)[3]},error=function(cond){message("Pubchem CID is empty")})),ifelse(!sjmisc::is_empty(PC),webchem::cs_convert(PC,from="csid",to="smiles"),NA))))
      ##########################################
    }
    
    
  }else if(!sjmisc::is_empty(as.character(InMEDA[["SMILES"]]))){

    SM<-stringr::str_trim(as.character(InMEDA[["SMILES"]]))
    out<-c(out,SM)
  }else{
    out<-c(out,"NA")
  }
  ############
  return(out)
  ############
}
#########################################################################
#########################################################################
MaKE.ONT.REC<-function(InMEDA)
{
  out<-c()
  ################
  print("entering the ontology area")
  ################
  if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI="))
  {
    ###################################
    print("enter the line 6 ...........")
    ###################################
    IN<-stringr::str_trim(as.character(InMEDA[["InChI"]]))
    mol <-tryCatch({rinchi::parse.inchi(IN)},error=function(cond){message("parse inchi is empty")})
    SM<-ifelse(!sjmisc::is_empty(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)})),tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}),tryCatch({rcdk::get.smiles(mol[[1]])},error=function(cond){message("get smiles is empty")}))
    IK<-ifelse(!sjmisc::is_empty(tryCatch({rinchi::get.inchi.key(SM)},error=function(cond){return(NA)})),tryCatch({rinchi::get.inchi.key(SM)},error=function(cond){return(NA)}),tryCatch({gsub("InChIKey=","",stringr::str_trim(tryCatch({getCactus(IN, "stdinchikey")},error=function(cond){return(NA)})))},error=function(cond){return(NA)}))
    IK1<-paste("INCHIKEY:",IK,sep=" ")
    ############################
    if(!sjmisc::is_empty(IK)){
      ##########################
      print("enter the line 8")
      ##########################
      IKCRV<-tryCatch({classyfireR::get_classification(IK)},error=function(cond){message("Classifier could not fetch the information")})
      ########################################
      ########################################
      if(!sjmisc::is_empty(IKCRV)){
        #############################
        ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
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
        ### adding this new part###
	if(!sjmisc::is_empty(SM)){
		####################################
	        tes<-tryCatch({classyfireR::submit_query(label = 'query_test', input = SM, type = 'STRUCTURE')},error=function(cond){message("Classyfire is empty")})
      		tes1<-tryCatch({ClassSmilesToOntolgy(SM)},error=function(cond){message("Classyfire smiles to Ontology empty")})
      		tes2<-getOntoSM(SM)
                ###########################################################
      		###########################################################
      		IK1<-paste("INCHIKEY:",IK,sep=" ")
      		##tes2<-paste("Ontology:",tes1,sep=" ")
      		FINCH<-paste("INCHI:",IN,sep=" ")
      		####################################
      		out<-c(out,tes2)
      		out<-c(out,IK1)
      		out<-c(out,FINCH)
		####################################

        }else{

        SM<-tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)})
        IN<-tryCatch({rinchi::get.inchi(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
	IK<-paste("INCHIKEY:",tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)}) ,sep=" ")
	tes1<-tryCatch({ClassSmilesToOntolgy(SM)},error=function(cond){message("Classyfire smiles to Ontology empty")})
	#####################################
	GVOF=getOntFromSMI(IK,InMEDA)
	######################################
        out<-c(out,GVOF[1])
        out<-c(out,GVOF[2])
        out<-c(out,GVOF[3])
        ########################################
	}
      #######################################################
      #######################################################
      }
   #######################################################################
    }else if(!sjmisc::is_empty(SM)){
      ##################################################################
      ##################################################################
      tes<-tryCatch({classyfireR::submit_query(label = 'query_test', input = SM, type = 'STRUCTURE')},error=function(cond){message("Classyfire is empty")})
    #############################################################
    #############################################################
      tes1<-tryCatch({ClassSmilesToOntolgy(SM)},error=function(cond){message("Classyfire smiles to ontology is empty")})
      GVOF=getOntFromSMI(IK,InMEDA)
      ###########################################################
      ###########################################################
      out<-c(out,GVOF[1])
      out<-c(out,GVOF[2])
      out<-c(out,GVOF[3])
      ############################################################
    }else{
  ################################################################
      F1ONT<-paste("Ontology:","",sep=" ")
      #############################################
      ##############################################
      if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){
        ############################################
        FINCH<-paste("INCHI:",stringr::str_trim(as.character(InMEDA[["InChI"]])),sep=" ")
        IK<-paste("INCHIKEY:",tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)}) ,sep=" ")
        GVOF=getOntFromSMI(IK,InMEDA)
	##################
        out<-c(out,GVOF[1])
        out<-c(out,GVOF[2])
        out<-c(out,GVOF[3])
        #################
      }else{
        ###################################
	FINCH<-tryCatch({rinchi::get.inchi(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})      
        IK<-paste("INCHIKEY:",stringr::str_trim(as.character(InMEDA[["InChI"]])),sep=" ")
	GVOF=getOntFromSMI(IK,InMEDA)
        ####################
        out<-c(out,GVOF[1])
        out<-c(out,GVOF[2])
        out<-c(out,GVOF[3])
        ####################
      }
      ######################
    }## end of else..else if ...if
  }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["InChI"]])))){
 ########################################################################
    ############################
    if(tryCatch({webchem::is.inchikey(stringr::str_trim(as.character(InMEDA[["InChI"]])))},error=function(cond){message("inchikey validation failed")}))
    {
      #########################################
      print("entering the inchikey area")
      ########################################
      tes<-tryCatch({webchem::get_cid(stringr::str_trim(as.character(InMEDA[["InChI"]])), from = "inchikey")},error=function(cond){message("webchem not able to get cid from Inchikey")})
      tes1<-tryCatch({tes$cid},error=function(cond){message("Inchikey to CID did not convert")})
      #####################################################################
      ######################################################################
      tes2<-tryCatch({webchem::pc_prop(as.numeric(stringr::str_trim(gsub("[[:punct:]]", "",tes1[1]))), properties = c("MolecularFormula", "MolecularWeight","CanonicalSMILES","InChI","InChIKey"))},error=function(cond){message("Inchikey to CID did not convert so did not get properties")})
      SM<-ifelse(!sjmisc::is_empty(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)})),tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}),tryCatch({tes2$CanonicalSMILES},error=function(cond){message("Inchi to Inchikey failed because of CID not converting")}))
      IN<-ifelse(!sjmisc::is_empty(tryCatch({tes2$InChI},error=function(cond){message("Inchi to Inchikey failed because of CID not converting")})),tryCatch({tes2$InChI},error=function(cond){message("Inchi to Inchikey failed because of CID not converting")}),tryCatch({rinchi::get.inchi(SM)},error=function(cond){message("smiles to Inchikey conversion ")}))
      IK<-ifelse(!sjmisc::is_empty(tryCatch({stringr::str_trim(as.character(InMEDA[["InChI"]]))},error=function(cond){message("Inchi to Inchikey failed because of CID not converting")})),tryCatch({stringr::str_trim(as.character(InMEDA[["InChI"]]))},error=function(cond){message("Inchi to Inchikey failed because of CID not converting")}),tryCatch({rinchi::get.inchi.key(SM)},error=function(cond){message("smiles to Inchikey conversion ")}))
      #######################################
      #######################################
      FINCH<-paste("INCHI:",IN,sep=" ")
      IK1<-paste("INCHIKEY:",IK,sep=" ")
      ############################
      if(!sjmisc::is_empty(IK)){
        ##########################
        print("enter the line 8")
        #############################
        IKCRV<-tryCatch({classyfireR::get_classification(IK)},error=function(cond){message("Classifier could not fetch the information")})
        #####################################
        #####################################
        if(!sjmisc::is_empty(IKCRV)){

          ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
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
          FINCH<-paste("INCHI:",IN,sep=" ")
	  GVOF=getOntFromSMI(IK,InMEDA)
          ################################
	  out<-c(out,GVOF[1])
	  out<-c(out,GVOF[2])
	  out<-c(out,GVOF[3])
          ################################
        }## end of else
      }else if(!sjmisc::is_empty(SM)){
        #############################################
        tes<-tryCatch({classyfireR::submit_query(label = 'query_test', input = SM, type = 'STRUCTURE')},error=function(cond){message("Classyfire is empty")})
        ##################################################
        IK<-tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})  
        #########################################################
        GVOF=getOntFromSMI(IK,InMEDA)
	################################
        out<-c(out,GVOF[1])
        out<-c(out,GVOF[2])
        out<-c(out,GVOF[3])
	#################################

      }else{
        F1ONT<-paste("Ontology:","",sep=" ")
        ######################################################
        if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){
          ####################################################
          FINCH<-paste("INCHI:",stringr::str_trim(as.character(InMEDA[["InChI"]])),sep=" ")
	  IK<-tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
          GVOF=getOntFromSMI(IK,InMEDA)
	  ##############################
	  out<-c(out,GVOF[1])
          out<-c(out,GVOF[2])
          out<-c(out,GVOF[3])
          ##############################
        }else{
          ###################################
	  FINCH<-tryCatch({rinchi::get.inchi(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})	
          IK<-stringr::str_trim(as.character(InMEDA[["InChI"]]))
          GVOF=getOntFromSMI(IK,InMEDA)
          ##########################
	  out<-c(out,GVOF[1])
          out<-c(out,GVOF[2])
          out<-c(out,GVOF[3])
          ##########################
        }

      }## end of else ...else if ..else
    }else if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){

	    print("enter the else if loop ...where it has the INCHI information")

	    FINCH<-paste("INCHI:",stringr::str_trim(as.character(InMEDA[["InChI"]])),sep=" ")
	    IK<-tryCatch({webchem::cs_convert(as.character(InMEDA[["InChI"]]),from = "inchi", to = "inchikey")},error=function(cond){message("webchecm could not fetch the info inchi")})
  	    SM<-tryCatch({webchem::cs_convert(as.character(InMEDA[["InChI"]]),from = "inchi", to = "smiles")},error=function(cond){message("webchecm could not fetch the info inchi")})
  	    tes1<-tryCatch({ClassSmilesToOntolgy(SM)},error=function(cond){message("smiles to ontology is failing")})
 	    GVOF=getOntFromSMI(IK,InMEDA)
	    ##########################
	    out<-c(out,GVOF[1])
            out<-c(out,GVOF[2])
            out<-c(out,GVOF[3])
	    ##########################
    
    }else{

      print("checking if entering the else loop in the inch or inchikey function")
      ##inchikey validation failed ..so must be CAS or try to get ontology from smiles ...
      ########################################################################################
      if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["InChI"]]))) & startsWith(as.character(InMEDA[["InChI"]]),'CAS:')){
        ###############################
        print("enter the CAS area ONtology")
        ##############################
        CV<-stringr::str_trim(as.character(InMEDA[["InChI"]]))
        CV1<-stringr::str_replace(CV,pattern='CAS:',replacement ="")
        CV2<-stringr::str_trim(as.character(CV1))
        ###########################################
        PCID<-tryCatch({webchem::get_cid(CV2, from = "xref/rn",match="first")},error=function(cond){message("Pubchem Id is empty")})
	##################################################################################
	##################################################################################
        PCID1<-tryCatch({webchem::pc_prop(as.numeric(stringr::str_trim(gsub("[[:punct:]]", "",PCID[[2]][1]))), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES","InChI","InChIKey"))},error=function(cond){message("Pubchem Id is empty")})
        ###############################
        gSMI<-tryCatch({PCID1$CanonicalSMILES},error=function(cond){message("smiles is not found")})
	IK<-ifelse(!sjmisc::is_empty(tryCatch({PCID1$InChIKey},error=function(cond){message("CID to inchikey conversion failed")})),tryCatch({PCID1$InChIKey},error=function(cond){message("CID to inchikey conversion failed")}),tryCatch({getCactus(CV2,"stdinchikey")},error=function(cond){message("CID to inchikey conversion is failed")}))
	IN<-ifelse(!sjmisc::is_empty(tryCatch({PCID1$InChI},error=function(cond){message("CID to inchi conversion failed")})),tryCatch({PCID1$InChI},error=function(cond){message("CID to inchi conversion failed")}),tryCatch({getCactus(CV2,"stdinchi")},error=function(cond){message("CID to inchikey conversion failed")}))
        IK1<-paste("INCHIKEY:",IK,sep=" ")
        #############################
        if(!sjmisc::is_empty(IK)){
          ##########################
          IKCRV<-tryCatch({classyfireR::get_classification(IK)},error=function(cond){message("Classifier could not fetch the information")})
          ######################################
          #####################################
          if(!sjmisc::is_empty(IKCRV)){

            ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
            ##########################
	    IN<-tryCatch({rinchi::get.inchi(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
	    IK<-tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
	    ###############################
            IK1<-paste("INCHIKEY:",IK,sep=" ")
            tes2<-paste("Ontology:",ONTV,sep=" ")
            FINCH<-paste("INCHI:",IN,sep=" ")
            ###########################
            out<-c(out,tes2)
            out<-c(out,IK1)
            out<-c(out,FINCH)
            ########################
          }## ikcrv END
        }else{

          if(!sjmisc::is_empty(as.character(InMEDA[["SMILES"]])) & !startsWith(as.character(InMEDA[["SMILES"]]),'not available')){
            ####################################################################
	    F1SM<-stringr::str_trim(as.character(InMEDA[["SMILES"]]))
            F1SM1<-tryCatch({rinchi::get.inchi.key(F1SM)},error=function(cond){message("webchecm could not fetch the info")})
            ##################################
            if(!sjmisc::is_empty(F1SM1)){
              ################################
              IK<-F1SM1
              IKCRV<-tryCatch({classyfireR::get_classification(IK)},error=function(cond){message("Classifier could not fetch the information")})
              ############################################
              ############################################
              if(!sjmisc::is_empty(IKCRV)){
                ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
                ##########################
                IK1<-paste("INCHIKEY:",IK,sep=" ")
                tes2<-paste("Ontology:",ONTV,sep=" ")
                FINCH<-paste("INCHI:",IN,sep=" ")
                #########################
                out<-c(out,tes2)
                out<-c(out,IK1)
                out<-c(out,FINCH)
                ##############################
              }##IKCRV
            }else{

              ## not able to get inchikey from smiles too check pubchemID
              if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(InMEDA[["PubChem CID"]]),'not available')){
                ##############################################
                ##FPUCID<-stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))
                FPUCID<-stringr::str_trim(gsub("[[:punct:]]", "",as.character(InMEDA[["PubChem CID"]])))
                FPUCID1<-as.numeric(FPUCID)
                FINSM<-tryCatch({webchem::pc_prop(FPUCID1)},error=function(cond){message("webchecm could not fetch the information")})
                FIINK<-tryCatch({FINSM$InChIKey},error=function(cond){message("webchecm could not fetch the info inchikey")})
		SM<-tryCatch({FINSM$CanonicalSMILES},error=function(cond){message("webchecm could not fetch the info for smile")})
		IN<-tryCatch({FINSM$InChI},error=function(cond){message("webchecm could not fetch the info inchi")})
                #############################
                IK<-FIINK
                #############################
                if(!sjmisc::is_empty(IK)){
                  ##########################
                  IKCRV<-tryCatch({classyfireR::get_classification(IK)},error=function(cond){message("Classifier could not fetch the information")})
                  #####################################
                  ######################################
                  if(!sjmisc::is_empty(IKCRV)){

                    ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
                    ##########################
		    IN<-tryCatch({rinchi::get.inchi(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
		    ###############################
                    IK1<-paste("INCHIKEY:",IK,sep=" ")
                    tes2<-paste("Ontology:",ONTV,sep=" ")
                    FINCH<-paste("INCHI:",IN,sep=" ")
                    ########################
                    out<-c(out,tes2)
                    out<-c(out,IK1)
                    out<-c(out,FINCH)
                    ########################
                  }else{
			  ############################
			  ###tes1<-tryCatch({ClassSmilesToOntolgy(SM)},error=function(cond){message("smiles to ontology is failing")})
			  #############################
			  IN<-tryCatch({rinchi::get.inchi(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
			  IK<-tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
			  ##############################
			  GVOF=getOntFromSMI(IK,InMEDA)
			  ###############################
			  out<-c(out,GVOF[1])
			  out<-c(out,GVOF[2])
			  out<-c(out,GVOF[3])
			  ###############################
		  }## added the else loop
	      ######################################################################	  
                }else{

			print("enter the else loop approximately around 1594 line")
                  #############################################################################
                  if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){
                  #############################################################################
		    if(!tryCatch({webchem::is.inchikey(as.character(InMEDA[["InChI"]]))},error=function(cond){message("inchikey validation failed")})){
		    
	            ############################	    
                    FINCH<-paste("INCHI:",as.character(InMEDA[["InChI"]]),sep=" ")
                    IK1<-paste("INCHIKEY:",tryCatch({webchem::cs_convert(as.character(InMEDA[["InChI"]]),from = "inchi", to = "inchikey")},error=function(cond){message("webchecm could not fetch the info inchi")}),sep=" ")
                    IK<-tryCatch({webchem::cs_convert(as.character(InMEDA[["InChI"]]),from = "inchi", to = "inchikey")},error=function(cond){message("webchecm could not fetch the info inchi")})
		    SM<-tryCatch({webchem::cs_convert(as.character(InMEDA[["InChI"]]),from = "inchi", to = "smiles")},error=function(cond){message("webchecm could not fetch the info inchi")})
		    ############################
		    GVOF=getOntFromSMI(IK,InMEDA)
		    ############################
		    out<-c(out,GVOF[1])
		    out<-c(out,GVOF[2])
		    out<-c(out,GVOF[3])
                    ############################

                  }else{

                    ###################################
		    FINCH<-tryCatch({rinchi::get.inchi(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})	  
		    IK=stringr::str_trim(as.character(InMEDA[["InChI"]])) 
		    IK1<-paste("INCHIKEY:",stringr::str_trim(as.character(InMEDA[["InChI"]])),sep=" ")
		    SM<-stringr::str_trim(as.character(InMEDA[["SMILES"]]))
                    ####################################
		    GVOF=getOntFromSMI(IK,InMEDA)
		    ####################################
		    out<-c(out,GVOF[1])
                    out<-c(out,GVOF[2])
                    out<-c(out,GVOF[3])
                    ####################################
                    }
                ######################
                   }
               ######################
                }## end of else
            ########################################

              }else{

                ##############################################
                if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["InChI"]]))) & startsWith(stringr::str_trim(as.character(InMEDA[["InChI"]])),"InChI=")){
		  if(!tryCatch({webchem::is.inchikey(stringr::str_trim(as.character(InMEDA[["InChI"]])))},error=function(cond){message("inchikey validation failed")})){
                  #####################################
                  FINCH<-paste("INCHI:",stringr::str_trim(as.character(InMEDA[["InChI"]])),sep=" ")
		  IK=tryCatch({webchem::cs_convert(as.character(InMEDA[["InChI"]]),from = "inchi", to = "inchikey")},error=function(cond){message("webchecm could not fetch the info inchi")})
		  IK1<-paste("INCHIKEY:",tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)}) ,sep=" ")
                  ###################################
		  GVOF=getOntFromSMI(IK,InMEDA)
		  ####################################
                  out<-c(out,GVOF[1])
                  out<-c(out,GVOF[2])
                  out<-c(out,GVOF[3])
		  ###################################
                }else{
                  ###################################
		  FINCH<-paste("INCHI:",tryCatch({rinchi::get.inchi(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)}),sep=" ")	
                  IK<-paste("INCHIKEY:",as.character(InMEDA[["InChI"]]),sep=" ")
                  ###################################
		  GVOF=getOntFromSMI(IK,InMEDA)
		  ###############################
		  out<-c(out,GVOF[1])
                  out<-c(out,GVOF[2])
                  out<-c(out,GVOF[3])
                  ###################################
                }### else loop
	      ###########################	
		}## end of if loop
	####################################	
              }## end of else## CID
        ########################################
            }## check else..smiles
        #############################################   
          }else{
            #############################################
		  print("enter the line areound 178 lines")
            F1ONT<-paste("Ontology:","",sep=" ")
            ##############################################
            if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["InChI"]]))) & startsWith(stringr::str_trim(as.character(InMEDA[["InChI"]])) ,"InChI=")){
		    if(!tryCatch({webchem::is.inchikey(stringr::str_trim(as.character(InMEDA[["InChI"]])))},error=function(cond){message("inchikey validation failed")})){
              ##################################
              FINCH<-paste("INCHI:",stringr::str_trim(as.character(InMEDA[["InChI"]])),sep=" ")
	      IK1<-paste("INCHIKEY:",tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)}) ,sep=" ")
              IK=tryCatch({webchem::cs_convert(as.character(InMEDA[["InChI"]]),from = "inchi", to = "inchikey")},error=function(cond){message("webchecm could not fetch the info inchi")})
	      ##################################
	      GVOF=getOntFromSMI(IK,InMEDA)
	      ######################################
              out<-c(out,GVOF[1])
              out<-c(out,GVOF[2])
              out<-c(out,GVOF[3])
	      ##################################
            }else{
              ###################################
	      IK<-tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
              GVOF=getOntFromSMI(IK,InMEDA)
              ######################
	      out<-c(out,GVOF[1])
	      out<-c(out,GVOF[2])
              out<-c(out,GVOF[3])
              ######################
            }## else loop
         ##########################		    
           }## end of if loop 
        ##########################
          }## inner smiles ..end ..else
       #############################	
        }## end of smiles## CAS
     ################################################
      }else{
        #############################################
	print("enter the if condition")	
        ##############################################
        if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){
		if(!tryCatch({webchem::is.inchikey(as.character(InMEDA[["InChI"]]))},error=function(cond){message("inchikey validation failed")})){
         ##########################################
            print("enter the if loop around 1722")		       	
	 ######################################
          FINCH<-paste("INCHI:",stringr::str_trim(as.character(InMEDA[["InChI"]])),sep=" ")
	  IK<-tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
	  #################################
	  GVOF=getOntFromSMI(IK,InMEDA)
	  #################################
	  out<-c(out,GVOF[1])
          out<-c(out,GVOF[2])
          out<-c(out,GVOF[3])
	  #################################
          #################################
        }else{
		print("enter the else loop around 1722")
          ###################################
	  FINCH<-tryCatch({rinchi::get.inchi(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)}) 	
          IK=tryCatch({webchem::cs_convert(as.character(InMEDA[["InChI"]]),from = "inchi", to = "inchikey")},error=function(cond){message("webchecm could not fetch the info inchi")})
	  ################################
          GVOF=getOntFromSMI(IK,InMEDA)
	  ################################
	  out<-c(out,GVOF[1])
          out<-c(out,GVOF[2])
          out<-c(out,GVOF[3])
	  ################################
        } ## else loop
     #####################################		
       }else{
	       ########################################
	       print("enter the checking if it is entering this vallue")
	       SM<-tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)})
	       IK<-tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
               ############################################
	       GVOF=getOntFromSMI(IK,InMEDA)
	       ###########################################
	       out<-c(out,GVOF[1])
               out<-c(out,GVOF[2])
               out<-c(out,GVOF[3])
	       ###########################################
	} 
      }#######
      ##################################
    }### end of else ..so starting checking CAS..smiles ..so ..on
  }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["InChI"]]))) & startsWith(as.character(InMEDA[["InChI"]]),'CAS:')){
  #################################################
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
    IK<-ifelse(!sjmisc::is_empty(tryCatch({PCID1$InChIKey},error=function(cond){message("CID to inchikey conversion failed")})),tryCatch({PCID1$InChIKey},error=function(cond){message("CID to inchikey conversion failed")}),tryCatch({getCactus(CV2,"stdinchikey")},error=function(cond){message("CID to inchikey conversion failed")}))
    IN<-ifelse(!sjmisc::is_empty(tryCatch({PCID1$InChI},error=function(cond){message("CID to inchi conversion failed")})),tryCatch({PCID1$InChI},error=function(cond){message("CID to inchi conversion failed")}),tryCatch({getCactus(CV2,"stdinchi")},error=function(cond){message("CID to inchikey conversion failed")}))
    IK1<-paste("INCHIKEY:",IK,sep=" ")
    #############################
    if(!sjmisc::is_empty(IK)){
      ##########################
      IKCRV<-tryCatch({classyfireR::get_classification(IK)},error=function(cond){message("Classifier could not fetch the information")})
      ######################################
      if(!sjmisc::is_empty(IKCRV)){

        ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
        ##########################
        IK1<-paste("INCHIKEY:",IK,sep=" ")
        tes2<-paste("Ontology:",ONTV,sep=" ")
        FINCH<-paste("INCHI:",IN,sep=" ")
        ##########################
        out<-c(out,tes2)
        out<-c(out,IK1)
        out<-c(out,FINCH)
        ###########################
      }## ikcrv END
    }else{
	    print("enter the line around 1799")
      if(!sjmisc::is_empty(as.character(InMEDA[["SMILES"]])) & !startsWith(as.character(InMEDA[["SMILES"]]),'not available')){
        #############################################
        SM<-stringr::str_trim(as.character(InMEDA[["SMILES"]]))
        F1SM<-stringr::str_trim(as.character(InMEDA[["SMILES"]]))
        F1SM1<-tryCatch({rinchi::get.inchi.key(F1SM)},error=function(cond){message("webchecm could not fetch the info")})
        ##############################
        if(!sjmisc::is_empty(F1SM1)){
          ############################
          IK<-F1SM1
          IKCRV<-tryCatch({classyfireR::get_classification(IK)},error=function(cond){message("Classifier could not fetch the information")})
          ###################################
          ############################################
          if(!sjmisc::is_empty(IKCRV)){
            ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
            ##########################
            IK1<-paste("INCHIKEY:",IK,sep=" ")
            tes2<-paste("Ontology:",ONTV,sep=" ")
            FINCH<-paste("INCHI:",IN,sep=" ")
            #####################
            out<-c(out,tes2)
            out<-c(out,IK1)
            out<-c(out,FINCH)
            ####################
          }##IKCRV
        }else{
		
		print("enter the line around 1826")
          if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(InMEDA[["PubChem CID"]]),'not available')){
            ############################################
            FPUCID<-stringr::str_trim(gsub("[[:punct:]]", "",as.character(InMEDA[["PubChem CID"]])))
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
              #######################################
              if(!sjmisc::is_empty(IKCRV)){

                ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
                ##########################
                IK1<-paste("INCHIKEY:",IK,sep=" ")
                tes2<-paste("Ontology:",ONTV,sep=" ")
                FINCH<-paste("INCHI:",IN,sep=" ")
                ##########################
                out<-c(out,tes2)
                out<-c(out,IK1)
                out<-c(out,FINCH)
                ##########################
              }## ikcrv END
            }else{

	      print("enter the line around 1858")
              F1ONT<-paste("Ontology:","",sep=" ")
              ##############################################
              if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){
		      if(!tryCatch({webchem::is.inchikey(as.character(InMEDA[["InChI"]]))},error=function(cond){message("inchikey validation failed")})){
                ###############################
                FINCH<-paste("INCHI:",stringr::str_trim(as.character(InMEDA[["InChI"]])),sep=" ")
		IK=tryCatch({webchem::cs_convert(as.character(InMEDA[["InChI"]]),from = "inchi", to = "inchikey")},error=function(cond){message("webchecm could not fetch the info inchi")})
		#######################################
		GVOF=getOntFromSMI(IK,InMEDA)
		######################################
		out<-c(out,GVOF[1])
		out<-c(out,GVOF[2])
		out<-c(out,GVOF[3])
		#######################################
              }else{
                ###################################
		FINCH<-tryCatch({rinchi::get.inchi(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})      
	        IK<-tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
		#######################################
                GVOF=getOntFromSMI(IK,InMEDA)
		#######################################
		out<-c(out,GVOF[1])
                out<-c(out,GVOF[2])
                out<-c(out,GVOF[3])
                #######################################
              }
           #########################
            }## end of else
           #########################
            }
	  ###########################  
          }else{
		  print("enter the line around 1898")
            ##############################################
            F1ONT<-paste("Ontology:","",sep=" ")
            ##############################################
            if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){
		    if(!tryCatch({webchem::is.inchikey(as.character(InMEDA[["InChI"]]))},error=function(cond){message("inchikey validation failed")})){
              #################################################
	      IK<-tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
	      ##########################
	      GVOF=getOntFromSMI(IK,InMEDA) 
              ##################
	      out<-c(out,GVOF[1])
              out<-c(out,GVOF[2])
              out<-c(out,GVOF[3])
              #################
            }else{
              ###################################
	      FINCH<-tryCatch({rinchi::get.inchi(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})	    
              IK<-tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
              #######################
	      GVOF=getOntFromSMI(IK,InMEDA)
	      ######################################
	      out<-c(out,GVOF[1])
	      out<-c(out,GVOF[2])
	      out<-c(out,GVOF[3])
	      #######################
              ####################
            }
         ################		    
          }## end of else## CID
        ##################
         }
       #######################	  
        }## check else..smiles
      }else{
        #############################################
        ##############################################
        if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){
		if(!tryCatch({webchem::is.inchikey(as.character(InMEDA[["InChI"]]))},error=function(cond){message("inchikey validation failed")})){
          ##########################################
          FINCH<-paste("INCHI:",stringr::str_trim(as.character(InMEDA[["InChI"]])),sep=" ")
	  IK<-paste("INCHIKEY:",tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)}) ,sep=" ")
	  GVOF=getOntFromSMI(IK,InMEDA)
	  ############################
	  out<-c(out,GVOF[1])
          out<-c(out,GVOF[2])
          out<-c(out,GVOF[3])
	  ##############################
        }else{
          ###################################
	  FINCH<-tryCatch({rinchi::get.inchi(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})	
          IK=stringr::str_trim(as.character(InMEDA[["InChI"]]))
	  #################################
	  GVOF=getOntFromSMI(IK,InMEDA)
	  ####################################
	  out<-c(out,GVOF[1])
          out<-c(out,GVOF[2])
          out<-c(out,GVOF[3])
       #######################################
        }
       }
      }## inner smiles ..end ..else
    }## end of smiles## CAS

  }else{
    ##########################################
    ###################################

    PCID<-stringr::str_trim(gsub("[[:punct:]]", "",as.character(InMEDA[["PubChem CID"]])))
    PCSM<- stringr::str_trim(as.character(InMEDA[["SMILES"]]))
    
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
        ####################################
        if(!sjmisc::is_empty(IKCRV)){
          ####################
          ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
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
	  SM<-tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)})
          IK<-tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
	  ####################################
	  GVOF=getOntFromSMI(IK,InMEDA)
	  ######################################
	  out<-c(out,GVOF[1])
	  out<-c(out,GVOF[2])
          out<-c(out,GVOF[3])
          #################################
        }###else
      }else if(!sjmisc::is_empty(gSMI)){
        ##################################################
        print("enter the line 184")
        ##################################################
        tes<-tryCatch({classyfireR::submit_query(label = 'query_test', input = gSMI, type = 'STRUCTURE')},error=function(cond){message("adduct value is missing")})
        IK<-tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
	################################
	GVOF=getOntFromSMI(IK,InMEDA)
        ######################################
        out<-c(out,GVOF[1])
        out<-c(out,GVOF[2])
        out<-c(out,GVOF[3])
	#####################################
        #################
      }else{

        ###############################################################
	###F1ONT<-tryCatch({MaKE.ONT.REC(InMEDA)[1]},error=function(cond){message("Classifier could not fecth the information")})      
        F1ONT<-paste("Ontology:","",sep=" ")
        if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){
          #########################################
          FINCH<-paste("INCHI:",stringr::str_trim(as.character(InMEDA[["InChI"]])),sep=" ")
	  IK<-tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
	  ################################
          GVOF=getOntFromSMI(IK,InMEDA)
          ######################################
	  out<-c(out,GVOF[1])
          out<-c(out,GVOF[2])
          out<-c(out,GVOF[3])
	  #######################################
          ################
        }else{
          ##################################
          IK<-tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})		
          GVOF=getOntFromSMI(IK,InMEDA)
	  out<-c(out,GVOF[1])
          out<-c(out,GVOF[2])
          out<-c(out,GVOF[3])
	  
          ####################
        }##else..if
        ###############################################################
      }###else...else if ..if
    }else{
      ######check the smiles is not empty...PCID is empty
      ###PCSM<-SMV

      print("enter the else loop ...PCID")	    
      SMV<-PCSM
      ###############################
      if(!sjmisc::is_empty(SMV)){
        #############################
        IK<-tryCatch({rinchi::get.inchi.key(SMV)},error=function(cond){message("rinchi could not fetch inchikey missing")})
        SMV1<-tryCatch({rinchi::get.inchi(SMV)},error=function(cond){message("rinchi could not fetch inchi missing")})
        IK1<-paste("INCHIKEY:",IK,sep=" ")
        FINCH<-paste("INCHI:",SMV1,sep=" ")
        #########################
        if(!sjmisc::is_empty(IK)){
          #############################
          IKCRV<-tryCatch({classyfireR::get_classification(IK)},error=function(cond){message("Classifier could not fetch the information")})
          ##################################################
          if(!sjmisc::is_empty(IKCRV)){
            ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
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
            IK<-tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
            ########################
	    GVOF=getOntFromSMI(IK,InMEDA)
	    ##########################
	    out<-c(out,GVOF[1])
            out<-c(out,GVOF[2])
            out<-c(out,GVOF[3])
            ########################
          }## end of IKCRV

        }else if(!sjmisc::is_empty(SMV)){
          ##########################################
	  IK<-tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
          #######################################
          GVOF=getOntFromSMI(IK,InMEDA)
          ######################################
          out<-c(out,GVOF[1])
          out<-c(out,GVOF[2])
          out<-c(out,GVOF[3])
          ######################################	  
          ##########################
        }else{
          ##############################################
          F1ONT<-paste("Ontology:","",sep=" ")
          if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),"InChI=")){
            IK=tryCatch({webchem::cs_convert(as.character(InMEDA[["InChI"]]),from = "inchi", to = "inchikey")},error=function(cond){message("webchecm could not fetch the info inchi")})
            ###########################
	    GVOF=getOntFromSMI(IK,InMEDA)
	    ######################################
            out<-c(out,GVOF[1])
            out<-c(out,GVOF[2])
            out<-c(out,GVOF[3])
            ###################
          }else{
            ###FINCH<-paste("INCHI:","",sep=" ")
	    FINCH<-tryCatch({rinchi::get.inchi(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})	  
            IK=stringr::str_trim(as.character(InMEDA[["InChI"]]))
	    ###############################
	    GVOF=getOntFromSMI(IK,InMEDA)
	    ######################################
            out<-c(out,GVOF[1])
            out<-c(out,GVOF[2])
            out<-c(out,GVOF[3])
            ##################################
          }
        #################
        }## end of else if ...else ...if
      #########
      }
    ############
    }## entering the else
  #################
  }## end of else
  ##########################
  return(out)
#############################
}## end of ontology function

#########################################################################
#########################################################################
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
############################################################################
############################################################################

comwithNu <- function(elem)
{
  ##############################
  tes1 <- c()
  ##############################
  if (grepl("^[[:digit:]]+", elem) &
      !grepl("^[0-9]*[.])?[0-9]+", elem)) {
    elem1 <- stringi::stri_extract_first(elem, regex = "[0-9]+")
    elem2 <- substring(elem, 2)
    if (numbers_only(elem1) &!gdata::startsWith(elem2, pattern = "M") &!gdata::startsWith(elem2, pattern = "*")){
      nelem1 = paste(elem1,"*",enviGCMS::getmass(elem2),sep="")
      tes1 <- c(tes1,nelem1)
    }else if(numbers_only(elem1) & gdata::startsWith(elem2, pattern = "M") & !gdata::startsWith(elem2, pattern = "*"))
    {
      nelem1 = paste(elem1,"*",elem2,sep="")
      tes1 <- c(tes1,nelem1)
    }else{
      
      tes1 <- c(tes1,elem)
    }
    
  }else{
    tes1 <- c(tes1,elem)
  }
  #######################
  return(paste(tes1, collapse = ""))
  #######################
}


################################################################################
################################################################################

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
        ########################
        ##print("entering this part")
        ##print(tes[i])
        ###print(VL)
        #######################      
        VL<-comwithNu(tes[i])
        tes1<-c(tes1,VL)
        #######################
        
      }else{
        ##########################################      
        ##########################################      
        if(tes[i]=="ACN")
        {
          ###tes1<-c(tes1,enviGCMS::getmass(tes[i]))
          tes1<-c(tes1,41.02599796)
        }else{
          tes1<-c(tes1,enviGCMS::getmass(tes[i]))
          
        }###end of else
      ############################################
      ############################################
      }## end if the else
      
      
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
    addu2<-tryCatch({qdapRegex::ex_between(addu1, "[", "]")[[1]]},error=function(cond){message("error happened in precursortype extraction")})
    addu3<-tryCatch({AIN[AIN$V1==addu2 & AIN$V3==cha,]$V2},error=function(cond){message("adduct match is empty")})
    ##############################################################################
    ########################################################
    if(!sjmisc::is_empty(addu3))
    ##########################################################	    
      ###if(!sjmisc::is_empty(str_replace_all(addu3, "\"", "")))
    #########################################################	    
    {
      ##print(addu3)
      faddu<-c(faddu,addu3)
      
    }else{
      NAF<-unlist(strsplit(sub(".*\\[([^][]+)].*", "\\1",addu2),"(?=[+-])", perl = TRUE))
      ####################################
      #####################################
      NAF1<-tryCatch({comNFMDA(NAF)},error=function(cond){message("comNFMDA(NAF) is failing")})
      ##########################
      ###########################
      faddu<-c(faddu,NAF1)
      ######################################
      
    }
    
  }else if(cha=="+"){
    cha<-"1+"
    addu2<-tryCatch({qdapRegex::ex_between(addu1, "[", "]")[[1]]},error=function(cond){message("error happened in precursortype extraction")})
    addu3<-tryCatch({AIN[AIN$V1==addu2 & AIN$V3==cha,]$V2},error=function(cond){message("adduct match is empty")})
    ################################################################
    #################################################################
    if(!sjmisc::is_empty(addu3))
    {
      faddu<-c(faddu,addu3)
      
    }else{
      NAF<-unlist(strsplit(sub(".*\\[([^][]+)].*", "\\1",addu2),"(?=[+-])", perl = TRUE))
      #############################################################
      ############################################################
      NAF1<-tryCatch({comNFMDA(NAF)},error=function(cond){message("comNFMDA(NAF) is failing...2")})
      faddu<-c(faddu,NAF1)
      ################################################
      ################################################
    }
    
  }else if(cha=="2-"){
    
    ##x<-"pass1"
    dIV<-"2"
    addu2<-tryCatch({qdapRegex::ex_between(addu1, "[", "]")[[1]]},error=function(cond){message("error happened in precursortype extraction")})
    addu3<-tryCatch({AIN[AIN$V1==addu2 & AIN$V3==cha,]$V2},error=function(cond){message("adduct match is empty")})
    #################################################################
    ############################################################
    if(!sjmisc::is_empty(addu3))
    ###############################################################	    
    ######################################################	    
    {
      ##print(addu1)
      faddu<-c(faddu,addu3)
      
    }else{
      NAF<-unlist(strsplit(sub(".*\\[([^][]+)].*", "\\1",addu2),"(?=[+-])", perl = TRUE))
      ########################################################
      ########################################################
      NAF1<-tryCatch({comNFMDA(NAF)},error=function(cond){message("comNFMDA(NAF) is failing...3")})
      NAF2<-tryCatch({finADWC(NAF1,dIV)},error=function(cond){message("finADWC(NAF1,dIV) is failing...3")})
      ########################################################
      #########################################################
      faddu<-c(faddu,NAF2)
      ########################################################
      #########################################################
    }
    
  }else if(cha=="2+"){
    dIV<-"2"
    addu2<-tryCatch({qdapRegex::ex_between(addu1, "[", "]")[[1]]},warning=function(cond){message("error happened in precursortype extraction")})
    addu3<-tryCatch({AIN[AIN$V1==addu2 & AIN$V3==cha,]$V2},warning=function(cond){message("adduct match is empty")})
    ################################################################
    ################################################################
    if(!sjmisc::is_empty(addu3))
    {
      ##print(addu1)
      faddu<-c(faddu,addu3)
      
    }else{
      NAF<-unlist(strsplit(sub(".*\\[([^][]+)].*", "\\1",addu2),"(?=[+-])", perl = TRUE))
      ################################################################
      ###############################################################
      NAF1<-tryCatch({comNFMDA(NAF)},error=function(cond){message("comNFMDA(NAF) is failing...4")})
      NAF2<-tryCatch({finADWC(NAF1,dIV)},error=function(cond){message("finADWC(NAF1,dIV) is failing...4")})
      ####################################################################
      ######################################################################
      faddu<-c(faddu,NAF2)
      ####################################################################
    }
    
  }else if(cha=="3-"){
    ##x<-"pass1"
    dIV<-"3"
    addu2<-tryCatch({qdapRegex::ex_between(addu1, "[", "]")[[1]]},warning=function(cond){message("error happened in precursortype extraction")})
    addu3<-tryCatch({AIN[AIN$V1==addu2 & AIN$V3==cha,]$V2},warning=function(cond){message("adduct match is empty")})
    ######################################################################
    ##########################################################################
    if(!sjmisc::is_empty(addu3))
    {
      ##print(addu1)
      faddu<-c(faddu,addu3)
      
    }else{
      ##NAF<-unlist(strsplit(addu2, "(?=[+-])", perl = TRUE))
      NAF<-unlist(strsplit(sub(".*\\[([^][]+)].*", "\\1",addu2),"(?=[+-])", perl = TRUE))
      NAF1<-tryCatch({comNFMDA(NAF)},error=function(cond){message("comNFMDA(NAF) is failing...5")})
      ###################################################################
      #################################################################
      NAF2<-tryCatch({finADWC(NAF1,dIV)},error=function(cond){message("finADWC(NAF1,dIV) is failing...5")})
      faddu<-c(faddu,NAF2)
      ######################################################################
      #######################################################################
    }
  }else if(cha=="3+"){
    dIV<-"3"
    addu2<-tryCatch({qdapRegex::ex_between(addu1, "[", "]")[[1]]},warning=function(cond){message("error happened in precursortype extraction")})
    addu3<-tryCatch({AIN[AIN$V1==addu2 & AIN$V3==cha,]$V2},warning=function(cond){message("adduct match is empty")})
    #########################################################################
    ########################################################################
    if(!sjmisc::is_empty(str_replace_all(addu3, "\"", "")))
    {
      ##print(addu1)
      faddu<-c(faddu,addu3)
      
    }else{
      NAF<-unlist(strsplit(sub(".*\\[([^][]+)].*", "\\1",addu2),"(?=[+-])", perl = TRUE))
      NAF1<-tryCatch({comNFMDA(NAF)},error=function(cond){message("comNFMDA(NAF) is failing...6")})
      ###########################################################
      #############################################################
      NAF2<-tryCatch({finADWC(NAF1,dIV)},error=function(cond){message("finADWC(NAF1,dIV) is failing...6")})
      faddu<-c(faddu,NAF2)
      
    }
  }else{
    x<-"pass"
    
  }
  return(faddu)
}
#######################################################################
#######################################################################

#############################################################################################
#############################################################################################
#############################################################################################

MaKlist<-function(gFile)
{
  #############################################
  lines <- readLines(gFile)
  ## This is the original one commenting
  lst <-split(lines, cumsum(lines==""))
  ####lst <-split(lines, cummean(lines==""))
  ####lst <-split(lines, GMCM:::cummean(lines==""))
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
  ### adding this new code here
  NPTYL<-lapply(PTYL,function(x) length(x>0))
  NPTYL1<-which(unlist(NPTYL)!=0)
  PTYL1<-NPTYL1
  ## commenting this original code PTYL1<-Filter(length,PTYL)
  ##PTYL1<-Filter(length,PTYL)
  ##############################################
  ###### Adduct list ###########################
  ##############################################
  ##############################################
  faddu<-c()
  if(length(PTYL1)>0){
    for(i in c(1:length(PTYL1))){
     ##############################
	    ## added this new code
	    FLV<-PTYL1[i]
	    ## commenting the original addu and addu1
	    addu<-grep("PRECURSORTYPE:|ADDUCTIONNAME:",lst2[[FLV]])
	    addu1<-stringr::str_trim(gsub('PRECURSORTYPE:|ADDUCTIONNAME:','',lst2[[FLV]][addu]))

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
          NAF<-unlist(strsplit(sub(".*\\[([^][]+)].*", "\\1",addu2),"(?=[+-])", perl = TRUE))
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
          NAF<-unlist(strsplit(sub(".*\\[([^][]+)].*", "\\1",addu2),"(?=[+-])", perl = TRUE))	
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
          NAF<-unlist(strsplit(sub(".*\\[([^][]+)].*", "\\1",addu2),"(?=[+-])", perl = TRUE))
          ####################################################	
	  ####################################################
          NAF1<-tryCatch({comNFMDA(NAF)},error=function(cond){message("comNFMDA(NAF) failing check the charge 2")})
          NAF2<-tryCatch({finADWC(NAF1,dIV)},error=function(cond){message("finADWC(NAF1,dIV) failing check the charge 2")})
          faddu<-c(faddu,NAF2)
	  ###################################################
        }
        
      }else if(cha=="2+"){
        dIV<-"2"
        addu2<-tryCatch({qdapRegex::ex_between(addu1, "[", "]")[[1]]},warning=function(cond){message("error happened in precursortype extraction")})
        addu3<-tryCatch({AIN[AIN$V1==addu2 & AIN$V3==cha,]$V2},warning=function(cond){message("adduct match is empty")})
        if(!sjmisc::is_empty(addu3))
        {
          faddu<-c(faddu,addu3)
          
        }else{
          NAF<-unlist(strsplit(sub(".*\\[([^][]+)].*", "\\1",addu2),"(?=[+-])", perl = TRUE))
          ####################################################################	
	  ###################################################################
          NAF1<-tryCatch({comNFMDA(NAF)},error=function(cond){message("comNFMDA(NAF) failing check the charge 2 ..pos")})
          NAF2<-tryCatch({finADWC(NAF1,dIV)},error=function(cond){message("finADWC(NAF1,dIV) failing check the charge 2..pos")})
          faddu<-c(faddu,NAF2)
	  ####################################
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
          NAF<-unlist(strsplit(sub(".*\\[([^][]+)].*", "\\1",addu2),"(?=[+-])", perl = TRUE))
          ##############################################################	
	  ###############################################################
          NAF1<-tryCatch({comNFMDA(NAF)},error=function(cond){message("comNFMDA(NAF) failing check the charge 3 ..neg")})
          NAF2<-tryCatch({finADWC(NAF1,dIV)},error=function(cond){message("finADWC(NAF1,dIV) failing check the charge 3..neg")})
          faddu<-c(faddu,NAF2)
	  ##############################################################
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
          NAF<-unlist(strsplit(sub(".*\\[([^][]+)].*", "\\1",addu2),"(?=[+-])", perl = TRUE))
          ########################################################	
	  #########################################################
          NAF1<-tryCatch({comNFMDA(NAF)},error=function(cond){message("comNFMDA(NAF) failing check the charge 3 ..pos")})
          NAF2<-tryCatch({finADWC(NAF1,dIV)},error=function(cond){message("finADWC(NAF1,dIV) failing check the charge 3..pos")})
          faddu<-c(faddu,NAF2)
	  ###########################################################
          
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
  ######### Coping to new list ##########################################################################
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
  #######################################################################################################
  
}
##########################################################################################################
##########################################################################################################
##########################################################################################################

##########################################################################################################
##########################################################################################################
##########################################################################################################
NFFilter1<-function(InMEDA,InAdVA,InMSPL,InPMZ,InRTL)
{
  ###############################################
       FMa<-c()
  ###############################################
  print("entering the NFFilter1 function")
  #############################################
  if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & !startsWith(as.character(InMEDA[["InChI"]]),'not available') & !startsWith(as.character(InMEDA[["InChI"]]),'CAS:') & !startsWith(as.character(InMEDA[["InChI"]]),'InChI=')){
    if(tryCatch({webchem::is.inchikey(stringr::str_trim(as.character(InMEDA[["InChI"]])))},error=function(cond){message("inchikey..file must be empty")})){
      ########################################################################
      print("enter the Inchi key AREA..inchikey is not empty")
      ########################################################################
      IK<-as.character(InMEDA[["InChI"]])
      IK1<-tryCatch({webchem::get_cid(IK, from = "inchikey")},error=function(cond){message("Inchi name must be empty or rinchi not abe to fetch")})
      PCID<-tryCatch({IK1[[2]][1]},error=function(cond){message("Pubchem CID is empty")})
      FPUCID1<-as.numeric(PCID)
      ########################################################################
      PCID1<-tryCatch({webchem::pc_prop(as.numeric(stringr::str_trim(gsub("[[:punct:]]", "",PCID))), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem Id is empty")})
      ###############################################################################
      ################################################################################
      PCID2<-ifelse(!sjmisc::is_empty(as.numeric(PCID1$ExactMass)),as.numeric(PCID1$ExactMass),ifelse(!sjmisc::is_empty(tryCatch({as.numeric(ConvCIDtoOID1(FPUCID1)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")})),tryCatch({as.numeric(ConvCIDtoOID1(FPUCID1)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}),tryCatch({PuCIDtoEM(FPUCID1)},error=function(cond){message("Pubchem CID is empty..did not get exact mass1")})))
      ##########################################################################
      ##########################################################################
      if(!sjmisc::is_empty(PCID2) & (PCID2 !=0))
      {
        #####################################      
        print("enter the if loop...inkikey")
        #####################################
        FMa<-c(FMa,PCID2)
        #####################################
      }else{
        ###########################################################      
        print("enter the else part")      
        ###########################################################
        if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["InChI"]]))) & startsWith(stringr::str_trim(as.character(InMEDA[["InChI"]])),'InChI=')){
          ##############################################
          print("enter the Inchi part ...in else loop...Inchikey")
          ##############################################
	  IK<-stringr::str_trim(as.character(InMEDA[["InChI"]]))
          IK1<-tryCatch({webchem::get_cid(IK, from = "inchi")},error=function(cond){message("Inchi name must be empty or rinchi not abe to fetch")})
          PCID<-tryCatch({IK1[[2]][1]},error=function(cond){message("Pubchem Id is empty")})
          #######################################
          #########################################
          PCID1<-tryCatch({webchem::pc_prop(as.numeric(stringr::str_trim(gsub("[[:punct:]]", "",PCID))), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem Id is empty")})
          PCID2<-tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){message("Pubchem Id is empty")})
          ###################################
          if(!sjmisc::is_empty(PCID2)){
            ###########################################	  
            print("enter the if loop ...inchi area")
            ###########################################	  
            FMa<-c(FMa,PCID2)
            ############################################
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
            PCID1<-tryCatch({webchem::pc_prop(as.numeric(PCID[[2]][1]), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES","InChI","InChIKey"))},error=function(cond){message("Pubchem CId is empty")})
	    FPUCID1<-tryCatch({as.numeric(PCID[[2]][1])},error=function(cond){message("Pubchem CID is empty")})
	    PCID2<-ifelse(!sjmisc::is_empty(as.numeric(PCID1$ExactMass)),as.numeric(PCID1$ExactMass),ifelse(!sjmisc::is_empty(tryCatch({as.numeric(ConvCIDtoOID1(FPUCID1)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")})),tryCatch({as.numeric(ConvCIDtoOID1(FPUCID1)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}),tryCatch({PuCIDtoEM(FPUCID1)},error=function(cond){message("Pubchem CID is empty..did not get exact mass1")})))
	    #############################################
            #############################################
            if(!sjmisc::is_empty(PCID2) & (PCID2 !=0)){
              ###########################################	    
              print("enter the if loop ...cas area")
              ###########################################	    
              FMa<-c(FMa,PCID2)
              ###########################################
            }else if(!sjmisc::is_empty(as.character(InMEDA[["SMILES"]])) & !startsWith(as.character(InMEDA[["SMILES"]]),'not available')){
              ######################################
              print("entering smiles area in InchiKey")
              ########################################
              IK<-stringr::str_trim(as.character(InMEDA[["SMILES"]]))
              ########################################
	      ########################################
              tes<-tryCatch({rcdk::parse.smiles(IK)},error=function(cond){message("smiles not abe to parse")})
              tes1<-tryCatch({tes[[1]]},error=function(cond){message("smiles parse information is empty")})
              tes2<-tryCatch({rcdk::get.exact.mass(tes1)},error=function(cond){message("smiles not abe to fetch")})
              PCID2<-tryCatch({tes2},error=function(cond){message("smiles not abe to fetch")})
              #######################################
              #######################################
              if(!sjmisc::is_empty(PCID2)){
                ############################      
                FMa<-c(FMa,PCID2)
                #############################
              }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(InMEDA[["PubChem CID"]]),'not available')){
                ######################################
		#########################################
                FPUCID<-stringr::str_trim(gsub("[[:punct:]]", "",as.character(InMEDA[["PubChem CID"]])))
                FPUCID1<-tryCatch({as.numeric(FPUCID)},error=function(cond){message("FPUCID must be empty")})
	        ######################################
		######################################
                PCID1<-tryCatch({webchem::pc_prop(as.numeric(FPUCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem CId is empty..did not get exact mass")})
                PCID2<-ifelse(!sjmisc::is_empty(as.numeric(PCID1$ExactMass)),as.numeric(PCID1$ExactMass),ifelse(!sjmisc::is_empty(tryCatch({as.numeric(ConvCIDtoOID1(FPUCID1)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")})),tryCatch({as.numeric(ConvCIDtoOID1(FPUCID1)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}),tryCatch({PuCIDtoEM(FPUCID1)},error=function(cond){message("Pubchem CID is empty..did not get exact mass1")})))
		######################################
                ######################################
                if(!sjmisc::is_empty(PCID2) & (PCID2 !=0)){
                  ############################################	
                  print("entering the if loop..Pubchem CID")
                  ###########################################
                  FMa<-c(FMa,PCID2)
                  ###########################################
                }else if(tryCatch({as.numeric(ConvCIDtoOID1(FPUCID)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}) != 0)
                   {
                  FMa<-c(FMa,as.numeric(ConvCIDtoOID1(FPUCID)[5]))

               }else if(tryCatch({as.numeric(PuCIDtoEM(FPUCID))},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}) != 0)
                   {  
                 FMa<-c(FMa,as.numeric(PuCIDtoEM(FPUCID)))
               }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["Formula"]]))) & !startsWith(stringr::str_trim(as.character(InMEDA[["Formula"]])),'not available')){
                  ########################################################
                  ########################################################	
                  EM<-tryCatch({Rdisop::getMolecule(stringr::str_trim(as.character(InMEDA[["Formula"]])))},error=function(cond){message("Pubchem Id is empty")})
                  EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem CID is empty")})
                  #########################################################
                  #########################################################
                  if(!sjmisc::is_empty(EM1)){
                    FMa<-c(FMa,EM1)
                  }else{
                    
                    ########################################
                    ########################################	  
                    PMA<-tryCatch({ConvPCIDtoOCN(FPUCID)},error=function(cond){message("Pubchem value is empty")})
                    FMa<-c(FMa,ifelse(!sjmisc::is_empty(tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")})),tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")}),FuFtoRe1(InMEDA)))
		    #########################################
                    ########################################
                  }
                  ######################################
                }else{
                  #####################################
                  #####################################	
                  PMA<-tryCatch({FuFtoRe1(InMEDA)},error=function(cond){message("FuFtoRe1 is empty")})	
                  FMa<-c(FMa,PMA)
                  ############################
                  #############################
                }
              }else{
                ### PubchemId is not found ..entering the else loop
                if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["Formula"]]))) & !startsWith(stringr::str_trim(as.character(InMEDA[["Formula"]])),'not available')){
                  ############################################
                  ############################################	
                  EM<-tryCatch({Rdisop::getMolecule(stringr::str_trim(as.character(InMEDA[["Formula"]])))},error=function(cond){message("Pubchem Id is empty")})
                  EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
                  ############################################
                  ############################################
                  if(!sjmisc::is_empty(EM1)){
                    ##print("formula area.....")
                    ##print(EM1)
                    FMa<-c(FMa,EM1)
                  }else{
                    #############################
                    #############################		  
                    PMA<-tryCatch({FuFtoRe1(InMEDA)},error=function(cond){message("FuFtoRe1 is empty")})		  
                    FMa<-c(FMa,PMA)
                    #############################
                    #############################
                  }
                }
             ###############################
              }
            ####################  PCID2 ...end ###
            }else{
              ##############################################
              ## smiles end ... not found
              ###############################################
              if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(InMEDA[["PubChem CID"]]),'not available')){
                ###############################################
                print("enter the CID ...area")
                ###############################################
                ##FPUCID<-stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))
	        ################################################
                FPUCID<-stringr::str_trim(gsub("[[:punct:]]", "",as.character(InMEDA[["PubChem CID"]])))      
                FPUCID1<-tryCatch({as.numeric(FPUCID)},error=function(cond){message("FPUCID value is empty")})
		##########################################
                ##########################################
                PCID1<-tryCatch({webchem::pc_prop(as.numeric(PCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem CID is empty")})
		################################################
		################################################
                PCID2<-ifelse(!sjmisc::is_empty(as.numeric(PCID1$ExactMass)),as.numeric(PCID1$ExactMass),ifelse(!sjmisc::is_empty(tryCatch({as.numeric(ConvCIDtoOID1(FPUCID1)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")})),tryCatch({as.numeric(ConvCIDtoOID1(FPUCID1)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}),tryCatch({PuCIDtoEM(FPUCID1)},error=function(cond){message("Pubchem CID is empty..did not get exact mass1")})))
		##########################################
		##########################################
                if(!sjmisc::is_empty(PCID2) & (PCID2 !=0)){
                  
                  FMa<-c(FMa,PCID2)
                }else if(tryCatch({as.numeric(ConvCIDtoOID1(FPUCID)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}) != 0){
                 FMa<-c(FMa,tryCatch({as.numeric(ConvCIDtoOID1(FPUCID)[5])},error=function(cond){message("ConvCIDtoOID1 value is empty")}))

                }else if(tryCatch({as.numeric(PuCIDtoEM(FPUCID))},error=function(cond){message("Pubchem CID is empty.did not get exact mass")}) !=0){
                 FMa<-c(FMa,as.numeric(PuCIDtoEM(FPUCID)))

                 }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["Formula"]]))) & !startsWith(stringr::str_trim(as.character(InMEDA[["Formula"]])),'not available')){
                  ##############################################
                  ##############################################	
                  EM<-tryCatch({Rdisop::getMolecule(stringr::str_trim(as.character(InMEDA[["Formula"]])))},error=function(cond){message("Pubchem CID is empty")})
                  EM1<-tryCatch({EM$exactmass},error=function(cond){message("Rdisop::getMolecule exact mass is empty")})
                  ###############################################
                  ###############################################
                  if(!sjmisc::is_empty(EM1)){
                    
                    FMa<-c(FMa,EM1)
                  }else{
                    
                    ####################################
		    ####################################	  
                    PMA<-tryCatch({ConvPCIDtoOCN(FPUCID)},error=function(cond){message("Pubchem value is empty")})
		    #####################################
		    ########################################
                    FMa<-c(FMa,ifelse(!sjmisc::is_empty(tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")})),tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")}), tryCatch({FuFtoRe1(InMEDA)},error=function(cond){message("FuFtoRe1 conversion is empty")})))
                    ####################################
                    ####################################
                  }
                  ######################################
                }else{
                  
	          #################################		
                  PMA<-tryCatch({FuFtoRe1(InMEDA)},error=function(cond){message("FuFtoRe1 conversion is empty")})	
                  FMa<-c(FMa,PMA)
                  #################################
                }
              #########################################################
              }else{
                ###### PubchemID is not found... so getting exact mass from
                if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["Formula"]]))) & !startsWith(stringr::str_trim(as.character(InMEDA[["Formula"]])),'not available')){
                  ############################################
		  ############################################	
                  EM<-tryCatch({stringr::str_trim(Rdisop::getMolecule(as.character(InMEDA[["Formula"]])))},error=function(cond){message("formula to mass conversion is empty")})
                  EM1<-tryCatch({EM$exactmass},error=function(cond){message("Rdisop::getMolecule is empty")})
                  ##########################################
                  ##########################################
                  if(!sjmisc::is_empty(EM1)){
                    FMa<-c(FMa,EM1)
                  }else{
                    #############################
                    PMA<-tryCatch({FuFtoRe1(InMEDA)},error=function(cond){message("FuFtoRe1 conversion is empty")})	  
                    FMa<-c(FMa,PMA)
                    ##############################
                  }
                }
              #######################################
              }
            #############################
            } ##end if loop Pubchem CID
          ###############################
          } ### end of else ...smiles end ... not found
       #################################
        }else{
       ### Inchikey end start of smiles##
          if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["SMILES"]]))) & !startsWith(stringr::str_trim(as.character(InMEDA[["SMILES"]])),'not available')){
            ###########################################
            print("enter the smiles area ..smiles")
            ###########################################
            IK<-stringr::str_trim(as.character(InMEDA[["SMILES"]]))
            ###########################################
            tes<-tryCatch({rcdk::parse.smiles(IK)},error=function(cond){message("smiles not abe to parse")})
            tes1<-tryCatch({tes[[1]]},error=function(cond){message("smiles parse information is empty")})
            tes2<-tryCatch({rcdk::get.exact.mass(tes1)},error=function(cond){message("smiles not abe to fetch")})
            PCID2<-tryCatch({tes2},error=function(cond){message("smiles not abe to fetch")})
            Fval<-tryCatch({PuSmilesToEM(InMEDA[["SMILES"]])},error=function(cond){return(0)})
            ################################################
            ################################################
            if(!sjmisc::is_empty(PCID2)){
              FMa<-c(FMa,PCID2)
            }else if(!sjmisc::is_empty(Fval) & Fval != 0){
              
              FMa<-c(FMa,Fval)	    
              
            }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(InMEDA[["PubChem CID"]]),'not available')){
              ######################################
              print("enter the CID pass area ...")
              ######################################
              ######################################	    
              FPUCID<-stringr::str_trim(gsub("[[:punct:]]", "",as.character(InMEDA[["PubChem CID"]])))
              FPUCID1<-tryCatch({as.numeric(FPUCID)},error=function(cond){return(NA)})
	      ######################################
              ######################################
              PCID1<-tryCatch({webchem::pc_prop(as.numeric(FPUCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem CID is empty")})
	      ########################################
	      #########################################
              PCID2<-ifelse(!sjmisc::is_empty(as.numeric(PCID1$ExactMass)),as.numeric(PCID1$ExactMass),ifelse(!sjmisc::is_empty(tryCatch({as.numeric(ConvCIDtoOID1(FPUCID1)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")})),tryCatch({as.numeric(ConvCIDtoOID1(FPUCID1)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}),tryCatch({PuCIDtoEM(FPUCID1)},error=function(cond){message("Pubchem CID is empty..did not get exact mass1")})))
	      #######################################
              ########################################
              if(!sjmisc::is_empty(PCID2) & (PCID2 !=0)){
                
                FMa<-c(FMa,PCID2)
                
              }else if(tryCatch({as.numeric(ConvCIDtoOID1(FPUCID)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}) != 0){
                FMa<-c(FMa,tryCatch({as.numeric(ConvCIDtoOID1(FPUCID)[5])},error=function(cond){message("Pubchem CID to other value is empty")}))
              }else if(tryCatch({as.numeric(PuCIDtoEM(FPUCID))},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}) !=0){
                 FMa<-c(FMa,as.numeric(PuCIDtoEM(FPUCID)))
              }else if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){
                ##################################
                EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem CID is empty")})
                EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem CID is empty")})
                ##################################
                if(!sjmisc::is_empty(EM1)){
                  FMa<-c(FMa,EM1)
                }else{
                  ############################
                  ############################
		  ############################	
                  PMA<-tryCatch({ConvPCIDtoOCN(FPUCID)},error=function(cond){message("Pubchem to other conversion is empty")})
		 ###################################################################
		 ##########################################################
                  FMa<-c(FMa,ifelse(!sjmisc::is_empty(tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")})),tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")}),tryCatch({FuFtoRe1(InMEDA)},error=function(cond){message("FuFtoRe1 is empty")})))
                  #############################
                  #############################
                }
                ################################
              }else{
                ## formula not found and --smiles and pubchem failed to get PubchemID
                ##FMa<-c(FMa,PCID2)
                ##############################
                PMA<-tryCatch({FuFtoRe1(InMEDA)},error=function(cond){message("FuFtoRe1 is empty")})      
                FMa<-c(FMa,PMA)
                ###############################
              }
            }else{
              ### PubchemId is not found ..entering the else loop
              if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["Formula"]]))) & !startsWith(stringr::str_trim(as.character(InMEDA[["Formula"]])),'not available')){
                ###################################
                ###################################      
                EM<-tryCatch({Rdisop::getMolecule(stringr::str_trim(as.character(InMEDA[["Formula"]])))},error=function(cond){message("Pubchem CID is empty")})
                EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem CID is empty")})
                ###################################		
                ###################################
                if(!sjmisc::is_empty(EM1)){
                  FMa<-c(FMa,EM1)
                }else{
                  ##############################
                  PMA<-tryCatch({FuFtoRe1(InMEDA)},error=function(cond){message("FuFtoRe1 is empty")})	
                  FMa<-c(FMa,PMA)
                  #############################
                }
              }
            #######################################
            }
          ####################  PCID2 ...end ######
          }else{
            ## smiles end ... not found
            ############################################################
            if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(stringr::str_trim(InMEDA[["PubChem CID"]])),'not available')){
              ##########################################################
              print("enter the CID ...area ...")
              ##########################################################
              ##FPUCID<-stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))
              FPUCID<-stringr::str_trim(gsub("[[:punct:]]", "",as.character(InMEDA[["PubChem CID"]])))
              FPUCID1<-tryCatch({as.numeric(FPUCID)},error=function(cond){message("Convert to numeric Pubchem CID is empty")})
	      ##########################################################
	      ##########################################################
              PCID1<-tryCatch({webchem::pc_prop(as.numeric(FPUCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem CID is empty")})
              ##PCID2<-tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){message("Pubchem CID is empty")})
              PCID2<-ifelse(!sjmisc::is_empty(as.numeric(PCID1$ExactMass)),as.numeric(PCID1$ExactMass),ifelse(!sjmisc::is_empty(tryCatch({as.numeric(ConvCIDtoOID1(FPUCID1)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")})),tryCatch({as.numeric(ConvCIDtoOID1(FPUCID1)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}),tryCatch({PuCIDtoEM(FPUCID1)},error=function(cond){message("Pubchem CID is empty..did not get exact mass1")})))
	      #####################################################
              #####################################################
              if(!sjmisc::is_empty(PCID2)){
                
                FMa<-c(FMa,PCID2)
                
              }else if(tryCatch({as.numeric(ConvCIDtoOID1(FPUCID)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}) != 0){
                FMa<-c(FMa, tryCatch({as.numeric(ConvCIDtoOID1(FPUCID)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}))
              }else if(tryCatch({as.numeric(PuCIDtoEM(FPUCID))},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}) !=0){
                 FMa<-c(FMa,as.numeric(PuCIDtoEM(FPUCID)))
              }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["Formula"]]))) & !startsWith(stringr::str_trim(as.character(InMEDA[["Formula"]])),'not available')){
                
                ###############################################
                EM<-tryCatch({Rdisop::getMolecule(stringr::str_trim(as.character(InMEDA[["Formula"]])))},error=function(cond){message("Rdisop::getMolecule conversion is empty")})
                EM1<-tryCatch({EM$exactmass},error=function(cond){message("Exact mass convesrsion failed")})
                ################################################
                if(!sjmisc::is_empty(EM1)){
                  
                  FMa<-c(FMa,EM1)
                }else{
                  #################################
                  #################################	
                  PMA<-tryCatch({ConvPCIDtoOCN(FPUCID)},error=function(cond){message("Pubchem value is empty")})
                  FMa<-c(FMa,ifelse(!sjmisc::is_empty(tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")})),tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")}),tryCatch({FuFtoRe1(InMEDA)},error=function(cond){message("FuFtoRe1 is empty")})))
                  ##################################
                  ##################################
                }
                #######################################
              }else{
                ### Trying to get the exact mass from name provided
                ########################################
                ##PMA<-FuFtoRe(InMEDA)
		#########################################      
                PMA<-tryCatch({FuFtoRe1(InMEDA)},error=function(cond){message("FuFtoRe1 is empty")})      
                FMa<-c(FMa,PMA)
                ########################################                
                 }
              ##########################################
            }else{
              ###### PubchemID is not found... so getting exact mass from
              if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["Formula"]]))) & !startsWith(stringr::str_trim(as.character(InMEDA[["Formula"]])),'not available')){
                #######################################
                EM<-tryCatch({Rdisop::getMolecule(stringr::str_trim(as.character(InMEDA[["Formula"]])))},error=function(cond){message("formula to exact mass conversion failed")})
                EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem CID is empty")})
                ########################################
                if(!sjmisc::is_empty(EM1)){
		 ######################################	
                  ##print("formula area.....")
                  ##print(EM1)
		  ##################################	
                  FMa<-c(FMa,EM1)
                }else{
		  ####################################	
                  #########################
                  ##PMA<-FuFtoRe(InMEDA)
                  PMA<-tryCatch({FuFtoRe1(InMEDA)},error=function(cond){message("FuFtoRe1 is empty")})	
                  FMa<-c(FMa,PMA)
                  #########################
                }
             ###############################
              } ## formula
            } ## end of else loop
          } ## end of else
        }### InCHi.. end
      } ## end of main else
    }else{
   ###################################################	    
      if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["InChI"]]))) & startsWith(stringr::str_trim(as.character(InMEDA[["InChI"]])),'InChI=')){
        ##############################################
        print("enter the Inchi part ...")
        IK<-stringr::str_trim(as.character(InMEDA[["InChI"]]))
        ##############################################
        IK1<-tryCatch({webchem::get_cid(IK, from = "inchi")},error=function(cond){message("Inchi name must be empty or rinchi not abe to fetch")})
        PCID<-tryCatch({IK1[[2]][1]},error=function(cond){message("Pubchem Id is empty")})
        #########################################
        #########################################
        PCID1<-tryCatch({webchem::pc_prop(as.numeric(stringr::str_trim(gsub("[[:punct:]]", "",PCID))), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem Id is empty")})
        FPUCID1<-as.numeric(as.numeric(stringr::str_trim(gsub("[[:punct:]]", "",PCID))))
	PCID2<-ifelse(!sjmisc::is_empty(as.numeric(PCID1$ExactMass)),as.numeric(PCID1$ExactMass),ifelse(!sjmisc::is_empty(tryCatch({as.numeric(ConvCIDtoOID1(FPUCID1)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")})),tryCatch({as.numeric(ConvCIDtoOID1(FPUCID1)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}),tryCatch({PuCIDtoEM(FPUCID1)},error=function(cond){message("Pubchem CID is empty..did not get exact mass1")})))
	###########################################
	#####################################################################################################
        if(!sjmisc::is_empty(PCID2) & (PCID2 !=0)){
          ###########################################	  
          print("enter the if loop ...inchi area")
          ###########################################	  
          FMa<-c(FMa,PCID2)
          ############################################
        }else if(tryCatch({as.numeric(ConvCIDtoOID1(FPUCID)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}) != 0){
           FMa<-c(FMa,tryCatch({as.numeric(ConvCIDtoOID1(FPUCID)[5])},error=function(cond){return(0)}))
        }else if(tryCatch({as.numeric(PuCIDtoEM(FPUCID))},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}) !=0){
           FMa<-c(FMa, tryCatch({as.numeric(PuCIDtoEM(FPUCID))},error=function(cond){return(0)})) 
        }else if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & startsWith(as.character(InMEDA[["InChI"]]),'InChI=') & !sjmisc::is_empty(tryCatch({CONinctoEM(InMEDA)},error=function(cond){message("Inchi to exact mass is empty")}))){
          
          FMa<-c(FMa,tryCatch({CONinctoEM(InMEDA)},error=function(cond){message("Inchi to exact mass is empty")}))
        }else{
          PMA<-tryCatch({FuFtoRe1(InMEDA)},error=function(cond){message("FuFtoRe1 to exact mass is empty")})	
          FMa<-c(FMa,PMA)
        }
      }else{
        PMA<-tryCatch({FuFtoRe1(InMEDA)},error=function(cond){message("FuFtoRe1 to exact mass is empty")})
        FMa<-c(FMa,PMA)
      }
    }###end of else
  }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["InChI"]]))) & startsWith(stringr::str_trim(as.character(InMEDA[["InChI"]])) ,'CAS:')){
  ###############################################
    ##browser()	  
    print("enter the cas function area.... cas is avilable")
    ############################################
    CV<-stringr::str_trim(as.character(InMEDA[["InChI"]]))
    CV1<-stringr::str_replace(CV,pattern='CAS:',replacement ="")
    CV2<-stringr::str_trim(as.character(CV1))
    ############################################
    PCID<-tryCatch({webchem::get_cid(CV2, from = "xref/rn",match="first")},error=function(cond){message("Pubchem CID is empty")})
    #############################################
    PCID1<-tryCatch({webchem::pc_prop(as.numeric(stringr::str_trim(gsub("[[:punct:]]", "",PCID$cid))), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem Id is empty")})
    ##########################################################
    ##PCID1<-tryCatch({webchem::pc_prop(as.numeric(PCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem Id is empty")})
    ###########################################################
    FPUCID1<-as.numeric(as.numeric(stringr::str_trim(gsub("[[:punct:]]", "",PCID$cid))))
    PCID2<-ifelse(!sjmisc::is_empty(tryCatch({as.numeric(PCID1$ExactMass)},error=function(cond){message("ExactMass in CAS empty")})),as.numeric(PCID1$ExactMass),ifelse(!sjmisc::is_empty(tryCatch({as.numeric(ConvCIDtoOID1(FPUCID1)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")})),tryCatch({as.numeric(ConvCIDtoOID1(FPUCID1)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}),tryCatch({PuCIDtoEM(FPUCID1)},error=function(cond){message("Pubchem CID is empty..did not get exact mass1")})))
    ##############################################################
    print("entering the PCID2 value")
    print(PCID2)
    print(!sjmisc::is_empty(PCID2))
    #############################################################
    if(!sjmisc::is_empty(PCID2) & (PCID2 != 0)){
      
      FMa<-c(FMa,PCID2)
      
    }else if(!sjmisc::is_empty(as.character(InMEDA[["SMILES"]])) & !startsWith(stringr::str_trim(as.character(InMEDA[["SMILES"]])),'not available')){
      ##########################################
      print("enter the smiles..in cas area")
      ##########################################
      ##########################################	    
      IK<-stringr::str_trim(as.character(InMEDA[["SMILES"]]))
      tes<-tryCatch({rcdk::parse.smiles(IK)},error=function(cond){message("smiles not abe to parse")})
      tes1<-tryCatch({tes[[1]]},error=function(cond){message("smiles parse information is empty")})
      tes2<-tryCatch({rcdk::get.exact.mass(tes1)},error=function(cond){message("smiles not abe to fetch")})
      PCID2<-tryCatch({tes2},error=function(cond){message("smiles not abe to fetch")})
      ##########################################
      ##########################################
      if(!sjmisc::is_empty(PCID2)){
        
        FMa<-c(FMa,PCID2)
        
      }else if(!sjmisc::is_empty(as.character(InMEDA[["SMILES"]])) & !startsWith(as.character(InMEDA[["SMILES"]]),'not available')){
	      PCID2<-ifelse(!sjmisc::is_empty(tryCatch({ConvSMItoOID1(InMEDA[["SMILES"]])[5]}, error = function(x) {return(0)})),tryCatch({ConvSMItoOID1(InMEDA[["SMILES"]])[5]}, error = function(x) {return(0)}),tryCatch({PuSmilesToEM(InMEDA[["SMILES"]])}, error = function(x) {return(0)}))
	      ##FMa<-c(FMa,PCID2)
	      if(!sjmisc::is_empty(PCID2) & PCID2 != 0)
	      {
		      FMa<-c(FMa,PCID2)

	      }else{
		      PCID2<-ifelse(!sjmisc::is_empty(tryCatch({CONSMItoEM(InMEDA)},error=function(cond){message("SMILES to exact mass is empty")})),tryCatch({CONSMItoEM(InMEDA)},error=function(cond){message("SMILES to exact mass is empty")}),ifelse(!sjmisc::is_empty(tryCatch({mzAnnotation::smileToAccurateMass(as.character(InMEDA[["SMILES"]]))}, error = function(x) {return(0)})), tryCatch({mzAnnotation::smileToAccurateMass(as.character(InMEDA[["SMILES"]]))}, error = function(x) {return(0)}), tryCatch({PuSmilesToEM(InMEDA[["SMILES"]])},error = function(x) {return(0)})))
		      FMa<-c(FMa,PCID2)

	      }

      }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(stringr::str_trim(as.character(InMEDA[["PubChem CID"]])),'not available')){
        ###############################################################
	################################################################      
        PCID<-stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))      
        FPUCID<-stringr::str_trim(gsub("[[:punct:]]", "",PCID))      
        FPUCID1<-as.numeric(FPUCID)
	############################################
        ############################################
        PCID1<-tryCatch({webchem::pc_prop(as.numeric(FPUCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem Id is empty")})
        PCID2<-ifelse(!sjmisc::is_empty(as.numeric(PCID1$ExactMass)),as.numeric(PCID1$ExactMass),ifelse(!sjmisc::is_empty(tryCatch({as.numeric(ConvCIDtoOID1(FPUCID1)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")})),tryCatch({as.numeric(ConvCIDtoOID1(FPUCID1)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}),tryCatch({PuCIDtoEM(FPUCID1)},error=function(cond){message("Pubchem CID is empty..did not get exact mass1")})))
	#############################################
	###########################################
        if(!sjmisc::is_empty(PCID2) & (PCID2 !=0)){
          
          FMa<-c(FMa,PCID2)
          
        }else if(tryCatch({as.numeric(ConvCIDtoOID1(FPUCID)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}) != 0){
          FMa<-c(FMa,tryCatch({as.numeric(ConvCIDtoOID1(FPUCID)[5])},error=function(cond){return(0)}))

        }else if(tryCatch({as.numeric(PuCIDtoEM(FPUCID))},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}) !=0){
          FMa<-c(FMa,tryCatch({as.numeric(PuCIDtoEM(FPUCID))},error=function(cond){return(0)}))

        }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["Formula"]]))) & !startsWith(stringr::str_trim(as.character(InMEDA[["Formula"]])),'not available')){
          ##########################################
          ##########################################
          EM<-tryCatch({stringr::str_trim(Rdisop::getMolecule(as.character(InMEDA[["Formula"]])))},error=function(cond){message("Formula to exact mass is empty")})
          EM1<-tryCatch({EM$exactmass},error=function(cond){message("Exact mass is empty")})
          #########################################
          #########################################
          if(!sjmisc::is_empty(EM1)){
            
            FMa<-c(FMa,EM1)
          }else{

            ###################################
            ###################################		  
            PMA<-tryCatch({ConvPCIDtoOCN(FPUCID)},error=function(cond){message("Pubchem value is empty")})
	    ##########################################################
	    ##########################################################
            FMa<-c(FMa,ifelse(!sjmisc::is_empty(tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")})),tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")}), tryCatch({FuFtoRe1(InMEDA)},error=function(cond){message("FuFtoRe1 value is empty")})))
            ###################################
            ###################################
          }
          
        }else{

          ###################################
	  print("else loop ..NAME to exact mass..")	
          PMA<-tryCatch({FuFtoRe1(InMEDA)},error=function(cond){message("FuFtoRe1 is empty")})	
          FMa<-c(FMa,PMA)
          #################################
        }
        #######################################
      }else{
        ##print("enter the else area ..2")
        #### PubchemID is not found
        ###################################################
        if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["Formula"]]))) & !startsWith(stringr::str_trim(as.character(InMEDA[["Formula"]])),'not available')){
          ################################################
          EM<-tryCatch({stringr::str_trim(Rdisop::getMolecule(as.character(InMEDA[["Formula"]])))},error=function(cond){message("Pubchem Id is empty")})
          EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
          ###########################################
          ###########################################
          if(!sjmisc::is_empty(EM1)){
            ####################################
            FMa<-c(FMa,EM1)
            ############################
          }else{
            #############################
            PMA<tryCatch({FuFtoRe1(InMEDA)},error=function(cond){message("FuFtoRe1 is empty")})
            #########################	  
            ##PMA<-FuFtoRe(InMEDA)
	    ##########################
            FMa<-c(FMa,PMA)
            ###########################
          }
        }
        #######################################
      }
    ############## Smiles is empty and Pubchem CID
    }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(InMEDA[["PubChem CID"]]),'not available')){
      ########################################################
      print("entering the pubchemCID area")
      #########################################################
      PCID<-as.character(InMEDA[["PubChem CID"]])	    
      FPUCID<-stringr::str_trim(gsub("[[:punct:]]", "",PCID))
      FPUCID1<-tryCatch({as.numeric(FPUCID)},error=function(cond){return(0)})      
      #########################################################
      PCID1<-tryCatch({webchem::pc_prop(as.numeric(FPUCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem Id is empty")})
      ###############################################
      ###############################################
      PCID2<-ifelse(!sjmisc::is_empty(as.numeric(PCID1$ExactMass)),as.numeric(PCID1$ExactMass),ifelse(!sjmisc::is_empty(tryCatch({as.numeric(ConvCIDtoOID1(FPUCID1)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")})),tryCatch({as.numeric(ConvCIDtoOID1(FPUCID1)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}),tryCatch({PuCIDtoEM(FPUCID1)},error=function(cond){message("Pubchem CID is empty..did not get exact mass1")})))
      ##################################
      ##################################
      if(!sjmisc::is_empty(PCID2) & (PCID2 !=0)){
        
        FMa<-c(FMa,PCID2)
        
      }else if(tryCatch({as.numeric(ConvCIDtoOID1(FPUCID)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}) != 0){
        FMa<-c(FMa,tryCatch({as.numeric(ConvCIDtoOID1(FPUCID)[5])},error=function(cond){return(0)}))

      }else if(tryCatch({as.numeric(PuCIDtoEM(FPUCID))},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}) !=0){
        FMa<-c(FMa,tryCatch({as.numeric(PuCIDtoEM(FPUCID))},error=function(cond){return(0)}))

      }else if(!sjmisc::is_empty(as.character(InMEDA[["Formula"]])) & !startsWith(as.character(InMEDA[["Formula"]]),'not available')){
        #####################################
        #####################################      
        EM<-tryCatch({Rdisop::getMolecule(as.character(InMEDA[["Formula"]]))},error=function(cond){message("Pubchem Id is empty")})
        EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
        ####################################
        ####################################
        if(!sjmisc::is_empty(EM1)){
          
          FMa<-c(FMa,EM1)
        }else{
          
          ####################################
          PMA<-tryCatch({ConvPCIDtoOCN(FPUCID)},error=function(cond){message("Pubchem value is empty")})
	  ####################################
          ###FMa<-c(FMa,ifelse(!sjmisc::is_empty(tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")})),tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")}),0))
	  ##################################################
          FMa<-c(FMa,ifelse(!sjmisc::is_empty(tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")})),tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")}), tryCatch({FuFtoRe1(InMEDA)},error=function(cond){message("FuFtoRe1 is empty")})))
          ########################################
          ########################################
        }
        ####################################
        ####################################
        
      }else{
        print("entering Pubchem CID else loop---not got excat mass from pubchem CID")
        #######################
        PMA<-tryCatch({FuFtoRe1(InMEDA)},error=function(cond){message("FuFtoRe1 is empty")})
        FMa<-c(FMa,PMA)
        ######################
        
      }###end of else
      #######################################
    }else{
      ##print("enter the else area ..pubcehm CID not found")
      #### PubchemID is not found
      if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["Formula"]])))  & !startsWith(stringr::str_trim(as.character(InMEDA[["Formula"]])),'not available')){
        #########################################
        #########################################      
        EM<-tryCatch({stringr::str_trim(Rdisop::getMolecule(as.character(InMEDA[["Formula"]])))},error=function(cond){message("Pubchem Id is empty")})
        EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
        ########################################
        ########################################
        if(!sjmisc::is_empty(EM1)){
          
          FMa<-c(FMa,EM1)
          
        }else{
          ############################
          ############################		
          PMA<-tryCatch({FuFtoRe1(InMEDA)},error=function(cond){message("FuFtoRe1 is empty")})	
          FMa<-c(FMa,PMA)
          #############################
	  #############################
          
        }
        ################################
      }else{
        ###############################
	###############################      
        PMA<-tryCatch({FuFtoRe1(InMEDA)},error=function(cond){message("FuFtoRe1 value is empty")})      
        FMa<-c(FMa,PMA)
        ###############################
	###############################
        
      }
    }###end of else loop 
  }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["SMILES"]]))) & !startsWith(stringr::str_trim(as.character(InMEDA[["SMILES"]])),'not available')){
  ##########################################
    print("enter smiles area --2")
    #######################################
    IK<-stringr::str_trim(as.character(InMEDA[["SMILES"]]))
    IK1<-tryCatch({webchem::get_cid(IK, from = "smiles")},error=function(cond){message("smiles not abe to fetch")})
    PCID<-tryCatch({IK1[[2]][1]},error=function(cond){message("Pubchem Id is empty")})
    ########################### changing the smiles part to get exact mass
    tes<-tryCatch({rcdk::parse.smiles(IK)},error=function(cond){message("smiles not abe to parse")})
    tes1<-tryCatch({tes[[1]]},error=function(cond){message("smiles parse information is empty")})
    tes2<-tryCatch({rcdk::get.exact.mass(tes1)},error=function(cond){message("smiles not abe to fetch")})
    PCID2<-tryCatch({tes2},error=function(cond){message("smiles not abe to fetch")})
    Fval<-tryCatch({PuSmilesToEM(InMEDA[["SMILES"]])},error=function(cond){return(0)})
    #######################################
    #######################################
    if(!sjmisc::is_empty(PCID2)){
      print("enter the if loop ..smiles")
      
      FMa<-c(FMa,PCID2)
      
    }else if(!sjmisc::is_empty(Fval) & Fval != 0){
      
      FMa<-c(FMa,Fval)
      
    }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(InMEDA[["PubChem CID"]]),'not available')){
      ############################################
      print("enter the pubchem CID ..smiles area ..meaning smiles are there and not able to get exact mass..pubchem CID is avilable")
      #############################################
      FPUCID<-stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))
      FPUCID1<-tryCatch({as.numeric(FPUCID)},error=function(cond){return(0)})
      ###########################################
      ###########################################
      PCID1<-tryCatch({webchem::pc_prop(as.numeric(FPUCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem Id is empty")})
      PCID2<-ifelse(!sjmisc::is_empty(as.numeric(PCID1$ExactMass)),as.numeric(PCID1$ExactMass),ifelse(!sjmisc::is_empty(tryCatch({as.numeric(ConvCIDtoOID1(FPUCID1)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")})),tryCatch({as.numeric(ConvCIDtoOID1(FPUCID1)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}),tryCatch({PuCIDtoEM(FPUCID1)},error=function(cond){message("Pubchem CID is empty..did not get exact mass1")})))
      ###########################################
      ########################################
      if(!sjmisc::is_empty(PCID2) & (PCID2 !=0)){
        
        FMa<-c(FMa,PCID2)
        
      }else if(tryCatch({as.numeric(ConvCIDtoOID1(FPUCID)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}) != 0){
        FMa<-c(FMa,as.numeric(ConvCIDtoOID1(FPUCID)[5]))
       }else if(tryCatch({as.numeric(PuCIDtoEM(FPUCID))},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}) !=0){
         FMa<-c(FMa,as.numeric(PuCIDtoEM(FPUCID)))
       }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["Formula"]]))) & !startsWith(stringr::str_trim(as.character(InMEDA[["Formula"]])),'not available')){
        print("enter the if pubchem CID ..else if..")
        ##################################
        ##################################
        EM<-tryCatch({Rdisop::getMolecule(stringr::str_trim(as.character(InMEDA[["Formula"]])))},error=function(cond){message("Pubchem Id is empty")})
        EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
        ##################################
        ##################################
        if(!sjmisc::is_empty(EM1)){
          
          FMa<-c(FMa,EM1)
          
        }else{
          print("enter the pubchem CID area ..else part in else if..that means formula is empty")
          
          ###################################
          ###################################
          PMA<-tryCatch({ConvPCIDtoOCN(FPUCID)},error=function(cond){message("Pubchem value is empty")})
          FMa<-c(FMa,ifelse(!sjmisc::is_empty(tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")})),tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")}),FuFtoRe1(InMEDA)))
          ###################################
          ###################################
        }
        #####################################
        #####################################
      }else{
        ### formula not aviable ...else for ..else if else
        print("enter the pubchem CID area...else ...that means neither formula ...nothing is avilable for Pubchem CID")
        ######################################
        #####################################
        PMA<-tryCatch({ConvPCIDtoOCN(FPUCID)},error=function(cond){message("Pubchem value is empty")})
	#####################################
	######################################
        FMa<-c(FMa,ifelse(!sjmisc::is_empty(tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")})),tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")}),FuFtoRe1(InMEDA)))
        ######################################
        ######################################
      }
      #########################################
      
    }else{
      #### Pubchem CID #########################
      print("entering the else part---smiles ... will check formula...for exact mass")
      ##########################################
      if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["Formula"]]))) & !startsWith(stringr::str_trim(as.character(InMEDA[["Formula"]])),'not available')){
        ##################################
        EM<-tryCatch({Rdisop::getMolecule(stringr::str_trim(as.character(InMEDA[["Formula"]])))},error=function(cond){message("Pubchem Id is empty")})
        EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
        ##################################
        
        if(!sjmisc::is_empty(EM1)){

          FMa<-c(FMa,EM1)

        }else{
          
          ##########################
          ##PMA<-FuFtoRe(InMEDA)
          PMA<-tryCatch({FuFtoRe1(InMEDA)},error=function(cond){message("FuFtoRe1 value is empty")})	
          FMa<-c(FMa,PMA)
          ##########################
        }
        
        
      }else{
        ##############################
	##############################      
        PMA<-tryCatch({FuFtoRe1(InMEDA)},error=function(cond){message("FuFtoRe1 value is empty")})      
        FMa<-c(FMa,PMA)
        ##############################
	##############################
      }
      ###########################################
    }
 ################################################
  }else if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))) & !startsWith(as.character(InMEDA[["PubChem CID"]]),'not available')){
    ########################################
    print("entering the Pubchem CID area")
    #######################################
    ##FPUCID<-stringr::str_trim(as.character(InMEDA[["PubChem CID"]]))
    #######################################
    FPUCID<-stringr::str_trim(gsub("[[:punct:]]", "",as.character(InMEDA[["PubChem CID"]])))
    FPUCID1<-tryCatch({as.numeric(FPUCID)},error=function(cond){return(0)})
    #########################################
    #######################################
    PCID1<-tryCatch({webchem::pc_prop(as.numeric(FPUCID), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem CID is empty")})
    ########################################
    ########################################
    PCID2<-ifelse(!sjmisc::is_empty(as.numeric(PCID1$ExactMass)),as.numeric(PCID1$ExactMass),ifelse(!sjmisc::is_empty(tryCatch({as.numeric(ConvCIDtoOID1(FPUCID1)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")})),tryCatch({as.numeric(ConvCIDtoOID1(FPUCID1)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}),tryCatch({PuCIDtoEM(FPUCID1)},error=function(cond){message("Pubchem CID is empty..did not get exact mass1")})))
    ########################################
    ########################################
    if(!sjmisc::is_empty(PCID2) & (PCID2 !=0)){
      print("entering the if loop...PCID2")
      FMa<-c(FMa,PCID2)
    }else if(tryCatch({as.numeric(ConvCIDtoOID1(FPUCID)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}) != 0){
        FMa<-c(FMa,as.numeric(ConvCIDtoOID1(FPUCID)[5]))
      }else if(tryCatch({as.numeric(PuCIDtoEM(FPUCID))},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}) !=0){
          FMa<-c(FMa,as.numeric(PuCIDtoEM(FPUCID)))
        }else if(!sjmisc::is_empty(PCID2) & !sjmisc::is_empty(tryCatch({CONcidtoEM(InMEDA)},error=function(cond){message("some mistake happened in indexes")}))){
      ####FMa<-c(FMa,as.numeric(as.character(InMEDA[["Exact mass"]])))
          FMa<-c(FMa,tryCatch({CONcidtoEM(InMEDA)},error=function(cond){return(0)}))
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
        
        #####################################
        #####################################      
        PMA<-tryCatch({ConvPCIDtoOCN(FPUCID)},error=function(cond){message("Pubchem value is empty")})
        #########################################
        #########################################
        FMa<-c(FMa,ifelse(!sjmisc::is_empty(tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")})),tryCatch({PMA[4]},error=function(cond){message("Pubchem value is empty")}),FuFtoRe1(InMEDA)))
        #######################################
        #######################################
      }
      
    }else{
      print("entering the else of Pubchem CID")
      print("entering this new....in the else part")
      ##############################
      ##############################
      PMA<-tryCatch({FuFtoRe1(InMEDA)},error=function(cond){message("FuFtoRe1 value is empty")})
      FMa<-c(FMa,PMA)
      ##############################
      ##############################
    }
    ################################
    
  }else{
    ###################################
    print("pubchem CID is  not available ...so checking formula..adding new code")
   #######################################################
   #######################################################   
      if(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["Formula"]]))) & !startsWith(stringr::str_trim(as.character(InMEDA[["Formula"]])),'not available')){
      #######################################################
      EM<-tryCatch({Rdisop::getMolecule(stringr::str_trim(as.character(InMEDA[["Formula"]])))},error=function(cond){message("Pubchem Id is empty")})
      EM1<-tryCatch({EM$exactmass},error=function(cond){message("Pubchem Id is empty")})
      ########################################################
      if(!sjmisc::is_empty(EM1)){
        print("this is entering the formula if loop...")      
        FMa<-c(FMa,EM1)
      }else{
        print("adding this new code here")
        ###############################################################
        ###############################################################
        if(!sjmisc::is_empty(as.character(InMEDA[["Name"]])) & !startsWith(as.character(InMEDA[["Name"]]),'not available')){
		print("entering the name pass in the if loop")
          GCID<-tryCatch({webchem::get_cid(stringr::str_trim(InMEDA[["Name"]]))},error=function(cond){message("name to CID value empty")})
          GCID1<-tryCatch({GCID[[2]][1]},error=function(cond){message("Get the CID value is empty")})
          G1CID1<-tryCatch({stringr::str_trim(gsub("[[:punct:]]", "",GCID1))},error=function(cond){message("Get the G1CID1 value is empty")})
	  FPUCID1<-tryCatch({as.numeric(G1CID1)},error=function(cond){message("Get the G1CID1 value is empty")})
          PCID1<-tryCatch({webchem::pc_prop(as.numeric(G1CID1), properties = c("MolecularFormula", "ExactMass","CanonicalSMILES"))},error=function(cond){message("Pubchem CId is empty..did not get exact mass")})
	  PCID2<-ifelse(!sjmisc::is_empty(as.numeric(PCID1$ExactMass)),as.numeric(PCID1$ExactMass),ifelse(!sjmisc::is_empty(tryCatch({as.numeric(ConvCIDtoOID1(FPUCID1)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")})),tryCatch({as.numeric(ConvCIDtoOID1(FPUCID1)[5])},error=function(cond){message("Pubchem CID is empty..did not get exact mass")}),tryCatch({PuCIDtoEM(FPUCID1)},error=function(cond){message("Pubchem CID is empty..did not get exact mass1")})))
	  #############################################################
	  #############################################################
          ########################
          FMa<-c(FMa,PCID2)
          ########################
        }else{
          ### adding this new part here 
          ###########################
          PMA<-tryCatch({FuFtoRe1(InMEDA)},error=function(cond){message("FuFtoRe1 ..something breaks")})	
          FMa<-c(FMa,PMA)
          ##############################
          ##############################
        }### this is the end of else
        
      }###end of else
      print("this is entering the formula else...that means formula not avilable... in final else loop")
      
      ########################################
      ########################################
       PMA<-tryCatch({FuFtoRe1(InMEDA)},error=function(cond){message("FuFtoRe1 ..something breaks1")})
       FMa<-c(FMa,PMA)
      ###########################################
      ###########################################
      
      
      }else{
      ## Formula is not avilable
        #######################################
        #######################################
        PMA<-tryCatch({FuFtoRe1(InMEDA)},error=function(cond){message("FuFtoRe1 ..something breaks2")})
        FMa<-c(FMa,PMA)
        ##########################################
        ##########################################
    }### end of else
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
##########################################################################################################
Ikfilter <- function(InKeyVal,InMSPL,InMEDA,InAdVA,InPMZ,InRTL){
  ######################################
  ##print("checking if problem is in Inchikey area")
  #######################################
  out<-c()
  #######################################
  IINF<-tryCatch({webchem::cts_compinfo(InKeyVal)},error=function(cond){message("webchecm cts could not fetch the info")})
  ######################################
  ######################################
  if(length(IINF)> 0 & !sjmisc::is_empty(IINF)){	  
    #####################################
    #####################################
    IK<-tryCatch({IINF[[1]][[1]]},error=function(cond){message("Inchikey value is empty")})
    IN<-tryCatch({IINF[[1]][[2]]},error=function(cond){message("InchI value is empty")})
    PMZ<-tryCatch({IINF[[1]][[4]]},error=function(cond){message("PrecursorMZ value is empty")})
    FM<-tryCatch({IINF[[1]][[5]]},error=function(cond){message("Formula value is empty")})
    ####################################
    ####################################
    CID<-tryCatch({webchem::get_cid(tryCatch({IINF[[1]][[1]]},error=function(cond){message("Inchikey value is empty..CID value")}), from = "inchikey")},error=function(cond){message("webchecm could not fetch the info from the inchikey")})
    CID1<-tryCatch({CID[[2]]},error=function(cond){return(NA)})
    CID2<-tryCatch({as.character(CID1)},error=function(cond){return(NA)})
    CID3<-tryCatch({gsub("[[:punct:]]", "",CID2 )},error=function(cond){return(NA)})
    CID4<-tryCatch({unlist(strsplit(CID3, " "))},error=function(cond){return(NA)})
    TCID5<-tryCatch({paste(CID4, collapse = ';')},error=function(cond){return(NA)})
    CID5<-ifelse(!sjmisc::is_empty(TCID5),TCID5,tryCatch({webchem::cs_convert(IK,from = "inchikey", to = "csid")},error=function(cond){return(NA)}))
    #######################################
    ########################################
    SM<-ifelse(!sjmisc::is_empty(tryCatch({gETSmiles(InMEDA)},error=function(cond){message("smiles values is empty")})),tryCatch({gETSmiles(InMEDA)},error=function(cond){message("smiles values is empty")}),ifelse(!sjmisc::is_empty(tryCatch({webchem::cs_convert(IN,from="inchi",to="smiles")},error=function(cond){return(NA)})),tryCatch({webchem::cs_convert(IN,from="inchi",to="smiles")},error=function(cond){return(NA)}),tryCatch({getCactus(IN,"smiles")},error=function(cond){return(NA)})))
    ########################################
    ########################################
    SM1<-ifelse(!sjmisc::is_empty(tryCatch({SM},error=function(cond){message("smiles value is empty")})),tryCatch({SM},error=function(cond){message("smile value is empty")}),ifelse(!sjmisc::is_empty(tryCatch({gETSmiles(InMEDA)},error=function(cond){message("smiles values is empty")})),tryCatch({gETSmiles(InMEDA)},error=function(cond){message("smiles values is empty")}),tryCatch({ConvINKtoOID1(IK)[3]},error=function(cond){return(NA)})))
    ##########################################
    ##########################################
    InchiV<-ifelse(!sjmisc::is_empty(IN),tryCatch({IN},error=function(cond){message("Inchi value is empty")}),tryCatch({rinchi::get.inchi(SM1)},error=function(cond){return(NA)}))
    ##########################################
    ##########################################
    InchiV<-ifelse(!sjmisc::is_empty(IN),tryCatch({IN},error=function(cond){message("Inchi value is empty")}),tryCatch({rinchi::get.inchi(SM1)},error=function(cond){return(NA)}))
    ##########################################
    AUIN<-as.character(InMEDA[["Adduct"]])
    ##########################################
    InKeyVal<-ifelse(!sjmisc::is_empty(IK),tryCatch({IK},error=function(cond){message("Inchi value is empty")}),tryCatch({rinchi::get.inchi.key(SM1)},error=function(cond){return(NA)}))
    ###########################################
    IN<-InchiV
    ###########################################
    if(!sjmisc::is_empty(AUIN) & !sjmisc::is_empty(PMZ)){
      ########################################
      print("enter the are ...PMZ(Exact mass) and AUIN(Adduct) is not empty ...in the function Inchikey function")
      ##########################
      AUIN1<-tryCatch({qdapRegex::ex_between(AUIN, "[", "]")[[1]]},error=function(cond){message("Adduct value is missing")})
      AUIN2<-tryCatch({FADINF(AUIN)},error=function(cond){message("adduct value matching is not found")})      
      AAMS<-tryCatch({stringr::str_replace(AUIN2, "M",as.character(PMZ$exactmass))},error=function(cond){message("Missing adduct replacement")})
      AAMS1<-tryCatch({as.numeric(pander::evals(AAMS)[[1]]$result)},error=function(cond){message("Error in adduct replacement step")})
      #########################
      #########################
      PPm=AAMS1*(mz_Tol/(1000000))
      #########################
      MPPmL=AAMS1-PPm
      MPPmU=AAMS1+PPm
      ###########################
      Tmass<-InPMZ[InPMZ >= MPPmL & InPMZ <= MPPmU]
      ITmass<-which(InPMZ %in% Tmass)
      ############################
      VRT<-as.numeric(as.character(InMEDA[["RT (min)"]]))
      ###################################
      VRTL<-VRT-RT_Tol
      VRTU<-VRT+RT_Tol
      ###################################
      TRTL<-InRTL[InRTL >= VRTL & InRTL <= VRTU]
      ITRTL<-which(InRTL %in% TRTL)
      ###################################
      ##ITRTL<-match(TRTL,InRTL)
      ##################################
      ##################################
      ##print(ITmass)
      ###############################################
      ###############################################
      print("enter the line ... 1195")
      if(length(ITRTL) >= 1){
        print("enter the line ... 1197")
        if(length(ITmass) >= 1){
          print("enter the line ... 1198")
          ################################
          INLL<-intersect(ITmass,ITRTL)
          ################################
          if(length(INLL) == 1){
            #####################
            print("enter the line ... 1199")
            ####################
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
            PT1<-c("Num Peaks:","NumPeaks:")
            FNPA<-which(stringi::stri_detect_fixed(FNA,"Num Peaks:"))
            F1NPA<-FNA[FNPA]
            F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:","",F1NPA)))
            #################################
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
            ############################################################
            ############################################################
            if((!sjmisc::is_empty(P1TV2) || !sjmisc::is_empty(as.character(InMEDA[["Adduct"]]))) & (F2NPA != 0)){
              ############################################		  
              ############################################
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
              #################################
              NAT<-as.character(InMEDA[["Adduct"]])
              NAT1<-paste("PRECURSORTYPE:",NAT,sep=" ")
              out<-c(out,NAT1)
              ##################################
              ###################################
              FIN1<-InMEDA[["Ionization mode"]]
              F1IN1<-as.character(FIN1)
              F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
              out<-c(out,F2IN1)
              ##################################################
              ##################################################
              IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
              ###################################################
              ###################################################
              if(!sjmisc::is_empty(IKCRV)){
                ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
                F1ONT<-paste("Ontology:",ONTV,sep=" ")
                out<-c(out,F1ONT)
              }else{
                F1ONT<-tryCatch({MaKE.ONT.REC(InMEDA)[1]},error=function(cond){message("Classifier could not fecth the information")})      
                out<-c(out,F1ONT)
              }              
              ###############################################
              ###############################################
              ### Changing this part inchi ..inchikey and smiles
              if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & !startsWith(as.character(InMEDA[["InChI"]]),'not available') & !startsWith(as.character(InMEDA[["InChI"]]),'CAS:') & !startsWith(as.character(InMEDA[["InChI"]]),'InChI=')){
                if(tryCatch({webchem::is.inchikey(stringr::str_trim(as.character(InMEDA[["InChI"]])))},error=function(cond){message("inchikey..file must be empty")})){
                  FINK<-ifelse(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])),paste("INCHIKEY:",as.character(InMEDA[["InChI"]])),paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" "))
                  FINCH<-paste("INCHI:",InchiV,sep=" ")
                  out<-c(out,FINK)
                  out<-c(out,FINCH)
                }else{
                  FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                  FINCH<-ifelse(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])),paste("INCHI:",as.character(InMEDA[["InChI"]])),paste("INCHI:",InchiV,sep=" "))
                  out<-c(out,FINK)
                  out<-c(out,FINCH)
                }
              }else{
                FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                out<-c(out,FINK)
                FINCH<-paste("INCHI:",InchiV,sep=" ")
                out<-c(out,FINCH)
                
              }
	      ###############################################
              FSIM<-ifelse(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["SMILES"]]))),paste("SMILES:",stringr::str_trim(as.character(InMEDA[["SMILES"]])),sep=" "),paste("SMILES:",SM1,sep=" "))
              ##############################################
              ##############################################
              out<-c(out,FSIM)
              #################################################
              #################################################
              FFOR=ifelse(!sjmisc::is_empty(gETSmiles(InMEDA)),ifelse(!sjmisc::is_empty(tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)})),tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})),ifelse(!sjmisc::is_empty(InMEDA[["InChI"]]),tryCatch({getCactus(InMEDA[["InChI"]], "formula")},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})))
              ########################################
              FFOR1<-paste("FORMULA:",FFOR,sep=" ")
              out<-c(out,FFOR1)
              ########################################
              FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
              FINS1<-FNA[FINS]
              FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
              out<-c(out,FINS2)
              ############################
              FAUT<-as.character(InMEDA[["Authors"]])
              FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
              out<-c(out,FAUT1)
              ##########################
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
              ########################
              FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
              out<-c(out,FCOM)
              ########################
              FNPA<-which(stringi::stri_detect_fixed(FNA,"Num Peaks:"))
              F1NPA<-FNA[FNPA]
              ########################
              PT1<-c("Num Peaks:","NumPeaks:")
              FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
              F1NPA<-FNA[FNPA]
              F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
              Np<-F2NPA
              ###################
              Fpea<-FNA[(FNPA+1):Find]
              Fpea1<-Fpea[!is.na(Fpea)]
              Fpea2<-PeakVali(Fpea1)
              #################################################
              if(F2NPA != 0)
              {
                
                MV=AAMS1
                
                tes1<-unlist(strsplit(Fpea, "\t|\t\t"))
                tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
                tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
                tes4<-which(tes2 > (3+MV))
                NTES4<-which(tes3 <= 100)
                
		NTES5<-which(tes3 <= 0.005*(mean(tes3[order(tes3, decreasing =TRUE)][1:5])))
                
                NTES<-c(NTES4,tes4)
                
                if(Np >= 60)
                {
                  
                  TES<-Centroid(Np,tes2,tes3,tes4,NTES,NTES4,NTES5,DEFAULT_MZ_TOLERANCE,Fpea2)
                  
                  F1NPA<-TES[[1]]
                  tes7<-TES[[2]]       
                  
                  out<-c(out,F1NPA)
                  out<-c(out,tes7)
                }else{
                  #################################
                  F1NPA<-paste0("Num Peaks: ",tryCatch({length(Fpea2)},error=function(cond){message("Fpea is empty")}))
                  out<-c(out,F1NPA)
                  out<-c(out,Fpea2)      
                  #################################      
                  
                }
                
              } ### end of else Fea
              ###########################################################
              ### adding this }
           ############################################################### 
            }else{
              #####################################################
              print("enter the part 1..... P1TV2...Adduct")
              #####################################################
              #########################################
              if((as.character(InMEDA[["Adduct"]]) == "[M]+" || as.character(InMEDA[["Adduct"]]) == "[M]-") & (F2NPA != 0)){
                ######################################	    
                ######################################
                print("enter the part 1..if loop")
                ###############################
                FNA1<-which(stringi::stri_detect_fixed(FNA,"NAME:"))
                FNA2<-InMEDA[["Name"]]
                FNA3<-as.character(FNA2)
                ###########################
                FNAM<-paste("NAME:",FNA3,sep=" ")
                out<-c(out,FNAM)
                ########################
                FRA1<-which(stringi::stri_detect_fixed(FNA,"RETENTIONTIME:"))
                F1RA1<-FNA[FRA1]
                out<-c(out,F1RA1)
                ########################
                FMZ1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORMZ:"))
                F1MZ1<-FNA[FMZ1]
                out<-c(out,F1MZ1)
                ##############################
                NAT<-as.character(InMEDA[["Adduct"]])
                NAT1<-paste("PRECURSORTYPE:",NAT,sep=" ")
                out<-c(out,NAT1)
                #################################
                FIN1<-InMEDA[["Ionization mode"]]
                F1IN1<-as.character(FIN1)
                F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                out<-c(out,F2IN1)
                ###################################################
                ###################################################
                IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
                ###################################################
                ###################################################
                if(!sjmisc::is_empty(IKCRV)){
                  ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
                  F1ONT<-paste("Ontology:",ONTV,sep=" ")
                  out<-c(out,F1ONT)
                }else{
                  F1ONT<-tryCatch({MaKE.ONT.REC(InMEDA)[1]},error=function(cond){message("Classifier could not fecth the information")})
                  out<-c(out,F1ONT)
                }
                
                #############################################
                #############################################
                if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & !startsWith(as.character(InMEDA[["InChI"]]),'not available') & !startsWith(as.character(InMEDA[["InChI"]]),'CAS:') & !startsWith(as.character(InMEDA[["InChI"]]),'InChI=')){
                  if(tryCatch({webchem::is.inchikey(stringr::str_trim(as.character(InMEDA[["InChI"]])))},error=function(cond){message("inchikey..file must be empty")})){
                    FINK<-ifelse(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])),paste("INCHIKEY:",as.character(InMEDA[["InChI"]])),paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" "))
                    FINCH<-paste("INCHI:",InchiV,sep=" ")
                    out<-c(out,FINK)
                    out<-c(out,FINCH)
                  }else{
                    FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                    FINCH<-ifelse(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])),paste("INCHI:",as.character(InMEDA[["InChI"]])),paste("INCHI:",InchiV,sep=" "))
                    out<-c(out,FINK)
                    out<-c(out,FINCH)
                  }
                }else{
                  FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                  out<-c(out,FINK)
                  FINCH<-paste("INCHI:",InchiV,sep=" ")
                  out<-c(out,FINCH)
                }
		##############################################
                FSIM<-ifelse(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["SMILES"]]))),stringr::str_trim(as.character(InMEDA[["SMILES"]])),SM1)
                #############################################
                #############################################
                out<-c(out,FSIM)
                ############################################
                ############################################
                FFOR=ifelse(!sjmisc::is_empty(gETSmiles(InMEDA)),ifelse(!sjmisc::is_empty(tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)})),tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})),ifelse(!sjmisc::is_empty(InMEDA[["InChI"]]),tryCatch({getCactus(InMEDA[["InChI"]], "formula")},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})))
                ###########################################
                ###########################################
                FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                out<-c(out,FFOR1)
                #######################################
                FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                FINS1<-FNA[FINS]
                FINS2<-ifelse(!sjmisc::is_empty(stringr::str_trim(gsub("INTENSITY:","",FINS1))),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
                out<-c(out,FINS2)
                ###############################
                FAUT<-as.character(InMEDA[["Authors"]])
                FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
                out<-c(out,FAUT1)
                #######################
                FLIC<-paste("LICENSE:","CC BY",sep=" ")
                out<-c(out,FLIC)
                #####################
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
                FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
                out<-c(out,FCOM)
                ##############################
                PT1<-c("Num Peaks:","NumPeaks:")
                FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
                F1NPA<-FNA[FNPA]
                F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
                ############################
                ############################
                Np<-F2NPA
                Fpea<-FNA[(FNPA+1):Find]
                Fpea1<-Fpea[!is.na(Fpea)]	      
                Fpea2<-PeakVali(Fpea1)
                ############################################
                if(F2NPA != 0)
                {
                  
                  MV=AAMS1
                  
                  tes1<-unlist(strsplit(Fpea, "\t|\t\t"))
                  tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
                  tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
                  tes4<-which(tes2 > (3+MV))
                  NTES4<-which(tes3 <= 100)
                  NTES<-c(NTES4,tes4)
                  
                  if(Np >= 60)
                  {
                    
                    TES<-Centroid(Np,tes2,tes3,tes4,NTES,NTES4,NTES5,DEFAULT_MZ_TOLERANCE,Fpea2)
                    
                    F1NPA<-TES[[1]]
                    tes7<-TES[[2]]       
                    
                    out<-c(out,F1NPA)
                    out<-c(out,tes7)
                  }else{
                    #################################
                    F1NPA<-paste0("Num Peaks: ",tryCatch({length(Fpea2)},error=function(cond){message("Fpea is empty")}))
                    out<-c(out,F1NPA)
                    out<-c(out,Fpea2)      
                    #################################      
                    
                  }
              ##################################################    
                }### end of the else ...F2NPA    
             ###end of Fea ################################################       
             ##############################################################
              }  ### end of "[M]+"
            ################################################################
           #### end of else ###############################################
            }  
         ###################################################################  
          }else if(length(INLL) > 1){
            #####################################		  
            print("entering the 1569")            
            #####################################
            MONMS=InMSPL[INLL]
            #####################################
            #####################################
            TRA<-unname(rapply(MONMS, function(x) grep("RETENTIONTIME:",x, value=TRUE)))
            TRA1<-as.numeric(stringr::str_trim(stringr::str_replace(TRA, "RETENTIONTIME:", "")))
            TRA2<-abs(VRT-TRA1)
            ##############################
            NPV1=unname(rapply(MONMS, function(x) grep("Num Peaks:",x, value=TRUE)))
            NPV=as.numeric(stringr::str_trim(stringr::str_replace(NPV1, "Num Peaks:", "")))
            ##############################
            NPNZ=which(NPV!=0)
            ##############################
            TRA3<-which(TRA2==min(TRA2[TRA2 > 0]))
            NTRA3<-intersect(TRA3,NPNZ)
            TRA4<-INLL[NTRA3]
            ######################
            TRA5<-InMSPL[TRA4]
            #######################
            F1FPL<-TRA5
            #######################
            ######################
            F2FPL<-F1FPL
            #######################
            #######################
            Find<-tryCatch({length(F2FPL[[1]])},error=function(cond){message("List value is empty")})
            ########################
            FNA<-tryCatch({F1FPL[[1]]},error=function(cond){message("List value is empty")})
            ########################### adding this new
            #############################
            PT1<-c("Num Peaks:","NumPeaks:")
            FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
            F1NPA<-FNA[FNPA]
            F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
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
            ################################################
            if((!sjmisc::is_empty(P1TV2) || !sjmisc::is_empty(as.character(InMEDA[["Adduct"]]))) & (F2NPA != 0)){
              ####################################################
              print("entering the if loop")	    
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
              ##############################
              NAT<-as.character(InMEDA[["Adduct"]])
              NAT1<-paste("PRECURSORTYPE:",NAT,sep=" ")
              out<-c(out,NAT1)
              #############################
              FIN1<-InMEDA[["Ionization mode"]]
              F1IN1<-as.character(FIN1)
              F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
              out<-c(out,F2IN1)
              #############################################
              #############################################	      
              IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
              #############################################
              #############################################
              if(!sjmisc::is_empty(IKCRV)){
                ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
                F1ONT<-paste("Ontology:",ONTV,sep=" ")
                out<-c(out,F1ONT)
              }else{
                F1ONT<-tryCatch({MaKE.ONT.REC(InMEDA)[1]},error=function(cond){message("Classifier could not fecth the information")})
                out<-c(out,F1ONT)
              }
              ##############################################
              ##############################################
              if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & !startsWith(as.character(InMEDA[["InChI"]]),'not available') & !startsWith(as.character(InMEDA[["InChI"]]),'CAS:') & !startsWith(as.character(InMEDA[["InChI"]]),'InChI=')){
                if(tryCatch({webchem::is.inchikey(stringr::str_trim(as.character(InMEDA[["InChI"]])))},error=function(cond){message("inchikey..file must be empty")})){
                  FINK<-ifelse(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])),paste("INCHIKEY:",as.character(InMEDA[["InChI"]])),paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" "))
                  FINCH<-paste("INCHI:",InchiV,sep=" ")
                  out<-c(out,FINK)
                  out<-c(out,FINCH)
                }else{
                  FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                  FINCH<-ifelse(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])),paste("INCHI:",as.character(InMEDA[["InChI"]])),paste("INCHI:",InchiV,sep=" "))
                  out<-c(out,FINK)
                  out<-c(out,FINCH)
                }
              }else{
                FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                out<-c(out,FINK)
                FINCH<-paste("INCHI:",InchiV,sep=" ")
                out<-c(out,FINCH)
              }
              ###############################################
	      ###############################################
	      FSIM<-ifelse(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["SMILES"]]))),paste("SMILES:",stringr::str_trim(as.character(InMEDA[["SMILES"]])),sep=" "),paste("SMILES:",SM1,sep=" ")) 
              ##############################################
              #############################################
              out<-c(out,FSIM)
              #######################################
              #######################################
              FFOR=ifelse(!sjmisc::is_empty(gETSmiles(InMEDA)),ifelse(!sjmisc::is_empty(tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)})),tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})),ifelse(!sjmisc::is_empty(InMEDA[["InChI"]]),tryCatch({getCactus(InMEDA[["InChI"]], "formula")},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})))
              ######################################
              FFOR1<-paste("FORMULA:",FFOR,sep=" ")
              out<-c(out,FFOR1)
              ###############################
              FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
              FINS1<-FNA[FINS]
              FINS2<-ifelse(!sjmisc::is_empty(stringr::str_trim(gsub("INTENSITY:","",FINS1))),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
              ##FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
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
              PT1<-c("Num Peaks:","NumPeaks:")
              FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
              F1NPA<-FNA[FNPA]
              F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
              ############################
              ############################
              Np<-F2NPA
              Fpea<-FNA[(FNPA+1):Find]
              Fpea1<-Fpea[!is.na(Fpea)]	      
              Fpea2<-PeakVali(Fpea1)
              #########################
              if(F2NPA != 0)
              {
                ###############################################
                print("Enter the if loop ...Fpea")
                MV=AAMS1
               
                tes1<-unlist(strsplit(Fpea, "\t|\t\t"))
                tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
                tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
                tes4<-which(tes2 > (3+MV))
                NTES4<-which(tes3 <= 100)
                
                NTES5<-which(tes3 <= 0.005*(mean(tes3[order(tes3, decreasing =TRUE)][1:5])))
                NTES<-c(NTES4,tes4)
                
                if(Np >= 60)
                {
                  
                  TES<-Centroid(Np,tes2,tes3,tes4,NTES,NTES4,NTES5,DEFAULT_MZ_TOLERANCE,Fpea2)
                  
                  F1NPA<-TES[[1]]
                  tes7<-TES[[2]]       
                  
                  out<-c(out,F1NPA)
                  out<-c(out,tes7)
                }else{
                  #################################
                  F1NPA<-paste0("Num Peaks: ",tryCatch({length(Fpea2)},error=function(cond){message("Fpea is empty")}))
                  out<-c(out,F1NPA)
                  out<-c(out,Fpea2)      
                  #################################      
                  
                }
                
            ####################################
              }### End of the else,F2NPA
          ##############################################################
          ##############################################################
            }else{
              ###############################################################
              print("enter the part 2")
              ##############################
              
              #####################################
              if((as.character(InMEDA[["Adduct"]]) == "[M]+" || as.character(InMEDA[["Adduct"]]) == "[M]-") & (F2NPA != 0)){
                
                ###############################
                print("enter the part 2...if loop")
                ##################################
                FNA1<-which(stringi::stri_detect_fixed(FNA,"NAME:"))
                FNA2<-InMEDA[["Name"]]
                FNA3<-as.character(FNA2)
                #######################
                FNAM<-paste("NAME:",FNA3,sep=" ")
                out<-c(out,FNAM)
                #########################
                FRA1<-which(stringi::stri_detect_fixed(FNA,"RETENTIONTIME:"))
                F1RA1<-FNA[FRA1]
                out<-c(out,F1RA1)
                ################################
                FMZ1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORMZ:"))
                F1MZ1<-FNA[FMZ1]
                out<-c(out,F1MZ1)
                ################################
                NAT<-as.character(InMEDA[["Adduct"]])
                NAT1<-paste("PRECURSORTYPE:",NAT,sep=" ")
                out<-c(out,NAT1)
                #################################
                FIN1<-InMEDA[["Ionization mode"]]
                F1IN1<-as.character(FIN1)
                F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                out<-c(out,F2IN1)
                ###################################################
                ###################################################
                IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
                ###################################################
                ###################################################
                if(!sjmisc::is_empty(IKCRV)){
                  ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
                  F1ONT<-paste("Ontology:",ONTV,sep=" ")
                  out<-c(out,F1ONT)
                }else{
                  F1ONT<-tryCatch({MaKE.ONT.REC(InMEDA)[1]},error=function(cond){message("Classifier could not fecth the information")})
                  out<-c(out,F1ONT)
                }
                #####################################################
                #####################################################
                if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & !startsWith(as.character(InMEDA[["InChI"]]),'not available') & !startsWith(as.character(InMEDA[["InChI"]]),'CAS:') & !startsWith(as.character(InMEDA[["InChI"]]),'InChI=')){
                  if(tryCatch({webchem::is.inchikey(stringr::str_trim(as.character(InMEDA[["InChI"]])))},error=function(cond){message("inchikey..file must be empty")})){
                    FINK<-ifelse(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])),paste("INCHIKEY:",as.character(InMEDA[["InChI"]])),paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" "))
                    FINCH<-paste("INCHI:",InchiV,sep=" ")
                    out<-c(out,FINK)
                    out<-c(out,FINCH)
                  }else{
                    FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                    FINCH<-ifelse(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])),paste("INCHI:",as.character(InMEDA[["InChI"]])),paste("INCHI:",InchiV,sep=" "))
                    out<-c(out,FINK)
                    out<-c(out,FINCH)
                  }
                }else{
                  FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                  out<-c(out,FINK)
                  FINCH<-paste("INCHI:",InchiV,sep=" ")
                  out<-c(out,FINCH)
                }
	       #######################################	
                FSIM<-ifelse(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["SMILES"]]))),paste("SMILES:",stringr::str_trim(as.character(InMEDA[["SMILES"]])),sep=" "),paste("SMILES:",SM1,sep=" "))
                out<-c(out,FSIM)
                ####################################
                ####################################
		 FFOR=ifelse(!sjmisc::is_empty(gETSmiles(InMEDA)),ifelse(!sjmisc::is_empty(tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)})),tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})),ifelse(!sjmisc::is_empty(InMEDA[["InChI"]]),tryCatch({getCactus(InMEDA[["InChI"]], "formula")},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})))
                FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                out<-c(out,FFOR1)
                #############################
                FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                FINS1<-FNA[FINS]
                FINS2<-ifelse(!sjmisc::is_empty(stringr::str_trim(gsub("INTENSITY:","",FINS1))),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
                ##FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
                out<-c(out,FINS2)
                ############################
                FAUT<-as.character(InMEDA[["Authors"]])
                FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
                out<-c(out,FAUT1)
                ##########################
                FLIC<-paste("LICENSE:","CC BY",sep=" ")
                out<-c(out,FLIC)
                ###########################
                FCIE<-paste("COLLISIONENERGY:",as.character(InMEDA[["Collision energy"]]),sep=" ")
                out<-c(out,FCIE)
                ############################
                FINST<-as.character(InMEDA[["INSTRUMENT_TYPE"]])
                FINST1<-paste("INSTRUMENTTYPE:",FINST,sep=" ")
                out<-c(out,FINST1)
                ########################
                FINS<-as.character(InMEDA[["INSTRUMENT"]])
                FINS1<-paste("INSTRUMENT:",FINS,sep=" ")
                out<-c(out,FINS1)
                ########################
                PT1<-c("Num Peaks:","NumPeaks:")
                FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
                F1NPA<-FNA[FNPA]
                F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
                ############################
                ############################
                Np<-F2NPA
                Fpea<-FNA[(FNPA+1):Find]
                Fpea1<-Fpea[!is.na(Fpea)]	      
                Fpea2<-PeakVali(Fpea1)
              ######################################
                if(F2NPA != 0)
              ######################################
                {
                  
                  MV=AAMS1
                  tes1<-unlist(strsplit(Fpea, "\t|\t\t"))
                  tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
                  tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
                  tes4<-which(tes2 > (3+MV))
                  NTES4<-which(tes3 <= 100)
                  NTES<-c(NTES4,tes4)
                  
		  NTES5<-which(tes3 <= 0.005*(mean(tes3[order(tes3, decreasing =TRUE)][1:5])))
                  
                  if(Np >= 60)
                  {
                    
                    TES<-Centroid(Np,tes2,tes3,tes4,NTES,NTES4,NTES5,DEFAULT_MZ_TOLERANCE,Fpea2)
                    
                    F1NPA<-TES[[1]]
                    tes7<-TES[[2]]       
                    
                    out<-c(out,F1NPA)
                    out<-c(out,tes7)
                  }else{
                    #################################
                    F1NPA<-paste0("Num Peaks: ",tryCatch({length(Fpea2)},error=function(cond){message("Fpea is empty")}))
                    out<-c(out,F1NPA)
                    out<-c(out,Fpea2)      
                    #################################      
                  }
              ########################################    
                }## End of the main else
                ##########################################################
                ##}
          #############################################################
              }### end of [M]+ and [M]-
          ############################################################
            } 
            ###############################################################
            ### entering the
            ## }    
            ##################################################################  
          }else{
            #print("entering the line 750")
            PASS<-RRV
          }
          ####### This is testing , if this works
          return(out)
          ######################
        } ### this is mz closing brace
      } ## this is RT closing braces
    }else{
      
      #####################################################
      ########################################
      NCIDV<-tryCatch({webchem::get_cid(InKeyVal, from = "inchikey")},error=function(cond){message("webchecm could not fetch the info")})
      NCIDV1<-tryCatch({as.numeric(NCIDV[[2]][1])},error=function(cond){message("webchecm could not fetch the info")})
      ########################################
      ########################################
      FMWFS<-tryCatch({webchem::pc_prop(NCIDV1, properties = c("MolecularFormula", "ExactMass","CanonicalSMILES","InChI","InChIKey"))},error=function(cond){message("webchem conversion from CID is not sucessful")})
      FMWFS1<-tryCatch({FMWFS$ExactMass},error=function(cond){message("CID to exact mass is empty")})
      SM<-tryCatch({FMWFS$CanonicalSMILES},error=function(cond){message("CID to canonical smiles values is empty")})
      IN<-tryCatch({FMWFS$InChI},error=function(cond){message("CID to inchi value is empty")})
      IK<-tryCatch({FMWFS$InChIKey},error=function(cond){message("CID to inchikey is empty")})
      #######################################
      ######################################
      if(!sjmisc::is_empty(AUIN) & !sjmisc::is_empty(FMWFS1)){
        ############################################################
        print("enter the line 325")
        ############################################################
        AUIN1<-tryCatch({qdapRegex::ex_between(AUIN, "[", "]")[[1]]},error=function(cond){message("Adduct value is missing")})
        AUIN2<-tryCatch({FADINF(AUIN)},error=function(cond){message("adduct value matching is not found")})
        ##AUIN2<-tryCatch({InAdVA[InAdVA$V1==AUIN1,]$V8},warning=function(cond){message("Adduct value is missing")})
        AAMS<-tryCatch({stringr::str_replace(AUIN2, "M",as.character(FMWFS1))},error=function(cond){message("Missing adduct replacement")})
        AAMS1<-tryCatch({as.numeric(pander::evals(AAMS)[[1]]$result)},error=function(cond){message("Error in adduct replacement step")})
        ###################################
        ##PPm=AAMS1*(25/(1000000))
        PPm=AAMS1*(mz_Tol/(1000000))
        #########################
        MPPmL=AAMS1-PPm
        MPPmU=AAMS1+PPm
        ###########################
        Tmass<-InPMZ[InPMZ >= MPPmL & InPMZ <= MPPmU]
        ITmass<-which(InPMZ %in% Tmass)
        #############################
        ##ITmass<-match(Tmass,InPMZ)
        ###########################
        VRT<-as.numeric(as.character(InMEDA[["RT (min)"]]))
        VRTL<-VRT-RT_Tol
        VRTU<-VRT+RT_Tol
        #############################
        #############################
        TRTL<-InRTL[InRTL >= VRTL & InRTL <= VRTU]
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
              SM<-as.character(InMEDA[["SMILES"]])
              ##################################
              ##################################
              SM1<-ifelse(!sjmisc::is_empty(tryCatch({SM},error=function(cond){message("smiles value is empty")})),tryCatch({SM},error=function(cond){message("smile value is empty")}),ifelse(!sjmisc::is_empty(tryCatch({gETSmiles(InMEDA)},error=function(cond){message("smiles values is empty")})), trycatch({gETSmiles(InMEDA)},error=function(cond){message("smiles values is empty")}),"NA"))
              ##############################
              ##print(SM1)
              #############################
              F2FPL<-tryCatch({F1FPL},error=function(cond){message("List value is empty")})
              ######################
              Find<-tryCatch({length(F2FPL[[1]])},error=function(cond){message("List value is empty")})
              ########################
              FNA<-tryCatch({F1FPL[[1]]},error=function(cond){message("there is an error in list")})
              ######################## adding this new ################
              PT1<-c("Num Peaks:","NumPeaks:")
              FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
              F1NPA<-FNA[FNPA]
              F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
              ##########################################################
              PT <- c("PRECURSORTYPE:", "ADDUCTIONNAME:")
              NPT<-tryCatch({grep(paste(PT,collapse="|"), FNA, value=TRUE)},error=function(cond){message("List value is empty")})
              NPT1<-tryCatch({match(NPT,FNA)},error=function(cond){message("List value is empty")})
              PTV <- tryCatch({stringr::str_remove(FNA[NPT1],c("PRECURSORTYPE:","ADDUCTIONNAME:"))},error=function(cond){message("List value is empty")})
              PTV1<- tryCatch({PTV[1]},error=function(cond){message("List value is empty")})
              PTV2<-tryCatch({stringr::str_trim(PTV1)},error=function(cond){message("List value is empty")})
              P1TV2<-tryCatch({stringr::str_trim(gsub("ADDUCTIONNAME:","",PTV2))},error=function(cond){message("List value is empty")})
              PTV3<-tryCatch({paste("PRECURSORTYPE:",P1TV2)},error=function(cond){message("List value is empty")})
              ###########################################################
              ###########################################################
              print("enter my test...3")
              print(P1TV2)
              print(as.character(InMEDA[["Adduct"]]))
              ##################################################
              if((!sjmisc::is_empty(P1TV2) || !sjmisc::is_empty(as.character(InMEDA[["Adduct"]]))) & (F2NPA != 0)){
                ################################################
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
                ##########################
                NAT<-as.character(InMEDA[["Adduct"]])
                NAT1<-paste("PRECURSORTYPE:",NAT,sep=" ")
                out<-c(out,NAT1)
                ###########################
                FIN1<-InMEDA[["Ionization mode"]]
                F1IN1<-as.character(FIN1)
                F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                out<-c(out,F2IN1)
                ##################################################
                ##################################################
                IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
                ##################################################
                ##################################################
                if(!sjmisc::is_empty(IKCRV)){
                  ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
                  F1ONT<-paste("Ontology:",ONTV,sep=" ")
                  out<-c(out,F1ONT)
                }else{
                  ##F1ONT<-paste("Ontology:","",sep=" ")
                  F1ONT<-tryCatch({MaKE.ONT.REC(InMEDA)[1]},error=function(cond){message("Classifier could not fecth the information")})	
                  out<-c(out,F1ONT)
                }
                #################################################
                #################################################
                ##################################################
                if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & !startsWith(as.character(InMEDA[["InChI"]]),'not available') & !startsWith(as.character(InMEDA[["InChI"]]),'CAS:') & !startsWith(as.character(InMEDA[["InChI"]]),'InChI=')){
                  if(tryCatch({webchem::is.inchikey(stringr::str_trim(as.character(InMEDA[["InChI"]])))},error=function(cond){message("inchikey..file must be empty")})){
                    FINK<-ifelse(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])),paste("INCHIKEY:",as.character(InMEDA[["InChI"]])),paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" "))
                    FINCH<-paste("INCHI:",IN,sep=" ")
                    out<-c(out,FINK)
                    out<-c(out,FINCH)
                  }else{
                    FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                    FINCH<-ifelse(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])),paste("INCHI:",as.character(InMEDA[["InChI"]])),paste("INCHI:",IN,sep=" "))
                  }
                }else{
                  FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                  out<-c(out,FINK)
                  FINCH<-paste("INCHI:",IN,sep=" ")
                  out<-c(out,FINCH)
                }
                ####################################################
		FSIM<-ifelse(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["SMILES"]]))),paste("SMILES:",stringr::str_trim(as.character(InMEDA[["SMILES"]])),sep=" "),paste("SMILES:",SM1,sep=" "))
                out<-c(out,FSIM)
                ##################################################
                ###################################################
                FFOR=ifelse(!sjmisc::is_empty(gETSmiles(InMEDA)),ifelse(!sjmisc::is_empty(tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)})),tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})),ifelse(!sjmisc::is_empty(InMEDA[["InChI"]]),tryCatch({getCactus(InMEDA[["InChI"]], "formula")},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})))
                ################################
                ################################
                FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                out<-c(out,FFOR1)
                #######################
                FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                FINS1<-FNA[FINS]
                FINS2<-ifelse(!sjmisc::is_empty(stringr::str_trim(gsub("INTENSITY:","",FINS1))),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
                ###FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
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
                PT1<-c("Num Peaks:","NumPeaks:")
                FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
                F1NPA<-FNA[FNPA]
                F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
                ############################
                ############################
                Np<-F2NPA
                Fpea<-FNA[(FNPA+1):Find]
                Fpea1<-Fpea[!is.na(Fpea)]	      
                Fpea2<-PeakVali(Fpea1)
                ##############################################
                if(F2NPA != 0)
                  ###############################
                {
                  
                  MV=AAMS1
                  tes1<-unlist(strsplit(Fpea, "\t|\t\t"))
                  tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
                  tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
                  tes4<-which(tes2 > (3+MV))
                  NTES4<-which(tes3 <= 100)
                  NTES<-c(NTES4,tes4)
                  
		  NTES5<-which(tes3 <= 0.005*(mean(tes3[order(tes3, decreasing =TRUE)][1:5])))
                  
                  if(Np >= 60)
                  {
                    
                    TES<-Centroid(Np,tes2,tes3,tes4,NTES,NTES4,NTES5,DEFAULT_MZ_TOLERANCE,Fpea2)
                    
                    F1NPA<-TES[[1]]
                    tes7<-TES[[2]]       
                    
                    out<-c(out,F1NPA)
                    out<-c(out,tes7)
                  }else{
                    #################################
                    F1NPA<-paste0("Num Peaks: ",tryCatch({length(Fpea2)},error=function(cond){message("Fpea is empty")}))
                    out<-c(out,F1NPA)
                    out<-c(out,Fpea2)      
                    #################################      
                    
                  }
                  
                } ## end of the else
            #############################################################################    
            #############################################################################
              }else{
                #########################################
                print("enter the part 3")
                #########################################
                if((as.character(InMEDA[["Adduct"]]) == "[M]+" || as.character(InMEDA[["Adduct"]]) == "[M]-") & (F2NPA != 0)){
                  #######################################
                  print("enter the part 3..if loop")
                  #######################################
                  FNA1<-which(stringi::stri_detect_fixed(FNA,"NAME:"))
                  FNA2<-InMEDA[["Name"]]
                  FNA3<-as.character(FNA2)
                  #######################
                  FNAM<-paste("NAME:",FNA3,sep=" ")
                  out<-c(out,FNAM)
                  #######################
                  FRA1<-which(stringi::stri_detect_fixed(FNA,"RETENTIONTIME:"))
                  F1RA1<-FNA[FRA1]
                  out<-c(out,F1RA1)
                  ################################
                  FMZ1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORMZ:"))
                  F1MZ1<-FNA[FMZ1]
                  out<-c(out,F1MZ1)
                  #################################
                  NAT<-as.character(InMEDA[["Adduct"]])
                  NAT1<-paste("PRECURSORTYPE:",NAT,sep=" ")
                  out<-c(out,NAT1)
                  #################################
                  FIN1<-InMEDA[["Ionization mode"]]
                  F1IN1<-as.character(FIN1)
                  F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                  out<-c(out,F2IN1)
                  #################################################
                  #################################################
                  IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
                  #################################################
                  #################################################
                  if(!sjmisc::is_empty(IKCRV)){
                    ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
                    F1ONT<-paste("Ontology:",ONTV,sep=" ")
                    out<-c(out,F1ONT)
                  }else{
                    F1ONT<-tryCatch({MaKE.ONT.REC(InMEDA)[1]},error=function(cond){message("Classifier could not fecth the information")})	  
                    ##F1ONT<-paste("Ontology:","",sep=" ")
                    out<-c(out,F1ONT)
                  }
                  
                  #############################################
                  #############################################
                  ##print("enter the line ...1878")
                  ##############################################
                  ###############################################
                  if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & !startsWith(as.character(InMEDA[["InChI"]]),'not available') & !startsWith(as.character(InMEDA[["InChI"]]),'CAS:') & !startsWith(as.character(InMEDA[["InChI"]]),'InChI=')){
                    if(tryCatch({webchem::is.inchikey(stringr::str_trim(as.character(InMEDA[["InChI"]])))},error=function(cond){message("inchikey..file must be empty")})){
                      FINK<-ifelse(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])),paste("INCHIKEY:",as.character(InMEDA[["InChI"]])),paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" "))
                      FINCH<-paste("INCHI:",IN,sep=" ")
                      out<-c(out,FINK)
                      out<-c(out,FINCH)
                    }else{
                      FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                      FINCH<-ifelse(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])),paste("INCHI:",as.character(InMEDA[["InChI"]])),paste("INCHI:",IN,sep=" "))
                      out<-c(out,FINK)
                      out<-c(out,FINCH)
                    }
                  }else{
                    FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                    out<-c(out,FINK)
                    FINCH<-paste("INCHI:",IN,sep=" ")
                    out<-c(out,FINCH)
                  }
                 ##################################################
                  FSIM<-ifelse(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["SMILES"]]))),paste("SMILES:",stringr::str_trim(as.character(InMEDA[["SMILES"]])),sep=" "),paste("SMILES:",SM1,sep=" "))
                  
                  out<-c(out,FSIM)
                  #############################################
                  #############################################
                  FFOR=ifelse(!sjmisc::is_empty(gETSmiles(InMEDA)),ifelse(!sjmisc::is_empty(tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)})),tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})),ifelse(!sjmisc::is_empty(InMEDA[["InChI"]]),tryCatch({getCactus(InMEDA[["InChI"]], "formula")},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})))
                  ######################################
                  FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                  out<-c(out,FFOR1)
                  #############################
                  FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                  FINS1<-FNA[FINS]
                  FINS2<-ifelse(!sjmisc::is_empty(stringr::str_trim(gsub("INTENSITY:","",FINS1))),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
                  ###FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
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
                  FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
                  out<-c(out,FCOM)
                  ##################
                  PT1<-c("Num Peaks:","NumPeaks:")
                  FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
                  F1NPA<-FNA[FNPA]
                  F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
                  ############################
                  ############################
                  Np<-F2NPA
                  Fpea<-FNA[(FNPA+1):Find]
                  Fpea1<-Fpea[!is.na(Fpea)]	      
                  Fpea2<-PeakVali(Fpea1)
                  ###########################
                  if(F2NPA != 0)
                  {
                    
                    MV=AAMS1
                    tes1<-unlist(strsplit(Fpea, "\t|\t\t"))
                    tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
                    tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
                    tes4<-which(tes2 > (3+MV))
                    NTES4<-which(tes3 <= 100)
                    NTES<-c(NTES4,tes4)
                    
		    NTES5<-which(tes3 <= 0.005*(mean(tes3[order(tes3, decreasing =TRUE)][1:5])))
                    
                    if(Np >= 60)
                    {
                      
                      TES<-Centroid(Np,tes2,tes3,tes4,NTES,NTES4,NTES5,DEFAULT_MZ_TOLERANCE,Fpea2)
                      
                      F1NPA<-TES[[1]]
                      tes7<-TES[[2]]       
                      
                      out<-c(out,F1NPA)
                      out<-c(out,tes7)
                    }else{
                      #################################
                      F1NPA<-paste0("Num Peaks: ",tryCatch({length(Fpea2)},error=function(cond){message("Fpea is empty")}))
                      out<-c(out,F1NPA)
                      out<-c(out,Fpea2)      
                      #################################      
                      
                    }
                  ###############  
                  }### end of the else
              ############################################
                }
            ##############################################
            ##############################################
              }
          ############################################################################################  
            } else if(length(INLL) > 1){
              ###################################		    
              print("entering the 1569")
              ###################################
              MONMS=InMSPL[INLL]
              #####################################
              #########################
              TRA<-unname(rapply(MONMS, function(x) grep("RETENTIONTIME:",x, value=TRUE)))
              TRA1<-as.numeric(stringr::str_trim(stringr::str_replace(TRA, "RETENTIONTIME:", "")))
              TRA2<-abs(VRT-TRA1)
              #########TRA3<-which.min(TRA2)
              TRA3<-which(TRA2==min(TRA2[TRA2 > 0]))
              ###############################
              NPV1=unname(rapply(MONMS, function(x) grep("Num Peaks:",x, value=TRUE)))
              NPV=as.numeric(stringr::str_trim(stringr::str_replace(NPV1, "Num Peaks:", "")))
              ##################################
              NPNZ=which(NPV!=0)
              NTRA3<-intersect(TRA3,NPNZ)
              TRA4<-INLL[NTRA3]
              ##########################
              TRA5<-InMSPL[TRA4]
              #######################
              F1FPL<-TRA5
              #######################################
              #######################################
              SM<-as.character(InMEDA[["SMILES"]])
              #######################################
              #######################################
              SM1<-ifelse(!sjmisc::is_empty(tryCatch({SM},error=function(cond){message("smiles value is empty")})),tryCatch({SM},error=function(cond){message("smile value is empty")}),ifelse(!sjmisc::is_empty(tryCatch({gETSmiles(InMEDA)},error=function(cond){message("smiles values is empty")})),tryCatch({gETSmiles(InMEDA)},error=function(cond){message("smiles values is empty")}),"NA"))
              ###################################
              ###################################
              F2FPL<-F1FPL
              ###################################
              Find<-tryCatch({length(F2FPL[[1]])},error=function(cond){message("List value is empty")})
              #####################################
              FNA<-tryCatch({F1FPL[[1]]},error=function(cond){message("List value is empty")})
              ########################### adding this new
              PT1<-c("Num Peaks:","NumPeaks:")
              FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
              F1NPA<-FNA[FNPA]
              F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
              ###############################
              ################################
              Fpea<-FNA[(FNPA+1):Find]
              Fpea1<-Fpea[!is.na(Fpea)]
              Fpea2<-PeakVali(Fpea1)
              Np<-F2NPA
	      ############################################
              ############################################
              PT <- c("PRECURSORTYPE:", "ADDUCTIONNAME:")
              NPT<-tryCatch({grep(paste(PT,collapse="|"), FNA, value=TRUE)},error=function(cond){message("List value is empty")})
              NPT1<-tryCatch({match(NPT,FNA)},error=function(cond){message("List value is empty")})
              PTV <-tryCatch({stringr::str_remove(FNA[NPT1],c("PRECURSORTYPE:","ADDUCTIONNAME:"))},error=function(cond){message("List value is empty")})
              PTV1<-tryCatch({PTV[1]},error=function(cond){message("List value is empty")})
              PTV2<-tryCatch({stringr::str_trim(PTV1)},error=function(cond){message("List value is empty")})
              P1TV2<-tryCatch({stringr::str_trim(gsub("ADDUCTIONNAME:","",PTV2))},error=function(cond){message("List value is empty")})
              PTV3<-tryCatch({paste("PRECURSORTYPE:",P1TV2)},error=function(cond){message("List value is empty")})
              ##PTV3<-tryCatch({paste("PRECURSORTYPE:",PTV2)},error=function(cond){message("List value is empty")})
              ####################################
              print("enter my test...4")
              print(P1TV2)
              print(as.character(InMEDA[["Adduct"]]))
              ######################## adding this new
              if((!sjmisc::is_empty(P1TV2) || !sjmisc::is_empty(as.character(InMEDA[["Adduct"]]))) & (F2NPA != 0)){
                ###################################
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
                NAT<-as.character(InMEDA[["Adduct"]])
                NAT1<-paste("PRECURSORTYPE:",NAT,sep=" ")
                out<-c(out,NAT1)
                #############################
                FIN1<-InMEDA[["Ionization mode"]]
                F1IN1<-as.character(FIN1)
                F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                out<-c(out,F2IN1)
                ###################################
                IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
                ###################################
                ###################################
                if(!sjmisc::is_empty(IKCRV)){
                  ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
                  F1ONT<-paste("Ontology:",ONTV,sep=" ")
                  out<-c(out,F1ONT)
                }else{
                  F1ONT<-tryCatch({MaKE.ONT.REC(InMEDA)[1]},error=function(cond){message("Classifier could not fecth the information")})	
                  ##F1ONT<-paste("Ontology:","",sep=" ")
                  out<-c(out,F1ONT)
                }
                ################################################
                #####################################
                if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & !startsWith(as.character(InMEDA[["InChI"]]),'not available') & !startsWith(as.character(InMEDA[["InChI"]]),'CAS:') & !startsWith(as.character(InMEDA[["InChI"]]),'InChI=')){
                  if(tryCatch({webchem::is.inchikey(stringr::str_trim(as.character(InMEDA[["InChI"]])))},error=function(cond){message("inchikey..file must be empty")})){
                    FINK<-ifelse(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])),paste("INCHIKEY:",as.character(InMEDA[["InChI"]])),paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" "))
                    FINCH<-paste("INCHI:",IN,sep=" ")
                    out<-c(out,FINK)
                    out<-c(out,FINCH)
                  }else{
                    FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                    FINCH<-ifelse(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])),paste("INCHI:",as.character(InMEDA[["InChI"]])),paste("INCHI:",IN,sep=" "))
                    out<-c(out,FINK)
                    out<-c(out,FINCH)
                  }
                }else{
                  FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                  out<-c(out,FINK)
                  FINCH<-paste("INCHI:",IN,sep=" ")
                  out<-c(out,FINCH)
                }
                #####################################
		FSIM<-ifelse(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["SMILES"]]))),paste("SMILES:",stringr::str_trim(as.character(InMEDA[["SMILES"]])),sep=" "),paste("SMILES:",SM1,sep=" "))
                
                out<-c(out,FSIM)
                ##################################
                print("enter the smiles value")
                print(FSIM)
                print(stringr::str_trim(as.character(InMEDA[["SMILES"]])))
                ####################################
                FFOR=ifelse(!sjmisc::is_empty(gETSmiles(InMEDA)),ifelse(!sjmisc::is_empty(tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)})),tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})),ifelse(!sjmisc::is_empty(InMEDA[["InChI"]]),tryCatch({getCactus(InMEDA[["InChI"]], "formula")},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})))
                ##################################
                FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                out<-c(out,FFOR1)
                #################################
                FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                FINS1<-FNA[FINS]
                FINS2<-ifelse(!sjmisc::is_empty(stringr::str_trim(gsub("INTENSITY:","",FINS1))),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
                out<-c(out,FINS2)
                ##############################
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
                PT1<-c("Num Peaks:","NumPeaks:")
                FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
                F1NPA<-FNA[FNPA]
                F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
                ############################
                ############################
                Np<-F2NPA
                Fpea<-FNA[(FNPA+1):Find]
                Fpea1<-Fpea[!is.na(Fpea)]	      
                Fpea2<-PeakVali(Fpea1)
                #########################
                if(F2NPA != 0)
              #########################
                {
                  
                  MV=AAMS1
                  tes1<-unlist(strsplit(Fpea, "\t|\t\t"))
                  tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
                  tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
                  tes4<-which(tes2 > (3+MV))
                  NTES4<-which(tes3 <= 100)
                  
		  NTES5<-which(tes3 <= 0.005*(mean(tes3[order(tes3, decreasing =TRUE)][1:5])))
                  NTES<-c(NTES4,tes4)
                  
                  if(Np >= 60)
                    #####if(length(tes4)>1)
                  {
                    
                    TES<-Centroid(Np,tes2,tes3,tes4,NTES,NTES4,NTES5,DEFAULT_MZ_TOLERANCE,Fpea2)
                    
                    F1NPA<-TES[[1]]
                    tes7<-TES[[2]]       
                    
                    out<-c(out,F1NPA)
                    out<-c(out,tes7)
                  }else{
                    #################################
                    F1NPA<-paste0("Num Peaks: ",tryCatch({length(Fpea2)},error=function(cond){message("Fpea is empty")}))
                    out<-c(out,F1NPA)
                    out<-c(out,Fpea2)      
                    #################################      
                    
                  }
                 
                  
                  
                }#### end of the if and else 
            ###############################################
              }else{
            ########################################
                
                ####################################
                if((as.character(InMEDA[["Adduct"]]) == "[M]+" || as.character(InMEDA[["Adduct"]]) == "[M]-") & (F2NPA != 0)){
                  
                  ############################
                  print("enter the part 4..if loop")
                  #############################
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
		  #################################
		  NAT<-as.character(InMEDA[["Adduct"]])
                  NAT1<-paste("PRECURSORTYPE:",NAT,sep=" ")
                  out<-c(out,NAT1)
                  #################################
                  FIN1<-InMEDA[["Ionization mode"]]
                  F1IN1<-as.character(FIN1)
                  F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                  out<-c(out,F2IN1)
                  ##################################################
                  
                  IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
                  ###################################################
                  #################################################
                  if(!sjmisc::is_empty(IKCRV)){
                    ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
                    F1ONT<-paste("Ontology:",ONTV,sep=" ")
                    out<-c(out,F1ONT)
                  }else{
                    ##F1ONT<-paste("Ontology:","",sep=" ")
                    F1ONT<-tryCatch({MaKE.ONT.REC(InMEDA)[1]},error=function(cond){message("Classifier could not fecth the information")})	  
                    out<-c(out,F1ONT)
                  }
                  
                  #############################################
                  ################################################
                  if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & !startsWith(as.character(InMEDA[["InChI"]]),'not available') & !startsWith(as.character(InMEDA[["InChI"]]),'CAS:') & !startsWith(as.character(InMEDA[["InChI"]]),'InChI=')){
                    if(tryCatch({webchem::is.inchikey(stringr::str_trim(as.character(InMEDA[["InChI"]])))},error=function(cond){message("inchikey..file must be empty")})){
                      FINK<-ifelse(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])),paste("INCHIKEY:",as.character(InMEDA[["InChI"]])),paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" "))
                      FINCH<-paste("INCHI:",IN,sep=" ")
                      out<-c(out,FINK)
                      out<-c(out,FINCH)
                    }else{
                      FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                      FINCH<-ifelse(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])),paste("INCHI:",as.character(InMEDA[["InChI"]])),paste("INCHI:",IN,sep=" "))
                      out<-c(out,FINK)
                      out<-c(out,FINCH)
                    }
                  }else{
                    FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                    out<-c(out,FINK)
                    FINCH<-paste("INCHI:",IN,sep=" ")
                    out<-c(out,FINCH)
                  }
		  ###############################################
                  FSIM<-ifelse(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["SMILES"]]))),paste("SMILES:",stringr::str_trim(as.character(InMEDA[["SMILES"]])),sep=" "),paste("SMILES:",SM1,sep=" "))
                  
                  
                  out<-c(out,FSIM)
                  ##############################################
                  ##############################################
                  FFOR=ifelse(!sjmisc::is_empty(gETSmiles(InMEDA)),ifelse(!sjmisc::is_empty(tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)})),tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})),ifelse(!sjmisc::is_empty(InMEDA[["InChI"]]),tryCatch({getCactus(InMEDA[["InChI"]], "formula")},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})))
                  ##################################
                  FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                  out<-c(out,FFOR1)
                  ###############################
                  FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                  FINS1<-FNA[FINS]
                  FINS2<-ifelse(!sjmisc::is_empty(stringr::str_trim(gsub("INTENSITY:","",FINS1))),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
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
                  FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
                  out<-c(out,FCOM)
                  ##################
                  PT1<-c("Num Peaks:","NumPeaks:")
                  FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
                  F1NPA<-FNA[FNPA]
                  F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
                  ############################
                  ############################
                  Np<-F2NPA
                  Fpea<-FNA[(FNPA+1):Find]
                  Fpea1<-Fpea[!is.na(Fpea)]	      
                  Fpea2<-PeakVali(Fpea1)
                  ############################
                  if(F2NPA != 0)
                  {
                    
                    MV=AAMS1
                    ####tes1<-unlist(strsplit(Fpea, "\t"))
                    tes1<-unlist(strsplit(Fpea, "\t|\t\t"))
                    tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
                    tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
                    tes4<-which(tes2 > (3+MV))
                    NTES4<-which(tes3 <= 100)
                    
		    NTES5<-which(tes3 <= 0.005*(mean(tes3[order(tes3, decreasing =TRUE)][1:5])))
                    NTES<-c(NTES4,tes4)
                    if(Np >= 60)
                    {
                      
                      TES<-Centroid(Np,tes2,tes3,tes4,NTES,NTES4,NTES5,DEFAULT_MZ_TOLERANCE,Fpea2)
                      
                      F1NPA<-TES[[1]]
                      tes7<-TES[[2]]       
                      
                      out<-c(out,F1NPA)
                      out<-c(out,tes7)
                    }else{
                      #################################
                      F1NPA<-paste0("Num Peaks: ",tryCatch({length(Fpea2)},error=function(cond){message("Fpea is empty")}))
                      out<-c(out,F1NPA)
                      out<-c(out,Fpea2)      
                      #################################      
                      
                    }
                    
                  }### end of the else
                  ############################
                }             
                #############################
                ##############################
              }
              #### need to add } to  ######## 
              ################################  
            } else{
              #print("entering the line 750")
              PASS<-RRV
            }
            ####### This is testing , if this works
            return(out)
            #############################
          } ### this is mz closing brace
        } ## this is RT closing braces
      #################################
      }else{
        ### THis is the new function I am adding from REST API conversion ..Inchikey to exact value 
        ###########################################
        ###########################################
        FMWFS1<-ifelse(!sjmisc::is_empty(tryCatch({ConvINKtoOID(InKeyVal)[4]},error=function(cond){message("Inchikey conversion is empty")})),tryCatch({ConvINKtoOID(InKeyVal)[4]},error=function(cond){message("iNCHIKEY conversion is empty")}),"NA")
        IN<-tryCatch({ConvINKtoOID1(InKeyVal)[1]},error=function(cond){message("conversion to inchi value is empty")})
        IK<-tryCatch({ConvINKtoOID1(InKeyVal)[2]},error=function(cond){message("conversion to inchikey value is empty")})
        ############################################
        ############################################
        if(!sjmisc::is_empty(AUIN) & !sjmisc::is_empty(FMWFS1)){
          ############################################################
          print("enter the line 325...enter the if loop")
          ############################################################
          AUIN1<-tryCatch({qdapRegex::ex_between(AUIN, "[", "]")[[1]]},error=function(cond){message("Adduct value is missing")})
          AUIN2<-tryCatch({FADINF(AUIN)},error=function(cond){message("adduct value matching is not found")})
          ##AUIN2<-tryCatch({InAdVA[InAdVA$V1==AUIN1,]$V8},warning=function(cond){message("Adduct value is missing")})
          AAMS<-tryCatch({stringr::str_replace(AUIN2, "M",as.character(FMWFS1))},error=function(cond){message("Missing adduct replacement")})
          AAMS1<-tryCatch({as.numeric(pander::evals(AAMS)[[1]]$result)},error=function(cond){message("Error in adduct replacement step")})
          ###################################
          ##PPm=AAMS1*(25/(1000000))
          PPm=AAMS1*(mz_Tol/(1000000))
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
          ##################################
          print("enter temp test ---function 1 ...else part")
          #############################
          TRTL<-InRTL[InRTL >= VRTL & InRTL <= VRTU]
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
                ##SM1<-as.character(InMEDA[["SMILES"]])
                #####################
                F2FPL<-tryCatch({F1FPL},error=function(cond){message("List value is empty")})
                ######################
                Find<-tryCatch({length(F2FPL[[1]])},error=function(cond){message("List value is empty")})
                ########################
                FNA<-tryCatch({F1FPL[[1]]},error=function(cond){message("there is an error in list")})
                ######################## adding this new ################
                PT1<-c("Num Peaks:","NumPeaks:")
                FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
                F1NPA<-FNA[FNPA]
                F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
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
                #########################################################
                if((!sjmisc::is_empty(P1TV2) || !sjmisc::is_empty(as.character(InMEDA[["Adduct"]]))) & (F2NPA != 0)){
                 
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
		  ###############################
		  NAT<-as.character(InMEDA[["Adduct"]])
                  NAT1<-paste("PRECURSORTYPE:",NAT,sep=" ")
                  out<-c(out,NAT1)
                  ##########################
                  ###########################
                  FIN1<-InMEDA[["Ionization mode"]]
                  F1IN1<-as.character(FIN1)
                  F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                  out<-c(out,F2IN1)
                  ##################################################
                  ##################################################
                  IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
                  ####################################################
                  ####################################################
                  if(!sjmisc::is_empty(IKCRV)){
                    ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
                    F1ONT<-paste("Ontology:",ONTV,sep=" ")
                    out<-c(out,F1ONT)
                  }else{
                    F1ONT<-tryCatch({MaKE.ONT.REC(InMEDA)[1]},error=function(cond){message("Classifier could not fecth the information")})
                    out<-c(out,F1ONT)
                  }
                  #################################################
                  ##################################################
                  if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & !startsWith(as.character(InMEDA[["InChI"]]),'not available') & !startsWith(as.character(InMEDA[["InChI"]]),'CAS:') & !startsWith(as.character(InMEDA[["InChI"]]),'InChI=')){
                    if(tryCatch({webchem::is.inchikey(stringr::str_trim(as.character(InMEDA[["InChI"]])))},error=function(cond){message("inchikey..file must be empty")})){
                      FINK<-ifelse(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])),paste("INCHIKEY:",as.character(InMEDA[["InChI"]])),paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" "))
                      FINCH<-paste("INCHI:",IN,sep=" ")
                      out<-c(out,FINK)
                      out<-c(out,FINCH)
                    }else{
                      FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                      FINCH<-ifelse(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])),paste("INCHI:",as.character(InMEDA[["InChI"]])),paste("INCHI:",IN,sep=" "))
                      out<-c(out,FINK)
                      out<-c(out,FINCH)
                    }
                  }else{
                    FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                    out<-c(out,FINK)
                    FINCH<-paste("INCHI:",IN,sep=" ")
                    out<-c(out,FINCH)
                  }
		  ##########################################
                  FSIM<-ifelse(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["SMILES"]]))),paste("SMILES:",stringr::str_trim(as.character(InMEDA[["SMILES"]])),sep=" "),paste("SMILES:",SM1,sep=" "))
                  
                  out<-c(out,FSIM)
                  ##################################################
                  ##################################################
                  FFOR=ifelse(!sjmisc::is_empty(gETSmiles(InMEDA)),ifelse(!sjmisc::is_empty(tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)})),tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})),ifelse(!sjmisc::is_empty(InMEDA[["InChI"]]),tryCatch({getCactus(InMEDA[["InChI"]], "formula")},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})))
                  #####################################
                  FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                  out<-c(out,FFOR1)
                  #############################
                  FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                  FINS1<-FNA[FINS]
                  FINS2<-ifelse(!sjmisc::is_empty(stringr::str_trim(gsub("INTENSITY:","",FINS1))),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
                  out<-c(out,FINS2)
                  ###########################
                  FAUT<-as.character(InMEDA[["Authors"]])
                  FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
                  out<-c(out,FAUT1)
                  ##########################
                  ##FLIC<-paste("LICENSE:",sep=" ")
                  FLIC<-paste("LICENSE:","CC BY",sep=" ")
                  out<-c(out,FLIC)
                  ####################
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
                  PT1<-c("Num Peaks:","NumPeaks:")
                  FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
                  F1NPA<-FNA[FNPA]
                  F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
                  ############################
                  ############################
                  Np<-F2NPA
                  Fpea<-FNA[(FNPA+1):Find]
                  Fpea1<-Fpea[!is.na(Fpea)]	      
                  Fpea2<-PeakVali(Fpea1)
                  ##########################
                  if(F2NPA != 0)
                  ###########################
                  {
                    
                    MV=AAMS1
                    tes1<-unlist(strsplit(Fpea, "\t|\t\t"))
                    tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
                    tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
                    tes4<-which(tes2 > (3+MV))
                    NTES4<-which(tes3 <= 100)
                    
		    NTES5<-which(tes3 <= 0.005*(mean(tes3[order(tes3, decreasing =TRUE)][1:5])))
                    NTES<-c(NTES4,tes4)
                   
                    
                    if(Np >= 60)
                    {
                      
                      TES<-Centroid(Np,tes2,tes3,tes4,NTES,NTES4,NTES5,DEFAULT_MZ_TOLERANCE,Fpea2)
                      
                      F1NPA<-TES[[1]]
                      tes7<-TES[[2]]       
                      
                      out<-c(out,F1NPA)
                      out<-c(out,tes7)
                    }else{
                      #################################
                      F1NPA<-paste0("Num Peaks: ",tryCatch({length(Fpea2)},error=function(cond){message("Fpea is empty")}))
                      out<-c(out,F1NPA)
                      out<-c(out,Fpea2)      
                      #################################      
                      
                    }
                    
                    
                    
                  }
               #############################################################################    
               #############################################################################
                }else{
                  #########################################
                  print("enter the part 3")
                  
                  #########################################
                  if((as.character(InMEDA[["Adduct"]]) == "[M]+" || as.character(InMEDA[["Adduct"]]) == "[M]-") & (F2NPA != 0)){
                    
                    #######################################
                    print("enter the part 3..if loop")
                    #######################################
                    FNA1<-which(stringi::stri_detect_fixed(FNA,"NAME:"))
                    FNA2<-InMEDA[["Name"]]
                    FNA3<-as.character(FNA2)
                    #######################
                    FNAM<-paste("NAME:",FNA3,sep=" ")
                    out<-c(out,FNAM)
                    #######################
                    FRA1<-which(stringi::stri_detect_fixed(FNA,"RETENTIONTIME:"))
                    F1RA1<-FNA[FRA1]
                    out<-c(out,F1RA1)
                    ################################
                    FMZ1<-which(stringi::stri_detect_fixed(FNA,"PRECURSORMZ:"))
                    F1MZ1<-FNA[FMZ1]
                    out<-c(out,F1MZ1)
                    ################################
                    #################################
                    FIN1<-InMEDA[["Ionization mode"]]
                    F1IN1<-as.character(FIN1)
                    F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                    out<-c(out,F2IN1)
                    #################################################
                    #################################################
                    IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
                    #################################################
                    #################################################
                    if(!sjmisc::is_empty(IKCRV)){
                      ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
                      F1ONT<-paste("Ontology:",ONTV,sep=" ")
                      out<-c(out,F1ONT)
                    }else{
                      F1ONT<-tryCatch({MaKE.ONT.REC(InMEDA)[1]},error=function(cond){message("Classifier could not fecth the information")})
                      ##F1ONT<-paste("Ontology:","",sep=" ")
                      out<-c(out,F1ONT)
                    }
                    
                    #############################################
                    ###########################################
                    if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & !startsWith(as.character(InMEDA[["InChI"]]),'not available') & !startsWith(as.character(InMEDA[["InChI"]]),'CAS:') & !startsWith(as.character(InMEDA[["InChI"]]),'InChI=')){
                      if(tryCatch({webchem::is.inchikey(stringr::str_trim(as.character(InMEDA[["InChI"]])))},error=function(cond){message("inchikey..file must be empty")})){
                        FINK<-ifelse(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])),paste("INCHIKEY:",as.character(InMEDA[["InChI"]])),paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" "))
                        FINCH<-paste("INCHI:",IN,sep=" ")
                        out<-c(out,FINK)
                        out<-c(out,FINCH)
                      }else{
                        FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                        FINCH<-ifelse(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])),paste("INCHI:",as.character(InMEDA[["InChI"]])),paste("INCHI:",IN,sep=" "))
                        out<-c(out,FINK)
                        out<-c(out,FINCH)
                      }
                    }else{
                      FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                      out<-c(out,FINK)
                      FINCH<-paste("INCHI:",IN,sep=" ")
                      out<-c(out,FINCH)
                    }
		   ############################################### 
                    FSIM<-ifelse(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["SMILES"]]))),paste("SMILES:",stringr::str_trim(as.character(InMEDA[["SMILES"]])),sep=" "),paste("SMILES:",SM1,sep=" "))
                    
                    out<-c(out,FSIM)
                    ###################################
                    print("enter the smiles value")
                    print(FSIM)
                    print(stringr::str_trim(as.character(InMEDA[["SMILES"]])))
                    #####################################
                    FFOR=ifelse(!sjmisc::is_empty(gETSmiles(InMEDA)),ifelse(!sjmisc::is_empty(tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)})),tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})),ifelse(!sjmisc::is_empty(InMEDA[["InChI"]]),tryCatch({getCactus(InMEDA[["InChI"]], "formula")},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})))
                    FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                    out<-c(out,FFOR1)
                    #################################
                    FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                    FINS1<-FNA[FINS]
                    FINS2<-ifelse(!sjmisc::is_empty(stringr::str_trim(gsub("INTENSITY:","",FINS1))),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
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
                    FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
                    out<-c(out,FCOM)
                    ##################
                    PT1<-c("Num Peaks:","NumPeaks:")
                    FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
                    F1NPA<-FNA[FNPA]
                    F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
                    ############################
                    ############################
                    Np<-F2NPA
                    Fpea<-FNA[(FNPA+1):Find]
                    Fpea1<-Fpea[!is.na(Fpea)]	      
                    Fpea2<-PeakVali(Fpea1)
                    #################################
                    if(F2NPA != 0)
                      ##############################	    
                    {
                      
                      MV=AAMS1
                      tes1<-unlist(strsplit(Fpea, "\t|\t\t"))
                      tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
                      tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
                      tes4<-which(tes2 > (3+MV))
                      NTES4<-which(tes3 <= 100)
                      NTES<-c(NTES4,tes4)
                      
		      NTES5<-which(tes3 <= 0.005*(mean(tes3[order(tes3, decreasing =TRUE)][1:5])))
                      
                      if(Np >= 60)
                      {
                        
                        TES<-Centroid(Np,tes2,tes3,tes4,NTES,NTES4,NTES5,DEFAULT_MZ_TOLERANCE,Fpea2)
                        
                        F1NPA<-TES[[1]]
                        tes7<-TES[[2]]       
                        
                        out<-c(out,F1NPA)
                        out<-c(out,tes7)
                      }else{
                        #################################
                        F1NPA<-paste0("Num Peaks: ",tryCatch({length(Fpea2)},error=function(cond){message("Fpea is empty")}))
                        out<-c(out,F1NPA)
                        out<-c(out,Fpea2)      
                        #################################      
                        
                      }
                      
                    } ## end of the else
                ###########################################
                  }
              ###########################################
              ###########################################
                }
            ############################################################################################  
              }else if(length(INLL) > 1){
                ###################################		    
                print("entering the 1569")
                ###################################
                MONMS=InMSPL[INLL]
                #########################
                TRA<-unname(rapply(MONMS, function(x) grep("RETENTIONTIME:",x, value=TRUE)))
                TRA1<-as.numeric(stringr::str_trim(stringr::str_replace(TRA, "RETENTIONTIME:", "")))
                TRA2<-abs(VRT-TRA1)
                ###TRA3<-which.min(TRA2)
                TRA3<-which(TRA2==min(TRA2[TRA2 > 0]))
                ######################################
                NPV1=unname(rapply(MONMS, function(x) grep("Num Peaks:",x, value=TRUE)))
                NPV=as.numeric(stringr::str_trim(stringr::str_replace(NPV1, "Num Peaks:", "")))
                ######################################
                NPNZ=which(NPV!=0)
                NTRA3<-intersect(TRA3,NPNZ)
                TRA4<-INLL[NTRA3]
                ###############################
                TRA5<-InMSPL[TRA4]
                #######################
                F1FPL<-TRA5
                ###################################
                ##################################
                SM<-as.character(InMEDA[["SMILES"]])
                SM1<-ifelse(!sjmisc::is_empty(tryCatch({SM},error=function(cond){message("smiles value is empty")})),tryCatch({SM},error=function(cond){message("smile value is empty")}),ifelse(sjmisc::is_empty(tryCatch({gETSmiles(InMEDA)},error=function(cond){message("smiles values is empty")})),!sjmisc::is_empty(tryCatch({gETSmiles(InMEDA)},error=function(cond){message("smiles values is empty")})),"NA"))
                #################################
                print("entering the test area ....2")
                print(SM1)
                ####################################
                ##############################
                F2FPL<-F1FPL
                #################################
                Find<-tryCatch({length(F2FPL[[1]])},error=function(cond){message("List value is empty")})
                ########################
                FNA<-tryCatch({F1FPL[[1]]},error=function(cond){message("List value is empty")})
                ########################### adding this new
                PT1<-c("Num Peaks:","NumPeaks:")
                FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
                F1NPA<-FNA[FNPA]
                F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
                ############################################
                PT <- c("PRECURSORTYPE:", "ADDUCTIONNAME:")
                NPT<-tryCatch({grep(paste(PT,collapse="|"), FNA, value=TRUE)},error=function(cond){message("List value is empty")})
                NPT1<-tryCatch({match(NPT,FNA)},error=function(cond){message("List value is empty")})
                PTV <-tryCatch({stringr::str_remove(FNA[NPT1],c("PRECURSORTYPE:","ADDUCTIONNAME:"))},error=function(cond){message("List value is empty")})
                PTV1<-tryCatch({PTV[1]},error=function(cond){message("List value is empty")})
                PTV2<-tryCatch({stringr::str_trim(PTV1)},error=function(cond){message("List value is empty")})
                P1TV2<-tryCatch({stringr::str_trim(gsub("ADDUCTIONNAME:","",PTV2))},error=function(cond){message("List value is empty")})
                PTV3<-tryCatch({paste("PRECURSORTYPE:",P1TV2)},error=function(cond){message("List value is empty")})
                ################################################
                ################################################
                print("enter my test...4")
                print(P1TV2)
                print(as.character(InMEDA[["Adduct"]]))
                ######################## adding this new
                #################################################
                if((!sjmisc::is_empty(P1TV2) || !sjmisc::is_empty(as.character(InMEDA[["Adduct"]]))) & (F2NPA != 0)){
                  
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
		  #################################
		  NAT<-as.character(InMEDA[["Adduct"]])
                  NAT1<-paste("PRECURSORTYPE:",NAT,sep=" ")
                  out<-c(out,NAT1)
                  ################################
                  #############################
                  FIN1<-InMEDA[["Ionization mode"]]
                  F1IN1<-as.character(FIN1)
                  F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                  out<-c(out,F2IN1)
                  ###################################
                  IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
                  ###################################
                  ###################################
                  if(!sjmisc::is_empty(IKCRV)){
                    ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
                    F1ONT<-paste("Ontology:",ONTV,sep=" ")
                    out<-c(out,F1ONT)
                  }else{
                    F1ONT<-tryCatch({MaKE.ONT.REC(InMEDA)[1]},error=function(cond){message("Classifier could not fecth the information")})	  
                    out<-c(out,F1ONT)
                  }
                  ###################################
                  ##################################
                  if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & !startsWith(as.character(InMEDA[["InChI"]]),'not available') & !startsWith(as.character(InMEDA[["InChI"]]),'CAS:') & !startsWith(as.character(InMEDA[["InChI"]]),'InChI=')){
                    if(tryCatch({webchem::is.inchikey(stringr::str_trim(as.character(InMEDA[["InChI"]])))},error=function(cond){message("inchikey..file must be empty")})){
                      FINK<-ifelse(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])),paste("INCHIKEY:",as.character(InMEDA[["InChI"]])),paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" "))
                      FINCH<-paste("INCHI:",IN,sep=" ")
                      out<-c(out,FINK)
                      out<-c(out,FINCH)
                    }else{
                      FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                      FINCH<-ifelse(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])),paste("INCHI:",as.character(InMEDA[["InChI"]])),paste("INCHI:",IN,sep=" "))
                      out<-c(out,FINK)
                      out<-c(out,FINCH)
                    }
                  }else{
                    FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                    out<-c(out,FINK)
                    FINCH<-paste("INCHI:",IN,sep=" ")
                    out<-c(out,FINCH)
                  }
                  #######################################
                  FSIM<-ifelse(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["SMILES"]]))),paste("SMILES:",stringr::str_trim(as.character(InMEDA[["SMILES"]])),sep=" "),paste("SMILES:",SM1,sep=" "))
                  
                  out<-c(out,FSIM)
                  ########################################
                  print("enter the smiles value")
                  print(FSIM)
                  print(stringr::str_trim(as.character(InMEDA[["SMILES"]])))
                  ########################################
                  ########################################
                  FFOR=ifelse(!sjmisc::is_empty(gETSmiles(InMEDA)),ifelse(!sjmisc::is_empty(tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)})),tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})),ifelse(!sjmisc::is_empty(InMEDA[["InChI"]]),tryCatch({getCactus(InMEDA[["InChI"]], "formula")},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})))
                  #########################################
                  #########################################
                  FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                  out<-c(out,FFOR1)
                  ###############################
                  FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                  FINS1<-FNA[FINS]
                  FINS2<-ifelse(!sjmisc::is_empty(stringr::str_trim(gsub("INTENSITY:","",FINS1))),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
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
                  PT1<-c("Num Peaks:","NumPeaks:")
                  FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
                  F1NPA<-FNA[FNPA]
                  F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
                  ############################
                  ############################
                  Np<-F2NPA
                  Fpea<-FNA[(FNPA+1):Find]
                  Fpea1<-Fpea[!is.na(Fpea)]	      
                  Fpea2<-PeakVali(Fpea1)
                  #########################
                  if(F2NPA != 0)
                #############################
                  {
                    MV=AAMS1
                    tes1<-unlist(strsplit(Fpea, "\t|\t\t"))
                    tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
                    tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
                    tes4<-which(tes2 > (3+MV))
                    NTES4<-which(tes3 <= 100)
                    
		    NTES5<-which(tes3 <= 0.005*(mean(tes3[order(tes3, decreasing =TRUE)][1:5])))
                    NTES<-c(NTES4,tes4)
                    
                    if(Np >= 60)
                    {
                      
                      TES<-Centroid(Np,tes2,tes3,tes4,NTES,NTES4,NTES5,DEFAULT_MZ_TOLERANCE,Fpea2)
                      
                      F1NPA<-TES[[1]]
                      tes7<-TES[[2]]       
                      
                      out<-c(out,F1NPA)
                      out<-c(out,tes7)
                    }else{
                      #################################
                      F1NPA<-paste0("Num Peaks: ",tryCatch({length(Fpea2)},error=function(cond){message("Fpea is empty")}))
                      out<-c(out,F1NPA)
                      out<-c(out,Fpea2)      
                      #################################      
                      
                    }
                 ################################ 
                  }
              ####################################################
                }else{
                  ########################################
                  print("enter the part 4")
                  ####################################
                  if((as.character(InMEDA[["Adduct"]]) == "[M]+" || as.character(InMEDA[["Adduct"]]) == "[M]-") & (F2NPA != 0)){
                   ####################
                    print("enter the part 4..if loop")
                    ##########################
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
                    NAT<-as.character(InMEDA[["Adduct"]])
                    NAT1<-paste("PRECURSORTYPE:",NAT,sep=" ")
                    out<-c(out,NAT1)
                    ################################
                    #################################
                    FIN1<-InMEDA[["Ionization mode"]]
                    F1IN1<-as.character(FIN1)
                    F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                    out<-c(out,F2IN1)
                    ##################################################
                    IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
                    ###################################################
                    #################################################
                    if(!sjmisc::is_empty(IKCRV)){
                      ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
                      F1ONT<-paste("Ontology:",ONTV,sep=" ")
                      out<-c(out,F1ONT)
                    }else{
                      ##F1ONT<-paste("Ontology:","",sep=" ")
                      F1ONT<-tryCatch({MaKE.ONT.REC(InMEDA)[1]},error=function(cond){message("Classifier could not fecth the information")})	    
                      out<-c(out,F1ONT)
                    }
                    
                    ############################################
                    ############################################
                    if(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])) & !startsWith(as.character(InMEDA[["InChI"]]),'not available') & !startsWith(as.character(InMEDA[["InChI"]]),'CAS:') & !startsWith(as.character(InMEDA[["InChI"]]),'InChI=')){
                      if(tryCatch({webchem::is.inchikey(stringr::str_trim(as.character(InMEDA[["InChI"]])))},error=function(cond){message("inchikey..file must be empty")})){
                        FINK<-ifelse(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])),paste("INCHIKEY:",as.character(InMEDA[["InChI"]])),paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" "))
                        FINCH<-paste("INCHI:",IN,sep=" ")
                        out<-c(out,FINK)
                        out<-c(out,FINCH)
                      }else{
                        FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                        FINCH<-ifelse(!sjmisc::is_empty(as.character(InMEDA[["InChI"]])),paste("INCHI:",as.character(InMEDA[["InChI"]])),paste("INCHI:",IN,sep=" "))
                        out<-c(out,FINK)
                        out<-c(out,FINCH)
                      }
                    }else{
                      FINK<-paste("INCHIKEY:",tryCatch({IK},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                      out<-c(out,FINK)
                      FINCH<-paste("INCHI:",IN,sep=" ")
                      out<-c(out,FINCH)
                    }
                   ##############################################
                    FSIM<-ifelse(!sjmisc::is_empty(stringr::str_trim(as.character(InMEDA[["SMILES"]]))),paste("SMILES:",stringr::str_trim(as.character(InMEDA[["SMILES"]])),sep=" "),paste("SMILES:",SM1,sep=" "))
                    
                    out<-c(out,FSIM)
                    ###########################################
                    ########################################
                    FFOR=ifelse(!sjmisc::is_empty(gETSmiles(InMEDA)),ifelse(!sjmisc::is_empty(tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)})),tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})),ifelse(!sjmisc::is_empty(InMEDA[["InChI"]]),tryCatch({getCactus(InMEDA[["InChI"]], "formula")},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})))
                    FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                    out<-c(out,FFOR1)
                    ########################################
                    FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                    FINS1<-FNA[FINS]
                    FINS2<-ifelse(!sjmisc::is_empty(stringr::str_trim(gsub("INTENSITY:","",FINS1))),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
                    ####FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
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
                    FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
                    out<-c(out,FCOM)
                    ##################
                    PT1<-c("Num Peaks:","NumPeaks:")
                    FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
                    F1NPA<-FNA[FNPA]
                    F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
                    ############################
                    ############################
                    Np<-F2NPA
                    Fpea<-FNA[(FNPA+1):Find]
                    Fpea1<-Fpea[!is.na(Fpea)]	      
                    Fpea2<-PeakVali(Fpea1)
                    #########################
                    if(F2NPA != 0)
                    ################################
                    {
                      
                      MV=AAMS1
                      tes1<-unlist(strsplit(Fpea, "\t|\t\t"))
                      tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
                      tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
                      tes4<-which(tes2 > (3+MV))
                      NTES4<-which(tes3 <= 100)
                      
		      NTES5<-which(tes3 <= 0.005*(mean(tes3[order(tes3, decreasing =TRUE)][1:5])))
                      NTES<-c(NTES4,tes4)
                      
                      if(Np >= 60)
                      {
                        
                        TES<-Centroid(Np,tes2,tes3,tes4,NTES,NTES4,NTES5,DEFAULT_MZ_TOLERANCE,Fpea2)
                        
                        F1NPA<-TES[[1]]
                        tes7<-TES[[2]]       
                        
                        out<-c(out,F1NPA)
                        out<-c(out,tes7)
                      }else{
                        #################################
                        F1NPA<-paste0("Num Peaks: ",tryCatch({length(Fpea2)},error=function(cond){message("Fpea is empty")}))
                        out<-c(out,F1NPA)
                        out<-c(out,Fpea2)      
                        #################################      
                        
                      }
                      
                      
                    }## end of the else
                    ############################
                    ############################
                  }             
                  #####################
                  #####################
                }
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
          #
        }else{
          PASS1 <-RRV
        } 
        ########################################
        ########################################
      }
    ########################################
      ##If I have to add further function 
    ########################################
    }
    
  } 
  
}



##############################################################################################################################################
##############################################################################################################################################
##print("enter the area of the function NFfilter")
##############################################################################################################################################
##############################################################################################################################################
NFFilter<-function(InMEDA,InAdVA,InMSPL,InPMZ,InRTL)
{
  ###########################
  out<-c()
  ###########################
  print("entering the function NFFilter")
  ###########################
  IV<-as.character(InMEDA[["InChI"]])
  EM<-as.character(InMEDA[["Exact mass"]])
  #############################
  EM1<-as.numeric(EM)
  #############################
  NEM1<-tryCatch({NFFilter1(InMEDA,InAdVA,InMSPL,InPMZ,InRTL)[[1]]},error=function(cond){message("NFFilter1 value is empty")})
  #############################
  print("entering my test area")
  print(NEM1)
  print(NEM1[[1]])
  ###############################
  if(!sjmisc::is_empty(NEM1)){
    #############################
    print("enter the line ..196")
    ####################################
    ####################################
    ADV<-as.character(InMEDA[["Adduct"]])
    ADV1<-tryCatch({qdapRegex::ex_between(ADV, "[", "]")[[1]]},error=function(cond){message("adduct value is missing")})
    ADV2<-tryCatch({FADINF(ADV)},error=function(cond){message("adduct value matching is not found")})
    ADV3<-tryCatch({stringr::str_replace(ADV2, "M",as.character(NEM1))},error=function(cond){message("adduct value replacement is not found")})
    ADV4<-tryCatch({as.numeric(pander::evals(ADV3)[[1]]$result)},error=function(cond){message("getting the result")})
    #####################
    #####################
    PPm=ADV4*(mz_Tol/(1000000))
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
    ####################
    #######################
    VRTL<-VRT-RT_Tol
    VRTU<-VRT+RT_Tol
    ##############################
    ################################
    TRTL<-InRTL[InRTL >= VRTL & InRTL <= VRTU]
    ITRTL<-which(InRTL %in% TRTL)
    ######################################
    INLL<-intersect(ITmass,ITRTL)
    ########################################
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
	  #################################################################
	  PT1<-c("Num Peaks:","NumPeaks:")
          FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
          F1NPA<-FNA[FNPA]
          F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
          Np<-F2NPA
          #######################################################
          PT <- c("PRECURSORTYPE:", "ADDUCTIONNAME:")
          NPT<-tryCatch({grep(paste(PT,collapse="|"), FNA, value=TRUE)},error=function(cond){message("List value is empty")})
          NPT1<-tryCatch({match(NPT,FNA)},error=function(cond){message("List value is empty")})
          PTV <-tryCatch({stringr::str_remove(FNA[NPT1],c("PRECURSORTYPE:","ADDUCTIONNAME:"))},error=function(cond){message("List value is empty")})
          PTV1<-tryCatch({PTV[1]},error=function(cond){message("List value is empty")})
          PTV2<-tryCatch({stringr::str_trim(PTV1)},error=function(cond){message("List value is empty")})
          P1TV2<-tryCatch({stringr::str_trim(gsub("ADDUCTIONNAME:","",PTV2))},error=function(cond){message("List value is empty")})
          PTV3<-tryCatch({paste("PRECURSORTYPE:",P1TV2)},error=function(cond){message("List value is empty")})
          ##PTV3<-tryCatch({paste("PRECURSORTYPE:",PTV2)},error=function(cond){message("List value is empty")})
          ####################################################
          print("enter the test area...1")
	  #########################################################
          ######################## adding this new#########
	  if((!sjmisc::is_empty(P1TV2) || !sjmisc::is_empty(as.character(InMEDA[["Adduct"]]))) & (F2NPA != 0)){
            ##############################################
            print("entering the line 234.............")
            ##############################################
            FNA1<-which(stringi::stri_detect_fixed(FNA,"NAME:"))
            FNA2<-tryCatch({InMEDA[["Name"]]},error=function(cond){message("NAME value is empty")})
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
	    ###########################
	    NAT<-as.character(InMEDA[["Adduct"]])
            NAT1<-paste("PRECURSORTYPE:",NAT,sep=" ")
            out<-c(out,NAT1)
            ###########################
            FIN1<-tryCatch({InMEDA[["Ionization mode"]]},error=function(cond){message("Ionization value is empty")})
            F1IN1<-as.character(FIN1)
            F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
            out<-c(out,F2IN1)
            ########################### need to add this new ## this is the ontology part
            ##print("check if problem in ontology part")
            ########### Adding this new ####################
            ################################################
            FINIKOT<-tryCatch({MaKE.ONT.REC(InMEDA)},error=function(cond){message("convert MaKE.ONT.REC is empty")})
            ##############################
	    ##############################
            if(!sjmisc::is_empty(FINIKOT)){
              #####################
              out<-c(out,FINIKOT)
              ####################
            }else{
		    print("enter the else loop ...")
	      ############################################	    
	      InKeyVal<-tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})	    
              ################################
              IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},error=function(cond){message("Classifier could not fecth the information")})
              ################################
              if(!sjmisc::is_empty(IKCRV)){
                ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
                F1ONT<-paste("Ontology:",ONTV,sep=" ")
                out<-c(out,F1ONT)
              }else{

		F1ONT<-tryCatch({MaKE.ONT.REC(InMEDA)[1]},error=function(cond){return(NA)})
                SM<-tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)})
	        tes1<-tryCatch({ClassSmilesToOntolgy(SM)},error=function(cond){message("smiles to ontology is failing")})
	        tes2<-paste("Ontology:",tes1,sep=" ")	
	       
	        F1ONT<-ifelse(!sjmisc::is_empty(tryCatch({MaKE.ONT.REC(InMEDA)[1]},error=function(cond){return(NA)})),tryCatch({MaKE.ONT.REC(InMEDA)[1]},error=function(cond){return(NA)}),tes2)
	        
                out<-c(out,F1ONT)
              }
              ###############################
              FINK<-paste("INCHIKEY:",tryCatch({InKeyVal},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
              ######out<-c(out,FINK)
	      out<-c(out,tryCatch({SelectInchiorInchiKey(InMEDA,IK,InchiV)[1]},error=function(cond){return(NA)}))
	      ################################
	      INV3<-tryCatch({rinchi::get.inchi(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
              FINCH<-paste("INCHI:",INV3,sep=" ")
	      out<-c(out,tryCatch({SelectInchiorInchiKey(InMEDA,IK,InchiV)[2]},error=function(cond){return(NA)}))
              ######################################
              
            }
            ################################################
	    #################################################
	    FSIM<-paste("SMILES:",tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}),sep=" ")
            out<-c(out,FSIM)
	    ##############################################
	    ##############################################
	    ##############################################
	    FFOR=ifelse(!sjmisc::is_empty(gETSmiles(InMEDA)),ifelse(!sjmisc::is_empty(tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)})),tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})),ifelse(!sjmisc::is_empty(InMEDA[["InChI"]]),tryCatch({getCactus(InMEDA[["InChI"]], "formula")},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})))
	    #################################################
	    #########################################
            FFOR1<-paste("FORMULA:",FFOR,sep=" ")
            out<-c(out,FFOR1)
            ################################
            ################################
            FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
            FINS1<-FNA[FINS]
	    FINS2<-ifelse(!sjmisc::is_empty(stringr::str_trim(gsub("INTENSITY:","",FINS1))),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
            ###FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
            out<-c(out,FINS2)
            ###############################
	    ###############################
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
            ################################
	    ################################
	    PT1<-c("Num Peaks:","NumPeaks:")
            FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
            F1NPA<-FNA[FNPA]
            F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
            Np<-F2NPA
            #######################
            Fpea<-FNA[(FNPA+1):Find]
            Fpea1<-Fpea[!is.na(Fpea)]
	    Fpea2<-PeakVali(Fpea1)
	    #######################################
	    #######################################
            if(F2NPA != 0)	    
	    #############################
            {
	      print("this is entering inside the if loop")
              ############################
              MV=ADV4
	      ###########################
	      tes1<-unlist(strsplit(Fpea2, "\t|\t\t|;| +"))
              tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
              tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
              tes4<-which(tes2 > (3+MV))
	      NTES4<-which(tes3 <= 100)
              NTES5<-which(tes3 <= 0.005*(mean(tes3[order(tes3, decreasing =TRUE)][1:10])))	      
	      NTES<-c(NTES4,tes4)



              if((Np >= 60) & (length(NTES) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))
  	      {


	        NTES1<-c(NTES4,tes4)


		if((Np >= 60) & (length(NTES1) > 1))
		{

			

			tes5<-tes2[-NTES1]
			tes6<-tes3[-NTES1]

                        mzV<-tes5
			intenV<-tes6


                        
			if(length(mzV) > 1)
  			{

    				nMZV=tryCatch({cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]},error=function(cond){return(0)})
    				nINV=tryCatch({cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]},error=function(cond){return(0)})
    				tes7<-paste(nMZV,nINV,sep="\t")
    				F1NPA<-paste0("Num Peaks: ",length(tes7))

    				out<-c(out,F1NPA)
    				out<-c(out,tes7)

  			}else{


    				tes5<-tes2
    				tes6<-tes3
                               



    				nMZV=tryCatch({cenPeaks(tes2,tes3, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]},error=function(cond){return(0)})
    				nINV=tryCatch({cenPeaks(tes2,tes3, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]},error=function(cond){return(0)})



    				if(length(nMZV) > 1 & length(nINV) > 1)
    				{
      					tes7<-paste(nMZV,nINV,sep="\t")

      					F1NPA<-paste0("Num Peaks: ",length(tes7))
      					out<-c(out,F1NPA)
      					out<-c(out,tes7)

    				}else{

    					tes7<-paste(tes5,tes6,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))
    					out<-c(out,F1NPA)
    					out<-c(out,tes7)
    				}## end of if else

  			}


		}else{



			if((length(tes2[-NTES]) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))
  			{
    				tes5<-tes2[-NTES]
    				tes6<-tes3[-NTES]
    				tes7<-paste(tes5,tes6,sep="\t")
    				F1NPA<-paste0("Num Peaks: ",length(tes7))
    				out<-c(out,F1NPA)
    				out<-c(out,tes7)

  			}else{

    				tes7<-paste(tes2,tes3,sep="\t")
    				F1NPA<-paste0("Num Peaks: ",length(tes7))
    				out<-c(out,F1NPA)
    				out<-c(out,tes7)

  			}## end of if else loop
		
		}#### end of the else

  	      }else{


	      if((length(NTES) > 1) & (length(NTES) != Np) & (length(Fpea2) > length(NTES)))

              {

                   if(length(tes2[-NTES]) > 1)
    		   {

      			tes5<-tes2[-NTES]
      			tes6<-tes3[-NTES]


			tes7<-paste(tes5,tes6,sep="\t")
      			### Adding this new here
      			F1NPA<-paste0("Num Peaks: ",length(tes7))
      			out<-c(out,F1NPA)
      			#############################
      			out<-c(out,tes7)
      			##############################
    		   }else{


      			tes5<-tes2
      			tes6<-tes3
      			tes7<-paste(tes5,tes6,sep="\t")
      			F1NPA<-paste0("Num Peaks: ",length(tes7))
      			out<-c(out,F1NPA)
      			out<-c(out,tes7)

    			}### if,else loop


              }else{


		      if((Np >= 60) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))
		      {


			      mzV<-tes2
			      intenV<-tes3

			      nMZV=cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]
			      nINV=cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]

			      tes7<-paste(nMZV,nINV,sep="\t")
			      F1NPA<-paste0("Num Peaks: ",length(tes7))

			      out<-c(out,F1NPA)
			      out<-c(out,tes7)


		      }else{



			      ######## Adding this new ####
			      F1NPA<-paste0("Num Peaks: ",tryCatch({length(Fpea2)},error=function(cond){message("Fpea is empty")}))
                	      out<-c(out,F1NPA)
			      out<-c(out,Fpea2)
	        	      #############################

	          }### end of else

	 #########################################################
         #########################################################      
              }
	 ########################################################
	      }### end of the else and Np >= 120                  
	 ##############################################################     
	 ##############################################################
	 ##############################################################  
            } ## end of else
       #########################################################
       #########################################################
          }else{
            ####################################	  
            print("enter the else part...1")
            ####################################
            ####################################
	    if((as.character(InMEDA[["Adduct"]]) == "[M]+" || as.character(InMEDA[["Adduct"]]) == "[M]-") & (F2NPA != 0)){
              #####################################
              #####################################
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
	      #################################
	      NAT<-as.character(InMEDA[["Adduct"]])
              NAT1<-paste("PRECURSORTYPE:",NAT,sep=" ")
              out<-c(out,NAT1)
              #################################
              FIN1<-InMEDA[["Ionization mode"]]
              F1IN1<-as.character(FIN1)
              F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
              out<-c(out,F2IN1)
              ##################################################
	      ##################################################
	      FINIKOT<-tryCatch({MaKE.ONT.REC(InMEDA)},error=function(cond){message("convert MaKE.ONT.REC is empty")})
              ##################################################
              ##################################################
              if(!sjmisc::is_empty(FINIKOT)){
                #############################
                out<-c(out,FINIKOT)
                #############################
              }else{
                ###################################
		###################################      
		InKeyVal<-tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})      
                ################################
                IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
                #############################
                if(!sjmisc::is_empty(IKCRV)){
                  ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
                  F1ONT<-paste("Ontology:",ONTV,sep=" ")
                  out<-c(out,F1ONT)
                }else{
                  ##F1ONT<-paste("Ontology:","",sep=" ")
		  F1ONT<-tryCatch({MaKE.ONT.REC(InMEDA)[1]},error=function(cond){message("Classifier could not fecth the information")})    	
                  out<-c(out,F1ONT)
                }
                ##########################
                FINK<-paste("INCHIKEY:",tryCatch({InKeyVal},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
		out<-c(out,tryCatch({SelectInchiorInchiKey(InMEDA,IK,InchiV)[1]},error=function(cond){return(NA)}))
                ####out<-c(out,FINK)
                ################################
		###################################
		INV3<-tryCatch({rinchi::get.inchi(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
                FINCH<-paste("INCHI:",INV3,sep=" ")
                ####out<-c(out,FINCH)
		out<-c(out,tryCatch({SelectInchiorInchiKey(InMEDA,IK,InchiV)[2]},error=function(cond){return(NA)}))
                #############################
                
                
              }
              ###########################
              ###########################
	      FSIM<-paste("SMILES:",tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}),sep=" ")
              out<-c(out,FSIM)
              ############################
	      ############################
	      FFOR=ifelse(!sjmisc::is_empty(gETSmiles(InMEDA)),ifelse(!sjmisc::is_empty(tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)})),tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})),ifelse(!sjmisc::is_empty(InMEDA[["InChI"]]),tryCatch({getCactus(InMEDA[["InChI"]], "formula")},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})))
	      ############################
	      #############################
              FFOR1<-paste("FORMULA:",FFOR,sep=" ")
              out<-c(out,FFOR1)
              #############################
              #############################
              FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
              FINS1<-FNA[FINS]
	      FINS2<-ifelse(!sjmisc::is_empty(stringr::str_trim(gsub("INTENSITY:","",FINS1))),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
              ###FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
              out<-c(out,FINS2)
              ############################
              FAUT<-as.character(InMEDA[["Authors"]])
              FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
              out<-c(out,FAUT1)
              #############################
              #############################
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
              FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
              out<-c(out,FCOM)
              #########################
	      #########################
	      PT1<-c("Num Peaks:","NumPeaks:")
              FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
              F1NPA<-FNA[FNPA]
              F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
              Np<-F2NPA
              ###################
              Fpea<-FNA[(FNPA+1):Find]
	      Fpea1<-Fpea[!is.na(Fpea)]
	      Fpea2<-PeakVali(Fpea1)
	      #####################################
	      #####################################
	      if(F2NPA != 0)
	      #####################################
              {
                
                MV=AAMS1
	        
	        tes1<-unlist(strsplit(Fpea2, "\t|\t\t|;| +"))
                tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
                tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
                tes4<-which(tes2 > (3+MV))
		NTES4<-which(tes3 <= 100)
		
		NTES5<-which(tes3 <= 0.005*(mean(tes3[order(tes3, decreasing =TRUE)][1:10])))
		NTES<-c(NTES4,tes4)

		###if(Np >= 60)
		if((Np >= 60) & (length(NTES) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))
  		{
			NTES1<-c(NTES4,tes4)


			if((Np >= 60) & (length(NTES1) > 1))
			{
			

				tes5<-tes2[-NTES1]
				tes6<-tes3[-NTES1]

				mzV<-tes5
				intenV<-tes6

				 if(length(mzV) > 1)
  				 {

    					nMZV=tryCatch({cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]},error=function(cond){return(0)})
    					nINV=tryCatch({cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]},error=function(cond){return(0)})
    					tes7<-paste(nMZV,nINV,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))

    					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  				  } else{

    					tes5<-tes2
    					tes6<-tes3

    					nMZV=tryCatch({cenPeaks(tes2,tes3, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]},error=function(cond){return(0)})
    					nINV=tryCatch({cenPeaks(tes2,tes3, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]},error=function(cond){return(0)})

    					if(length(nMZV) > 1 & length(nINV) > 1)
    					{

      						tes7<-paste(nMZV,nINV,sep="\t")

      						F1NPA<-paste0("Num Peaks: ",length(tes7))
      						out<-c(out,F1NPA)
      						out<-c(out,tes7)

    					}else{

    						tes7<-paste(tes5,tes6,sep="\t")
    						F1NPA<-paste0("Num Peaks: ",length(tes7))
    						out<-c(out,F1NPA)
    						out<-c(out,tes7)
    					}## end of if else

  				}


			}else{


				if((length(tes2[-NTES]) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES))) 
  				{
    					tes5<-tes2[-NTES]
    					tes6<-tes3[-NTES]
    					tes7<-paste(tes5,tes6,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))
    					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  				}else{

    					tes7<-paste(tes2,tes3,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))
    					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  				}## end of if else loop
			}### end of the else
  		}else{

	        if((length(NTES) > 1) & (length(NTES) != Np) & (length(Fpea2) > length(NTES)))
                {

		 ##########################################	
                  if(length(tes2[-NTES]) > 1)
    		  {

      			tes5<-tes2[-NTES]
      			tes6<-tes3[-NTES]
      			tes7<-paste(tes5,tes6,sep="\t")
      			### Adding this new here
      			F1NPA<-paste0("Num Peaks: ",length(tes7))
      			out<-c(out,F1NPA)
      			#############################
      			out<-c(out,tes7)
      			##############################
    		}else{

      			tes5<-tes2
      			tes6<-tes3
      			tes7<-paste(tes5,tes6,sep="\t")
      			F1NPA<-paste0("Num Peaks: ",length(tes7))
      			out<-c(out,F1NPA)
      			out<-c(out,tes7)

    			}### if,else loop


                }else{

			if((Np >= 60))
    			{
      				mzV<-tes2
      				intenV<-tes3

      				nMZV=cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]
      				nINV=cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]

      				tes7<-paste(nMZV,nINV,sep="\t")
      				F1NPA<-paste0("Num Peaks: ",length(tes7))

      				out<-c(out,F1NPA)
      				out<-c(out,tes7)


    			}else{

	          ###########################
		  F1NPA<-paste0("Num Peaks: ",tryCatch({length(Fpea2)},error=function(cond){message("Fpea is empty")}))		
                  out<-c(out,F1NPA)
		  out<-c(out,Fpea2)
	          ###########################
                      }
                 ###################################
               ###############################
                }
            ##########################################################
		}### end of else and Np >= 120
	    ##########################################################	

              }### end of the if loop and Num Peaks 
           ############################################
            }
          #################################################
          }
       #########################################################  
       #########################################################
        }else if(length(INLL) > 1){
          ##############################
          print("enter the line ...511...........")
          ##############################
          MONMS=InMSPL[INLL]
          #########################
          TRA<-unname(rapply(MONMS, function(x) grep("RETENTIONTIME:",x, value=TRUE)))
          ##########print("enter the value of TRA")
	  ##########print(TRA)
	  TRA1<-as.numeric(stringr::str_trim(stringr::str_replace(TRA, "RETENTIONTIME:", "")))
          TRA2<-abs(VRT-TRA1)
	  ###########print("enter the value of TRA2")
	  TRA3<-which(TRA2==min(TRA2[TRA2 > 0]))
	  #################################################
	  PT1<-c("Num Peaks:","NumPeaks:")
	  ####################################################
	  NPV1=unname(rapply(MONMS, function(x) grep(paste(PT1,collapse="|"), x, value=TRUE)))
          NPV=as.numeric(stringr::str_trim(stringr::str_replace(NPV1, "Num Peaks:|NumPeaks:", "")))
          ##############################################################################
	  ####################################################
	  NPNZ=which(NPV!=0)
	  NTRA3<-intersect(TRA3,NPNZ)
	  TRA4<-INLL[NTRA3]
          ##########################
	  ##########################
          TRA5<-InMSPL[TRA4]
          ###########################
          F1FPL<-TRA5
          ###########################
          ###########################
          F1FPL<-TRA5
          F2FPL<-F1FPL
          ##########################
          Find<-tryCatch({length(F2FPL[[1]])},error=function(cond){message("there is an error in list")})
          ###########################
          FNA<-tryCatch({F1FPL[[1]]},error=function(cond){message("there is an error in list")})
          ####### adding this new ###################################################
	  ###################################################################
	  PT1<-c("Num Peaks:","NumPeaks:")
          FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
          F1NPA<-FNA[FNPA]
          F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
	  Np<-F2NPA
          ##############################
          ##Fpea<-FNA[(FNPA+1):Find]
          ###########################################################################
          PT <- c("PRECURSORTYPE:", "ADDUCTIONNAME:")
          NPT<-tryCatch({grep(paste(PT,collapse="|"), FNA, value=TRUE)},error=function(cond){message("List value is empty")})
          NPT1<-tryCatch({match(NPT,FNA)},error=function(cond){message("List value is empty")})
          PTV <-tryCatch({stringr::str_remove(FNA[NPT1],c("PRECURSORTYPE:","ADDUCTIONNAME:"))},error=function(cond){message("List value is empty")})
          PTV1<-tryCatch({PTV[1]},error=function(cond){message("List value is empty")})
          PTV2<-tryCatch({stringr::str_trim(PTV1)},error=function(cond){message("List value is empty")})
          P1TV2<-tryCatch({stringr::str_trim(gsub("ADDUCTIONNAME:","",PTV2))},error=function(cond){message("List value is empty")})
          PTV3<-tryCatch({paste("PRECURSORTYPE:",P1TV2)},error=function(cond){message("List value is empty")})
          ##PTV3<-tryCatch({paste("PRECURSORTYPE:",PTV2)},error=function(cond){message("List value is empty")})
          ##########################################################################
          ##################################################
	  if((!sjmisc::is_empty(P1TV2) || !sjmisc::is_empty(as.character(InMEDA[["Adduct"]]))) & (F2NPA != 0)){
          ##if(!sjmisc::is_empty(P1TV2) || !sjmisc::is_empty(as.character(InMEDA[["Adduct"]]))){
            ##################################################		  
            ##################################################
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
	    ####################################################
            NAT<-as.character(InMEDA[["Adduct"]])
            NAT1<-paste("PRECURSORTYPE:",NAT,sep=" ")
            out<-c(out,NAT1)
            ############################################################
            FIN1<-InMEDA[["Ionization mode"]]
            F1IN1<-as.character(FIN1)
            F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
            out<-c(out,F2IN1)
            ##########################################################
            ##################### adding this new ####################
            FINIKOT<-MaKE.ONT.REC(InMEDA)
	    ##########################################################
            #####################################################
            #####################################################
            if(!sjmisc::is_empty(FINIKOT)){
              #####################
              out<-c(out,FINIKOT)
              ####################
            }else{
              
	      ##############################################	    
	      InKeyVal<-tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})	    
              ################################
              ################################
              IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
              #################################################
              if(!sjmisc::is_empty(IKCRV)){
                ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
                F1ONT<-paste("Ontology:",ONTV,sep=" ")
                out<-c(out,F1ONT)
              }else{
                ##F1ONT<-paste("Ontology:","",sep=" ")
		F1ONT<-tryCatch({MaKE.ONT.REC(InMEDA)[1]},error=function(cond){message("Classifier could not fecth the information")})      
                out<-c(out,F1ONT)
              }
              ##################################
              FINK<-paste("INCHIKEY:",tryCatch({InKeyVal},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
              ###out<-c(out,FINK)
	      out<-c(out,tryCatch({SelectInchiorInchiKey(InMEDA,IK,InchiV)[1]},error=function(cond){return(NA)}))
              ################################
	      INV3<-tryCatch({rinchi::get.inchi(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
	      ###################################
              FINCH<-paste("INCHI:",INV3,sep=" ")
              ####out<-c(out,FINCH)
	      out<-c(out,tryCatch({SelectInchiorInchiKey(InMEDA,IK,InchiV)[2]},error=function(cond){return(NA)}))
              ###################################
            }
            ##########################################################
            ##########################################################
            FSIM<-paste("SMILES:",tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}),sep=" ")
            out<-c(out,FSIM)
	    #########################################################
            #############################################################
	    FFOR=ifelse(!sjmisc::is_empty(gETSmiles(InMEDA)),ifelse(!sjmisc::is_empty(tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)})),tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})),ifelse(!sjmisc::is_empty(InMEDA[["InChI"]]),tryCatch({getCactus(InMEDA[["InChI"]], "formula")},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})))
	    ##############################################################
	    #############################################################
            FFOR1<-paste("FORMULA:",FFOR,sep=" ")
            out<-c(out,FFOR1)
            ##########################################################
            ##########################################################
            FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
            FINS1<-FNA[FINS]
	    FINS2<-ifelse(!sjmisc::is_empty(stringr::str_trim(gsub("INTENSITY:","",FINS1))),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
            ###########################################
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
            FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
            out<-c(out,FCOM)
            ##################################################
	    ##################################################
	    PT1<-c("Num Peaks:","NumPeaks:")
            FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
            F1NPA<-FNA[FNPA]
            F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
	    Np<-F2NPA
            ###################
            Fpea<-FNA[(FNPA+1):Find]
	    Fpea1<-Fpea[!is.na(Fpea)]
	    Fpea2<-PeakVali(Fpea1)
	    ########################
	    ######################################
	    if(F2NPA != 0)
	    #####################################
            {
	       
	      MV=ADV4
	      tes1<-unlist(strsplit(Fpea2, "\t|\t\t|;| +"))
              tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
              tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
              tes4<-which(tes2 > (3+MV))
	      NTES4<-which(tes3 <= 100)
	      NTES5<-which(tes3 <= 0.005*(mean(tes3[order(tes3, decreasing =TRUE)][1:10])))
	      NTES<-c(NTES4,tes4)
              
	      if((Np >= 60) & (length(NTES) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))
  	      {
		NTES1<-c(NTES4,tes4)      

	        if((Np >= 60) & (length(NTES1) > 1))
		{
			tes5<-tes2[-NTES1]
			tes6<-tes3[-NTES1]

			mzV<-tes5
			intenV<-tes6


			if(length(mzV) > 1)
  			{

    				nMZV=tryCatch({cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]},error=function(cond){return(0)})
    				nINV=tryCatch({cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]},error=function(cond){return(0)})
    				tes7<-paste(nMZV,nINV,sep="\t")
    				F1NPA<-paste0("Num Peaks: ",length(tes7))

    				out<-c(out,F1NPA)
    				out<-c(out,tes7)

  			}else{

    				tes5<-tes2
    				tes6<-tes3

    				nMZV=tryCatch({cenPeaks(tes2,tes3, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]},error=function(cond){return(0)})
   				nINV=tryCatch({cenPeaks(tes2,tes3, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]},error=function(cond){return(0)})

    				if(length(nMZV) > 1 & length(nINV) > 1)
    				{

      					tes7<-paste(nMZV,nINV,sep="\t")

      					F1NPA<-paste0("Num Peaks: ",length(tes7))
      					out<-c(out,F1NPA)
      					out<-c(out,tes7)

    				}else{

    					tes7<-paste(tes5,tes6,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))
    					out<-c(out,F1NPA)
    					out<-c(out,tes7)
    					}## end of if else

  			}


		}else{


			if((length(tes2[-NTES]) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))
  			{
    				tes5<-tes2[-NTES]
    				tes6<-tes3[-NTES]
    				tes7<-paste(tes5,tes6,sep="\t")
    				F1NPA<-paste0("Num Peaks: ",length(tes7))
    				out<-c(out,F1NPA)
    				out<-c(out,tes7)

  			}else{

    				tes7<-paste(tes2,tes3,sep="\t")
    				F1NPA<-paste0("Num Peaks: ",length(tes7))
    				out<-c(out,F1NPA)
    				out<-c(out,tes7)

  			}## end of if else loop

		}### end of the else

  		}else{
	
	      if((length(NTES) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))		
              {
                
		if(length(tes2[-NTES]) > 1)
    		{

      			tes5<-tes2[-NTES]
      			tes6<-tes3[-NTES]
      			tes7<-paste(tes5,tes6,sep="\t")
      			### Adding this new here
      			F1NPA<-paste0("Num Peaks: ",length(tes7))
      			out<-c(out,F1NPA)
      			#############################
      			out<-c(out,tes7)
      			##############################
    		}else{

      			tes5<-tes2
      			tes6<-tes3
      			tes7<-paste(tes5,tes6,sep="\t")
      			F1NPA<-paste0("Num Peaks: ",length(tes7))
      			out<-c(out,F1NPA)
      			out<-c(out,tes7)

    			}### if,else loop

              }else{

		      if((Np >= 60))
    		      {
      			mzV<-tes2
      			intenV<-tes3

      			nMZV=cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]
      			nINV=cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]

      			tes7<-paste(nMZV,nINV,sep="\t")
      			F1NPA<-paste0("Num Peaks: ",length(tes7))

      			out<-c(out,F1NPA)
      			out<-c(out,tes7)


    		    }else{


		######################
		F1NPA<-paste0("Num Peaks: ",tryCatch({length(Fpea2)},error=function(cond){message("Fpea is empty")}))	    
                out<-c(out,F1NPA)
		out<-c(out,Fpea2)
		######### Adding this new

	       }###end of else amd Np >= 80
           #################################
              }
	    #######################################################
	      }### end of the if loop and Np >= 120
	   ########################################################   
             } ### this is the else MV= ADV4...closing and end of if loop
           ############################## adding this new
	   ##### print test the value
	    ##print("entering the out value")
	    ##print(out)
	    #############################################
            #############################################

          }else{
            #############################################
            print("enter the else part...2")
            #############################################
            #############################################
	    if((as.character(InMEDA[["Adduct"]]) == "[M]+" || as.character(InMEDA[["Adduct"]]) == "[M]-") & (F2NPA != 0)){
              ############################################		    
              ############################################
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
	      #################################
	      NAT<-as.character(InMEDA[["Adduct"]])
              NAT1<-paste("PRECURSORTYPE:",NAT,sep=" ")
              out<-c(out,NAT1)
              #################################
              FIN1<-InMEDA[["Ionization mode"]]
              F1IN1<-as.character(FIN1)
              F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
              out<-c(out,F2IN1)
              ###########################################################
              ############### Adding this new############################
              FINIKOT<-MaKE.ONT.REC(InMEDA)
              ######################################################
              ######################################################
              if(!sjmisc::is_empty(FINIKOT)){
                #####################
                out<-c(out,FINIKOT)
                ####################
              }else{
                
                ###############################################
		InKeyVal<-tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})      
                ################################
                ################################
                IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
                #################################################
                if(!sjmisc::is_empty(IKCRV)){
                  ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
                  F1ONT<-paste("Ontology:",ONTV,sep=" ")
                  out<-c(out,F1ONT)
                }else{
                  ##F1ONT<-paste("Ontology:","",sep=" ")
		  F1ONT<-tryCatch({MaKE.ONT.REC(InMEDA)[1]},error=function(cond){message("Classifier could not fecth the information")})	
                  out<-c(out,F1ONT)
                }
                ##################################
                FINK<-paste("INCHIKEY:",tryCatch({InKeyVal},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                out<-c(out,tryCatch({SelectInchiorInchiKey(InMEDA,IK,InchiV)[1]},error=function(cond){return(NA)}))
		###out<-c(out,FINK)
                ################################
                INV3<-tryCatch({rinchi::get.inchi(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
		FINCH<-paste("INCHI:",INV3,sep=" ")
		out<-c(out,tryCatch({SelectInchiorInchiKey(InMEDA,IK,InchiV)[2]},error=function(cond){return(NA)}))
                ####out<-c(out,FINCH)
                #################################
              }
              ##############################
              ##############################
              ###FSIM<-paste("SMILES:",gETSmiles(InMEDA),sep=" ")
	      FSIM<-paste("SMILES:",tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}),sep=" ")
              out<-c(out,FSIM)
              ####################################
	      print("check the smiles value in 2nd function")
	      print(FSIM)
	      ####################################
	      FFOR=ifelse(!sjmisc::is_empty(gETSmiles(InMEDA)),ifelse(!sjmisc::is_empty(tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)})),tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})),ifelse(!sjmisc::is_empty(InMEDA[["InChI"]]),tryCatch({getCactus(InMEDA[["InChI"]], "formula")},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})))
              ###########################################
	      FFOR1<-paste("FORMULA:",FFOR,sep=" ")
              out<-c(out,FFOR1)             
              ####################################################
              ###################################################
              FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
              FINS1<-FNA[FINS]
	      FINS2<-ifelse(!sjmisc::is_empty(stringr::str_trim(gsub("INTENSITY:","",FINS1))),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
              ####FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
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
              FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
              out<-c(out,FCOM)
              ##################
	      ###############################################
              PT1<-c("Num Peaks:","NumPeaks:")
              FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
              F1NPA<-FNA[FNPA]
              F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
              Np<-F2NPA
              Fpea<-FNA[(FNPA+1):Find]
	      Fpea1<-Fpea[!is.na(Fpea)]
	      Fpea2<-PeakVali(Fpea1)
	      ###################################################
              ###################################################
	      if(F2NPA != 0)
	      ###############################
              {
                
                MV=AAMS1
	        tes1<-unlist(strsplit(Fpea2, "\t|\t\t|;| +"))
                tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
                tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
                tes4<-which(tes2 > (3+MV))
		NTES4<-which(tes3 <= 100)
		NTES<-c(NTES4,tes4)


		if((Np >= 60) & (length(NTES) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))
  		{
			NTES1<-c(NTES4,tes4)


			if((Np >= 60) & (length(NTES1) > 1))
			{
				tes5<-tes2[-NTES1]
				tes6<-tes3[-NTES1]

				mzV<-tes5
				intenV<-tes6

				if(length(mzV) > 1)
  				{

    					nMZV=tryCatch({cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]},error=function(cond){return(0)})
    					nINV=tryCatch({cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]},error=function(cond){return(0)})
    					tes7<-paste(nMZV,nINV,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))

    					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  				}else{

    					tes5<-tes2
    					tes6<-tes3

    					nMZV=tryCatch({cenPeaks(tes2,tes3, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]},error=function(cond){return(0)})
    					nINV=tryCatch({cenPeaks(tes2,tes3, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]},error=function(cond){return(0)})

    					if(length(nMZV) > 1 & length(nINV) > 1)
    					{

      						tes7<-paste(nMZV,nINV,sep="\t")

      						F1NPA<-paste0("Num Peaks: ",length(tes7))
      						out<-c(out,F1NPA)
      						out<-c(out,tes7)

    					}else{

    						tes7<-paste(tes5,tes6,sep="\t")
    						F1NPA<-paste0("Num Peaks: ",length(tes7))
    						out<-c(out,F1NPA)
    						out<-c(out,tes7)

    					}## end of if else

  				}


			}else{


				if((length(tes2[-NTES]) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))
  				{
    					tes5<-tes2[-NTES]
   					tes6<-tes3[-NTES]
    					tes7<-paste(tes5,tes6,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))
    					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  				}else{

    					tes7<-paste(tes2,tes3,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))
    					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  				}## end of if else loop

			}### End of the else
  		}else{

		if((length(NTES) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))	
                {
		  ########################

	          if(length(tes2[-NTES]) > 1)
    		  {

      			tes5<-tes2[-NTES]
      			tes6<-tes3[-NTES]
      			tes7<-paste(tes5,tes6,sep="\t")
      			### Adding this new here
      			F1NPA<-paste0("Num Peaks: ",length(tes7))
      			out<-c(out,F1NPA)
      			#############################
      			out<-c(out,tes7)
      			##############################
    		}else{

      			tes5<-tes2
      			tes6<-tes3
      			tes7<-paste(tes5,tes6,sep="\t")
      			F1NPA<-paste0("Num Peaks: ",length(tes7))
      			out<-c(out,F1NPA)
      			out<-c(out,tes7)

    			}### if,else loop		

                }else{


                  
		  if((Np >= 60))	
    		   {
      			mzV<-tes2
      			intenV<-tes3
      
      			nMZV=cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]
      			nINV=cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]
      			tes7<-paste(nMZV,nINV,sep="\t")
      			F1NPA<-paste0("Num Peaks: ",length(tes7))
      			out<-c(out,F1NPA)
      			out<-c(out,tes7)
      
      
                   }else{	


		  ################################
		  F1NPA<-paste0("Num Peaks: ",tryCatch({length(Fpea2)},error=function(cond){message("Fpea is empty")}))	   
                  out<-c(out,F1NPA)
		  out<-c(out,Fpea2)
		  ################################
              #####################################
		   }### end of else and (Np >= 80)
              #################################
                }
            ##########################################
		}### end of the if loop and Np >= 120
	    ##########################################	
              }### end of the if loop and F1Npa 
           ############################
            }
          #########################
          #########################
          }
         ###############################
         }else{
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
    ############################################
    print("entering this line..line 771")
    ###########################################
    FM<-as.character(InMEDA[["Exact mass"]])
    ###########################################
    if(!sjmisc::is_empty(FM))
    {
      ###################################################################
      print("entering the line ... 772")
      ###################################################################
      ####################################################################
      if(!sjmisc::is_empty(FM) & !sjmisc::is_empty(as.character(InMEDA[["Adduct"]]))){
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
        ##PPm=ADV4*(25/(1000000))
	PPm=ADV4*(mz_Tol/(1000000))
        #######################
        MPPmL=ADV4-PPm
        MPPmU=ADV4+PPm
        #######################
        Tmass<-InPMZ[InPMZ >= MPPmL & InPMZ <= MPPmU]
        ITmass<-which(InPMZ %in% Tmass)
        ##ITmass<-match(Tmass,InPMZ)
        #####################
        VRT<-as.numeric(as.character(InMEDA[["RT (min)"]]))
	#####################
        ##VRTL<-VRT-0.20
        ##VRTU<-VRT+0.20
	#####################
	VRTL<-VRT-RT_Tol
	VRTU<-VRT+RT_Tol
        #####################
        print("enter the function 2 ..else part")
        ############################
        TRTL<-InRTL[InRTL >= VRTL & InRTL <= VRTU]
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
              FNA<-tryCatch({F1FPL[[1]]},error=function(cond){message("there is an error in list value")})
              ######################## adding this new ###########
	      ####################################################
	      PT1<-c("Num Peaks:","NumPeaks:")
              FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
              F1NPA<-FNA[FNPA]
              F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
	      Np<-F2NPA
              ####################################################
              PT <- c("PRECURSORTYPE:", "ADDUCTIONNAME:")
              NPT<-tryCatch({grep(paste(PT,collapse="|"), FNA, value=TRUE)},error=function(cond){message("List value is empty")})
              NPT1<-tryCatch({match(NPT,FNA)},error=function(cond){message("List value is empty")})
              PTV <-tryCatch({stringr::str_remove(FNA[NPT1],c("PRECURSORTYPE:","ADDUCTIONNAME:"))},error=function(cond){message("List value is empty")})
              PTV1<-tryCatch({PTV[1]},error=function(cond){message("List value is empty")})
              PTV2<-tryCatch({stringr::str_trim(PTV1)},error=function(cond){message("List value is empty")})
              P1TV2<-tryCatch({stringr::str_trim(gsub("ADDUCTIONNAME:","",PTV2))},error=function(cond){message("List value is empty")})
              PTV3<-tryCatch({paste("PRECURSORTYPE:",P1TV2)},error=function(cond){message("List value is empty")})
	      ####################################################
              #####################################################
              print("enter the test area..3")
              ####################################################
              ################################ adding this new ###
              ###############################################
	      if((!sjmisc::is_empty(P1TV2) || !sjmisc::is_empty(as.character(InMEDA[["Adduct"]]))) & (F2NPA != 0)){
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
		#################################
		NAT<-as.character(InMEDA[["Adduct"]])
                NAT1<-paste("PRECURSORTYPE:",NAT,sep=" ")
                out<-c(out,NAT1)
                #################################
                FIN1<-InMEDA[["Ionization mode"]]
                F1IN1<-as.character(FIN1)
                F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                out<-c(out,F2IN1)
                ###############################################
                ###############################################
		InKeyVal<-tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
                ################################
                ###############################
                IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
                #################################################
                if(!sjmisc::is_empty(IKCRV)){
                  ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
                  F1ONT<-paste("Ontology:",ONTV,sep=" ")
                  out<-c(out,F1ONT)
                }else{
			F1ONT<-tryCatch({MaKE.ONT.REC(InMEDA)[1]},error=function(cond){message("Classifier could not fecth the information")})
                        ##F1ONT<-paste("Ontology:","",sep=" ")
                        out<-c(out,F1ONT)
                }
                ##################################
                FINK<-paste("INCHIKEY:",tryCatch({InKeyVal},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                ###out<-c(out,FINK)
		out<-c(out,tryCatch({SelectInchiorInchiKey(InMEDA,IK,InchiV)[1]},error=function(cond){return(NA)}))
                ###################################
		###################################
		INV3<-tryCatch({rinchi::get.inchi(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
                FINCH<-paste("INCHI:",INV3,sep=" ")
		out<-c(out,tryCatch({SelectInchiorInchiKey(InMEDA,IK,InchiV)[2]},error=function(cond){return(NA)}))
                ####out<-c(out,FINCH)
                #######################################
		#######################################
		FSIM<-paste("SMILES:",tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}),sep=" ")
                out<-c(out,FSIM)
		##################################################
                ################################
		FFOR=ifelse(!sjmisc::is_empty(gETSmiles(InMEDA)),ifelse(!sjmisc::is_empty(tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)})),tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})),ifelse(!sjmisc::is_empty(InMEDA[["InChI"]]),tryCatch({getCactus(InMEDA[["InChI"]], "formula")},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})))
		####################################
		#####################################
                FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                out<-c(out,FFOR1)
                ###############################
                ##############################
                FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                FINS1<-FNA[FINS]
		FINS2<-ifelse(!sjmisc::is_empty(stringr::str_trim(gsub("INTENSITY:","",FINS1))),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
                ####FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
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
                PT1<-c("Num Peaks:","NumPeaks:")
                FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
                F1NPA<-FNA[FNPA]
                F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
                Np<-F2NPA
                Fpea<-FNA[(FNPA+1):Find]
		Fpea1<-Fpea[!is.na(Fpea)]
		Fpea2<-PeakVali(Fpea1)
		################################################
		################################################
		if(F2NPA != 0)
		###########################
                {
                ############}else{
                  MV=ADV4
		  #################################################
		  tes1<-unlist(strsplit(Fpea2, "\t|\t\t|;| +"))
                  tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
                  tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
                  tes4<-which(tes2 > (3+MV))
		  NTES4<-which(tes3 <= 100)
		  NTES<-c(NTES4,tes4)


		  if((Np >= 60) & (length(NTES) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))
  		   {
    			###NTES1<-c(NTES4,tes4,NTES5)
			NTES1<-c(NTES4,tes4)   
    			
		         

		        if((Np >= 60))
			{
				tes5<-tes2[-NTES1]
				tes6<-tes3[-NTES1]

				mzV<-tes5
				intenV<-tes6

				if(length(mzV) > 1)
  				{

    					print("enter inside if loop")
    					nMZV=tryCatch({cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]},error=function(cond){return(0)})
    					nINV=tryCatch({cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]},error=function(cond){return(0)})
    					tes7<-paste(nMZV,nINV,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))

    					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  				}else{

    					tes5<-tes2
    					tes6<-tes3

    					nMZV=tryCatch({cenPeaks(tes2,tes3, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]},error=function(cond){return(0)})
    					nINV=tryCatch({cenPeaks(tes2,tes3, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]},error=function(cond){return(0)})

   					if(length(nMZV) > 1 & length(nINV) > 1)
    					{

      						print("enter the value of nMZV and nINV")
      						tes7<-paste(nMZV,nINV,sep="\t")

      						F1NPA<-paste0("Num Peaks: ",length(tes7))
      						out<-c(out,F1NPA)
      						out<-c(out,tes7)

    					}else{

    						tes7<-paste(tes5,tes6,sep="\t")
    						F1NPA<-paste0("Num Peaks: ",length(tes7))
    						out<-c(out,F1NPA)
    						out<-c(out,tes7)

    						}## end of if else

  					}


			}else{


				if((length(tes2[-NTES]) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))  
  				{
    					tes5<-tes2[-NTES]
    					tes6<-tes3[-NTES]
    					tes7<-paste(tes5,tes6,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))
    					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  				}else{

    					tes7<-paste(tes2,tes3,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))
    					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  				}## end of if else loop

			}### end of the else

  		   }else{
		
	          if((length(NTES) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))		   
                  {

		    if(length(tes2[-NTES]) > 1)
    		    {

      			tes5<-tes2[-NTES]
      			tes6<-tes3[-NTES]
      			tes7<-paste(tes5,tes6,sep="\t")
     			### Adding this new here
      			F1NPA<-paste0("Num Peaks: ",length(tes7))
      			out<-c(out,F1NPA)
      			#############################
      			out<-c(out,tes7)
      			##############################
    		    }else{

      			tes5<-tes2
      			tes6<-tes3
      			tes7<-paste(tes5,tes6,sep="\t")
      			F1NPA<-paste0("Num Peaks: ",length(tes7))
      			out<-c(out,F1NPA)
      			out<-c(out,tes7)

    			}### if,else loop

                  }else{


                        if((Np >= 60))
			{
 				mzV<-tes2
  				intenV<-tes3
  
  				nMZV=cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]
 				nINV=cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]
  				tes7<-paste(nMZV,nINV,sep="\t")
  				F1NPA<-paste0("Num Peaks: ",length(tes7))
  				out<-c(out,F1NPA)
  				out<-c(out,tes7)
  
			}else{
  
  
  				################################
  				F1NPA<-paste0("Num Peaks: ",tryCatch({length(Fpea2)},error=function(cond){message("Fpea is empty")}))
  				out<-c(out,F1NPA)
				out<-c(out,Fpea2)
  				##out<-c(out,Fpea)
  				################################
  
  
			}

                  }

                #################################################
	        #################################################
		    } ### end of if loop and else Np >= 120
	        #################################################	  
                } ## end of else
             #################################### adding this } new
             ###################################################
              }else{
                ###########################################    
                print("enter the else part...3")
                ###########################################
                ###########################################
	        if((as.character(InMEDA[["Adduct"]]) == "[M]+" || as.character(InMEDA[["Adduct"]]) == "[M]-") & (F2NPA != 0)){      
                  ###########################################
                  ###########################################
                  print("enter the else part...3..if loop")
                  ###########################
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
		  NAT<-as.character(InMEDA[["Adduct"]])
		  NAT1<-paste("PRECURSORTYPE:",NAT,sep=" ")
		  out<-c(out,NAT1)
		  ###############################
                  #################################
                  FIN1<-InMEDA[["Ionization mode"]]
                  F1IN1<-as.character(FIN1)
                  F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                  out<-c(out,F2IN1)
                  #################Adding this new##################
		  ##################################################
		  InKeyVal<-tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
                  ########################
                  IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
                  ######################
                  if(!sjmisc::is_empty(IKCRV)){
                    ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
                    F1ONT<-paste("Ontology:",ONTV,sep=" ")
                    out<-c(out,F1ONT)
                  }else{
                    ##F1ONT<-paste("Ontology:","",sep=" ")
		    F1ONT<-tryCatch({MaKE.ONT.REC(InMEDA)[1]},error=function(cond){message("Classifier could not fecth the information")})
                    out<-c(out,F1ONT)
                  }
                  ########################
                  FINK<-paste("INCHIKEY:",tryCatch({InKeyVal},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                  out<-c(out,tryCatch({SelectInchiorInchiKey(InMEDA,IK,InchiV)[1]},error=function(cond){return(NA)}))
		  ###out<-c(out,FINK)
                  ############################## 
		  #################################
		  INV3<-tryCatch({rinchi::get.inchi(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
                  FINCH<-paste("INCHI:",INV3,sep=" ")
                  ####out<-c(out,FINCH)
		  out<-c(out,tryCatch({SelectInchiorInchiKey(InMEDA,IK,InchiV)[2]},error=function(cond){return(NA)}))
                  ###########################
                  ###########################
                  ###FSIM<-paste("SMILES:",gETSmiles(InMEDA),sep=" ")
		  FSIM<-paste("SMILES:",tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}),sep=" ")
                  out<-c(out,FSIM)
		  ##########################
		  print("check the smiles value in second function")
		  print(FSIM)
                  ############################
		  FFOR=ifelse(!sjmisc::is_empty(gETSmiles(InMEDA)),ifelse(!sjmisc::is_empty(tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)})),tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})),ifelse(!sjmisc::is_empty(InMEDA[["InChI"]]),tryCatch({getCactus(InMEDA[["InChI"]], "formula")},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})))
		  ##################################
		  ###################################################
                  FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                  out<-c(out,FFOR1) 
                  ###################################################		
                  ###################################################
                  FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                  FINS1<-FNA[FINS]
		  FINS2<-ifelse(!sjmisc::is_empty(stringr::str_trim(gsub("INTENSITY:","",FINS1))),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
                  ###FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
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
                  FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
                  out<-c(out,FCOM)
                  ###############################
		  ###############################
		  PT1<-c("Num Peaks:","NumPeaks:")
                  FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
                  F1NPA<-FNA[FNPA]
                  F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
                  Np<-F2NPA
                  Fpea<-FNA[(FNPA+1):Find]
                  Fpea1<-Fpea[!is.na(Fpea)]
                  Fpea2<-PeakVali(Fpea1)
		  ############################
                  ############################
		  if(F2NPA != 0)
		  ##############################
                  {
                    ################################
                    MV=AAMS1
		    ################################# 
		    tes1<-unlist(strsplit(Fpea, "\t|\t\t"))
                    tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
                    tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
                    tes4<-which(tes2 > (3+MV))
		    NTES4<-which(tes3 <= 100)
		    NTES<-c(NTES4,tes4)
		    NTES5<-which(tes3 <= 0.005*(mean(tes3[order(tes3, decreasing =TRUE)][1:10])))

		    if((Np >= 60) & (length(NTES) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))
  		     {
			NTES1<-c(NTES4,tes4)     

		        if((Np >= 60) & (length(NTES1) > 1))
			{
				tes5<-tes2[-NTES1]
				tes6<-tes3[-NTES1]

				mzV<-tes5
				intenV<-tes6

				if(length(mzV) > 1)
  				{

    					print("enter inside if loop")
    					nMZV=tryCatch({cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]},error=function(cond){return(0)})
    					nINV=tryCatch({cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]},error=function(cond){return(0)})
    					tes7<-paste(nMZV,nINV,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))

    					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  				}else if(length(tes2[-NTES]) > 1){



    					tes5<-tes2[-NTES]
    					tes6<-tes3[-NTES]
    					tes7<-paste(tes5,tes6,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))
    					
					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  				}else{

    					tes5<-tes2
    					tes6<-tes3

    					nMZV=tryCatch({cenPeaks(tes2,tes3, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]},error=function(cond){return(0)})
    					nINV=tryCatch({cenPeaks(tes2,tes3, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]},error=function(cond){return(0)})

    					if(length(nMZV) > 1 & length(nINV) > 1)
    					{

      						print("enter the value of nMZV and nINV")
      						tes7<-paste(nMZV,nINV,sep="\t")

      						F1NPA<-paste0("Num Peaks: ",length(tes7))
      						out<-c(out,F1NPA)
      						out<-c(out,tes7)

    					}else{

    						tes7<-paste(tes5,tes6,sep="\t")
    						F1NPA<-paste0("Num Peaks: ",length(tes7))
    						out<-c(out,F1NPA)
    						out<-c(out,tes7)
    					}## end of if else

  					}


			}else{
                        

				if((length(tes2[-NTES]) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))
  				{
    					tes5<-tes2[-NTES]
    					tes6<-tes3[-NTES]
    					tes7<-paste(tes5,tes6,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))
    					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  				}else{

    					tes7<-paste(tes2,tes3,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))
    					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  				}## end of if else loop


			}### end of the else


  		    }else{

		    if((length(NTES) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))	    
                    {

			    if(length(tes2[-NTES]) > 1)
    			    {

      				tes5<-tes2[-NTES]
      				tes6<-tes3[-NTES]
      				tes7<-paste(tes5,tes6,sep="\t")
     				### Adding this new here
      				F1NPA<-paste0("Num Peaks: ",length(tes7))
      				out<-c(out,F1NPA)
      				#############################
      				out<-c(out,tes7)
      				##############################
    			   }else{

      					tes5<-tes2
      					tes6<-tes3
      					tes7<-paste(tes5,tes6,sep="\t")
      					F1NPA<-paste0("Num Peaks: ",length(tes7))
      					out<-c(out,F1NPA)
      					out<-c(out,tes7)

   				 }### if,else loop


                    }else{

			    ###if((Np >= 80))
			    if((Np >= 60))
    			    {
      				mzV<-tes2
      				intenV<-tes3

      				nMZV=cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]
      				nINV=cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]

      				tes7<-paste(nMZV,nINV,sep="\t")
      				F1NPA<-paste0("Num Peaks: ",length(tes7))

      				out<-c(out,F1NPA)
      				out<-c(out,tes7)


    			 }else{


		      ################################
		      F1NPA<-paste0("Num Peaks: ",tryCatch({length(Fpea2)},error=function(cond){message("Fpea is empty")}))		 
		      ###F1NPA<-paste0("Num Peaks: ",tryCatch({length(Fpea)},error=function(cond){message("Fpea is empty")}))
                      out<-c(out,F1NPA)
		      out<-c(out,Fpea2)
                      ###out<-c(out,Fpea)	    
		      ################################

			    }## end of else and Np >= 80
                #########################################
                    }
	        ########################################
		     }### end of if loop and Np >= 120 
	        ########################################	    
                  }
                ############################
                }
              #############################
              #############################
              }
        #########################################################
            }else if(length(INLL) > 1){
              ##############################################
              print("enter the line ...1020")
              ##############################################
              MONMS=InMSPL[INLL]
	      ##########################################
              #########################
              TRA<-unname(rapply(MONMS, function(x) grep("RETENTIONTIME:",x, value=TRUE)))
              TRA1<-as.numeric(stringr::str_trim(stringr::str_replace(TRA, "RETENTIONTIME:", "")))
              TRA2<-abs(VRT-TRA1)
              ##TRA3<-which.min(TRA2)
	      TRA3<-which(TRA2==min(TRA2[TRA2 > 0]))
	      ###################################
	      ##########################################
	      PT1<-c("Num Peaks:","NumPeaks:")
	      NPV1=unname(rapply(MONMS, function(x) grep(paste(PT1,collapse="|"), x, value=TRUE)))
              NPV=as.numeric(stringr::str_trim(stringr::str_replace(NPV1, "Num Peaks:|NumPeaks:", "")))
              #####NPV=as.numeric(stringr::str_trim(gsub("Num Peaks:","",unname(sapply(MONMS, `[[`, unique(MONMS))))))
	      NPNZ=which(NPV!=0)
	      NTRA3<-intersect(TRA3,NPNZ)
	      TRA4<-INLL[NTRA3]
              ##########################
              TRA5<-InMSPL[TRA4]
              ###########################
              ###########################
              F1FPL<-TRA5
              F2FPL<-F1FPL
              ##########################
              Find<-tryCatch({length(F2FPL[[1]])},error=function(cond){message("there is an error in list")})
              ##########################
              FNA<-tryCatch({F1FPL[[1]]},error=function(cond){message("there is an error in list")})
              ######################## adding this new
	      #################################
	      FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
              F1NPA<-FNA[FNPA]
              F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
              ##############################
              ###########################################
              PT <- c("PRECURSORTYPE:", "ADDUCTIONNAME:")
              NPT<-tryCatch({grep(paste(PT,collapse="|"), FNA, value=TRUE)},error=function(cond){message("List value is empty")})
              NPT1<-tryCatch({match(NPT,FNA)},error=function(cond){message("List value is empty")})
              PTV <-tryCatch({stringr::str_remove(FNA[NPT1],c("PRECURSORTYPE:","ADDUCTIONNAME:"))},error=function(cond){message("List value is empty")})
              PTV1<-tryCatch({PTV[1]},error=function(cond){message("List value is empty")})
              PTV2<-tryCatch({stringr::str_trim(PTV1)},error=function(cond){message("List value is empty")})
              P1TV2<-tryCatch({stringr::str_trim(gsub("ADDUCTIONNAME:","",PTV2))},error=function(cond){message("List value is empty")})
              PTV3<-tryCatch({paste("PRECURSORTYPE:",P1TV2)},error=function(cond){message("List value is empty")})
              ##PTV3<-tryCatch({paste("PRECURSORTYPE:",P1TV2)},error=function(cond){message("List value is empty")})
              ########################################
              print("enter the test area..4")
              ########################################
              print(P1TV2)
              print(as.character(InMEDA[["Adduct"]]))
              ################################## adding this new
	      if((!sjmisc::is_empty(P1TV2) || !sjmisc::is_empty(as.character(InMEDA[["Adduct"]]))) & (F2NPA != 0)){
                ###################################################	      
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
		###########################
		NAT<-as.character(InMEDA[["Adduct"]])
                NAT1<-paste("PRECURSORTYPE:",NAT,sep=" ")
                out<-c(out,NAT1)
                ###########################
                FIN1<-InMEDA[["Ionization mode"]]
                F1IN1<-as.character(FIN1)
                F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                out<-c(out,F2IN1)
                #####################################################
		#####################################################
		InKeyVal<-tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
                ##############################
                IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
                ##################################
                if(!sjmisc::is_empty(IKCRV)){
                  ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
                  F1ONT<-paste("Ontology:",ONTV,sep=" ")
                  out<-c(out,F1ONT)
                }else{
                  ##F1ONT<-paste("Ontology:","",sep=" ")
		  F1ONT<-tryCatch({MaKE.ONT.REC(InMEDA)[1]},error=function(cond){message("Classifier could not fecth the information")})	
                  out<-c(out,F1ONT)
                }
                ###############################
               
             
                #############################
		out<-c(out,tryCatch({SelectInchiorInchiKey(InMEDA,IK,InchiV)[1]},error=function(cond){return(NA)}))
		###out<-c(out,paste("INCHIKEY:",InKeyVal,sep=" "))
		###########################
                INV3<-tryCatch({rinchi::get.inchi(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
                FINCH<-paste("INCHI:",INV3,sep=" ")
                ####out<-c(out,FINCH)
		out<-c(out,tryCatch({SelectInchiorInchiKey(InMEDA,IK,InchiV)[2]},error=function(cond){return(NA)}))
                ######################################
                ######################################
		FSIM<-paste("SMILES:",tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}),sep=" ")
                out<-c(out,FSIM)
		###################################
		print("enter the smiles value in second function")
		print(FSIM)
                ############################################
		FFOR<-ifelse(!sjmisc::is_empty(gETSmiles(InMEDA)),tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){message("smile value is empty")}),tryCatch({InMEDA[["Formula"]]},error=function(cond){message("Formula is empty")}))
		##########################################
                FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                out<-c(out,FFOR1)
                ############################################
                ############################################
                FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                FINS1<-FNA[FINS]
		FINS2<-ifelse(!sjmisc::is_empty(stringr::str_trim(gsub("INTENSITY:","",FINS1))),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
                ###FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
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
                #####################################
		#####################################
                FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
                F1NPA<-FNA[FNPA]
                F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
                Np<-F2NPA
                Fpea<-FNA[(FNPA+1):Find]
                Fpea1<-Fpea[!is.na(Fpea)]
                ##Fpea2<-PeakVali(Fpea)
                Fpea2<-PeakVali(Fpea1)
                ########################################
		########################################
		if(F2NPA != 0)
		################################
                {
		  ##############################
		  ##############################################	
		  tes1<-unlist(strsplit(Fpea, "\t|\t\t"))
                  tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
                  tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
                  tes4<-which(tes2 > (3+MV))
		  NTES4<-which(tes3 <= 100)
		  NTES<-c(NTES4,tes4)


		  if((Np >= 60) & (length(NTES) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))
  		  {
    			NTES1<-c(NTES4,tes4,NTES5)
                        

                        #############################################
		        if((Np >= 60) & (length(NTES1) > 1))
			{
				tes5<-tes2[-NTES1]
				tes6<-tes3[-NTES1]

				mzV<-tes5
				intenV<-tes6

				if(length(mzV) > 1)
  				{

    					print("enter inside if loop")
    					nMZV=tryCatch({cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]},error=function(cond){return(0)})
    					nINV=tryCatch({cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]},error=function(cond){return(0)})
    					tes7<-paste(nMZV,nINV,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))

    					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  				}else{

    					tes5<-tes2
   					tes6<-tes3

    					nMZV=tryCatch({cenPeaks(tes2,tes3, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]},error=function(cond){return(0)})
   					nINV=tryCatch({cenPeaks(tes2,tes3, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]},error=function(cond){return(0)})

    					if(length(nMZV) > 1 & length(nINV) > 1)
    					{

      						print("enter the value of nMZV and nINV")
      						tes7<-paste(nMZV,nINV,sep="\t")

      						F1NPA<-paste0("Num Peaks: ",length(tes7))
      						out<-c(out,F1NPA)
      						out<-c(out,tes7)

    					}else{

    						tes7<-paste(tes5,tes6,sep="\t")
    						F1NPA<-paste0("Num Peaks: ",length(tes7))
    						out<-c(out,F1NPA)
    						out<-c(out,tes7)
    					}## end of if else

  				}



			}else{


				if(length(tes2[-NTES]) > 1)
  				{
    					tes5<-tes2[-NTES]
    					tes6<-tes3[-NTES]
    					tes7<-paste(tes5,tes6,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))
    					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  				}else{

    					tes7<-paste(tes2,tes3,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))
    					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  				}## end of if else loop

			}### end of else


  		   }else{

		  if((length(NTES) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))	   
                  {

			  if(length(tes2[-NTES]) > 1)
    			  {

      				tes5<-tes2[-NTES]
      				tes6<-tes3[-NTES]
      				tes7<-paste(tes5,tes6,sep="\t")
     				### Adding this new here
      				F1NPA<-paste0("Num Peaks: ",length(tes7))
      				out<-c(out,F1NPA)
      				#############################
      				out<-c(out,tes7)
      				##############################
    			}else{

      				tes5<-tes2
      				tes6<-tes3
      				tes7<-paste(tes5,tes6,sep="\t")
      				F1NPA<-paste0("Num Peaks: ",length(tes7))
      				out<-c(out,F1NPA)
      				out<-c(out,tes7)

    			}### if,else loop
		#############################
                  }else{
	        #############################

	            if((Np >= 60))		  
    		    {
      			mzV<-tes2
      			intenV<-tes3
      
      			nMZV=cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]
      			nINV=cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]
      
      			tes7<-paste(nMZV,nINV,sep="\t")
      			F1NPA<-paste0("Num Peaks: ",length(tes7))
      
      			out<-c(out,F1NPA)
      			out<-c(out,tes7)
      
      
    		   }else{	  


		    ####################################
		    F1NPA<-paste0("Num Peaks: ",tryCatch({length(Fpea2)},error=function(cond){message("Fpea is empty")}))	   
	            ###F1NPA<-paste0("Num Peaks: ",tryCatch({length(Fpea)},error=function(cond){message("Fpea is empty")}))
		    out<-c(out,F1NPA)
		    out<-c(out,Fpea2)
                    ##out<-c(out,Fpea)
		    ####################################
		    }## end of else and Np >= 80
		 ######################################## 
                  }
	       #########################################
		  }### end of if loop and Np >= 120
	       #########################################	  
                } ### this is the else MV= ADV4...closing
             #############################################
             #############################################
              }else{
                #####################################     
                print("enter the else part...4")
                #####################################
                #####################################
	        if((as.character(InMEDA[["Adduct"]]) == "[M]+" || as.character(InMEDA[["Adduct"]]) == "[M]-") & (F2NPA != 0)){
                ################################
                  #####################################
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
		  NAT<-as.character(InMEDA[["Adduct"]])
                  NAT1<-paste("PRECURSORTYPE:",NAT,sep=" ")
                  out<-c(out,NAT1)
		  ####################################
                  #################################
                  FIN1<-InMEDA[["Ionization mode"]]
                  F1IN1<-as.character(FIN1)
                  F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                  out<-c(out,F2IN1)
                  ########################################
		  ########################################
		  InKeyVal<-tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
                  ###########################################
                  IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
                  #########################
                  if(!sjmisc::is_empty(IKCRV)){
                    ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
                    F1ONT<-paste("Ontology:",ONTV,sep=" ")
                    out<-c(out,F1ONT)
                  }else{
			  F1ONT<-tryCatch({MaKE.ONT.REC(InMEDA)[1]},error=function(cond){message("Classifier could not fecth the information")})
                    ##F1ONT<-paste("Ontology:","",sep=" ")
                    out<-c(out,F1ONT)
                  }
                  ####################
                  FINK<-paste("INCHIKEY:",tryCatch({InKeyVal},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                  out<-c(out,tryCatch({SelectInchiorInchiKey(InMEDA,IK,InchiV)[1]},error=function(cond){return(NA)}))
		  ####out<-c(out,FINK)
                  ##################################
		  ##################################
		  INV3<-tryCatch({rinchi::get.inchi(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
                  FINCH<-paste("INCHI:",INV3,sep=" ")
		  out<-c(out,tryCatch({SelectInchiorInchiKey(InMEDA,IK,InchiV)[2]},error=function(cond){return(NA)}))
                  ###out<-c(out,FINCH)
                  ########################
                  ########################
                  FSIM<-paste("SMILES:",tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}),sep=" ")
                  out<-c(out,FSIM)
                  ##############################
                  #######################################################
		  FFOR=ifelse(!sjmisc::is_empty(gETSmiles(InMEDA)),ifelse(!sjmisc::is_empty(tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)})),tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})),ifelse(!sjmisc::is_empty(InMEDA[["InChI"]]),tryCatch({getCactus(InMEDA[["InChI"]], "formula")},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})))
		  FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                  out<-c(out,FFOR1)
                  ###################################################		
                  ###################################################
                  FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                  FINS1<-FNA[FINS]
		  FINS2<-ifelse(!sjmisc::is_empty(stringr::str_trim(gsub("INTENSITY:","",FINS1))),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
                  ######FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
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
                  FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
                  out<-c(out,FCOM)
                  ##################
                  ##################################
		  PT1<-c("Num Peaks:","NumPeaks:")
                  FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
                  F1NPA<-FNA[FNPA]
                  F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
                  Np<-F2NPA
                  Fpea<-FNA[(FNPA+1):Find]
                  Fpea1<-Fpea[!is.na(Fpea)]
                  Fpea2<-PeakVali(Fpea1)
                  #########################
		  if(F2NPA != 0)
		  ###########################
                  {
                    
                    MV=AAMS1
		    tes1<-unlist(strsplit(Fpea, "\t|\t\t"))
                    tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
                    tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
                    tes4<-which(tes2 > (3+MV))
		    NTES4<-which(tes3 <= 100)
		    ###NTES5<-which(tes3 <= 0.05*max(tes3))
		    NTES5<-which(tes3 <= 0.005*(mean(tes3[order(tes3, decreasing =TRUE)][1:10])))
		    NTES<-c(NTES4,tes4)

                    ###if(Np >= 60)
		    if((Np >= 60) & (length(NTES) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))
  		    {
    			###NTES1<-c(NTES4,tes4,NTES5)
			NTES1<-c(NTES4,tes4)    


		        if((Np >= 60) & (length(NTES1) > 1))
			{
				tes5<-tes2[-NTES1]
				tes6<-tes3[-NTES1]

				mzV<-tes5
				intenV<-tes6

				if(length(mzV) > 1)
  				{

    					nMZV=tryCatch({cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]},error=function(cond){return(0)})
    					nINV=tryCatch({cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]},error=function(cond){return(0)})
    					tes7<-paste(nMZV,nINV,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))

    					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  				}else{

    					tes5<-tes2
    					tes6<-tes3

    					nMZV=tryCatch({cenPeaks(tes2,tes3, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]},error=function(cond){return(0)})
    					nINV=tryCatch({cenPeaks(tes2,tes3, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]},error=function(cond){return(0)})

    					if(length(nMZV) > 1 & length(nINV) > 1)
    					{

      						tes7<-paste(nMZV,nINV,sep="\t")

      						F1NPA<-paste0("Num Peaks: ",length(tes7))
      						out<-c(out,F1NPA)
      						out<-c(out,tes7)

    					}else{

    						tes7<-paste(tes5,tes6,sep="\t")
    						F1NPA<-paste0("Num Peaks: ",length(tes7))
    						out<-c(out,F1NPA)
    						out<-c(out,tes7)
    					}## end of if else

  					}



			}else{


				if((length(tes2[-NTES]) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))
  				{
    					tes5<-tes2[-NTES]
    					tes6<-tes3[-NTES]
    					tes7<-paste(tes5,tes6,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))
    					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  				}else{

    					tes7<-paste(tes2,tes3,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))
    					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  				}## end of if else loop

			}## end of the else

  		   }else{

		    if((length(NTES) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))	   
                    {
		      ###############################################

			    if(length(tes2[-NTES]) > 1)
    			    {

      				tes5<-tes2[-NTES]
      				tes6<-tes3[-NTES]
      				tes7<-paste(tes5,tes6,sep="\t")
      				### Adding this new here
      				F1NPA<-paste0("Num Peaks: ",length(tes7))
      				out<-c(out,F1NPA)
      				#############################
      				out<-c(out,tes7)
      				##############################
    			}else{

      				tes5<-tes2
      				tes6<-tes3
      				tes7<-paste(tes5,tes6,sep="\t")
      				F1NPA<-paste0("Num Peaks: ",length(tes7))
      				out<-c(out,F1NPA)
      				out<-c(out,tes7)

    			  }### if,else loop
                    }else{

			    ###if((Np >= 80))
			    if((Np >= 60))
    			    {
      				mzV<-tes2
      				intenV<-tes3

      				nMZV=cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]
      				nINV=cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]

      				tes7<-paste(nMZV,nINV,sep="\t")
      				F1NPA<-paste0("Num Peaks: ",length(tes7))

      				out<-c(out,F1NPA)
      				out<-c(out,tes7)


    			  }else{

		       ##########################
		       F1NPA<-paste0("Num Peaks: ",tryCatch({length(Fpea2)},error=function(cond){message("Fpea is empty")}))		  
                       out<-c(out,F1NPA)
		       out<-c(out,Fpea2)
                       ###out<-c(out,Fpea)	    
		       ######## Adding this new
		        }## emd of else and Np >= 80

                   ###############################
                    }
		 ################################
		     }### end of the if loop and Np >= 120
		 ################################   
                  }
               ###############################
                }
              ##########################
              ##########################
              }
           ##########################
	   ########################## 
            } else {
              #print("entering the line 764")
              PASS<-InMEDA
            }  ## main else part
            ###testing if this works
            return(out)
          #####################
          } #ITmass
        } #ITRTL
     ###################
      } ## adduct match and exacct mass present
    } else{
      
      ###############################################################
      #### this entered the FM does not exists 
      FM2<-tryCatch({FuFtoRe(InMEDA)},error=function(cond){message("FuFtoRe is missing")})
      ###############################################################
      ###############################################################
      if(!sjmisc::is_empty(FM2))
      {
        ###################################################################
        print("entering the line ... 772")
        ###################################################################
        ####################################################################
        if(!sjmisc::is_empty(FM2) & !sjmisc::is_empty(as.character(InMEDA[["Adduct"]]))){
          #############################################
          print("entering the line ... 779")
          #############################################
          ADV<-as.character(InMEDA[["Adduct"]])
          ADV1<-tryCatch({qdapRegex::ex_between(ADV, "[", "]")[[1]]},error=function(cond){message("adduct value is missing")})
          ADV2<-tryCatch({FADINF(ADV)},error=function(cond){message("adduct value matching is not found")})
          ADV3<-tryCatch({stringr::str_replace(ADV2, "M",as.character(FM))},error=function(cond){message("adduct value replacement is not found")})
          ADV4<-tryCatch({as.numeric(pander::evals(ADV3)[[1]]$result)},error=function(cond){message("getting the result")})
          ######################
          ##PPm=ADV4*(25/(1000000))
	  PPm=ADV4*(mz_Tol/(1000000))
          ######################
          MPPmL=ADV4-PPm
          MPPmU=ADV4+PPm
          #######################
          Tmass<-InPMZ[InPMZ >= MPPmL & InPMZ <= MPPmU]
          ITmass<-which(InPMZ %in% Tmass)
          ##ITmass<-match(Tmass,InPMZ)
          #####################
          VRT<-as.numeric(as.character(InMEDA[["RT (min)"]]))
	  VRTL<-VRT-RT_Tol
	  VRTU<-VRT+RT_Tol
	  ######################
          #############################
          TRTL<-InRTL[InRTL >= VRTL & InRTL <= VRTU]
          ITRTL<-which(InRTL %in% TRTL)
          #############################
          ##ITRTL<-match(TRTL,InRTL)
          ##############################
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
		###############################################
		PT1<-c("Num Peaks:","NumPeaks:")
		FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
                F1NPA<-FNA[FNPA]
                F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
                ####################################################
                PT <- c("PRECURSORTYPE:", "ADDUCTIONNAME:")
                NPT<-tryCatch({grep(paste(PT,collapse="|"), FNA, value=TRUE)},error=function(cond){message("List value is empty")})
                NPT1<-tryCatch({match(NPT,FNA)},error=function(cond){message("List value is empty")})
                PTV <-tryCatch({stringr::str_remove(FNA[NPT1],c("PRECURSORTYPE:","ADDUCTIONNAME:"))},error=function(cond){message("List value is empty")})
                PTV1<-tryCatch({PTV[1]},error=function(cond){message("List value is empty")})
                PTV2<-tryCatch({stringr::str_trim(PTV1)},error=function(cond){message("List value is empty")})
                P1TV2<-tryCatch({stringr::str_trim(gsub("ADDUCTIONNAME:","",PTV2))},error=function(cond){message("List value is empty")})
                PTV3<-tryCatch({paste("PRECURSORTYPE:",P1TV2)},error=function(cond){message("List value is empty")})
                ##PTV3<-tryCatch({paste("PRECURSORTYPE:",PTV2)},error=function(cond){message("List value is empty")})
                #####################################################
                print("enter the test area..3")
                ####################################################
                print(P1TV2)
                print(as.character(InMEDA[["Adduct"]]))
                ################################ adding this new ###
                ##if(identical(P1TV2,as.character(InMEDA[["Adduct"]]))){
                ####################################################
		if((!sjmisc::is_empty(P1TV2) || !sjmisc::is_empty(as.character(InMEDA[["Adduct"]]))) & (F2NPA != 0)){
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
		  #################################
		  NAT<-as.character(InMEDA[["Adduct"]])
                  NAT1<-paste("PRECURSORTYPE:",NAT,sep=" ")
                  out<-c(out,NAT1)
                  #################################
                  FIN1<-InMEDA[["Ionization mode"]]
                  F1IN1<-as.character(FIN1)
                  F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                  out<-c(out,F2IN1)
                  ###############################################
                  ###############################################
		  InKeyVal<-tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
                  #####################################
                  #####################################
                  IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
                  #################################################
		  #################################################
                  if(!sjmisc::is_empty(IKCRV)){
                    ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
                    F1ONT<-paste("Ontology:",ONTV,sep=" ")
                    out<-c(out,F1ONT)
                  }else{
                    ##F1ONT<-paste("Ontology:","",sep=" ")
			  F1ONT<-tryCatch({MaKE.ONT.REC(InMEDA)[1]},error=function(cond){message("Classifier could not fecth the information")})
                    out<-c(out,F1ONT)
                  }
                  ######################################################
		  ######################################################
                  FINK<-paste("INCHIKEY:",tryCatch({InKeyVal},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                  out<-c(out,tryCatch({SelectInchiorInchiKey(InMEDA,IK,InchiV)[1]},error=function(cond){return(NA)}))
		  ###out<-c(out,FINK)
                  ################################
		  ##################################
		  INV3<-tryCatch({rinchi::get.inchi(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
                  FINCH<-paste("INCHI:",INV3,sep=" ")
		  out<-c(out,tryCatch({SelectInchiorInchiKey(InMEDA,IK,InchiV)[2]},error=function(cond){return(NA)}))
                  ###out<-c(out,FINCH)
                  ################################
                  ################################
		  FSIM<-paste("SMILES:",tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}),sep=" ")
                  out<-c(out,FSIM)
                  ################################
		 ################################################
		  FFOR=ifelse(!sjmisc::is_empty(gETSmiles(InMEDA)),ifelse(!sjmisc::is_empty(tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)})),tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})),ifelse(!sjmisc::is_empty(InMEDA[["InChI"]]),tryCatch({getCactus(InMEDA[["InChI"]], "formula")},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})))
		  #################################
		  #####################################
                  FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                  out<-c(out,FFOR1)
                  ###############################
                  ###############################
                  FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                  FINS1<-FNA[FINS]
		  FINS2<-ifelse(!sjmisc::is_empty(stringr::str_trim(gsub("INTENSITY:","",FINS1))),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
                  ######FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
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
		  ################################
		  PT1<-c("Num Peaks:","NumPeaks:")
                  FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
                  F1NPA<-FNA[FNPA]
                  F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
                  Np<-F2NPA
                  Fpea<-FNA[(FNPA+1):Find]
                  Fpea1<-Fpea[!is.na(Fpea)]
                  Fpea2<-PeakVali(Fpea1)
                  ######################################
		  ######################################
		  if(F2NPA != 0)
	          ##########################################		  
                  {
		    ####################
                    MV=ADV4
		    tes1<-unlist(strsplit(Fpea2, "\t|\t\t|;| +"))
                    tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
                    tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
                    tes4<-which(tes2 > (3+MV))
		    NTES4<-which(tes3 <= 100)
		    NTES<-c(NTES4,tes4)
		    NTES5<-which(tes3 <= 0.005*(mean(tes3[order(tes3, decreasing =TRUE)][1:10])))

		    if((Np >= 60) & (length(NTES) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))
  		     {
			NTES1<-c(NTES4,tes4)     




		        if((Np >= 60) & (length(NTES1) > 1))
			{
				tes5<-tes2[-NTES1]
				tes6<-tes3[-NTES1]

				mzV<-tes5
				intenV<-tes6

				if(length(mzV) > 1)
  				{

    					print("enter inside if loop")
    					nMZV=tryCatch({cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]},error=function(cond){return(0)})
    					nINV=tryCatch({cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]},error=function(cond){return(0)})
    					tes7<-paste(nMZV,nINV,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))

    					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  				}else{

    					tes5<-tes2
    					tes6<-tes3

    					nMZV=tryCatch({cenPeaks(tes2,tes3, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]},error=function(cond){return(0)})
    					nINV=tryCatch({cenPeaks(tes2,tes3, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]},error=function(cond){return(0)})

    					if(length(nMZV) > 1 & length(nINV) > 1)
    					{

      						print("enter the value of nMZV and nINV")
      						tes7<-paste(nMZV,nINV,sep="\t")

      						F1NPA<-paste0("Num Peaks: ",length(tes7))
      						out<-c(out,F1NPA)
      						out<-c(out,tes7)

    					}else{

    						tes7<-paste(tes5,tes6,sep="\t")
    						F1NPA<-paste0("Num Peaks: ",length(tes7))
    						out<-c(out,F1NPA)
    						out<-c(out,tes7)

    					}## end of if else

 				 }


			}else{


				if((length(tes2[-NTES]) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))
  				{
    					tes5<-tes2[-NTES]
    					tes6<-tes3[-NTES]
    					tes7<-paste(tes5,tes6,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))
    					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  				}else{

    					tes7<-paste(tes2,tes3,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))
    					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  				}## end of if else loop

			}### end of the else


  		     }else{


		    if((length(NTES) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))	     
                    {

			    if(length(tes2[-NTES]) > 1)
    			    {

      				tes5<-tes2[-NTES]
      				tes6<-tes3[-NTES]
      				tes7<-paste(tes5,tes6,sep="\t")
      				### Adding this new here
      				F1NPA<-paste0("Num Peaks: ",length(tes7))
      				out<-c(out,F1NPA)
      				#############################
      				out<-c(out,tes7)
      				##############################
    			}else{

      				tes5<-tes2
      				tes6<-tes3
      				tes7<-paste(tes5,tes6,sep="\t")
      				F1NPA<-paste0("Num Peaks: ",length(tes7))
      				out<-c(out,F1NPA)
      				out<-c(out,tes7)

    				}### if,else loop
                    }else{

		       ##########################	    

		       if((Np >= 60))
	    
    			{
      			 	mzV<-tes2
      				intenV<-tes3
      
      				nMZV=cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]
      				nINV=cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]
      
      				tes7<-paste(nMZV,nINV,sep="\t")
      				F1NPA<-paste0("Num Peaks: ",length(tes7))
      
      				out<-c(out,F1NPA)
      				out<-c(out,tes7)
      
      
    			}else{	    


		       #########################################
		       F1NPA<-paste0("Num Peaks: ",tryCatch({length(Fpea2)},error=function(cond){message("Fpea is empty")}))		
                       out<-c(out,F1NPA)
		       out<-c(out,Fpea2)
		       #########################################
			}## end of else and Np >= 80
	           ########################################	       
                    }
		 ##########################################
		     }### end of if loop and Np >= 120
		 #########################################   
                  } ## end of else
                 #################################### adding this } new
                 ##########################################
                }else{
                  #########################################    
                  print("enter the else part...3")
                  ##########################################
                  ###########################################
	          if((as.character(InMEDA[["Adduct"]]) == "[M]+" || as.character(InMEDA[["Adduct"]]) == "[M]-") & (F2NPA != 0)){	
                    ###########################################
                    ###########################################
                    print("enter the else part...3..if loop")
                    ###########################################
		    ###########################################
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
		    #################################
		    NAT<-as.character(InMEDA[["Adduct"]])
                    NAT1<-paste("PRECURSORTYPE:",NAT,sep=" ")
                    out<-c(out,NAT1)
                    #################################
                    FIN1<-InMEDA[["Ionization mode"]]
                    F1IN1<-as.character(FIN1)
                    F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                    out<-c(out,F2IN1)
                    #################Adding this new##################
                    ##################################################
		    InKeyVal<-tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
                    ##############################################
                    IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
                    ############################################
                    if(!sjmisc::is_empty(IKCRV)){
                      ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
                      F1ONT<-paste("Ontology:",ONTV,sep=" ")
                      out<-c(out,F1ONT)
                    }else{
                      F1ONT<-paste("Ontology:","",sep=" ")
                      out<-c(out,F1ONT)
                    }
                    ################################################
                    FINK<-paste("INCHIKEY:",tryCatch({InKeyVal},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                    out<-c(out,tryCatch({SelectInchiorInchiKey(InMEDA,IK,InchiV)[1]},error=function(cond){return(NA)}))
		    ###out<-c(out,FINK)
                    ###################################################
                    INV3<-tryCatch({rinchi::get.inchi(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
		    FINCH<-paste("INCHI:",INV3,sep=" ")
		    out<-c(out,tryCatch({SelectInchiorInchiKey(InMEDA,IK,InchiV)[2]},error=function(cond){return(NA)}))
                    ###out<-c(out,FINCH)
                    ####################################################
                    ####################################################
                    FSIM<-paste("SMILES:",tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}),sep=" ")
                    out<-c(out,FSIM)
                    ###################################################
		    ###################################################
		    FFOR=ifelse(!sjmisc::is_empty(gETSmiles(InMEDA)),ifelse(!sjmisc::is_empty(tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)})),tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})),ifelse(!sjmisc::is_empty(InMEDA[["InChI"]]),tryCatch({getCactus(InMEDA[["InChI"]], "formula")},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})))
		    ###################################################
		    #################################################
                    FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                    out<-c(out,FFOR1) 
                    ###################################################		
                    ###################################################
                    FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                    FINS1<-FNA[FINS]
		    FINS2<-ifelse(!sjmisc::is_empty(stringr::str_trim(gsub("INTENSITY:","",FINS1))),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
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
                    FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
                    out<-c(out,FCOM)
                    ##################
		    ################################
		    PT1<-c("Num Peaks:","NumPeaks:")
                    FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
                    F1NPA<-FNA[FNPA]
                    F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
                    Np<-F2NPA
                    Fpea<-FNA[(FNPA+1):Find]
                    Fpea1<-Fpea[!is.na(Fpea)]
                    Fpea2<-PeakVali(Fpea1)
                    #########################
		    if(F2NPA != 0)
		    #################################	    
                    {
                      
                      MV=AAMS1
		      tes1<-unlist(strsplit(Fpea2, "\t|\t\t|;| +"))
                      tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
                      tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
                      tes4<-which(tes2 > (3+MV))
		      NTES4<-which(tes3 <= 100)
		      NTES<-c(NTES4,tes4)
		      ###NTES5<-which(tes3 <= 0.05*max(tes3))
		      NTES5<-which(tes3 <= 0.005*(mean(tes3[order(tes3, decreasing =TRUE)][1:10])))

		      if((Np >= 60) & (length(NTES) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))
		      ###if((Np >= 80) & (length(NTES1) > 1))
  			{
				NTES1<-c(NTES4,tes4)

                        

			if((Np >= 60) & (length(NTES1) > 1))	
			{
				tes5<-tes2[-NTES1]
				tes6<-tes3[-NTES1]

				mzV<-tes5
				intenV<-tes6

				if(length(mzV) > 1)
  				{

    					print("enter inside if loop")
    					nMZV=tryCatch({cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]},error=function(cond){return(0)})
    					nINV=tryCatch({cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]},error=function(cond){return(0)})
    					tes7<-paste(nMZV,nINV,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))

    					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  				}else{

    					tes5<-tes2
    					tes6<-tes3

    					nMZV=tryCatch({cenPeaks(tes2,tes3, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]},error=function(cond){return(0)})
    					nINV=tryCatch({cenPeaks(tes2,tes3, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]},error=function(cond){return(0)})

    					if(length(nMZV) > 1 & length(nINV) > 1)
    					{

      						print("enter the value of nMZV and nINV")
      						tes7<-paste(nMZV,nINV,sep="\t")

      						F1NPA<-paste0("Num Peaks: ",length(tes7))
      						out<-c(out,F1NPA)
      						out<-c(out,tes7)

    					}else{

    						tes7<-paste(tes5,tes6,sep="\t")
    						F1NPA<-paste0("Num Peaks: ",length(tes7))
    						out<-c(out,F1NPA)
    						out<-c(out,tes7)
    					}## end of if else

  				}


			}else{
    				###out<-c(out,tes7)
				 if((length(tes2[-NTES]) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))
  				 {
    					tes5<-tes2[-NTES]
    					tes6<-tes3[-NTES]
    					tes7<-paste(tes5,tes6,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))
    					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  				}else{

    					tes7<-paste(tes2,tes3,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))
    					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  					}## end of if else loop
                         }### end of the else

  			}else{
		      
		      if((length(NTES) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))		
                      {

			          if(length(tes2[-NTES]) > 1)
    				  {

      					tes5<-tes2[-NTES]
      					tes6<-tes3[-NTES]
      					tes7<-paste(tes5,tes6,sep="\t")
     					### Adding this new here
      					F1NPA<-paste0("Num Peaks: ",length(tes7))
      					out<-c(out,F1NPA)
      					#############################
      					out<-c(out,tes7)
      					##############################
    				}else{

      					tes5<-tes2
      					tes6<-tes3
      					tes7<-paste(tes5,tes6,sep="\t")
      					F1NPA<-paste0("Num Peaks: ",length(tes7))
      					out<-c(out,F1NPA)
      					out<-c(out,tes7)

    				}### if,else loo


                      }else{

			if((Np >= 60))      
    			{
      				mzV<-tes2
      				intenV<-tes3
      
      				nMZV=cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]
      				nINV=cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]
      
      				tes7<-paste(nMZV,nINV,sep="\t")
      				F1NPA<-paste0("Num Peaks: ",length(tes7))
      
      				out<-c(out,F1NPA)
      				out<-c(out,tes7)
      
      
    			}else{      

			########################
			F1NPA<-paste0("Num Peaks: ",tryCatch({length(Fpea2)},error=function(cond){message("Fpea is empty")}))	
                        out<-c(out,F1NPA)
			out<-c(out,Fpea2)
                        ########################
			  }## end of else and Np >= 80
	           ###################################
                      }
		  ###################################
			}## end of if loop and Np >= 120
		  ###################################    
                    }                    		    
		  ##############################
		  ##############################  
                  }
                 #############################
                 #############################
                }
	      #########################################################	
              #########################################################
              }else if(length(INLL) > 1){
                ##############################################
                print("enter the line ...1020")
                ##############################################
                MONMS=InMSPL[INLL]
                ##############################
                TRA1<-as.numeric(stringr::str_trim(stringr::str_replace(TRA, "RETENTIONTIME:", "")))
                TRA2<-abs(VRT-TRA1)
                ##TRA3<-which.min(TRA2)
		TRA3<-which(TRA2==min(TRA2[TRA2 > 0]))
		###################################
		PT1<-c("Num Peaks:","NumPeaks:")
		NPV1=unname(rapply(MONMS, function(x) grep(paste(PT1,collapse="|"), x, value=TRUE)))
                NPV=as.numeric(stringr::str_trim(stringr::str_replace(NPV1, "Num Peaks:|NumPeaks:", "")))
		NPNZ=which(NPV!=0)
		NTRA3<-intersect(TRA3,NPNZ)
		TRA4<-INLL[NTRA3]
                ##########################
		##########################
                TRA4<-INLL[TRA3]
                ##TRA4<-INLL1[TRA3]
                TRA5<-InMSPL[TRA4]
                #############################
                ###########################
                F1FPL<-TRA5
                F2FPL<-F1FPL
                ##########################
                Find<-tryCatch({length(F2FPL[[1]])},error=function(cond){message("there is an error in list")})
                ##########################
                FNA<-tryCatch({F1FPL[[1]]},error=function(cond){message("there is an error in list")})
                ######################## adding this new##
		PT1<-c("Num Peaks:","NumPeaks:")
                FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
                F1NPA<-FNA[FNPA]
                F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
                ###########################################
                PT <- c("PRECURSORTYPE:", "ADDUCTIONNAME:")
                NPT<-tryCatch({grep(paste(PT,collapse="|"), FNA, value=TRUE)},error=function(cond){message("List value is empty")})
                NPT1<-tryCatch({match(NPT,FNA)},error=function(cond){message("List value is empty")})
                PTV <-tryCatch({stringr::str_remove(FNA[NPT1],c("PRECURSORTYPE:","ADDUCTIONNAME:"))},error=function(cond){message("List value is empty")})
                PTV1<-tryCatch({PTV[1]},error=function(cond){message("List value is empty")})
                PTV2<-tryCatch({stringr::str_trim(PTV1)},error=function(cond){message("List value is empty")})
                P1TV2<-tryCatch({stringr::str_trim(gsub("ADDUCTIONNAME:","",PTV2))},error=function(cond){message("List value is empty")})
                PTV3<-tryCatch({paste("PRECURSORTYPE:",P1TV2)},error=function(cond){message("List value is empty")})
                ##PTV3<-tryCatch({paste("PRECURSORTYPE:",P1TV2)},error=function(cond){message("List value is empty")})
                ########################################
                print("enter the test area..4")
                ########################################
                print(P1TV2)
                print(as.character(InMEDA[["Adduct"]]))
                ################################## adding this new
		if((!sjmisc::is_empty(P1TV2) || !sjmisc::is_empty(as.character(InMEDA[["Adduct"]]))) & (F2NPA != 0)){
		######################################	
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
		  ###########################
		  NAT<-as.character(InMEDA[["Adduct"]])
                  NAT1<-paste("PRECURSORTYPE:",NAT,sep=" ")
                  out<-c(out,NAT1)
                  ###########################
                  FIN1<-InMEDA[["Ionization mode"]]
                  F1IN1<-as.character(FIN1)
                  F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                  out<-c(out,F2IN1)
                  #################################################
                  ########################### adding this new #####
		  ########################################
		  InKeyVal<-tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
                  ##############################
                  IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
                  ##################################
                  if(!sjmisc::is_empty(IKCRV)){
                    ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
                    F1ONT<-paste("Ontology:",ONTV,sep=" ")
                    out<-c(out,F1ONT)
                  }else{
                    ##F1ONT<-paste("Ontology:","",sep=" ")
		    F1ONT<-tryCatch({MaKE.ONT.REC(InMEDA)[1]},error=function(cond){message("Classifier could not fecth the information")})	  
                    out<-c(out,F1ONT)
                  }
                  #################################
		  #################################
		  INV3<-tryCatch({rinchi::get.inchi(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
                  FINCH<-paste("INCHI:",INV3,sep=" ")
		  ##########################
                  out<-c(out,tryCatch({SelectInchiorInchiKey(InMEDA,IK,InchiV)[1]},error=function(cond){return(NA)}))
                  out<-c(out,tryCatch({SelectInchiorInchiKey(InMEDA,IK,InchiV)[2]},error=function(cond){return(NA)}))
		  ####out<-c(out,FINCH)
                  ######################################
		  FSIM<-paste("SMILES:",tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}),sep=" ")
                  out<-c(out,FSIM)
                  ##################################################
		  FFOR=ifelse(!sjmisc::is_empty(gETSmiles(InMEDA)),ifelse(!sjmisc::is_empty(tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)})),tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})),ifelse(!sjmisc::is_empty(InMEDA[["InChI"]]),tryCatch({getCactus(InMEDA[["InChI"]], "formula")},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})))
		  #################################################
		  ############################################
                  FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                  out<-c(out,FFOR1)
                  ############################################
                  ############################################
                  FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                  FINS1<-FNA[FINS]
		  FINS2<-ifelse(!sjmisc::is_empty(stringr::str_trim(gsub("INTENSITY:","",FINS1))),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
                  ####FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
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
                  ####################################
		  ####################################
		  PT1<-c("Num Peaks:","NumPeaks:")
                  FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
                  F1NPA<-FNA[FNPA]
                  F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
                  Np<-F2NPA
                  Fpea<-FNA[(FNPA+1):Find]
                  Fpea1<-Fpea[!is.na(Fpea)]
                  Fpea2<-PeakVali(Fpea1)
                  ######################################
		  ######################################
		  if(F2NPA != 0)
		  ###########################
                  {
                    MV=ADV4
		    ################################
		    tes1<-unlist(strsplit(Fpea2, "\t|\t\t|;| +"))
                    tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
                    tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
                    tes4<-which(tes2 > (3+MV))


		    if((Np >= 60) & (length(NTES) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))
  		     {
			     NTES1<-c(NTES4,tes4)
                        


			if((Np >= 60) & (length(NTES1) > 1))     
			{
				tes5<-tes2[-NTES1]
				tes6<-tes3[-NTES1]

				mzV<-tes5
				intenV<-tes6

				if(length(mzV) > 1)
  				{

    					print("entering inside the if loop")
    					nMZV=tryCatch({cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]},error=function(cond){return(0)})
    					nINV=tryCatch({cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]},error=function(cond){return(0)})
    					tes7<-paste(nMZV,nINV,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))

    					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  				}else{

    					tes5<-tes2
    					tes6<-tes3

    					nMZV=tryCatch({cenPeaks(tes2,tes3, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]},error=function(cond){return(0)})
    					nINV=tryCatch({cenPeaks(tes2,tes3, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]},error=function(cond){return(0)})

    					if(length(nMZV) > 1 & length(nINV) > 1)
    					{

      						print("checking the value of nMZV and nINV")
      						tes7<-paste(nMZV,nINV,sep="\t")

      						F1NPA<-paste0("Num Peaks: ",length(tes7))
      						out<-c(out,F1NPA)
      						out<-c(out,tes7)

    					}else{

    						tes7<-paste(tes5,tes6,sep="\t")
    						F1NPA<-paste0("Num Peaks: ",length(tes7))
    						out<-c(out,F1NPA)
    						out<-c(out,tes7)

    					}## end of if else

  				}


			}else{

				if((length(tes2[-NTES]) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))
  				{
    					tes5<-tes2[-NTES]
    					tes6<-tes3[-NTES]
    					tes7<-paste(tes5,tes6,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))
    					out<-c(out,F1NPA)
   					out<-c(out,tes7)

  				}else{

    					tes7<-paste(tes2,tes3,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))
    					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  				}## end of if else loop

			}## end of else


  		    }else{

		    NTES<-c(NTES4,tes4)
		    
		    if((length(NTES) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))
                    {

			    if(length(tes2[-NTES]) > 1)
    			    {

      				tes5<-tes2[-NTES]
      				tes6<-tes3[-NTES]
      				tes7<-paste(tes5,tes6,sep="\t")
      				### Adding this new here
      				F1NPA<-paste0("Num Peaks: ",length(tes7))
      				out<-c(out,F1NPA)
      				#############################
      				out<-c(out,tes7)
      				##############################
    			}else{

      				tes5<-tes2
      				tes6<-tes3
      				tes7<-paste(tes5,tes6,sep="\t")
      				F1NPA<-paste0("Num Peaks: ",length(tes7))
      				out<-c(out,F1NPA)
      				out<-c(out,tes7)

    			}### if,else loop
                   ########################################			    
                    }else{
		      #####################################	    
		      
		      ##if((Np >= 80))
	              if((Np >= 60))		    
    		       {
      			  mzV<-tes2
      			  intenV<-tes3
      
      			nMZV=cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]
      			nINV=cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]
      
      			tes7<-paste(nMZV,nINV,sep="\t")
      			F1NPA<-paste0("Num Peaks: ",length(tes7))
      
      			out<-c(out,F1NPA)
      			out<-c(out,tes7)
      
      
    			}else{	    


		      #############################
		      ###F1NPA<-paste0("Num Peaks: ",tryCatch({length(Fpea2)},error=function(cond){message("Fpea is empty")}))
		      F1NPA<-paste0("Num Peaks: ",tryCatch({length(Fpea2)},error=function(cond){message("Fpea is empty")}))		
		      out<-c(out,F1NPA)
		      out<-c(out,Fpea2)
                      ##out<-c(out,Fpea)
		      ################################
		       }## end of the else and Np >= 80
	           ###############################	      
                    }
		 #######################################
		     }## end of if loop and Np >= 120
		 #######################################   
                  } ### this is the else MV= ADV4...closing
                  ##########################################
                  ##########################################
                }else{
                  #####################################     
                  print("enter the else part...4")
                  #####################################
                  #####################################
		  if((as.character(InMEDA[["Adduct"]]) == "[M]+" || as.character(InMEDA[["Adduct"]]) == "[M]-") & (F2NPA != 0)){
                    #####################################
                    #####################################
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
		    #################################
		    NAT<-as.character(InMEDA[["Adduct"]])
                    NAT1<-paste("PRECURSORTYPE:",NAT,sep=" ")
                    out<-c(out,NAT1)
                    #################################
                    FIN1<-InMEDA[["Ionization mode"]]
                    F1IN1<-as.character(FIN1)
                    F2IN1<-paste("IONMODE:",F1IN1,sep=" ")
                    out<-c(out,F2IN1)
                    ########################################
		    ########################################
		    InKeyVal<-tryCatch({rinchi::get.inchi.key(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
                    ###########################################
                    IKCRV<-tryCatch({classyfireR::get_classification(InKeyVal)},warning=function(cond){message("Classifier could not fecth the information")})
                    ###########################################
                    if(!sjmisc::is_empty(IKCRV)){
                      ONTV<-do.call(paste, c(as.list(tryCatch({IKCRV@classification$Classification},warning=function(cond){message("Classifier could not fecth the information")})), sep = ";"))
                      F1ONT<-paste("Ontology:",ONTV,sep=" ")
                      out<-c(out,F1ONT)
                    }else{
		      F1ONT<-tryCatch({MaKE.ONT.REC(InMEDA)[1]},error=function(cond){message("Classifier could not fecth the information")})	    
                      ##F1ONT<-paste("Ontology:","",sep=" ")
                      out<-c(out,F1ONT)
                    }
                    ##########################
                    FINK<-paste("INCHIKEY:",tryCatch({InKeyVal},error=function(cond){message("Inchikey value is empty")}) ,sep=" ")
                    ###out<-c(out,FINK)
		    out<-c(out,tryCatch({SelectInchiorInchiKey(InMEDA,IK,InchiV)[1]},error=function(cond){return(NA)}))
                    #################################
		    ################################
		    INV3<-tryCatch({rinchi::get.inchi(tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}))},error=function(cond){return(NA)})
                    FINCH<-paste("INCHI:",INV3,sep=" ")
                    ###out<-c(out,FINCH)
		    out<-c(out,tryCatch({SelectInchiorInchiKey(InMEDA,IK,InchiV)[2]},error=function(cond){return(NA)}))
                    ########################
                    ########################
		    FSIM<-paste("SMILES:",tryCatch({gETSmiles(InMEDA)},error=function(cond){return(NA)}),sep=" ")
                    out<-c(out,FSIM)
                    #################################
		    FFOR=ifelse(!sjmisc::is_empty(gETSmiles(InMEDA)),ifelse(!sjmisc::is_empty(tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)})),tryCatch({RChemMass::MolFormFromSmiles.rcdk(gETSmiles(InMEDA))},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})),ifelse(!sjmisc::is_empty(InMEDA[["InChI"]]),tryCatch({getCactus(InMEDA[["InChI"]], "formula")},error=function(cond){return(NA)}),tryCatch({InMEDA[["Formula"]]},error=function(cond){return(NA)})))
		   
		    ###############################################
                    FFOR1<-paste("FORMULA:",FFOR,sep=" ")
                    out<-c(out,FFOR1)
                    ###################################################		
                    ###################################################
                    FINS<-which(stringi::stri_detect_fixed(FNA,"INTENSITY:"))
                    FINS1<-FNA[FINS]
		    FINS2<-ifelse(!sjmisc::is_empty(stringr::str_trim(gsub("INTENSITY:","",FINS1))),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
                    ###FINS2<-ifelse(!sjmisc::is_empty(FINS1),FINS1,paste("INTENSITY:",sample(100:200,1),sep=""))
                    out<-c(out,FINS2)
                    ##################################
                    FAUT<-as.character(InMEDA[["Authors"]])
                    FAUT1<-paste("AUTHORS:",FAUT,sep=" ")
                    out<-c(out,FAUT1)
                    ###############################
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
                    FCOM<-paste("COMMENT:",as.character(InMEDA[["Confidence"]]),sep=" ")
                    out<-c(out,FCOM)
                    ##################################
		    ##################################
		    PT1<-c("Num Peaks:","NumPeaks:")
                    FNPA<-tryCatch({grep(paste(PT1,collapse="|"), FNA)},error=function(cond){message("Num value is empty")})
                    F1NPA<-FNA[FNPA]
                    F2NPA<-as.numeric(stringr::str_trim(gsub("Num Peaks:|NumPeaks:","",F1NPA)))
                    Np<-F2NPA
                    Fpea<-FNA[(FNPA+1):Find]
                    Fpea1<-Fpea[!is.na(Fpea)]
                    Fpea2<-PeakVali(Fpea1)
                    #########################
		    ###############################
		    if(F2NPA != 0)
		    #############################
                    {
                    ###########}else{
                      
                      MV=AAMS1
		      tes1<-unlist(strsplit(Fpea2, "\t|\t\t|;| +"))
                      tes2<-as.numeric(tes1[schoolmath::is.odd(seq_along(tes1))])
                      tes3<-as.numeric(tes1[schoolmath::is.even(seq_along(tes1))])
                      tes4<-which(tes2 > (3+MV))
		      NTES4<-which(tes3 <= 100)
		      NTES<-c(NTES4,tes4)


		      if((Np >= 60) & (length(NTES) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))
  			{
				NTES1<-c(NTES4,tes4)


			if((Np >= 60) & (length(NTES1) > 1))	
			{
				tes5<-tes2[-NTES1]
				tes6<-tes3[-NTES1]

				mzV<-tes5
				intenV<-tes6

				 if(length(mzV) > 1)
  				 {

    					nMZV=tryCatch({cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]},error=function(cond){return(0)})
    					nINV=tryCatch({cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]},error=function(cond){return(0)})
    					tes7<-paste(nMZV,nINV,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))

    					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  				}else{

    					tes5<-tes2
    					tes6<-tes3

    					nMZV=tryCatch({cenPeaks(tes2,tes3, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]},error=function(cond){return(0)})
    					nINV=tryCatch({cenPeaks(tes2,tes3, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]},error=function(cond){return(0)})

    					if(length(nMZV) > 1 & length(nINV) > 1)
    					{

      						tes7<-paste(nMZV,nINV,sep="\t")

      						F1NPA<-paste0("Num Peaks: ",length(tes7))
      						out<-c(out,F1NPA)
      						out<-c(out,tes7)

    					}else{

    						tes7<-paste(tes5,tes6,sep="\t")
    						F1NPA<-paste0("Num Peaks: ",length(tes7))
    						out<-c(out,F1NPA)
    						out<-c(out,tes7)
    					}## end of if else

  					}


			}else{


				if((length(tes2[-NTES]) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))
  				{
    					tes5<-tes2[-NTES]
    					tes6<-tes3[-NTES]
    					tes7<-paste(tes5,tes6,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))
    					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  				}else{

    					tes7<-paste(tes2,tes3,sep="\t")
    					F1NPA<-paste0("Num Peaks: ",length(tes7))
    					out<-c(out,F1NPA)
    					out<-c(out,tes7)

  				}## end of if else loop
                        
			}### end of the else

                       }else{


		      if((length(NTES) > 1) & (length(Fpea2) != length(NTES)) & (length(Fpea2) > length(NTES)))	       
                      {



			if(length(tes2[-NTES]) > 1)
    			{

      				tes5<-tes2[-NTES]
      				tes6<-tes3[-NTES]
      				tes7<-paste(tes5,tes6,sep="\t")
      				### Adding this new here
      				F1NPA<-paste0("Num Peaks: ",length(tes7))
      				out<-c(out,F1NPA)
      				#############################
      				out<-c(out,tes7)
      				##############################
    			}else{

      				tes5<-tes2
      				tes6<-tes3
      				tes7<-paste(tes5,tes6,sep="\t")
      				F1NPA<-paste0("Num Peaks: ",length(tes7))
      				out<-c(out,F1NPA)
      				out<-c(out,tes7)

    			}### if,else loop      




                      }else{
			#######################

			
			if((Np >= 60))      
    			{
      				mzV<-tes2
      				intenV<-tes3
      
      				nMZV=cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[1]]
      				nINV=cenPeaks(mzV,intenV, MZ_TOLERANCE=DEFAULT_MZ_TOLERANCE)[[2]]
      
      				tes7<-paste(nMZV,nINV,sep="\t")
      				F1NPA<-paste0("Num Peaks: ",length(tes7))
      
      				out<-c(out,F1NPA)
      				out<-c(out,tes7)
      
      
    			}else{      


			#################################
	                F1NPA<-paste0("Num Peaks: ",tryCatch({length(Fpea2)},error=function(cond){message("Fpea is empty")}))			
			###F1NPA<-paste0("Num Peaks: ",tryCatch({length(Fpea)},error=function(cond){message("Fpea is empty")}))
                        out<-c(out,F1NPA)
			out<-c(out,Fpea2)
                        ###out<-c(out,Fpea)      
			#################################
			 }##end of else amd Np >= 80
		    ##################################### 	       
                      }
		  ######################################
			}## end of if loop and Np >= 120
		  ######################################    
                    }
                 ###############################
		 ##############################  
                  }
                ##############################
                ##############################
                }
              ###############################
	      ###############################	
              } else {
                #print("entering the line 764")
                PASS<-RRV
              }  ## main else part
              ###testing if this works
              return(out)
            ###############################
            } #ITmass
          } #ITRTL
          ###################
        } ## adduct match and exacct mass present
      }## FM2 is not empty 
    ###############################################################
    ###############################################################
    }###this is end of else ...where it is checking from NAME
  } ## end of else
###################################################################
####################################################################
}###NFFilter

###########################################################################################
###########################################################################################


###########################################################################################
###########################################################################################
for(i in 1:length(LmeCmu1))
####################################################################
####################################################################	
#################################	
{
  ########################
  print("entering the main function")
  ########################
  
  Val=LmeCmu1[i]
  Val1=LmeCmu1[i+1]
  nVal=Val+1
  #########################
  ################################
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
    BNFiNAEX26<-gsub("NRG[\\]Set01 NRG[\\]","",BNFiNAEX25)
    ########### changing this#############
    FINAMSP<-paste(BNFiNAEX26,"msp",sep=".")
    ########### changing this############# 
    ###FINAMSP<-paste(BNFiNAEX25,"msp",sep=".")
    ####################################
    ##print(FINAMSP)
    ##print(Fi6)
    ####################################
    OUNA<-tools::file_path_sans_ext(FiNA)
    OUNA1<-paste(OUNA,"passed","msp",sep = ".")
    OUNA2<-paste(Fi9,OUNA1,sep="")
    ##########mkdir mz.25ppm.20RT
    OUNA3<-paste(Fi11,OUNA1,sep="")
    OUNA4=OUNA3
    ###############FiNA1 ... INput ..msp file #############################
    #######################################################################
    NFiNA<-tryCatch({list.files(path =Fi6 , recursive = TRUE, full.names = TRUE)},warning=function(cond){message("some mistake happened in file search files")})
    NFiNA1<-tryCatch({match(FINAMSP,basename(NFiNA))},warning=function(cond){message("some mistake happened in file search files")})
    NFiNA2<-tryCatch({NFiNA[NFiNA1]},warning=function(cond){message("some mistake happened in file search files")})
    #####OUNA2 ... Output file ##########################################
    FiNA1<-NFiNA2
    ######################################################################
    ######################################################################
    ###print(FiNA1)
    ######################################################################
    ######################################################################
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
      ########################################
      if(file.exists(FiNA1[1]) & !sjmisc::is_empty(fmass) & !sjmisc::is_empty(FRTL1))	      
      ############################################	      
      {
        ########################################
        print("entering the line ...1732")
        ######################################
        if(dim(NRXF3)[1]>0){
	  for (i in 1:dim(NRXF3)[1]){
         ####################################		  
         #####################################
            
            RRV<-NRXF3[i,]
            ############################################
            ############################################
            PIN<-ifelse(is.na(RRV[["InChI"]]), 1, 0)
            PSM<-ifelse(is.na(RRV[["SMILES"]]), 1, 0)
            PPC<-ifelse(is.na(RRV[["PubChem CID"]]), 1, 0)
            ##############################################
	    ##############################################
            ###################################
            if(PIN == 1 & PSM == 1 & PPC == 1){
              ##################################
              print("enter the if loop as Inchi and smiles and pubchemid all are empty")
              #################### adding this new
              if(sjmisc::is_empty(RRV)){
                te<-"PASS"
              }else{
                #####################################
                #####################################
                RRV1<-NRXF3[i,]
                #######################################
                #########################################
                #########################################
                FINKE<-NFFilter(RRV1,AIN,lst2,fmass,FRTL1)
                if(!sjmisc::is_empty(FINKE) & length(FINKE) > 1)
                {
                  print("enter the line ...1680")
		  ifelse(!dir.exists(file.path(dirname(OUNA4))),dir.create(file.path(dirname(OUNA4))), FALSE)
                  cat(sapply(FINKE, toString), file=OUNA4, sep="\n",append=TRUE)
                  cat(sapply("", toString), file=OUNA4, sep="\n",append=TRUE)
                }else{
                  ################################################
                  print("enter the line ...1691")
                  print("enter the else part")
                  #################################################
                  DN<-dirname(dirname(File1))
		  DN1<-paste(DN,"Error-Report",sep="/")
		  ifelse(!dir.exists(DN1), dir.create(DN1), FALSE)
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
              ###########################################
              RRV2<-NRXF3[i,]
              ######################################################
              print("enter the line ...1749")
              #####################################################
              if(!sjmisc::is_empty(stringr::str_trim(as.character(RRV2[["InChI"]]))) & !startsWith(as.character(RRV2[["InChI"]]),'not available') & !startsWith(as.character(RRV2[["InChI"]]),'CAS:') & !startsWith(as.character(RRV2[["InChI"]]),'InChI=')){
                ##################################################
                print("enter the line ...1871")
                ##################################################
                if(tryCatch({webchem::is.inchikey(stringr::str_trim(as.character(RRV2[["InChI"]])))},error=function(cond){message("not able to pass is.inchiKey")})){
                  ###############################################
                  IV<-stringr::str_trim(as.character(RRV2[["InChI"]]))
                  ###########################################
                  print("entering the inchikey pass area")
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
		      ifelse(!dir.exists(file.path(dirname(OUNA4))),dir.create(file.path(dirname(OUNA4))), FALSE)
                      cat(sapply(outn, toString), file=OUNA4, sep="\n",append=TRUE)
                      cat(sapply("", toString), file=OUNA4, sep="\n",append=TRUE)
                      
                    }else{
                      print("enter the line ...1879")
                      FINKE<-NFFilter(RRV2,AIN,lst2,fmass,FRTL1)

		      print("enter the value of FINKE")
		      print(length(FINKE))
                      ### adding this code new 
                      len3<-length(FINKE)
                      if(len3 > 1){
                        print("enter the line ...1880")
		        ifelse(!dir.exists(file.path(dirname(OUNA4))),dir.create(file.path(dirname(OUNA4))), FALSE)
                        cat(sapply(FINKE, toString), file=OUNA4, sep="\n",append=TRUE)
                        cat(sapply("", toString), file=OUNA4, sep="\n",append=TRUE)
                      }else{
                        print("enter the line ...1891")
                        print("enter the else part...Error report")
                        ################################
                        DN<-dirname(dirname(File1))
                        DN1<-paste(DN,"Error-Report",sep="/")
			ifelse(!dir.exists(DN1), dir.create(DN1), FALSE)
                        FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
                        write.table(unname(as.data.frame(RRV2)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
                        #################################
                      }
                      ########################################################################## 
                    } # end of the else loop...that means it did not pass Inchikey
                    ##########################################################################                                
                  }### end of IV
                  #############################  
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
		    ifelse(!dir.exists(file.path(dirname(OUNA4))),dir.create(file.path(dirname(OUNA4))), FALSE)
                    cat(sapply(FINKE, toString), file=OUNA4, sep="\n",append=TRUE)
                    cat(sapply("", toString), file=OUNA4, sep="\n",append=TRUE)
                  }else{
                    print("enter the line ...1785")
                    print("enter the else part")
                    ############################################################
                    DN<-dirname(dirname(File1))
                    DN1<-paste(DN,"Error-Report",sep="/")
		    ifelse(!dir.exists(DN1), dir.create(DN1), FALSE)
                    FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
                    write.table(unname(as.data.frame(RRV2)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
                    ###########################################################
                  }### end of inner else
                  
                  ################################################################
                }### end of else ...where it did not pass Inchikey
                #####################################################################################################	
                ######################################################################################################################################
              }else if(!sjmisc::is_empty(as.character(RRV2[["InChI"]])) & startsWith(as.character(RRV2[["InChI"]]),'InChI=')){
                #################################################      
                print("enter the line ...861")
                print("entering the inchi...... pass area")
                ##IV1<-stringr::str_trim(as.character(RRV[["InChI"]]))
                #####################################################
                IV1<-stringr::str_trim(as.character(RRV2[["InChI"]]))
                #############################################################
                FSMV<-tryCatch({rinchi::parse.inchi(IV1)},error=function(cond){message("Inchi name must be empty or rinchi not abe to fetch")})
                FSMV1<-tryCatch({rcdk::get.smiles(FSMV[[1]])},error=function(cond){message("rcdk get smiles conversion is empty")})
                FSMV2<-tryCatch({rinchi::get.inchi.key(FSMV1)},error=function(cond){message("webchecm could not fetch the info")})
                #############################################################
		print("the value of FSMV2")
		print(FSMV2)
                #############################################################
                if(!sjmisc::is_empty(FSMV2)){
                  ##################################################		    
                  print("enter the line ...868")
                  ######################################
                  outn1<-Ikfilter(FSMV2,lst2,RRV2,AIN,fmass,FRTL1)
                  len2<-length(outn1)
                  ####################################################
                  if(len2 > 1)
                  {
                    print("enter the line ...878")
		    print("enter the directory path")
                    print(OUNA4)
		    ifelse(!dir.exists(file.path(dirname(OUNA4))),dir.create(file.path(dirname(OUNA4))), FALSE)
                    cat(sapply(outn1, toString), file=OUNA4, sep="\n",append=TRUE)
                    cat(sapply("", toString), file=OUNA4, sep="\n",append=TRUE)
                    
                    ##} # end of if loop
                  }else{
                    print("enter the line ...880")	
                    FINKE1<-NFFilter(RRV2,AIN,lst2,fmass,FRTL1)
                    ### adding this code new 
                    len3<-length(FINKE1)
                    if(len3 > 1){
                      print("enter the line ...888")
		      ifelse(!dir.exists(file.path(dirname(OUNA4))),dir.create(file.path(dirname(OUNA4))), FALSE)
                      cat(sapply(FINKE1, toString), file=OUNA4, sep="\n",append=TRUE)
                      cat(sapply("", toString), file=OUNA4, sep="\n",append=TRUE)
                    }else{
                      print("enter the line ...889")
                      print("enter the else part")
                      ###################################
                      DN<-dirname(dirname(File1))
                      DN1<-paste(DN,"Error-Report",sep="/")
		      ifelse(!dir.exists(DN1), dir.create(DN1), FALSE)
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
		    ifelse(!dir.exists(file.path(dirname(OUNA4))),dir.create(file.path(dirname(OUNA4))), FALSE)
                    cat(sapply(FINKE1, toString), file=OUNA4, sep="\n",append=TRUE)
                    cat(sapply("", toString), file=OUNA4, sep="\n",append=TRUE)
                  }else{
                    print("enter the line ...875")
                    print("enter the else part")
                    #####################################
                    DN<-dirname(dirname(File1))
                    DN1<-paste(DN,"Error-Report",sep="/")
		    ifelse(!dir.exists(DN1), dir.create(DN1), FALSE)
                    FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
                    write.table(unname(as.data.frame(RRV2)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
                    #######################################
                  }###end of inner else
                  ################################################
                  
                }#### end of else 
                
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
                ###################################################
                FIV=tryCatch({webchem::aw_query(CV2, from = 'cas')},error=function(cond){message("some mistake happened in cas finding")})
                FIV1=tryCatch({FIV[[1]]},error=function(cond){message("some mistake happened in file search files")})
                FIV2=tryCatch({FIV1$inchikey},error=function(cond){message("some mistake happened in file search files")})
                ###################################
                #################################
                if(!sjmisc::is_empty(FIV2)){
                  ##if(!is.na(FIV1)){
                  outn3<-Ikfilter(FIV2,lst2,RRV2,AIN,fmass,FRTL1)
                  len4<-length(outn3)
                  if(len4 > 1)
                  {
                    print("enter the line ...912")
		    ifelse(!dir.exists(file.path(dirname(OUNA4))),dir.create(file.path(dirname(OUNA4))), FALSE)
                    cat(sapply(outn3, toString), file=OUNA4, sep="\n",append=TRUE)
                    cat(sapply("", toString), file=OUNA4, sep="\n",append=TRUE)
                    
                    ##} # end of len4
                  }else{
                    FSMIL<-NFFilter(RRV2,AIN,lst2,fmass,FRTL1)
                    len3<-length(FSMIL)
                    if(len3 > 1){
                      print("enter the line ...921")
		      ifelse(!dir.exists(file.path(dirname(OUNA4))),dir.create(file.path(dirname(OUNA4))), FALSE)
                      cat(sapply(FSMIL, toString), file=OUNA4, sep="\n",append=TRUE)
                      cat(sapply("", toString), file=OUNA4, sep="\n",append=TRUE)
                    }else{
                      print("enter the line ...923")
                      print("enter the else part")
                      #####################################
                      DN<-dirname(dirname(File1))
                      DN1<-paste(DN,"Error-Report",sep="/")
		      ifelse(!dir.exists(DN1), dir.create(DN1), FALSE)
                      FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
                      write.table(unname(as.data.frame(RRV2)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
                      ########################################
                    }### end of inner else
                    ##########################
                  } # end of else loop
                  
                }else{
                  ######################################################################
		  ######################################################################	
                  print("CAS to inchikey conversion failed")
                  print("added the new code")
		  ###FuFtoRe1(InMEDA)
		  ###FuFtoRe1(RRV2)
		  ##################################################
                  IV<-ifelse(!sjmisc::is_empty(tryCatch({PuCAStoOI(CV2)[2]},error=function(cond){message("CAS is empty")})),tryCatch({PuCAStoOI(CV2)[2]},error=function(cond){message("CAS value is empty")}),FuFtoRe1(RRV2))
                  ###################################################################
                  ###################################################################
                  if(!sjmisc::is_empty(IV)){  
                    outn<-Ikfilter(IV,lst2,RRV2,AIN,fmass,FRTL1)
                    len1<-length(outn)
                    ##############################################
                    if(len1 > 1)
                    {
                      print("enter the line ...1878")
		      ifelse(!dir.exists(file.path(dirname(OUNA4))),dir.create(file.path(dirname(OUNA4))), FALSE)
                      cat(sapply(outn, toString), file=OUNA4, sep="\n",append=TRUE)
                      cat(sapply("", toString), file=OUNA4, sep="\n",append=TRUE)
                      
                    }else{
                      print("enter the line ...1879")
                      FINKE<-NFFilter(RRV2,AIN,lst2,fmass,FRTL1)
                      ### adding this code new 
                      len3<-length(FINKE)
                      if(len3 > 1){
                        print("enter the line ...1880")
		        ifelse(!dir.exists(file.path(dirname(OUNA4))),dir.create(file.path(dirname(OUNA4))), FALSE)
                        cat(sapply(FINKE, toString), file=OUNA4, sep="\n",append=TRUE)
                        cat(sapply("", toString), file=OUNA4, sep="\n",append=TRUE)
                      }else{
                        print("enter the line ...1891")
                        print("enter the else part...Error report")
                        ################################
                        DN<-dirname(dirname(File1))
                        DN1<-paste(DN,"Error-Report",sep="/")
			ifelse(!dir.exists(DN1), dir.create(DN1), FALSE)
                        FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
                        write.table(unname(as.data.frame(RRV2)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
                        #################################
                      }###end of inner else
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
		      ifelse(!dir.exists(file.path(dirname(OUNA4))),dir.create(file.path(dirname(OUNA4))), FALSE)
                      cat(sapply(FINKE, toString), file=OUNA4, sep="\n",append=TRUE)
                      cat(sapply("", toString), file=OUNA4, sep="\n",append=TRUE)
                    }else{
                      print("enter the line ...1785")
                      print("enter the else part")
                      ############################################################
                      DN<-dirname(dirname(File1))
                      DN1<-paste(DN,"Error-Report",sep="/")
		      ifelse(!dir.exists(DN1), dir.create(DN1), FALSE)
                      FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
                      write.table(unname(as.data.frame(RRV2)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
                      ###########################################################
                    }
                    ############################################################## 
                  } ## end of else
                  #######################################################
                  ######################################################
                  #############################################
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
		    ifelse(!dir.exists(file.path(dirname(OUNA4))),dir.create(file.path(dirname(OUNA4))), FALSE)
                    cat(sapply(outn4, toString), file=OUNA4, sep="\n",append=TRUE)
                    cat(sapply("", toString), file=OUNA4, sep="\n",append=TRUE)
                    
                  }else{
                    print("enter the line ...959")
                    FSMIL<-NFFilter(RRV2,AIN,lst2,fmass,FRTL1)
                    len3<-length(FSMIL)
                    if(len3 > 1){
                      print("enter the line ...963")
		      ifelse(!dir.exists(file.path(dirname(OUNA4))),dir.create(file.path(dirname(OUNA4))), FALSE)
                      cat(sapply(FSMIL, toString), file=OUNA4, sep="\n",append=TRUE)
                      cat(sapply("", toString), file=OUNA4, sep="\n",append=TRUE)
                    }else{
                      print("enter the line ...967")
                      print("enter the else part")
                      ##################################
                      DN<-dirname(dirname(File1))
                      DN1<-paste(DN,"Error-Report",sep="/")
		      ifelse(!dir.exists(DN1), dir.create(DN1), FALSE)
                      FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
                      write.table(unname(as.data.frame(RRV2)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
                      ##########################################
                    }##end of inner else
                    ################################################################## 
                  } ## end of else loop
                  
                }else{
                  ###############################################		 
                  #print("enter the ...not passing..F1SM1  ...failed")
                  print("smiles to inchikey conversion is empty so checking other code")
                  NF1SM1<-tryCatch({PuSMItoINKandInchi(F1SM)},error=function(cond){message("smiles to other conversion failed")})
                  NF1SM2<-tryCatch({NF1SM2[2]},error=function(cond){message("smiles to other conversion failed")})
                  #########################################
                  #########################################
                  if(!sjmisc::is_empty(NF1SM2)){
                    outn4<-Ikfilter(NF1SM2,lst2,RRV2,AIN,fmass,FRTL1)
                    loutn5<-length(outn4)
                    ###################################################
                    if(loutn5 > 1)
                    {
                      print("enter the line ...954")
		      ifelse(!dir.exists(file.path(dirname(OUNA4))),dir.create(file.path(dirname(OUNA4))), FALSE)
                      cat(sapply(outn4, toString), file=OUNA4, sep="\n",append=TRUE)
                      cat(sapply("", toString), file=OUNA4, sep="\n",append=TRUE)
                      
                    }else{
                      print("enter the line ...959")
                      FSMIL<-NFFilter(RRV2,AIN,lst2,fmass,FRTL1)
                      len3<-length(FSMIL)
                      if(len3 > 1){
                        print("enter the line ...963")
		        ifelse(!dir.exists(file.path(dirname(OUNA4))),dir.create(file.path(dirname(OUNA4))), FALSE)
                        cat(sapply(FSMIL, toString), file=OUNA4, sep="\n",append=TRUE)
                        cat(sapply("", toString), file=OUNA4, sep="\n",append=TRUE)
                      }else{
                        print("enter the line ...967")
                        print("enter the else part")
                        ##################################
                        DN<-dirname(dirname(File1))
                        DN1<-paste(DN,"Error-Report",sep="/")
			ifelse(!dir.exists(DN1), dir.create(DN1), FALSE)
                        FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
                        write.table(unname(as.data.frame(RRV2)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
                        ##########################################
                      }##end of inner else
                      ################################################################## 
                    } ## end of else loop
                    
                  }else{
                    ## adding this new 
                    FSMIL<-NFFilter(RRV2,AIN,lst2,fmass,FRTL1)
                    len3<-length(FSMIL)
                    if(len3 > 1){
                      print("enter the line ...963")
		      ifelse(!dir.exists(file.path(dirname(OUNA4))),dir.create(file.path(dirname(OUNA4))), FALSE)
                      cat(sapply(FSMIL, toString), file=OUNA4, sep="\n",append=TRUE)
                      cat(sapply("", toString), file=OUNA4, sep="\n",append=TRUE)
                    }else{
                      
                      ##NF1SM2<-tryCatch({rinchi::get.inchi.key()},error=function(cond){message("smiles to other conversion failed")})
                      print("enter the line ...967")
                      print("enter the else part")
                      ##################################
                      DN<-dirname(dirname(File1))
                      DN1<-paste(DN,"Error-Report",sep="/")
		      ifelse(!dir.exists(DN1), dir.create(DN1), FALSE)
                      FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
                      write.table(unname(as.data.frame(RRV2)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
                      ##########################################
                    }##end of inner else
                  }
              #####################################################
              #####################################################          
                } ## end of else...F1SM...else
                ###################################################### 
                ####}### end of else ..smiles to inchikey conversion failed
                
                #############################################################################################################################
              }else if(!sjmisc::is_empty(stringr::str_trim(as.character(RRV2[["PubChem CID"]]))) & !startsWith(as.character(RRV2[["PubChem CID"]]),'not available')){
                print("enter the line ...939")
                print("entering the pubchem pass area")
                #############################################
                FPUCID<-stringr::str_trim(as.character(RRV2[["PubChem CID"]]))
                FPUCID1<-as.numeric(FPUCID)
                ##############################################
                print("entering the PubchemID test area")
                ##print(FPUCID1)
                ######################################################################
                ######################################################################
                FINSM<-tryCatch({webchem::pc_prop(FPUCID1)},error=function(cond){return(NA)})
                FIINK<-tryCatch({FINSM$InChIKey},error=function(cond){return(NA)})
                #####################################################################
                #####################################################################
                if(!sjmisc::is_empty(FIINK)){
                  outn6<-Ikfilter(FIINK,lst2,RRV2,AIN,fmass,FRTL1)
                  len7<-length(outn6)
                  if(len7 > 1)
                  {
                    print("enter the line ...994")
                    ##################
                    print("PASS")
                    print(FPUCID1)
                    ##################
		    ifelse(!dir.exists(file.path(dirname(OUNA4))),dir.create(file.path(dirname(OUNA4))), FALSE)
                    cat(sapply(outn6, toString), file=OUNA4, sep="\n",append=TRUE)
                    cat(sapply("", toString), file=OUNA4, sep="\n",append=TRUE)
                    
                  }else{
                    print("enter the line ...999")
                    FCIDR<-NFFilter(RRV2,AIN,lst2,fmass,FRTL1)
                    len3<-length(FCIDR)
                    if(len3 > 1){
                      print("enter the line ...1003")
                      ##################
                      ###print("PASS")
                      #####print(FPUCID1)
                      ########################
		      ifelse(!dir.exists(file.path(dirname(OUNA4))),dir.create(file.path(dirname(OUNA4))), FALSE)
                      cat(sapply(FCIDR, toString), file=OUNA4, sep="\n",append=TRUE)
                      cat(sapply("", toString), file=OUNA4, sep="\n",append=TRUE)
                    }else{
                      print("enter the line ...1005")
                      print("enter the else part")
                      ###############
                      print("FAIL")
                      print(FPUCID1)
		      print(FCIDR)
                      ######################################
                      DN<-dirname(dirname(File1))
                      DN1<-paste(DN,"Error-Report",sep="/")
		      ifelse(!dir.exists(DN1), dir.create(DN1), FALSE)
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
                  ##################################################
                  #############entering new code here ##############
                  ##IV<-ifelse(!sjmisc::is_empty(tryCatch({ConvPCIDtoOCN(FPUCID1)[1]},error=function(cond){message("Pubchem is empty")})),tryCatch({ConvPCIDtoOCN(FPUCID1)[1]},error=function(cond){message("CAS value is empty")}),"NA")
                   IV<-ifelse(!sjmisc::is_empty(tryCatch({ConvCIDtoOID1(FPUCID1)[2]},error=function(cond){message("Pubchem CID is empty")})),tryCatch({ConvCIDtoOID1(FPUCID1)[2]},error=function(cond){message("Pubchem CID is empty")}),"NA")

		  if(!sjmisc::is_empty(IV)){  
                    outn<-Ikfilter(IV,lst2,RRV2,AIN,fmass,FRTL1)
                    len1<-length(outn)
                    ##############################################
                    if(len1 > 1)
                    {
                      print("enter the line ...1878")
                      ##########################
                      print("PASS")
                      print(FPUCID1)
                      ##########################
		      ifelse(!dir.exists(file.path(dirname(OUNA4))),dir.create(file.path(dirname(OUNA4))), FALSE)
                      cat(sapply(outn, toString), file=OUNA4, sep="\n",append=TRUE)
                      cat(sapply("", toString), file=OUNA4, sep="\n",append=TRUE)
                      
                    }else{
                      print("enter the line ...1879")
                      FINKE<-NFFilter(RRV2,AIN,lst2,fmass,FRTL1)
                      ### adding this code new 
                      len3<-length(FINKE)
                      if(len3 > 1){
                        print("enter the line ...1880")
                        ########################
                        print("PASS")
                        print(FPUCID1)
                        ############################
			ifelse(!dir.exists(file.path(dirname(OUNA4))),dir.create(file.path(dirname(OUNA4))), FALSE)
                        cat(sapply(FINKE, toString), file=OUNA4, sep="\n",append=TRUE)
                        cat(sapply("", toString), file=OUNA4, sep="\n",append=TRUE)
                      }else{
                        print("enter the line ...1891")
                        print("enter the else part...Error report")
                        ################################
                        ####print("FAIL")
                        ####print(FPUCID1)
                        #####################################
                        DN<-dirname(dirname(File1))
                        DN1<-paste(DN,"Error-Report",sep="/")
			ifelse(!dir.exists(DN1), dir.create(DN1), FALSE)
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
		      ifelse(!dir.exists(file.path(dirname(OUNA4))),dir.create(file.path(dirname(OUNA4))), FALSE)
                      cat(sapply(FINKE, toString), file=OUNA4, sep="\n",append=TRUE)
                      cat(sapply("", toString), file=OUNA4, sep="\n",append=TRUE)
                    }else{
                      print("enter the line ...1785")
                      print("enter the else part")
                      ############################################################
                      DN<-dirname(dirname(File1))
                      DN1<-paste(DN,"Error-Report",sep="/")
		      ifelse(!dir.exists(DN1), dir.create(DN1), FALSE)
                      FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
                      write.table(unname(as.data.frame(RRV2)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)
                      ###############################################################
                    }
                    ############################################################## 
                  } ## end of else
                  ###################################################
                  ###################################################
                }# end of else loop....FIINK
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
      }else{
	      print("something went wrong in MaKlist...entering the final else part")
	      print(FiNA1[1])
	      print("Filename that does not exists is .")
	      print(FiNA)
	      DN<-dirname(dirname(File1))
	      DN1<-paste(DN,"Error-Report",sep="/")
	      ifelse(!dir.exists(DN1), dir.create(DN1), FALSE)
	      FILE<-paste(DN1,"Error.Report.25ppm.20.txt",sep="/")
	      write.table(unname(as.data.frame(NRXF3)), file = FILE, sep = "\t",quote=F,row.names = F, col.names = F,append=T)

      }
      
    } ## FiNA is closing one
    ###################################
  } ## is_empty(Val and Val1)
}
#############################################################################################################################################
####################################s#########################################################################################################



