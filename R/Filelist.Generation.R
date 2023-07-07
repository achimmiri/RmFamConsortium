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
MaKlist<-function(gFile)
{
  #####################################
  #####################################
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
  #####################################
  return(lst2)
}
####################################
#####################################
args <- commandArgs(TRUE)
DiN<-args[1]
FiN<-args[2]
EXN<-args[3]
########################################
########################################
## 1) The argument 1 is the RMassBank files path location 
## 2) The argument 2 is the combined msp file location in RMassBank folder
## 3) The excel file location 
#########################################
#########################################
alis=MaKlist(FiN)
##############################################################################
##############################################################################
###########################################################################
###########################################################################
File1<-DiN
File2<-EXN
#############################################################################
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
RL<-dim(RXF3)[1]
####################################
####################################
CN<-c("Files", "ID\n")
cat(CN,file =paste(dirname(File2),"Filelist.2.csv",sep="/"),sep=",")
#####################################
#####################################
###########################################
for(i in 1:RL)
{
  
  FN<-paste0(tools::file_path_sans_ext(RXF3[i,][["File"]]),".passed.msp")
  
  ID<-i


  NFN<-paste(DiN,FN,sep="/")
  df<-data.frame(NFN,ID)
 
  
  write.table(unname(df), file = paste(dirname(File2),"Filelist.2.csv",sep="/"), sep = ",",quote=TRUE,row.names = F, col.names = F,append=T)
  
  
}

##################################################
##################################################
