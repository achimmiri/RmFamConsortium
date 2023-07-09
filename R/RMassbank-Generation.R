############################################
args <- commandArgs(TRUE)
DiN<-args[1]
setwd(DiN)
#############################################
### 1) Argument 1 is the location of the RMassBank Folder
#############################################
#############################################
callFC<-function()
{
## Install the required package

if(require(remotes)){
install.packages("remotes",dependencies=TRUE)
library(remotes)
remotes::install_github("achimmiri/RMassBankChemOnt")
####################
###library("RMassBankChemOnt")
library("RMassBank")
########################
w <- newMsmsWorkspace()
files <- list.files(DiN, pattern='*.passed.msp', full.names=TRUE)
w@files <- files
loadList(paste(DiN,"Compoundlist.2.csv",sep="/"))
########################
print("what is printing the DiN")
print(DiN)
print(paste(DiN,"RMB_options.test2.ini",sep="/"))
########################
loadRmbSettings(paste(DiN,"RMB_options.test2.ini",sep="/"))
##w <- msmsWorkflow(w, readMethod='msp', filetable='./Filelist.2.csv', mode='pH', steps=1, archivename='msp_archive')
w <- msmsWorkflow(w, readMethod='msp', filetable=paste(DiN,"Filelist.2.csv",sep="/"), mode='pH', steps=1, archivename='msp_archive')
mb <- newMbWorkspace(w)
#########################
#########################
#### this is gathering the smiles from data bases and generate infolist.csv #####
mb <- mbWorkflow(mb)
#########################
#########################
mb <- resetInfolists(mb)
###mb <- loadInfolist(mb, "./infolist.csv")
mb <- loadInfolist(mb,paste(DiN,"infolist.csv",sep="/"))
##browser()
mb <- mbWorkflow(mb, filter=FALSE)

##################
}###end of if loop remotes
#####################

}
##############################################
##############################################
callFC()
