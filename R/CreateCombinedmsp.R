args <- commandArgs(TRUE)
File1<-args[1]
############################
LF<-Sys.glob(file.path(File1,"*.passed.msp"))
myfiles = lapply(LF, readLines)
DESN=paste(tail(unlist(strsplit(File1, "/")),3)[1],"mz.25ppm.20RT","combined.msp",sep=".")
DESN1=paste(File1,DESN,sep="/")
lapply(myfiles, write, file=DESN1,append=T)
