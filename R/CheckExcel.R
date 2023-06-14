args <- commandArgs(TRUE)
File1<-args[1]
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
#####################################################################
##print(dim(RXF3)[1])
RV=dim(RXF3)[1]
RV1=RV-2
print(RV1)
##num1="$RV"
##num2=2
##echo $RV
##RV1=$((num1-num2))"
##echo "$RV1"
