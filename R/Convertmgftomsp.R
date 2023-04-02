get_ScanNumber<-function(LV)
{
  INT=grep("TITLE=",LV)
  INT1=tryCatch({LV[INT]},error=function(cond){message("Title value is empty")})
  INT2=tryCatch({sjmisc::str_end(INT1, "ScanNumber:|scan",regex = TRUE)},error=function(cond){message("Scan value is empty")})
  INT3=tryCatch({substring(INT1,(INT2+1),nchar(INT1))},error=function(cond){message("Substring scan value is empty")})
  INT4=tryCatch({ifelse(grepl("=",INT3),qdap::genXtract(INT3,"=","_")[[1]], substring(INT3,1,((sjmisc::str_end(INT3, "_"))-1)))},error=function(cond){message("Classifier could not fecth the information")})
  SCN=paste("SCANNUMBER:",ifelse(!sjmisc::is_empty(INT4),INT4,"NA"))
  return(SCN)
}

get_NumPeaks<-function(LV)
{
  if((length(grep("CHARGE=",LV)) == 0) & (grep("PEPMASS=",LV) != 0))
  {
     
     
     NPV <- ifelse(!sjmisc::is_empty(LV[as.numeric(grep("PEPMASS=",LV)):as.numeric(grep("END IONS",LV))]),paste("Num Peaks:",length(LV[(as.numeric(grep("PEPMASS=",LV))+1) :(as.numeric(grep("END IONS",LV))-1)])),paste("Num Peaks:",0))
     
     
     return(NPV)
    
  }else{
    
    NPV <- ifelse(!sjmisc::is_empty(LV[as.numeric(grep("CHARGE=",LV)):as.numeric(grep("END IONS",LV))]),paste("Num Peaks:",length(LV[(as.numeric(grep("CHARGE=",LV))+1) :(as.numeric(grep("END IONS",LV))-1)])),paste("Num Peaks:",0))
    
    return(NPV)
    
  }  
  
}


get_NumPeakVal<-function(LV)
{
  if(length((grep("CHARGE=",LV)) == 0) & (grep("PEPMASS=",LV) != 0))
  {
    
    NPV <- ifelse(!sjmisc::is_empty(LV[as.numeric(grep("PEPMASS=",LV)):as.numeric(grep("END IONS",LV))]),paste("Num Peaks:",length(LV[(as.numeric(grep("PEPMASS=",LV))+1) :(as.numeric(grep("END IONS",LV))-1)])),paste("Num Peaks:",0))
    
    if(paste("Num Peaks:",0) != NPV)
    {
      return(append(LV[(as.numeric(grep("PEPMASS=",LV))+1):(as.numeric(grep("END IONS",LV))-1)],""))
    }else{
      return("")
    }
    
  }else{
    
    NPV <- ifelse(!sjmisc::is_empty(LV[as.numeric(grep("CHARGE=",LV)):as.numeric(grep("END IONS",LV))]),paste("Num Peaks:",length(LV[(as.numeric(grep("CHARGE=",LV))+1) :(as.numeric(grep("END IONS",LV))-1)])),paste("Num Peaks:",0))
    if(paste("Num Peaks:",0) != NPV)
    {
      return(append(LV[(as.numeric(grep("CHARGE=",LV))+1):(as.numeric(grep("END IONS",LV))-1)],""))
    }else{
      return("")
    }
    
  } 
}
  

Make_msp_from_mgf<-function(ifile,ofile)
{
#############################################
lines <- readLines(ifile) #read msp as lines, and then manipulate it
lst <-split(lines, cumsum(lines=="BEGIN IONS"))
lst1 <-lapply(lst, function(x) if (x[1] == "") x[-1] else x)
LL<-sapply(lst, length)
IR<-which(unname(LL) == 0)
if(length(IR)>0){LL1<-LL[-(IR)]}else{LL1 <- LL}
###### MSP List#############################
lst2 <-lst1[names(LL1)]
lst2<-purrr::compact(lst2)
############################################

for(i in 1:length(lst2))
{
  
  if((length(which(stringr::str_starts(lst2[[i]], "BEGIN IONS"))) == 1) & (length(which(stringr::str_ends(lst2[[i]], "END IONS"))) == 1))
  {
   
    LV=lst2[[i]]
    
    gSN=get_ScanNumber(LV)
    
    ##1) Get PrecursorMz value
    premzI <- tryCatch({LV[grep('^PEPMASS=',LV, ignore.case=TRUE)]},error=function(cond){message("PEPMASS value is empty")})
    
    premzV <- tryCatch({gsub('PEPMASS=', '', premzI)},error=function(cond){message("PEPMASS value is not found in mgf")})
    
    FpremzV <-tryCatch({strsplit(premzV," +")[[1]][1]},error=function(cond){message("PrecursorMZ value is empty")})
    
    Finval <- tryCatch({strsplit(premzV," +")[[1]][2]},error=function(cond){message("PrecursorMZ value is empty")})
    ##2) Get the Retention time 
    RTVAI <- tryCatch({LV[grep('^RTINSECONDS=',LV, ignore.case=TRUE)]},error=function(cond){message("Retention time value is empty")})
    RTV<-tryCatch({gsub('RTINSECONDS=', '', RTVAI)},error=function(cond){message("Retention is not found in mgf")})
    FRTV<-tryCatch({strsplit(RTV," +")[[1]][1]},error=function(cond){message("RT value is empty")})
    
    ##3)  Get the NAME 
    NAMV=ifelse(!sjmisc::is_empty(LV[grep('NAME=',LV, ignore.case=TRUE)]),paste("NAME:",gsub('NAME=', '', LV[grep('NAME=',LV, ignore.case=TRUE)])),paste("NAME:","Unknown"))
    
    PTYVA <- ifelse(!sjmisc::is_empty(LV[!is.na(stringr::str_extract(LV, "\\[?[0-9]?M(-|\\+).*$"))]),LV[!is.na(stringr::str_extract(LV, "\\[?[0-9]?M(-|\\+).*$"))],paste("PRECURSORTYPE:","NA"))
    
    FINTV <- ifelse(!sjmisc::is_empty(Finval),paste("INTENSITY:",Finval),paste("INTENSITY:",runif(1, min=10, max=100)))
    
    NPV <- get_NumPeaks(LV)
    
    NPV1 <- get_NumPeakVal(LV)
    
   
    
    FMSP=c(NAMV,paste("PRECURSORMZ:",FpremzV),PTYVA,gSN,paste("RETENTIONTIME:",FRTV),paste("INTENSITY:",Finval),NPV,NPV1)
    
    
    cat(sapply(FMSP, toString), file=ofile,sep="\n",append=TRUE)
    
    
    
  }
    
}

}## Function end of mgf to msp
