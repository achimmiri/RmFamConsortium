Validate_MetaData<-function(File1)
{
################################
RXF<-readxl::read_excel(File1, sheet = 1, col_names = TRUE,skip=1,.name_repair="minimal")
RXF[] <- lapply(RXF, function(x) type.convert(as.character(x)))
RXF1 <- which(is.na(as.character(RXF[["File"]])))
#################################################################
RXF3=list()
if(length(RXF1) >= 1)
{
  RXF2<- RXF[-RXF1,]
  RXF3=RXF2
}else
{
  RXF3=RXF
}
##################################################################
###1) Remove the exact duplictates
RXF4 <- RXF3 %>% dplyr::distinct()
##2) Remove empty rows 
RXF5<- RXF4[!apply(RXF4 == "", 1, all),]
###########################################
DN=dirname(dirname(File1))
DRF=paste0(DN,"/raw data/exported as raw msp")
FNMD=lapply(RXF3[["File"]], FUN = function(x) basename(tools::file_path_sans_ext(x)))
FNMD1=unlist(FNMD)
LFRD=list.files(DRF,pattern="*.msp")
FLFRD=list.files(DRF,pattern="*.msp",full.names = TRUE)
LFRD1=lapply(FNMD, FUN = function(x) basename(tools::file_path_sans_ext(x)))
LFRD2=unlist(LFRD1)
#############################
FC=setequal(FNMD1,LFRD2)
SSUI=RXF3[,c("InChI","SMILES","PubChem CID")]
################################
for(i in 1:length(SSUI))
{
  if(all(is.na(SSUI[i,])) & rlang::is_empty(SSUI[i,]))
  {
    RXF6=RXF5[-i,]
  }else{
    RXF6=RXF5
  }
  
}
##################
if(all(is.na(RXF3[["RT (min)"]])) & all(is.na(RXF3[["Adduct"]])) & is.na(unique(RXF3[["INSTRUMENT"]])) & is.na(unique(RXF3[["INSTRUMENT_TYPE"]])) & is.na(unique(RXF3[["IONIZATION"]])))
{
  RXF7=tibble::tibble(0, 1*0)
}else{
  RXF7=RXF6
}
########################
MaKlist<-function(gFile)
{
  #############################################
  lines <- readLines(gFile)
  ## This is the original one commenting
  lst <-split(lines, cumsum(lines==""))
  lst1 <-lapply(lst, function(x) if (x[1] == "") x[-1] else x)
  LL<-sapply(lst, length)
  IR<-which(unname(LL) == 0)
  if(length(IR)>0){LL1<-LL[-(IR)]}else{LL1 <- LL}
  ###### MSP List#############################
  lst2 <-lst1[names(LL1)]
  lst2<-purrr::compact(lst2)
  #############################
  PMZL<-unname(rapply(lst2, function(x) grep("RETENTIONTIME:",x, value=TRUE)))
  PMZL1=tryCatch({as.numeric(stringr::str_trim(gsub("RETENTIONTIME:","",PMZL)))},error=function(cond){return(NA)})
  PTY=unname(rapply(lst2, function(x) grep("PRECURSORTYPE:",x, value=TRUE)))
  PTY1=tryCatch({stringr::str_trim(gsub("PRECURSORTYPE:","",PTY))},error=function(cond){return(NA)})
  
  return(list(PMZL1,PTY1))
}
#################################
RF=sample(1:length(FLFRD),1)
CRDF=MaKlist(FLFRD[RF])[[1]]
PTV=MaKlist(FLFRD[RF])[[2]]
###########################

if(dim(RXF7)[1] > 2 & !all(is.na(CRDF)) & all(as.numeric(CRDF)) & !all(is.na(PTV)))
{
  return(1)
  ###print("select the script 1 for processing the data... it has all necessry data")
}

}
