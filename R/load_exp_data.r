load_exp_data <- function(){
  ###################################
  ### LOAD THE EXPRESSION DATA ######
  tempLoc = tempdir()
  url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE18597&format=file"
  destfile <- sprintf("%s/GSE18597_RAW.tar",tempLoc)
  download.file(url, destfile, mode="wb")
  untar(destfile,exdir=tempLoc)
  
  library(affy)
  data <- ReadAffy(celfile.path = tempLoc) # Read in the CEL files in the directory, then normalize the data
  eset <- affy::rma(data)
  #write.exprs(eset,file="data.txt") # Save data to file (Data is log2 transformed and normalized)
  #exprs = exprs(eset)
  all_dat = exprs(eset)
  colnames(all_dat) = gsub("_.*","",colnames(all_dat))
  return(all_dat)
}