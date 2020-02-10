get_sample_annots <- function(){
  ###################################
  ### GET THE SAMPLE ANNOTATIONS ####
  load("Data/sample_data.Rda")
  sample_age     = sample_data[,"Age"]
  sample_genotype  = sample_data[,"Genotype"]
  sample_dat     = data.frame(age=sample_age,genotype=sample_genotype,names=sample_data$ID_REF)
  rownames(sample_dat) = sample_dat$names
  return(sample_dat)
}