drop_data_without_annotations <- function(all_dat,annot2){
  ###############################################
  ### DROP EXPRESSION DATA LACKING ANNOTATION ###
  
  good_dat  = all_dat[rownames(all_dat) %in% annot2$affy_mouse430_2,]
  # colnames(good_dat) = gsub(".CEL","",colnames(good_dat))
  pure_dat = good_dat[,as.character(sample_dat$names)]
  
  ### SAVE DATA TO GIVE TO SEB (OCT 2018)
  #rownames(annot3) = annot3$ID_REF
  #annot4 = annot3[rownames(pure_dat),]
  return(pure_dat)
}