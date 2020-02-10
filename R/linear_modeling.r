linear_modelling <- function(sample_dat,pure_dat,probe_annots){
  ###############################
  ### REGRESS OUT THE AGE EFFECT
  library(limma)
  age      = as.factor(gsub(" days","",sample_dat$age))
  mod  = model.matrix(~age)
  fit = lmFit(pure_dat,mod)
  resid_dat = residuals(fit,pure_dat)
  
  genotype      = sample_dat$genotype
  
  rownames(annot2)=annot2$cell_id
  annot4 = annot2[rownames(pure_dat),]
  #rownames(pure_dat) = annot4$external_gene_name
  
  design = model.matrix(~0+age)
  for(aaa in unique(age)){
    tmp_mut = as.numeric(genotype!="wild-type" & age==aaa)
    tmp = data.frame(tmp_mut)
    colnames(tmp) = c(sprintf("MUT_p%s",aaa))
    design = cbind(design,tmp)
  }
  
  fit1 = lmFit(pure_dat,design)
  eb = eBayes(fit1)
  tt_28 = topTable(eb, coef="MUT_p42", adjust="BH",number=1000000)
  tt_42 = topTable(eb, coef="MUT_p42", adjust="BH",number=1000000)
  tt_56 = topTable(eb, coef="MUT_p56", adjust="BH",number=1000000)
  tt_70 = topTable(eb, coef="MUT_p70", adjust="BH",number=1000000)
  tt_98 = topTable(eb, coef="MUT_p98", adjust="BH",number=1000000)
  tt_112 = topTable(eb, coef="MUT_p112", adjust="BH",number=1000000)
  tt_126 = topTable(eb, coef="MUT_p126", adjust="BH",number=1000000)
  
  allTT = list(tt_28=tt_28,tt_42=tt_42,tt_56=tt_56,tt_70=tt_70,tt_98=tt_98,tt_112=tt_112,tt_126=tt_126)
  
  # colnames(allTT$tt_28)[1]="MGI.symbol"
  # colnames(allTT$tt_42)[1]="MGI.symbol"
  # colnames(allTT$tt_56)[1]="MGI.symbol"
  # colnames(allTT$tt_70)[1]="MGI.symbol"
  # colnames(allTT$tt_98)[1]="MGI.symbol"
  # colnames(allTT$tt_112)[1]="MGI.symbol"
  # colnames(allTT$tt_126)[1]="MGI.symbol"
  add_mgi <- function(tt){
    tt$affy_mouse430_2 = rownames(tt)
    tt2 = merge(tt,probe_annots,by="affy_mouse430_2")
    tt2 = tt2[order(tt2$P.Value),]
    return(tt2)
  }
  allTT2 = lapply(allTT,FUN=add_mgi)
  
  return(allTT2)
}