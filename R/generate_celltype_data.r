generate_celltype_data <- function(){
  library(tidyverse)
  
  # PREP CORTEX SCT DATASET
  
  # Download the Zeisel expression data
  exp_file_CORT_ORIG = read.csv(url("https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt"),sep="\t",stringsAsFactors = FALSE)
  
  # Split it into expression and annotation info
  colnames      = exp_file_CORT_ORIG %>% dplyr::filter(tissue=="cell_id") %>% select(-c("X","tissue"))
  exp_file_CORT = exp_file_CORT_ORIG %>% dplyr::filter(X!="") %>% dplyr::filter(X!="(none)")  %>% `rownames<-`(.$X) %>% 
    dplyr::select(-"X") %>% `colnames<-`(colnames)
  
  annot_CORT    = exp_file_CORT_ORIG %>% dplyr::filter(X %in% c("","(none)")) %>% dplyr::filter(tissue!="") %>% dplyr::select(-X) %>% 
    `rownames<-`(.$tissue) %>% dplyr::select(-"tissue") %>% `colnames<-`(colnames) %>% t %>% as.data.frame(.,stringsAsFactors=FALSE)
  
  # PREP OLIGO SCT DATASET
  library(data.table)
  exp_file = as.data.frame(fread("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE75330&format=file&file=GSE75330_Marques_et_al_mol_counts2.tab.gz"),stringsAsFactors=FALSE)
  colnames(exp_file)[1] = "symbol"
  exp_file = exp_file[!duplicated(exp_file$symbol),]
  rownames(exp_file) = exp_file$symbol
  exp_file = exp_file[,-1]
  colnames(exp_file)=gsub("\\.","-",colnames(exp_file))
  
  # Get oligo annot info
  data <- GEOquery::getGEO('GSE75330', GSEMatrix=TRUE)
  annot_file=as.data.frame(data$GSE75330_series_matrix.txt.gz)
  rownames(annot_file) = annot_file$title
  annot_file = annot_file[colnames(exp_file),]
  annot_file = annot_file[!is.na(annot_file$title),]
  annot_OLIGO = data.frame(cell_id = annot_file$title,level1class=annot_file$inferred.cell.type.ch1,level2class=annot_file$inferred.cell.type.ch1)	
  annot_OLIGO = annot_OLIGO[annot_OLIGO$level2class!="(none)",]
  annot_OLIGO$level2class = as.character(annot_OLIGO$level2class)
  
  # Merge the two datasets
  exp1 = exp_file
  exp2 = exp_file_CORT
  annot1=annot_OLIGO
  annot2=annot_CORT
  merged=EWCE::merge_two_expfiles(exp1=exp1,exp2=exp2,annot1=annot1,annot2=annot_CORT)
  merged$annot$level1class=as.character(merged$annot$level1class)
  merged$annot$level2class=as.character(merged$annot$level2class)	
  
  
  # Drop and merge some annotations
  # - Drop "(none)" and "oligodendrocytes" from cortex dataset
  merged$exp = merged$exp[,merged$annot$level1class!="oligodendrocytes"]
  merged$annot = merged$annot[merged$annot$level1class!="oligodendrocytes",]
  
  # - Reword some of the names
  #merged$annot = gsub("endothelial-mural","Endothelial",merged$annot)	
  merged$annot$level1class[merged$annot$level1class=="endothelial-mural"] = merged$annot$level2class[merged$annot$level1class=="endothelial-mural"]
  merged$annot$level1class = gsub("Vsmc","Vascular Smooth Muscle Cell",merged$annot$level1class)
  merged$annot$level1class = gsub("Peric","Pericytes",merged$annot$level1class)
  merged$annot$level1class = gsub("Vend1","Vascular Endothelial",merged$annot$level1class)
  merged$annot$level1class = gsub("Vend2","Vascular Endothelial",merged$annot$level1class)
  #merged$annot$level1class = gsub("","",merged$annot$level1class)
  #merged$annot$level1class = gsub("","",merged$annot$level1class)
  merged$annot$level1class = gsub("astrocytes_ependymal","Astrocytes",merged$annot$level1class)
  merged$annot$level1class = gsub("pyramidal SS","Pyramidal Neurons",merged$annot$level1class)
  merged$annot$level1class = gsub("pyramidal CA1","Pyramidal Neurons",merged$annot$level1class)
  merged$annot$level1class = gsub("microglia","Microglia",merged$annot$level1class)
  merged$annot$level1class = gsub("interneurons","Interneurons",merged$annot$level1class)
  merged$annot$level1class = gsub("NFOL.*","Oligodendrocytes",merged$annot$level1class)
  merged$annot$level1class = gsub("MOL.*","Oligodendrocytes",merged$annot$level1class)
  merged$annot$level1class = gsub("OPC","Oligodendrocyte Precursor",merged$annot$level1class)
  merged$annot$level1class = gsub("COP","Oligodendrocyte Precursor",merged$annot$level1class)
  merged$annot$level1class = gsub("PPR","Vascular and Leptomeningeal Cells",merged$annot$level1class)
  #merged$annot = gsub("MFOL.*","Myelin-forming Oligodendrocytes",merged$annot)
  merged$annot$level1class = gsub("MFOL.*","Oligodendrocytes",merged$annot$level1class)
  
  # Split endothelial into Vascular Endothelial (Vend1+2) and Vascular Mural (Peric + Vsmc)
  
  annot=list()
  annot[[1]]=merged$annot$level1class
  
  fNames_OligosNCortex = EWCE::generate.celltype.data(exp=merged$exp,annotLevels=annot,groupName="OligosNCortex")
  load(fNames_OligosNCortex)
  
  ## GENERATE VIOLIN PLOTS
  #save(merged,file="mergedINPUT.rda")
  #for(i in 11095:dim(merged$exp)[1]){
  #  dat = data.frame(exp=as.vector(unlist(merged$exp[i,])),celltype=as.vector(unlist(merged$annot$level1class)))
  #  pdf(sprintf("GenePlots/%s.pdf",rownames(merged$exp)[i]),width=10,height=4.5)
  #  print(ggplot(dat)+geom_violin(aes(x=celltype,y=exp),fill="blue")+ylab("Molecules per cell")+xlab("Cell type")+ theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  #  dev.off()
  #}
  return(ctd)
}