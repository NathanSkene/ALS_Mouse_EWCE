get_probe_annots <- function(){
  ##################################
  ### GET THE PROBE ANNOTATIONS ####
  
  library("biomaRt")
  #human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  #mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  mouse = useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
  #attrib_mus = listAttributes(mouse)
  #attrib_mus[grep("Affy",attrib_mus[,2]),]
  #affy_hg_u133_plus_2
  #annot = getBM(attributes=c("affy_mouse430_2","external_gene_name"), filters="affy_mouse430_2", values=row.names(all_dat), mart=mouse)
  annot = getBM(attributes=c("affy_mouse430_2","external_gene_name"), mart=mouse)
  annot = annot[annot$affy_mouse430_2 %in% rownames(all_dat),]
  dup_probes = annot$affy_mouse430_2[duplicated(annot$affy_mouse430_2)]
  annot2 = annot[!(annot$affy_mouse430_2 %in% dup_probes),]
  annot3 = annot2
  colnames(annot3)[1] = "ID_REF"
  colnames(annot2)[2] = "MGI.symbol"
  return(annot2)
}