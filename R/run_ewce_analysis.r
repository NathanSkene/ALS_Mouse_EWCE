run_ewce_analysis <- function(thresh = 100,annotLevel=1,reps=1000,graph_theme){
  ewce_call_func <- function(x,sct_data,annotLevel,thresh,reps){
    full_res_OUT = EWCE::ewce_expression_data(sct_data=sct_data,annotLevel=annotLevel,tt=x,sortBy="t",thresh=thresh,reps=reps)
    return(full_res_OUT)
  }
  combined_ewce = lapply(allTT,FUN=ewce_call_func,sct_data=ctd,annotLevel=annotLevel,thresh=thresh,reps=reps)
  
  # Get a single results table
  for(i in i:length(combined_ewce)){combined_ewce[[i]]$joint_results$Age=as.numeric(gsub(".*_","",names(combined_ewce)[i]))}
  all_res = combined_ewce[[1]]$joint_results
  for(i in 2:length(combined_ewce)){all_res=rbind(all_res,combined_ewce[[i]]$joint_results)}
  all_res$q = p.adjust(all_res$p,method="BH")
  
  # Ignore all depleted cell types
  all_res2 = all_res
  all_res2$sd_from_mean[all_res2$sd_from_mean<0]=0
  
  all_res2$CellType = gsub("Vascular and Leptomeningeal Cells","VLMC",all_res2$CellType)
  all_res2$CellType = gsub("Oligodendrocyte Precursor","OPC",all_res2$CellType)
  library(ggplot2)
  ast = all_res2$q
  ast[all_res2$q>0.05]=""
  ast[all_res2$q<0.05]="*"
  all_res2$ast = ast
  all_res2$CellType = gsub("Vascular Endothelial","Vascular \nEndothelial",all_res2$CellType)
  all_res2$CellType = gsub("Vascular Endothelial","Vascular \nEndothelial",all_res2$CellType)
  all_res2$CellType = gsub("Vascular Smooth Muscle Cell","Vascular\nSmooth\nMuscle Cell",all_res2$CellType)
  all_res2$CellType = gsub("Oligodendrocytes","Oligo-\ndendrocytes",all_res2$CellType)
  all_res2$CellType = gsub("Pyramidal Neurons","Pyramidal\nNeurons",all_res2$CellType)
  
  pdf(sprintf("Results/Figures/SOD1_EWCE_Enrichments_Thresh%s.pdf",thresh),width=12,height=5)
  print(ggplot(all_res2) + geom_line(aes(x=as.factor(Age),y=sd_from_mean,group=CellType),stat="identity") + geom_point(aes(x=as.factor(Age),y=sd_from_mean),stat="identity") + facet_grid(Direction~CellType) + geom_text(aes(x=as.factor(Age),y=sd_from_mean+0.5,label=ast,color="red",size=5),stat="identity")+
          xlab("Age (in days)")+ylab("z-score")+graph_theme+ theme(axis.text.x = element_text(angle = 90, hjust = 1)))
  dev.off()
  
  write.csv(all_res2,file=sprintf("Results/Tables/SOD1_enrichment_results_Thresh%s.csv",thresh))
  
  ewce_res = list(combined_ewce=combined_ewce,ewce_table=all_res2)
  return(ewce_res)
}