# Set vars
reps = 1000
thresh = 100
annotLevel = 1

# Load libraries
if (!require("pacman")) install.packages("pacman")
if (!require("BiocManager")) install.packages("BiocManager")
pacman::p_load("reshape","boot","ggplots2","cowplot","affy","oligo","limma","mouse4302cdf","R.utils","biomaRt","drake",try.bioconductor=TRUE)

R.utils::sourceDirectory("R")

# Setup directory structure
dir.create("Results",showWarnings = FALSE)
dir.create("Results/Tables",showWarnings = FALSE)
dir.create("Results/Figures",showWarnings = FALSE)

all_dat = load_exp_data()

sample_dat = get_sample_annots()

probe_annots = get_probe_annots()

pure_dat = drop_data_without_annotations(all_dat,probe_annots)

allTT = linear_modelling(sample_dat,pure_dat,probe_annots)

ctd = generate_celltype_data() 

graph_theme = get_graph_theme()

ewce_res = run_ewce_analysis(thresh = 100,annotLevel=1,reps=1000,graph_theme)

# plot_bootstrap_plots = function(thresh = 100,ctd,ewce_res,annotLevel,reps,allTT){
#     for(ageI in 1:length(ages)){
#         print(sprintf("Age: %s",ages[ageI]))
#         tt=tt_list[[ageI]]
#         generate.bootstrap.plots.for.transcriptome(sct_data=celltype_data,tt=tt,thresh=thresh,annotLevel=1,reps=1000,full_results=combined_ewce[[ageI]],listFileName=sprintf("%sdays_MOUSE",ages[ageI]))
#     }
# }
# 
# 
# # Plot the genes (across ages) which contribute most to each celltype change, at each age
# for(thresh in c(100,150,200,250,300,350)){
#     load(file=sprintf("SOD1_EWCE_RES_Thresh%s.rda",thresh))
#     for(ageJJJ in 1:length(ages)){
#         # Which celltypes are significant at age ageJJJ with thresh=thresh?
#         for(dirS in c("Up","Down")){
#             res = combined_ewce[[ageI]]$joint_results
#             res = res[res$Direction==dirS,]
#             ccs = as.character(res$CellType[res$p<0.05])
#             for(cc in ccs){
#                 # Find the genes which could have contributed (because they are within the thresh of top table)
#                 tt = tt_list[[ageJJJ]]
#                 if(dirS=="Up"){tt=tt[order(tt$t,decreasing=TRUE),]}
#                 if(dirS=="Down"){tt=tt[order(tt$t,decreasing=FALSE),]}
#                 tt = tt[!duplicated(tt$MGI.symbol),]
#                 possGenes = tt[1:thresh,"MGI.symbol"]
#                 possGenes = possGenes[possGenes %in% rownames(celltype_data[[1]]$cell_dists)]
#                 
#                 # Find the five of those with highest exp in the celltype
#                 possExp = celltype_data[[1]]$cell_dists[possGenes,]
#                 
#                 genes = names(sort(possExp[,cc],decreasing=TRUE))[1:5]#c("Nupr1","Apod","Timp1")
#                 expr = matrix(0,nrow=length(ages),ncol=length(genes))
#                 rownames(expr) = ages
#                 colnames(expr) = genes
#                 
#                 abs_max <- function(mat){
#                     return(mat[which.max( abs(mat) )])
#                 }
#                 for(gene in genes){
#                     for(ageI in 1:length(ages)){
#                         expr[ageI,gene] = abs_max(tt_list[[ageI]][tt_list[[ageI]]$MGI.symbol==gene,]$logFC)
#                     }
#                 }
#                 expr2 = melt(expr)
#                 colnames(expr2) = c("Age","Gene","logFC")
#                 graph_theme = theme_bw(base_size = 12, base_family = "Helvetica") +
#                     theme(panel.grid.major = element_line(size = .5, color = "grey"),
#                           axis.line = element_line(size=.7, color = "black"), text = element_text(size=14),
#                           axis.title.x = element_text(vjust = -0.35), axis.title.y = element_text(vjust = 0.6))# + theme(legend.position="none")
#                 the_plot = ggplot(expr2)+geom_path(aes(x=factor(Age),y=logFC,color=Gene,group=Gene))+graph_theme+ggtitle(sprintf("%s genes",cc))+
#                         xlab("Age (in days)")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
#                 pdf(file=sprintf("TimeCoursePlots/%s_%s_age%s_%s",cc,thresh,ages[ageJJJ],dirS),width=4.5,height=5)
#                 print(the_plot)
#                 dev.off()
#             }
#         }
#     }
# }
# 
# #pdf(sprintf("Fig_EWCE_MBB_annot1_thresh(%s)_trim(%s)_SchizGWAS.pdf",thresh,trim),width=10,height=5)
# #print(ewce.plot(full_results_1$joint_results,mtc_method="BH"))
# #dev.off()
# 
# 
