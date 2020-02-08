source("/Users/natske/Google Drive/EWCE/R/generate.bootstrap.plots.r")
source("/Users/natske/Google Drive/EWCE/R/read_celltype_data.r")
source("/Users/natske/Google Drive/EWCE/R/get_annot_below.r")
source("/Users/natske/Google Drive/EWCE/R/bootstrap.enrichment.test.r")
source("/Users/natske/Google Drive/EWCE/R/ewce.plot.r")
source("/Users/natske/Google Drive/EWCE/R/merge_two_expfiles.r")
source("/Users/natske/Google Drive/EWCE/R/ewce_expression_data.r")
source("/Users/natske/Google Drive/EWCE/R/get_summed_proportions.r")
source("/Users/natske/Google Drive/EWCE/R/cell.list.dist.r")
library(reshape)
library(boot)

#source("http://bioconductor.org/biocLite.R")
#biocLite("affy")
#biocLite("oligo")
#biocLite("limma")

#mouse4302cdf

###################################
### LOAD THE EXPRESSION DATA ######
setwd("/Users/natske/Datasets that are too large to store elsewhere/SOD1 Spinal Cord (GSE18597)/GSE18597_RAW")
library(affy)
data <- ReadAffy() # Read in the CEL files in the directory, then normalize the data
eset <- rma(data)
write.exprs(eset,file="data.txt") # Save data to file (Data is log2 transformed and normalized)
#exprs = exprs(eset)
all_dat = exprs(eset)
colnames(all_dat) = gsub("_.*","",colnames(all_dat))

setwd("/Users/natske/Datasets that are too large to store elsewhere/SOD1 Spinal Cord (GSE18597)")

###################################
### GET THE SAMPLE ANNOTATIONS ####

sample_data = read.csv("/Users/natske/Datasets that are too large to store elsewhere/SOD1 Spinal Cord (GSE18597)/sample_annot.csv")
sample_age     = sample_data[,"Age"]
sample_genotype  = sample_data[,"Genotype"]
sample_dat     = data.frame(age=sample_age,genotype=sample_genotype,names=sample_data$ID_REF)
rownames(sample_dat) = sample_dat$names

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

###############################################
### DROP EXPRESSION DATA LACKING ANNOTATION ###

good_dat  = all_dat[rownames(all_dat) %in% annot2$affy_mouse430_2,]
# colnames(good_dat) = gsub(".CEL","",colnames(good_dat))
pure_dat = good_dat[,as.character(sample_dat$names)]

### SAVE DATA TO GIVE TO SEB (OCT 2018)
rownames(annot3) = annot3$ID_REF
annot4 = annot3[rownames(pure_dat),]


###############################
### REGRESS OUT THE AGE EFFECT
library(limma)
age      = as.factor(gsub(" days","",sample_dat$age))
mod  = model.matrix(~age)
fit = lmFit(pure_dat,mod)
resid_dat = residuals(fit,pure_dat)


genotype      = sample_dat$genotype

rownames(annot2)=annot2[,1]
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

#############################
### EWCE: LOAD THE OLIGO DATA

regenerate_celltype_data = FALSE
thresh=0.00
trim=0.0
wtLEVELS="woLev"
if(regenerate_celltype_data){
	# PREP CORTEX SCT DATASET
	data_path = "/Users/natske/Datasets that are too large to store elsewhere/Cell type dataset---Zeisel et al Cell Type dataset"
	exp_file_CORT = read.csv(sprintf("%s/expression_mRNA_17-Aug-2014_expressionOnly.txt",data_path),sep="\t")
	colnames(exp_file_CORT) = gsub("^X","",colnames(exp_file_CORT))
	exp_file_CORT = exp_file_CORT[!duplicated(exp_file_CORT$symbol),]
	rownames(exp_file_CORT) = exp_file_CORT$symbol
	exp_file_CORT = exp_file_CORT[,-1]
	annot_file_CORT = read.csv(sprintf("%s/expression_mRNA_17-Aug-2014_cellAnnot.txt",data_path),sep="\t")
	annot_CORT = data.frame(cell_id = annot_file_CORT$cell_id,level1class=annot_file_CORT$level1class,level2class=annot_file_CORT$level2class)
	annot_CORT = annot_CORT[annot_CORT$level2class!="(none)",]
	annot_CORT$level2class = as.character(annot_CORT$level2class)
	
	# PREP OLIGO SCT DATASET
	data_path = "/Users/natske/Datasets that are too large to store elsewhere/Cell type dataset---Oligodendrocytes (Marques)/"
	exp_file = read.csv(sprintf("%s/GSE75330_Marques_et_al_mol_counts2.tab",data_path),sep="\t")
	colnames(exp_file)[1] = "symbol"
	exp_file = exp_file[!duplicated(exp_file$symbol),]
	rownames(exp_file) = exp_file$symbol
	exp_file = exp_file[,-1]
	colnames(exp_file)=gsub("\\.","-",colnames(exp_file))
	annot_file = read.csv(sprintf("%sCellAnnot.csv",data_path),sep=",")
	#rownames(annot_file) = gsub("-",".",annot_file$CellID)
	rownames(annot_file) = annot_file$CellID
	annot_file = annot_file[colnames(exp_file),]
	annot_file = annot_file[!is.na(annot_file$CellType),]
	annot_OLIGO = data.frame(cell_id = annot_file$CellID,level1class=annot_file$CellType,level2class=annot_file$CellType)	
	annot_OLIGO = annot_OLIGO[annot_OLIGO$level2class!="(none)",]
	annot_OLIGO$level2class = as.character(annot_OLIGO$level2class)
	
	# Merge the two datasets
	exp1 = exp_file
	exp2 = exp_file_CORT
	annot1=annot_OLIGO
	annot2=annot_CORT
	merged=merge_two_expfiles(exp1=exp1,exp2=exp2,annot1=annot1,annot2=annot_CORT)
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
	
	celltype_data = read_celltype_data(exp=merged$exp,annot=annot,thresh=thresh,trim=trim)
	
	acs = colnames(celltype_data[[1]]$all_scts)
	celltype_data[[1]]$all_scts = celltype_data[[1]]$all_scts[,order(acs)]
	celltype_data[[1]]$cell_dists = celltype_data[[1]]$cell_dists[,order(acs)]
	#acs = colnames(celltype_data[[1]]$all_scts)
	#celltype_data[[1]]$all_scts = celltype_data[[1]]$all_scts[,!is.na(acs)]
	#celltype_data[[1]]$cell_dists = celltype_data[[1]]$cell_dists[,!is.na(acs)]
	
	# Drop the (none) cells
	#celltype_data[[1]]$all_scts = celltype_data[[1]]$all_scts[,colnames(celltype_data[[1]]$all_scts)!="(none)"]
	#celltype_data[[1]]$cell_dists = celltype_data[[1]]$cell_dists[,colnames(celltype_data[[1]]$cell_dists)!="(none)"]
	
	save(celltype_data,file=sprintf("celltype_data_OligosNCortex_(%s)_thresh(%s)_trim(%s).rda",wtLEVELS,thresh,trim))
	
	## GENERATE VIOLIN PLOTS
	save(merged,file="mergedINPUT.rda")
	for(i in 11095:dim(merged$exp)[1]){
	    dat = data.frame(exp=as.vector(unlist(merged$exp[i,])),celltype=as.vector(unlist(merged$annot$level1class)))
	    pdf(sprintf("GenePlots/%s.pdf",rownames(merged$exp)[i]),width=10,height=4.5)
	    print(ggplot(dat)+geom_violin(aes(x=celltype,y=exp),fill="blue")+ylab("Molecules per cell")+xlab("Cell type")+ theme(axis.text.x = element_text(angle = 45, hjust = 1)))
	    dev.off()
	}
}else{
	load(file=sprintf("celltype_data_OligosNCortex_(%s)_thresh(%s)_trim(%s).rda",wtLEVELS,thresh,trim))
}

#############################
### EWCE: RUN THE ANALYSIS

tt_28 = topTable(eb, coef="MUT_p42", adjust="BH",number=1000000)
tt_42 = topTable(eb, coef="MUT_p42", adjust="BH",number=1000000)
tt_56 = topTable(eb, coef="MUT_p56", adjust="BH",number=1000000)
tt_70 = topTable(eb, coef="MUT_p70", adjust="BH",number=1000000)
tt_98 = topTable(eb, coef="MUT_p98", adjust="BH",number=1000000)
tt_112 = topTable(eb, coef="MUT_p112", adjust="BH",number=1000000)
tt_126 = topTable(eb, coef="MUT_p126", adjust="BH",number=1000000)

colnames(tt_28)[1]="MGI.symbol"
colnames(tt_42)[1]="MGI.symbol"
colnames(tt_56)[1]="MGI.symbol"
colnames(tt_70)[1]="MGI.symbol"
colnames(tt_98)[1]="MGI.symbol"
colnames(tt_112)[1]="MGI.symbol"
colnames(tt_126)[1]="MGI.symbol"

graph_theme = theme_bw(base_size = 12, base_family = "Helvetica") +
     theme(panel.grid.major = element_line(size = .5, color = "grey"),
     axis.line = element_line(size=.7, color = "black"),legend.position = c(0.85, 0.7), text = element_text(size=14),
     axis.title.x = element_text(vjust = -0.35), axis.title.y = element_text(vjust = 0.6)) + theme(legend.position="none")

for(thresh in c(100)) #,250,350)){#c(100,150,200,250,300,350)){
	
	#full_res_28 = ewce_expression_data(celltype_data,annot=celltype_data[[1]]$annot,annotLevel=1,tt_28,sortBy="t",thresh=thresh,reps=10000,useHGNC=FALSE)
	#full_res_42 = ewce_expression_data(celltype_data,annot=celltype_data[[1]]$annot,annotLevel=1,tt_42,sortBy="t",thresh=thresh,reps=10000,useHGNC=FALSE)
	#full_res_56 = ewce_expression_data(celltype_data,annot=celltype_data[[1]]$annot,annotLevel=1,tt_56,sortBy="t",thresh=thresh,reps=10000,useHGNC=FALSE)
	#full_res_70 = ewce_expression_data(celltype_data,annot=celltype_data[[1]]$annot,annotLevel=1,tt_70,sortBy="t",thresh=thresh,reps=10000,useHGNC=FALSE)
	#full_res_98 = ewce_expression_data(celltype_data,annot=celltype_data[[1]]$annot,annotLevel=1,tt_98,sortBy="t",thresh=thresh,reps=10000,useHGNC=FALSE)
	#full_res_112 = ewce_expression_data(celltype_data,annot=celltype_data[[1]]$annot,annotLevel=1,tt_112,sortBy="t",thresh=thresh,reps=10000,useHGNC=FALSE)
	#full_res_126 = ewce_expression_data(celltype_data,annot=celltype_data[[1]]$annot,annotLevel=1,tt_126,sortBy="t",thresh=thresh,reps=10000,useHGNC=FALSE)
	
	#combined_ewce = list(full_res_28=full_res_28,full_res_42=full_res_42,full_res_56=full_res_56,full_res_70=full_res_70,full_res_98=full_res_98,full_res_112=full_res_112,full_res_126=full_res_126)
	#save(combined_ewce,file=sprintf("SOD1_EWCE_RES_Thresh%s.rda",thresh))
    load(file=sprintf("SOD1_EWCE_RES_Thresh%s.rda",thresh))
    
    full_res_28 = combined_ewce$full_res_28
    full_res_42 = combined_ewce$full_res_42
    full_res_56 = combined_ewce$full_res_56
    full_res_70 = combined_ewce$full_res_70
    full_res_98 = combined_ewce$full_res_98
    full_res_112 = combined_ewce$full_res_112
    full_res_126 = combined_ewce$full_res_126
    
    	
	all_res = rbind(cbind(full_res_28$joint_results,Age=28),cbind(full_res_42$joint_results,Age=42),cbind(full_res_56$joint_results,Age=56),cbind(full_res_70$joint_results,Age=70))
	all_res = rbind(all_res,cbind(full_res_98$joint_results,Age=98),cbind(full_res_112$joint_results,Age=112),cbind(full_res_126$joint_results,Age=126))
	all_res$q = p.adjust(all_res$p,method="BH")
	
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
	
	
	pdf(sprintf("Enrichment Plots/SOD1_EWCE_Enrichments_Thresh%s.pdf",thresh),width=12,height=5)
	print(ggplot(all_res2) + geom_line(aes(x=as.factor(Age),y=sd_from_mean,group=CellType),stat="identity") + geom_point(aes(x=as.factor(Age),y=sd_from_mean),stat="identity") + facet_grid(Direction~CellType) + geom_text(aes(x=as.factor(Age),y=sd_from_mean+0.5,label=ast,color="red",size=5),stat="identity")+
		xlab("Age (in days)")+ylab("z-score")+graph_theme+ theme(axis.text.x = element_text(angle = 90, hjust = 1)))
	dev.off()
	
	write.csv(all_res2,file=sprintf("Enrichment Results tables/SOD1_enrichment_results_Thresh%s.csv",thresh))
}

load("/Users/natske/Datasets that are too large to store elsewhere/SOD1 Spinal Cord (GSE18597)/tt.rda")
setwd("/Users/natske/Datasets that are too large to store elsewhere/SOD1 Spinal Cord (GSE18597)/")
ages = c(28,42,56,70,98,112,126)
tt_list = list(tt_28,tt_42,tt_56,tt_70,tt_98,tt_112,tt_126)
source("/Users/natske/Datasets that are too large to store elsewhere/SOD1 Spinal Cord (GSE18597)/generate.bootstrap.plots.for.transcriptome_SOD1.r")
load("/Users/natske/Datasets that are too large to store elsewhere/SOD1 Spinal Cord (GSE18597)/celltype_data_OligosNCortex_(woLev)_thresh(0)_trim(0).rda")


#for(thresh in c(100,150,200,250,300,350)){
for(thresh in c(100)){
    load(file=sprintf("SOD1_EWCE_RES_Thresh%s.rda",thresh))
    for(ageI in 1:length(ages)){
        print(sprintf("Age: %s",ages[ageI]))
        tt=tt_list[[ageI]]
        generate.bootstrap.plots.for.transcriptome(sct_data=celltype_data,tt=tt,thresh=thresh,annotLevel=1,reps=1000,full_results=combined_ewce[[ageI]],listFileName=sprintf("%sdays_MOUSE",ages[ageI]))
    }
}


# Plot the genes (across ages) which contribute most to each celltype change, at each age
for(thresh in c(100,150,200,250,300,350)){
    load(file=sprintf("SOD1_EWCE_RES_Thresh%s.rda",thresh))
    for(ageJJJ in 1:length(ages)){
        # Which celltypes are significant at age ageJJJ with thresh=thresh?
        for(dirS in c("Up","Down")){
            res = combined_ewce[[ageI]]$joint_results
            res = res[res$Direction==dirS,]
            ccs = as.character(res$CellType[res$p<0.05])
            for(cc in ccs){
                # Find the genes which could have contributed (because they are within the thresh of top table)
                tt = tt_list[[ageJJJ]]
                if(dirS=="Up"){tt=tt[order(tt$t,decreasing=TRUE),]}
                if(dirS=="Down"){tt=tt[order(tt$t,decreasing=FALSE),]}
                tt = tt[!duplicated(tt$MGI.symbol),]
                possGenes = tt[1:thresh,"MGI.symbol"]
                possGenes = possGenes[possGenes %in% rownames(celltype_data[[1]]$cell_dists)]
                
                # Find the five of those with highest exp in the celltype
                possExp = celltype_data[[1]]$cell_dists[possGenes,]
                
                genes = names(sort(possExp[,cc],decreasing=TRUE))[1:5]#c("Nupr1","Apod","Timp1")
                expr = matrix(0,nrow=length(ages),ncol=length(genes))
                rownames(expr) = ages
                colnames(expr) = genes
                
                abs_max <- function(mat){
                    return(mat[which.max( abs(mat) )])
                }
                for(gene in genes){
                    for(ageI in 1:length(ages)){
                        expr[ageI,gene] = abs_max(tt_list[[ageI]][tt_list[[ageI]]$MGI.symbol==gene,]$logFC)
                    }
                }
                expr2 = melt(expr)
                colnames(expr2) = c("Age","Gene","logFC")
                graph_theme = theme_bw(base_size = 12, base_family = "Helvetica") +
                    theme(panel.grid.major = element_line(size = .5, color = "grey"),
                          axis.line = element_line(size=.7, color = "black"), text = element_text(size=14),
                          axis.title.x = element_text(vjust = -0.35), axis.title.y = element_text(vjust = 0.6))# + theme(legend.position="none")
                the_plot = ggplot(expr2)+geom_path(aes(x=factor(Age),y=logFC,color=Gene,group=Gene))+graph_theme+ggtitle(sprintf("%s genes",cc))+
                        xlab("Age (in days)")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
                pdf(file=sprintf("TimeCoursePlots/%s_%s_age%s_%s",cc,thresh,ages[ageJJJ],dirS),width=4.5,height=5)
                print(the_plot)
                dev.off()
            }
        }
    }
}

#pdf(sprintf("Fig_EWCE_MBB_annot1_thresh(%s)_trim(%s)_SchizGWAS.pdf",thresh,trim),width=10,height=5)
#print(ewce.plot(full_results_1$joint_results,mtc_method="BH"))
#dev.off()


