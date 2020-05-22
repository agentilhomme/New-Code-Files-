library(xlsx)
library(ape)
library(phytools)
library(phylogram)
library(dendextend)
library(dplyr)
library(gplots)

#load files 
hits = read.table("uniqhits.txt")
hits$V1 = gsub(":","|",hits$V1)
feat = read.table("Final_Colwellia_Group_Features.tsv",fill=T,header=T, sep = '\t')
feathit = feat[match(hits$V1,feat$feature.patric_id),]

write.table(unique(feathit$feature.pgfam_id),file = "best_GPF_blast_match.txt")
saveRDS(unique(feathit$feature.pgfam_id),file = "best_GPF_blast_match.RDS")

Gene_Prescence = read.csv("Gene_Prescence.csv",fill = T,header = T, check.names = F)

#OGT.Data = read.table("./input/All_OGT_Data.tsv",header = T, fill = T,sep = '\t', colClasses = "character")
#OGT.Data = read.table("./input/All_OGT_Data.tsv",header = T, fill = T,sep = '\t', colClasses = "character")
#rownames(OGT.Data) = OGT.Data$patricID


treefile="My_Colwellia_Tree.nwk"


blast_match = as.character(readRDS("best_GPF_blast_match.RDS"))
rownames(Gene_Prescence) = Gene_Prescence[,1]; Gene_Prescence[,1] = NULL
GeneP1 = Gene_Prescence[blast_match,]
GeneP1 = na.omit(GeneP1)
GeneP1 = GeneP1[rowSums(GeneP1) > 0,colSums(GeneP1) > 0] 

tre.phylo = read.tree(treefile)
common_ids = intersect(tre.phylo$tip.label,colnames(GeneP1))
GeneP1 = GeneP1[,common_ids]
dend = keep.tip(tre.phylo,common_ids)

dend$node.label = NULL # Need to remove node labels
dend = midpoint.root(dend)
dend = ladderize(dend)
dend = drop.tip(dend,"1816218.31") #drop long tip (poor genome assembly)


# convert from 'phylo' to 'dendrogram'
dend1 = read.dendrogram(textConnection(write.tree(dend))) 
dend1 = as.dendrogram(dend1)
dend2 = as.hclust(dend1)
#Test_OGT = OGT.Data[ !(OGT.Data$patricID %in% "1816218.31"),]
#Test_OGT = Test_OGT$strain[dend2$order]
#labels(dend1) = Test_OGT
dend1 = set(dend1, "branches_lwd", 5)

#change colum names to math tree 
#colnames(GeneP1) = OGT.Data[colnames(GeneP1),"strain"]

# #change row global protein family id to description 
# #PFG_name = feathit[,c(3,5)]
# #PFG_name = na.omit(PFG_name)
# PFG_name = unique(PFG_name)
# # if protein family has same protein id but diferent discriptors then combine them
# setDT(PFG_name)[,feature.product := as.character(feature.product)][,feature.product := paste0(as.character(feature.product), collapse = ","), by = feature.pgfam_id]
# PFG_name = unique(PFG_name)
# rownames(PFG_name) = PFG_name$feature.pgfam_id
# common_PGF = intersect(PFG_name$feature.pgfam_id,rownames(GeneP1))
# PFG_name = as.data.frame(PFG_name)
# rownames(PFG_name) = PFG_name$feature.pgfam_id
# PFG_name = PFG_name[common_PGF,]



# dup_product = data.frame(PFG_name[duplicated(PFG_name$feature.product),])
# dup_product = as.vector(unique(dup_product$feature.product))

# for(i in 1:nrow(PFG_name)){
#   if(stri_duplicated(PFG_name$feature.product[[i]]) == T ){
#     PFG_name[i,2] = paste0(PFG_name$feature.pgfam_id[[i]]," (",PFG_name$feature.pgfam_id[[i]],")")
#   }
# }
# 
# PFG_name = make.unique(PFG_name$feature.product, sep = "_")
# 
# rownames(GeneP1) = PFG_name[rownames(GeneP1),"feature.product"]


GeneP1 = GeneP1[,labels(dend1)]

tiff("Usual_Suspects_Presence.tiff",width = 20,height = 8.5, units = "in", res = 500)
heatmap.2(as.matrix(t(GeneP1>0)+0),scale = "none",trace="none") # presence or abscence
dev.off()

tiff("Usual_Suspects_Count.tiff",width = 20,height = 8.5, units = "in", res = 500)
heatmap.2(as.matrix(t(GeneP1)^.25),scale = "none",trace="none",Rowv = dend1,margins=c(8,16),col=viridis) # count of proteins in protein families rescaled to the 4th root 
dev.off()


                  


# --------------------------------Make Usual Suspects dataframe with PGFams ID -----------------------------

# Try searching Usual Suspects acronyms and subsetting GF Patric ID first
Acr_genes = as.vector(Usual_Suspects$Gene_Acronym)
Acr_genes = na.omit(Acr_genes)
