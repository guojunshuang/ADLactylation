setwd("~/singlecell/Lactylation/AD/GSE157827")
fs=list.files('./GSE157827_RAW/','^GSM')
fs
library(tidyverse)
samples=str_split(fs,'_',simplify = T)[,1]
lapply(unique(samples),function(x){
  # x = unique(samples)[1]
  y=fs[grepl(x,fs)]
  folder=paste0("GSE157827_RAW/", paste(str_split(y[1],'_',simplify = T)[,1:2], collapse = "_"))
  dir.create(folder,recursive = T)

  file.rename(paste0("GSE157827_RAW/",y[1]),file.path(folder,"barcodes.tsv.gz"))

  file.rename(paste0("GSE157827_RAW/",y[2]),file.path(folder,"features.tsv.gz"))
  file.rename(paste0("GSE157827_RAW/",y[3]),file.path(folder,"matrix.mtx.gz"))
})
library(Seurat)
samples=list.files("GSE157827_RAW/")
samples
dir <- file.path('./GSE157827_RAW',samples)
names(dir) <- samples
counts <- Read10X(data.dir = dir)
scRNA = CreateSeuratObject(counts)
dim(scRNA)   
table(scRNA@meta.data$orig.ident)
head(scRNA@meta.data)
scRNA$orig.ident = stringr::str_split(rownames(scRNA@meta.data), "_",simplify = T)[,2]
scRNA$group = substr(scRNA$orig.ident,1,2)
head(scRNA@meta.data)
table(scRNA$group)
table(scRNA$orig.ident)
ranks = order(as.numeric(substr(unique(scRNA$orig.ident),3,4)))
scRNA$orig.ident = factor(scRNA$orig.ident, levels = unique(scRNA$orig.ident)[ranks])
table(scRNA$orig.ident, scRNA$group)
feats <- c("nFeature_RNA", "nCount_RNA")
library(patchwork)
p_filt_b_1=VlnPlot(scRNA, features = feats, pt.size = 0, ncol = 2, 
                   group.by = "group") + NoLegend()
p_filt_b_2=VlnPlot(scRNA, features = feats, pt.size = 0, ncol = 1, 
                   group.by = "orig.ident") + NoLegend()
mito_genes=rownames(scRNA)[grep("^MT-", rownames(scRNA))] ;mito_genes 
scRNA=PercentageFeatureSet(scRNA, "^MT-", col.name = "percent_mito")
fivenum(scRNA@meta.data$percent_mito)
ribo_genes=rownames(scRNA)[grep("^Rp[sl]", rownames(scRNA),ignore.case = T)];ribo_genes
scRNA=PercentageFeatureSet(scRNA, "^RP[SL]", col.name = "percent_ribo")
fivenum(scRNA@meta.data$percent_ribo)
feats <- c("percent_mito","percent_ribo")
p_filt_b_3 =VlnPlot(scRNA, features = feats, pt.size = 0, ncol = 2, 
                    group.by = "group") + NoLegend()
p_filt_b_4=VlnPlot(scRNA, features = feats, pt.size = 0, ncol = 1, 
                   group.by = "orig.ident") + NoLegend()
(p_filt_b_1 / p_filt_b_2) | (p_filt_b_3 / p_filt_b_4)
retained_c_umi_low <- scRNA$nFeature_RNA > 300
retained_c_umi_high <- scRNA$nFeature_RNA < 8000
retained_c_mito <- scRNA$percent_mito < 14
retained_c_ribo <- scRNA$percent_ribo < 3
table(retained_c_mito & retained_c_ribo & retained_c_umi_low & retained_c_umi_high)
table(scRNA$group[retained_c_mito & retained_c_ribo & retained_c_umi_low & retained_c_umi_high])
feature_rowsum = Matrix::rowSums(scRNA@assays$RNA@counts>0) 
head(feature_rowsum)
retained_f_low <- feature_rowsum > ncol(scRNA)*0.005
table(retained_f_low)
rankplot = data.frame(feature_count = sort(feature_rowsum),
                      gene = names(sort(retained_f_low)),
                      Rank = 1:length(feature_rowsum))
scRNA_filt = scRNA[retained_f_low, retained_c_mito & retained_c_ribo & 
                     retained_c_umi_low & retained_c_umi_high]
dim(scRNA_filt)
table(scRNA_filt$group)
feats <- c("nFeature_RNA", "nCount_RNA")
p_filt_a_1=VlnPlot(scRNA_filt, features = feats, pt.size = 0, ncol = 2, 
                   group.by = "group") + NoLegend()
p_filt_a_2=VlnPlot(scRNA_filt, features = feats, pt.size = 0, ncol = 1, 
                   group.by = "orig.ident") + NoLegend()
feats <- c("percent_mito","percent_ribo")
p_filt_a_3=VlnPlot(scRNA_filt, features = feats, pt.size = 0, ncol = 2, 
                   group.by = "group") + NoLegend()
p_filt_a_4=VlnPlot(scRNA_filt, features = feats, pt.size = 0, ncol = 1, 
                   group.by = "orig.ident") + NoLegend()
(p_filt_a_1 / p_filt_a_2) | (p_filt_a_3 / p_filt_a_4)

saveRDS(scRNA_filt, file = "sce_filt.rds")

scRNA = readRDS("sce_filt.rds")
options(future.globals.maxSize = 1000 * 1024^40) 
library(future)
plan()
plan("multiprocess", workers = 4)
plan(multisession, workers = 4) 
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000) 
scRNA <- ScaleData(scRNA, features = VariableFeatures(scRNA), 
                   vars.to.regress = c("nFeature_RNA","percent_mito"))

scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA)) 
scRNA = RunUMAP(scRNA, dims = 1:15)
p_dim = DimPlot(scRNA)
p_dim = DimPlot(scRNA, group.by = "group")
scRNA <- FindNeighbors(scRNA, dims = 1:15) 
scRNA <- FindClusters(scRNA, resolution = c(0.01,0.05,0.1,0.2,0.5))
Idents(scRNA) = scRNA$RNA_snn_res.0.1
table(scRNA@active.ident)
p_dim_1 = DimPlot(scRNA, label = T)
p_dim | p_dim_1
# astrocytes: AQP4, ADGRV1, GPC5, RYR3
# endothelial cells: CLDN5, ABCB1, EBF1
# excitatory neurons: CAMK2A, CBLN2, LDB2
# inhibitory neurons: GAD1, LHFPL3, PCDH15
# microglia: C3, LRMDA, DOCK8
# oligodendrocytes: MBP, PLP1, ST18
Cell types were annotated using the following canonical marker genes: 
astrocytes (AQP4), microglia (LRMDA, DOCK8),
 excitatory neurons (CAMK2A, CBLN2,), inhibitory neurons (GAD1, LHFPL3,), 
oligodendrocytes (MBP, PLP1), and endothelial cells (CLDN5).

astrocytes = c("AQP4", "ADGRV1", "GPC5", "RYR3")
DotPlot(scRNA, features = astrocytes, 
        assay = "RNA") # 2

endothelial = c("CLDN5", "ABCB1", "EBF1")
DotPlot(scRNA, features = endothelial, 
        assay = "RNA") # 11

excitatory = c("CAMK2A", "CBLN2", "LDB2")
DotPlot(scRNA, features = excitatory, 
        assay = "RNA") # 1,3,7,9,10

inhibitory = c("GAD1", "LHFPL3", "PCDH15")
DotPlot(scRNA, features = inhibitory, 
        assay = "RNA") # 4,5,6

microglia = c("C3", "LRMDA", "DOCK8")
DotPlot(scRNA, features = microglia, 
        assay = "RNA") # 8

oligodendrocytes = c("MBP", "PLP1", "ST18")
DotPlot(scRNA, features = oligodendrocytes, 
        assay = "RNA") # 0

gene_list = list(
  Astro = astrocytes,
  Endo = endothelial,
  Excit = excitatory,
  Inhib = inhibitory,
  Mic = microglia,
  Oligo = oligodendrocytes
)

scRNA@active.ident = factor(scRNA@active.ident, levels = c(2,11,1,3,7,9,10,4,5,6,8,0))
p_dot_1 = DotPlot(scRNA, features = gene_list, 
                  assay = "RNA") +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.8, hjust=0.8))

p_dim_1 | p_dot_1

scRNA$celltype = ifelse(scRNA@active.ident %in% c(2), "Astro",
                        ifelse(scRNA@active.ident %in% c(11), "Endo",
                               ifelse(scRNA@active.ident %in% c(8), "Mic",
                                      ifelse(scRNA@active.ident %in% c(0), "Oligo",
                                             ifelse(scRNA@active.ident %in% c(4,5,6), "Inhib",
                                                    "Excit")))))
table(scRNA$celltype)
Idents(scRNA) = scRNA$celltype
scRNA@active.ident = factor(scRNA@active.ident, levels = c("Astro","Endo","Excit","Inhib","Mic","Oligo"))
table(scRNA@active.ident)
p_dot_2= DotPlot(scRNA, features = gene_list, 
                 assay = "RNA") +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.8, hjust=0.8))
p_dim_2= DimPlot(scRNA, label = T)
p_dim_2 | p_dot_2 
saveRDS(scRNA, file = "scRNAzhushi.rds")
setwd("~/singlecell/Lactylation/AD/GSE157827")
scRNA <- readRDS("scRNAzhushi.rds")

#--------------------------------------------
library(data.table)
Gene <- data.table::fread("lactate_unique_genes.txt")
gene1 <- Gene$Symbol
geneset=list(gene1) 
scRNA=AddModuleScore(scRNA,features = geneset,name = 'Add')

FeaturePlot(scRNA,features = 'Add1', label = T,reduction = 'umap')
VlnPlot(scRNA,features = 'Add1', pt.size = 0)
scRNA@meta.data$group_lactate <- ifelse(scRNA@meta.data$Add1> median(scRNA@meta.data$Add1),"score_UP","score_DOWN")
Idents(scRNA)=scRNA@meta.data$group_lactate
DimPlot(scRNA, reduction = 'umap', label = F, pt.size = 0.5) 

library(monocle)
imm_all_EP <- subset(scRNA,celltype=="Inhib")
data=as.matrix(imm_all_EP@assays$RNA@counts)
data <- as(data, 'sparseMatrix')
library(Biobase)
pdata <- new('AnnotatedDataFrame', data = imm_all_EP@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fdata <- new('AnnotatedDataFrame', data = fData)
mndata <- newCellDataSet(data,
                         phenoData = pdata,
                         featureData = fdata,
                         expressionFamily = negbinomial.size())

mndata <- estimateSizeFactors(mndata)
mndata <- estimateDispersions(mndata, cores=4, relative_expr = TRUE)

disp_table <- dispersionTable(mndata)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mndata <- setOrderingFilter(mndata, disp.genes)
plot_ordering_genes(mndata)
mndata <- reduceDimension(mndata, max_components = 2, method = 'DDRTree')#
mndata <- orderCells(mndata)#
p1 <- plot_cell_trajectory(mndata, color_by = "Pseudotime")  
p2 <- plot_cell_trajectory(mndata, color_by = "State")  
p3 <- plot_cell_trajectory(mndata, color_by = "celltype") 
p4 <- plot_pseudotime_heatmap(mndata[gene1,],show_rownames = TRUE,
                              return_heatmap = TRUE)
gene1 <- c("ARF1","MECP2","ACSF2","LDHA") 
p5 <- plot_genes_in_pseudotime(mndata[gene1,],color_by = "Pseudotime",ncol=2)#
p6 <- plot_genes_in_pseudotime(mndata[gene1,],color_by = "State",ncol=2)
p7 <- plot_genes_in_pseudotime(mndata[gene1,],color_by = "celltype",ncol=2)

setwd("~/singlecell/Lactylation/AD/GSE157827")
scRNA <- readRDS("scRNAzhushi.rds")

fig_1b=DimPlot(scRNA,  cols = c("#E97171","#FBD46D","#F7D6BF","#4F8A8B","#B0CAC7","skyblue"),raster=FALSE
)

sce_TC <- scRNA
clus_freq <-data.frame(clus=sce_TC@meta.data[["celltype"]],sname=sce_TC$orig.ident,
                       group=sce_TC$group_lactate)
head(clus_freq)

table(clus_freq$clus)
df <- as.data.frame(xtabs(~clus+sname+group,data=clus_freq))
colnames(df)
head(df)
table(sce_TC$orig.ident)
df$group<-factor(df$group,levels = c("score_DOWN","score_UP"))
cplt <- c("#E97171","#FBD46D","#F7D6BF","#4F8A8B","#B0CAC7","skyblue")
f1c <- ggplot(df, aes(fill=clus, y=Freq, x=group)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values = cplt[1:10])+xlab("Samples")+ylab("Percentage")+
  theme_classic()+theme(axis.text=element_text(size=10,color = "black"))+coord_flip()
f1c


df$group<-factor(df$group,levels = c("score_DOWN","score_UP"))
cplt1<-c("#E97171","#FBD46D","#F7D6BF","#4F8A8B","#B0CAC7","skyblue")
f11 <- ggplot(df, aes(fill=clus, y=group, x=Freq)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values = cplt1[1:11])+xlab("Percentage")+ylab("Samples")+
  theme_classic()+theme(axis.text=element_text(size=10,color = "black"))+
  coord_flip()+RotatedAxis()
f11

#--------------------------------------------------------------------
library(tidyverse)
library(CellChat)
library(ggplot2)
library(patchwork)
library(ggalluvial)
library(igraph)
library(dplyr)
library(SeuratData)
library(ggplot2)
library(svglite)
library(ggalluvial)
library(NMF)
library(ComplexHeatmap)
library(patchwork)
scRNA@meta.data$group_lactate <- ifelse(scRNA@meta.data$Add1> median(scRNA@meta.data$Add1),"score_UP","score_DOWN")
score_UP <- subset(scRNA, subset = group_lactate == "score_UP")
score_DOWN <- subset(scRNA, subset = group_lactate == "score_DOWN")
#  score_UP  score_DOWN cellchat objects
# score_UP  score_DOWN
meta.data1 <- score_UP@meta.data
meta.data1$group <- meta.data1$celltype
meta.data1$celltype1 <- meta.data1$celltype
head(meta.data1)
meta <- as.data.frame(score_UP@active.ident)
colnames(meta) <- "meta" 
unique(meta$meta)
head(meta)
cellchat1 <- createCellChat(object = score_UP, 
                            meta =meta.data1, 
                            group.by = "celltype1")
levels(cellchat1@idents)
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
str(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB,search = "Secreted Signaling")
cellchat1@DB <- CellChatDB.use
options(future.globals.maxSize = 10 * 1024^3)  
future::plan("multicore", workers = 10)
cellchat1 <- subsetData(cellchat1) 
cellchat1 <- identifyOverExpressedGenes(cellchat1)
cellchat1 <- identifyOverExpressedInteractions(cellchat1)
cellchat1 <- projectData(cellchat1, PPI.human) 
cellchat1 <- computeCommunProb(cellchat1, raw.use = TRUE)
cellchat1 <- filterCommunication(cellchat1, min.cells = 10)
cellchat1 <- computeCommunProbPathway(cellchat1)
cellchat1 <- aggregateNet(cellchat1)
groupSize <- as.numeric(table(cellchat1@idents))
class(cellchat1)
score_UP <- cellchat1
class(score_UP)
#-------------------------------------
meta.data1 <- score_DOWN@meta.data
meta.data1$group <- meta.data1$celltype
meta.data1$celltype1 <- meta.data1$celltype
head(meta.data1)
meta <- as.data.frame(score_DOWN@active.ident)
colnames(meta) <- "meta" 
unique(meta$meta)
head(meta)
cellchat1 <- createCellChat(object = score_DOWN, 
                            meta =meta.data1, 
                            group.by = "celltype1")
levels(cellchat1@idents)
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
str(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB,search = "Secreted Signaling")
cellchat1@DB <- CellChatDB.use
options(future.globals.maxSize = 10 * 1024^3)  
future::plan("multicore", workers = 10)
cellchat1 <- subsetData(cellchat1) 
cellchat1 <- identifyOverExpressedGenes(cellchat1)
cellchat1 <- identifyOverExpressedInteractions(cellchat1)
cellchat1 <- projectData(cellchat1, PPI.human)
cellchat1 <- computeCommunProb(cellchat1, raw.use = TRUE)
cellchat1 <- filterCommunication(cellchat1, min.cells = 10)
cellchat1 <- computeCommunProbPathway(cellchat1)
cellchat1 <- aggregateNet(cellchat1)
groupSize <- as.numeric(table(cellchat1@idents))
class(cellchat1)
score_DOWN <- cellchat1
#----------------------------------------
score_UP <- updateCellChat(score_UP) 
score_DOWN <- updateCellChat(score_DOWN)

levels(score_UP@ident)
levels(score_DOWN@active.ident)
identical(levels(score_UP@active.ident),levels(score_DOWN@active.ident))

object.list <- list(score_DOWN = score_DOWN, 
                    score_UP = score_UP)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat

p1<- compareInteractions(cellchat, show.legend = F, group = c(1,2))
p2<- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p1+ p2

cellchat <- identifyOverExpressedGenes(cellchat)

cellchat <- computeNetPairwiseAssociation(cellchat)

cellchat <- computeNetSimilarityPairwise(cellchat)
netDiffs <- cellchat@net$prob  
netDiffs <- cellchat@net$pval
dim(netDiffs)  
par(mfrow = c(1,2), xpd = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count")
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

p3<- netVisual_heatmap(cellchat)
p4<- netVisual_heatmap(cellchat, measure = "weight")
p3+ p4

p5<- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE) 
p6<- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE) 
p5+ p6

levels(cellchat@idents$joint) 
p7<- netVisual_bubble(cellchat, 
                      sources.use = c(1:7), 
                      targets.use = c(1:7), 
                      comparison= c(1, 2), 
                      angle.x = 45)

