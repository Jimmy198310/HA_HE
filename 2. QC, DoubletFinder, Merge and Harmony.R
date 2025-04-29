Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

rm(list=ls())
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force = TRUE)

library(DoubletFinder)
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(viridis)
library(scCustomize)

getwd()        
data_directory=c("E18NCM/","E18NCF/","E18HEM/","E18HEF/","E18NAM/","E18NAF/","E18HAM/","E18HAF/")
project_name=c("E18_NCM","E18_NCF","E18_HEM","E18_HEF","E18_NAM","E18_NAF","E18_HAM","E18_HAF")
        
samples=project_name
        
sample1=make_seurat_object_and_doublet_removal(data_directory[1], samples[1])

###merge
seu_list=sample1
for (i in 2:8){
sc.i = make_seurat_object_and_doublet_removal(data_directory[i], samples[i])
seu_list=merge(seu_list,sc.i)
             }
        
table(seu_list$orig.ident)
##############################################################################
##############################################################################
##############################################################################
##############################################################################
###harmony
scRNA_harmony=seu_list
scRNA_harmony=NormalizeData(scRNA_harmony ) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
library(harmony)
scRNA_harmony=RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")

        
###Plot
scRNA_harmony=FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution =0.15)
scRNA_harmony=RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:20)
scRNA_harmony=RunTSNE(scRNA_harmony,reduction = "harmony",dims = 1:20)
DimPlot_scCustom(scRNA_harmony, reduction = "umap",label = T, raster=FALSE)
DimPlot_scCustom(scRNA_harmony, reduction = "tsne",label = T, raster=FALSE)
DimPlot(scRNA_harmony, reduction = "umap", split.by ='orig.ident', raster=FALSE,ncol = 8)
DimPlot_scCustom(scRNA_harmony, reduction = "umap", group.by='orig.ident',raster=FALSE)
table(scRNA_harmony$orig.ident)  

save(scRNA_harmony,file = "scRNA_harmony_0.15.rdata")

