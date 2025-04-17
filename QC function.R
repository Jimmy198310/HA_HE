cluster_resolution=1
seurat_standard_normalize_and_scale <- function(colon, cluster, cluster_resolution){
  # colon is seurat object, 
  colon <- NormalizeData(colon, normalization.method = "LogNormalize", scale.factor = 10000)
  colon <- FindVariableFeatures(colon, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(colon)
  colon <- ScaleData(colon, features = all.genes)
  colon <- RunPCA(colon, features = VariableFeatures(object = colon))
  if (cluster){
    colon <- FindNeighbors(colon, dims = 1:20)
    colon <- FindClusters(colon, resolution = cluster_resolution)
  }
  colon <- RunUMAP(colon, dims = 1:20)
  return(colon)
}
getwd()
make_seurat_object_and_doublet_removal <- function(data_directory, project_name){
  # function for basic seurat based qc and doubletfinder based doublet removal
  setwd("F:/SynologyDrive/Drive/高雌小鼠课题/HAvsHE/sc_data") 
  colon.data <- Read10X(data.dir = data_directory)
  currentSample <- CreateSeuratObject(counts = colon.data, project = project_name, min.cells = 3, min.features = 40)
  currentSample[["percent.mt"]] <- PercentageFeatureSet(currentSample, pattern = "^mt-")
  #去红细胞
  rownames(currentSample)[grep("^Hb[^(p)]",rownames(currentSample))]
  currentSample=PercentageFeatureSet(currentSample,"Hb[^(p)]",col.name="percent.hb")
  
  # qc plot-pre filtering
  setwd("F:/SynologyDrive/Drive/高雌小鼠课题/HAvsHE/sc_data/qcplot")
  pdf(paste0("./qc_plots_", project_name, "_prefiltered.pdf"))
  print(VlnPlot(currentSample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hb"), ncol = 3, pt.size = 0.05))
  dev.off()
  pdf(paste0("./qc_plots_", project_name, "_prefiltered_no_points.pdf"))
  print(VlnPlot(currentSample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hb"), ncol = 3, pt.size = 0))
  dev.off()
  
  # filter everything to 400 unique genes/cell
  currentSample=subset(currentSample, subset =  nFeature_RNA >300&nFeature_RNA<8000 & nCount_RNA >1000 & percent.mt<20 & percent.hb < 10)
  
  # Normalize and make UMAP
  currentSample <- seurat_standard_normalize_and_scale(currentSample, FALSE)
  
  # Run doublet finder
  nExp_poi <- round(0.08*length(currentSample@meta.data$orig.ident)*length(currentSample@meta.data$orig.ident)/10000)  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  seu_colon <- doubletFinder(currentSample, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = NULL, sct = FALSE)
  print(head(seu_colon@meta.data))
  
  # rename columns
  seu_colon$doublet.class <- seu_colon[[paste0("DF.classifications_0.25_0.09_",nExp_poi)]]
  seu_colon[[paste0("DF.classifications_0.25_0.09_",nExp_poi)]] <- NULL
  pann <- grep(pattern="^pANN", x=names(seu_colon@meta.data), value=TRUE)
  seu_colon$pANN <- seu_colon[[pann]]
  seu_colon[[pann]] <- NULL
  
  # plot pre and post doublet finder results
  pdf(paste0("./UMAP_pre_double_removal", project_name, ".pdf"))
  print(DimPlot(seu_colon, reduction = "umap", group.by = "doublet.class", cols = c("#D51F26", "#272E6A")))
  dev.off()
  seu_colon <- subset(seu_colon, subset = doublet.class != "Doublet")
  pdf(paste0("./UMAP_post_double_removal", project_name, ".pdf"))
  print(DimPlot(seu_colon, reduction = "umap", cols = c("#D51F26")))
  dev.off()
  
  # Remove extra stuff and return filtered Seurat object
  seu_colon <- DietSeurat(seu_colon, counts=TRUE, data=TRUE, scale.data=FALSE, assays="RNA")
  return(seu_colon)
}

seurat_qc_plots <- function(colon, sample_name){
  # Make some basic qc plots
  setwd("F:/SynologyDrive/Drive/高雌小鼠课题/HAvsHE/sc_data/qcplot")
  pdf(paste0("./seurat_nFeature_plots_", sample_name, ".pdf"), width = 40, height = 15)
  print(VlnPlot(colon, features = c("nFeature_RNA"), ncol = 1, pt.size = 0.2))
  dev.off()
  
  pdf(paste0("./seurat_nCount_plots_", sample_name, ".pdf"), width = 40, height = 15)
  print(VlnPlot(colon, features = c("nCount_RNA"), ncol = 1, pt.size = 0.2))
  dev.off()
  
  pdf(paste0("./seurat_pMT_plots_", sample_name, ".pdf"), width = 40, height = 15)
  print(VlnPlot(colon, features = c("percent.mt"), ncol = 1, pt.size = 0.2))
  dev.off()
  
  pdf(paste0("./seurat_pHb_plots_", sample_name, ".pdf"), width = 40, height = 15)
  print(VlnPlot(colon, features = c("percent.hb"), ncol = 1, pt.size = 0.2))
  dev.off()
}
