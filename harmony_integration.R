library(Seurat)
library(harmony)

seurat_files <- c(
  "vol/RU1065C_MET_LI/data_RU1065C_MET_LI_cb_DF.rds",
  "vol/RU1144_MET_LN/data_RU1144_MET_LN_cb_DF.rds",
  "vol/RU1181_PRI_LU/data_RU1181_PRI_LU_cb_DF.rds",
  "vol/RU1293A_MET_LN/data_RU1293A_MET_LN_cb_DF.rds",
  "vol/RU779D_MET_LI/data_RU779D_MET_LI_cb_DF.rds",
  "vol/RU1066_PRI_LU/data_RU1066_PRI_LU_cb_DF.rds",
  "vol/RU1144_REC_LU/data_RU1144_REC_LU_cb_DF.rds",
  "vol/RU1195A_REC_LU/data_RU1195A_REC_LU_cb_DF.rds",
  "vol/RU1311_REC_LU/data_RU1311_REC_LU_cb_DF.rds",
  "vol/RU1080C_MET_KI/data_RU1080C_MET_KI_cb_DF.rds",
  "vol/RU1145_PRI_LU/data_RU1145_PRI_LU_cb_DF.rds",
  "vol/RU1215_MET_PL/data_RU1215_MET_PL_cb_DF.rds",
  "vol/RU1322_MET_LN/data_RU1322_MET_LN_cb_DF.rds",
  "vol/RU1108a_REC_LU/data_RU1108a_REC_LU_cb_DF.rds",
  "vol/RU1152_MET_LN/data_RU1152_MET_LN_cb_DF.rds",
  "vol/RU1229A_PRI_LU/data_RU1229A_PRI_LU_cb_DF.rds",
  "vol/RU325_MET_PL/data_RU325_MET_PL_cb_DF.rds",
  "vol/RU1124A_MET_LN/data_RU1124A_MET_LN_cb_DF.rds",
  "vol/RU1181_MET_LN/data_RU1181_MET_LN_cb_DF.rds",
  "vol/RU1231A_MET_LN/data_RU1231A_MET_LN_cb_DF.rds",
  "vol/RU426B_PRI_LU/data_RU426B_PRI_LU_cb_DF.rds"
)

# Load Seurat objects
seurat_objects <- lapply(seurat_files, function(file_path) {
  readRDS(file_path)
})

# Set the default assay to 'RNA'
seurat_objects <- lapply(seurat_objects, function(seu) {
  DefaultAssay(seu) <- 'RNA'
  return(seu)
})

# Run NormalizeData, FindVariableFeatures, and ScaleData on each Seurat object
seurat_objects <- lapply(seurat_objects, function(seu) {
  seu <- NormalizeData(seu) %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(npcs = 50)
  return(seu)
})

integrated_data <- RunHarmony(object = seu, group.by.vars = 'patient', dims.use = 1:30,
                  assay.use = 'RNA', plot_convergence = TRUE)

integrated_data <- RunUMAP(integrated_data, dims = 1:30)

# Plot the UMAP
pdf("UMAP_after_Harmony.pdf")
DimPlot(integrated_data, label = TRUE, group.by = "celltype_bped_main", raster = TRUE, shuffle = TRUE)

# Use DimPlot with group.by set to 'sample'
DimPlot(integrated_data, label = TRUE, group.by = "sample", raster = TRUE, shuffle = TRUE)

dev.off()
