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

# Perform integration using Harmony
integrated_data <- HarmonySeurat(
  seurat_objects,
  group.by.vars = c("sample"), # You might need to adjust this based on your metadata
  max.iter.harmony = 20
)

# Run PCA on the integrated data (optional)
integrated_data <- RunPCA(integrated_data, npcs = 30, verbose = FALSE)

# Run UMAP on the integrated data (optional)
integrated_data <- RunUMAP(integrated_data, dims = 1:30)

# Visualize integrated data with DimPlot (optional)
pdf("Integrated_harmony_DimPlot.pdf")
DimPlot(integrated_data, group.by = "sample", label = TRUE)
dev.off()

# Save the integrated Seurat object
saveRDS(integrated_data, "harmony_integrated_seurat_object.rds")
