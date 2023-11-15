# File: integrate_samples.R

# Load required libraries
library(Seurat)

# List of file paths for the 22 Seurat objects
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

# Read Seurat objects into a list
seurat_list <- lapply(seurat_files, readRDS)

# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:30)

# Integrate data
seurat_integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

# Finalize integration (additional processing steps)
seurat_integrated <- ScaleData(seurat_integrated)
seurat_integrated <- RunPCA(seurat_integrated, npcs = 30)
seurat_integrated <- FindNeighbors(seurat_integrated)
seurat_integrated <- FindClusters(seurat_integrated)

# Save the integrated Seurat object
saveRDS(seurat_integrated, "path/to/integrated_seurat_object.rds")
