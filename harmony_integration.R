library(Seurat)
library(harmony)


#seurat_files <- c(
#  "vol/RU1065C_MET_LI/data_RU1065C_MET_LI_cb_DF.rds",
#  "vol/RU1144_MET_LN/data_RU1144_MET_LN_cb_DF.rds",
#  "vol/RU1181_PRI_LU/data_RU1181_PRI_LU_cb_DF.rds",
#  "vol/RU1293A_MET_LN/data_RU1293A_MET_LN_cb_DF.rds",
#  "vol/RU779D_MET_LI/data_RU779D_MET_LI_cb_DF.rds",
#  "vol/RU1066_PRI_LU/data_RU1066_PRI_LU_cb_DF.rds",
#  "vol/RU1144_REC_LU/data_RU1144_REC_LU_cb_DF.rds",
#  "vol/RU1195A_REC_LU/data_RU1195A_REC_LU_cb_DF.rds",
#  "vol/RU1311_REC_LU/data_RU1311_REC_LU_cb_DF.rds",
#  "vol/RU1080C_MET_KI/data_RU1080C_MET_KI_cb_DF.rds",
#  "vol/RU1145_PRI_LU/data_RU1145_PRI_LU_cb_DF.rds",
#  "vol/RU1215_MET_PL/data_RU1215_MET_PL_cb_DF.rds",
#  "vol/RU1322_MET_LN/data_RU1322_MET_LN_cb_DF.rds",
#  "vol/RU1108a_REC_LU/data_RU1108a_REC_LU_cb_DF.rds",
#  "vol/RU1152_MET_LN/data_RU1152_MET_LN_cb_DF.rds",
#  "vol/RU1229A_PRI_LU/data_RU1229A_PRI_LU_cb_DF.rds",
#  "vol/RU325_MET_PL/data_RU325_MET_PL_cb_DF.rds",
#  "vol/RU1124A_MET_LN/data_RU1124A_MET_LN_cb_DF.rds",
#  "vol/RU1181_MET_LN/data_RU1181_MET_LN_cb_DF.rds",
#  "vol/RU1231A_MET_LN/data_RU1231A_MET_LN_cb_DF.rds",
#  "vol/RU426B_PRI_LU/data_RU426B_PRI_LU_cb_DF.rds"
#)

obj1 <- readRDS("vol/RU1065C_MET_LI/data_RU1065C_MET_LI_cb_DF.rds")
obj2 <- readRDS("vol/RU1144_MET_LN/data_RU1144_MET_LN_cb_DF.rds")
obj3 <- readRDS("vol/RU1181_PRI_LU/data_RU1181_PRI_LU_cb_DF.rds")
obj4 <- readRDS("vol/RU1293A_MET_LN/data_RU1293A_MET_LN_cb_DF.rds")
obj5 <- readRDS("vol/RU779D_MET_LI/data_RU779D_MET_LI_cb_DF.rds")
obj6 <- readRDS("vol/RU1066_PRI_LU/data_RU1066_PRI_LU_cb_DF.rds")
obj7 <- readRDS("vol/RU1144_REC_LU/data_RU1144_REC_LU_cb_DF.rds")
obj8 <- readRDS("vol/RU1195A_REC_LU/data_RU1195A_REC_LU_cb_DF.rds")
obj9 <- readRDS("vol/RU1311_REC_LU/data_RU1311_REC_LU_cb_DF.rds")
obj10 <- readRDS("vol/RU1080C_MET_KI/data_RU1080C_MET_KI_cb_DF.rds")
obj11 <- readRDS("vol/RU1145_PRI_LU/data_RU1145_PRI_LU_cb_DF.rds")
obj12 <- readRDS("vol/RU1215_MET_PL/data_RU1215_MET_PL_cb_DF.rds")
obj13 <- readRDS("vol/RU1322_MET_LN/data_RU1322_MET_LN_cb_DF.rds")
obj14 <- readRDS("vol/RU1108a_REC_LU/data_RU1108a_REC_LU_cb_DF.rds")
obj15 <- readRDS("vol/RU1152_MET_LN/data_RU1152_MET_LN_cb_DF.rds")
obj16 <- readRDS("vol/RU1229A_PRI_LU/data_RU1229A_PRI_LU_cb_DF.rds")
obj17 <- readRDS("vol/RU325_MET_PL/data_RU325_MET_PL_cb_DF.rds")
obj18 <- readRDS("vol/RU1124A_MET_LN/data_RU1124A_MET_LN_cb_DF.rds")
obj19 <- readRDS("vol/RU1181_MET_LN/data_RU1181_MET_LN_cb_DF.rds")
obj20 <- readRDS("vol/RU1231A_MET_LN/data_RU1231A_MET_LN_cb_DF.rds")
obj21 <- readRDS("vol/RU426B_PRI_LU/data_RU426B_PRI_LU_cb_DF.rds")

merge_seurat <- merge(obj1, y = c(obj2, obj3, obj4, obj5, obj6, obj7, obj8, obj9, obj10,
                                       obj11, obj12, obj13, obj14, obj15, obj16, obj17, obj18, obj19,
                                       obj20, obj21), add.cell.ids = c("obj1", "obj2", "obj3", "obj4", "obj5", "obj6", "obj7", 
                                                 "obj8", "obj9", "obj10", "obj11", "obj12", "obj13", "obj14", 
                                                 "obj15", "obj16", "obj17", "obj18", "obj19", "obj20", "obj21"), project = "PBMC15K")

# Load Seurat objects
#seu <- lapply(seurat_files, function(file_path) {
#  readRDS(file_path)
#})
# Load Seurat objects
#seu_list <- lapply(seurat_files, readRDS)

#merged_seurat <- do.call(merge, c(seu_list, merge.vars = "all"))

DefaultAssay(merge_seurat) <- 'RNA'

# Run other pre-processing steps
merge_seurat <- DietSeurat(merge_seurat, assays = 'RNA')
merge_seurat <- NormalizeData(merge_seurat) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs = 50)

# Run Harmony on the Seurat object
merge_seurat <- RunHarmony(object = merge_seurat, group.by.vars = 'sample', dims.use = 1:30,
                             assay.use = 'RNA', plot_convergence = TRUE)

# Run UMAP on the integrated data
merge_seurat <- RunUMAP(merge_seurat, reduction = 'harmony', dims = 1:30)

# Create a PDF file to save the plots
pdf("UMAP_after_Harmony_1.pdf")

# Plot UMAP with labels
DimPlot(merge_seurat, group.by = "celltype_bped_main", label = TRUE)

# Add more plots if needed, e.g., by sample
#DimPlot(merge_seurat, group.by = "sample", label = TRUE)

DimPlot(merge_seurat, group.by = "celltype_bped_fine", label = TRUE)

DimPlot(merge_seurat, group.by = "celltype_iced_main", label = TRUE)

DimPlot(merge_seurat, group.by = "celltype_iced_fine", label = TRUE)

# Close the PDF file
dev.off()




# Set the default assay to 'RNA'
#seurat_objects <- lapply(seurat_objects, function(seu) {
#  DefaultAssay(seu) <- 'RNA'
#  return(seu)
#})

# Run NormalizeData, FindVariableFeatures, and ScaleData on each Seurat object
#seurat_objects <- lapply(seurat_objects, function(seu) {
#  seu <- NormalizeData(seu) %>%
#    FindVariableFeatures() %>%
#    ScaleData() %>%
#    RunPCA(npcs = 50)
#  return(seu)
#})

#pca_list <- lapply(seurat_objects, function(seu) {
#  return(seu[["pca"]])
#})

# Combine the PCA results into a single matrix
#combined_pca <- do.call(cbind, pca_list)

# Run Harmony on the combined PCA matrix
#integrated_data <- RunHarmony(data = combined_pca, meta_data = seurat_objects, vars.use = 'sample', max_iter = 20)

# Add the Harmony results back to each Seurat object
#seurat_objects <- Map(function(seu, harmony_result) {
#  seu[["harmony"]] <- harmony_result
#  return(seu)
#}, seurat_objects, harmony_result$corrected)

#assay_name <- 'RNA'

# Run UMAP on the integrated data
#integrated_data <- RunUMAP(integrated_data, dims = 1:30, assay = assay_name)

# Plot the UMAP
# pdf("UMAP_after_Harmony.pdf")

# DimPlot(integrated_data, label = TRUE, group.by = "celltype_bped_main", raster = TRUE, shuffle = TRUE)

# Use DimPlot with group.by set to 'sample'
# DimPlot(integrated_data, label = TRUE, group.by = "sample", raster = TRUE, shuffle = TRUE)

#dev.off()
