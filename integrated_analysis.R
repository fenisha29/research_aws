library(Seurat)
file_path_1 <- "vol/integrated_seurat_object.rds"
integrated_data <- readRDS(file_path_1)
head(integrated_data)
summary(integrated_data)
#names <- names(integrated_data@meta.data)

unique_cell_types <- unique(integrated_data$celltype_bped_main)


pdf("DimPlot_output.pdf")
# Use DimPlot with group.by set to 'celltype_bped_main' and labels from unique_cell_types
DimPlot(integrated_data, label = TRUE, group.by = "celltype_bped_main", raster = TRUE, shuffle = TRUE)
DimPlot(integrated_data, label = TRUE, group.by = "sample", raster = TRUE, shuffle = TRUE)
dev.off()
