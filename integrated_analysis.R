file_path_1 <- "vol/integrated_seurat_object.rds"
integrated_data <- readRDS(file_path_1)
head(integrated_data)
summary(integrated_data)
names <- names(integrated_data@meta.data)
DimPlot(integrated_data, label = T,group.by = names ,raster = T,shuffle = T)