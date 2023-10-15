#!/usr/bin/env Rscript

### Seurat QC plots

print(paste('Start:',Sys.time()))

library(dplyr)
library(Seurat)
library(gplots)
library(ggplot2)
library(ggrastr)
library(purrr)
library(cowplot)
library(scater)
library(patchwork)

#patients<-c('C51ctr_3p', 'C52ctr_3p', 'C53ctr_3p', 'C54ctr_3p', 'C55ctr_3p', 'C56ctr_3p', 'C57ctr_3p',
#            'sclc_n3525_3p', 'sclc_n3626_3p', 'sclc_n3628_3p', 'sclc_n3657_3p', 'sclc_ns005_3p', 'sclc_ns007_3p',
#            'sclc_pa013_5pv2', 'sclc_pa031_5pv2', 'sclc_pa096_3p', 'sclc_pa098_3p', 'sclc_pa121_5pv2', 
#            'sclc_pa124_3p', 'sclc_pa124_5pv2', 'sclc_pa146_3p', 'sclc_pa150_3p', 'sclc_pa151_3p', 'sclc_pa168_3p')

patients <- "RU1065C_MET_LI"

ifelse(!dir.exists(file.path("data/SCLC_Ctr")), dir.create(file.path("data/SCLC_Ctr")), FALSE)
pdf(file = paste0("data/SCLC_Ctr/plots_SCLC_Ctr_QC.pdf"),width = 15,height = 10)
for(pat in patients){
  seu<-readRDS(paste0("data/",pat,'/data_',pat,'_cb.rds'))
  
  ### stats
  stats<-as.data.frame(matrix(data=NA,nrow = 1, ncol = 8))
  colnames(stats)<-c('sample','n_features','n_cells','median_features','median_counts')
  rownames(stats)<-seu$patient[1]
  stats$sample<-seu$patient[1]
  stats$n_features<-dim(seu@assays$RNA@counts)[1]
  stats$n_cells<-dim(seu@assays$RNA@counts)[2]
  stats$median_features<-round(median(seu@meta.data$nFeature_RNA))
  stats$median_counts<-round(median(seu@meta.data$nCount_RNA))

  #p1<-textplot(t(stats),cex=1.2,halign='left')
  
  p2<-DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_bped_fine',repel = T,label.size = 2.5,raster = T,shuffle = T) + 
    guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
    theme(legend.text=element_text(size=6))+ ggtitle(pat)+
    theme(legend.position="none")
  p3<-DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_hpca_main',repel = T,label.size = 2.5,raster = T,shuffle = T) + 
    guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
    theme(legend.text=element_text(size=6))+ 
    theme(legend.position="none")
  
  p4<-DimPlot(seu, reduction = "pca",group.by = 'ident',raster = T,shuffle = T)+ 
    theme(legend.position="none")
  p5<-FeaturePlot(seu, features = c('ASCL1',"NEUROD1", 'POU2F3', "YAP1"), min.cutoff = "q05",
              max.cutoff = 'q95',order=T, raster = T)
  p10<-FeaturePlot(seu, features = c('PTPRC',"CD8A", 'CD68', "MKI67"), min.cutoff = "q05",
              max.cutoff = 'q95',order=T, raster = T)
  
  #print(p2+p3+p4+p5+p6+p10)
  print((p2 | p3 | p4) / (p5 | p10))
}
dev.off()


print(paste('End:',Sys.time()))
