#!/usr/bin/env Rscript

### Seurat analysis for cellbender output provided as an argument
## needs to run in /home/ubuntu

print(paste('Start:',Sys.time()))

library(dplyr)
library(Seurat)
library(gplots)
library(ggplot2)
library(ggrastr)
library(purrr)
library(cowplot)
library(DropletUtils)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(scater)
library(pheatmap)
library(Matrix)
library(DoubletFinder)

#pat <- commandArgs()[6]
pat <- "RU1065C_MET_LI"
doublet_rate <- 0.0644
print(pat)

##### Loading, merging, QC, dimension reduction #####
### Load dataset
system(paste0("aws s3 sync s3://sclc-seq/cellbender/v0.2.0_CellRanger6.1.1/", pat,"/ data/",pat,"/ ","--exclude '*' --include '",pat,"_filtered.h5' --quiet"))
seu.data <- Read10X_h5(paste0("data/",pat,"/",pat,'_filtered.h5'), use.names = TRUE, unique.features = TRUE)

### Initialize the Seurat object with the raw (non-normalized data)
seu_raw <- CreateSeuratObject(counts = seu.data, project = pat, min.cells = 1, min.features = 1)

# Annotate MT genes
seu_raw <- PercentageFeatureSet(seu_raw, pattern = "^MT-", col.name = "percent.mt")
seu_raw <- PercentageFeatureSet(seu_raw, pattern = "^RPS", col.name = "percent.rps")
seu_raw <- PercentageFeatureSet(seu_raw, pattern = "^RPL", col.name = "percent.rpl")
seu_raw$percent.rp=seu_raw$percent.rps + seu_raw$percent.rpl

# Annotate
seu_raw[["patient"]]<-pat
if(grepl('sclc',pat)){
  seu_raw[["ID"]]<-strsplit(pat,'_')[[1]][2]
  seu_raw[["sequencing"]]<-strsplit(pat,'_')[[1]][3]
  seu_raw[["cancer"]]<-'SCLC'
}else if(grepl('RU',pat)){
  seu_raw[["ID"]]<-strsplit(pat,'_')[[1]][1]
  seu_raw[["sequencing"]]<-'scRNA-seq'
  seu_raw[["cancer"]]<-'SCLC'
  seu_raw[["tissue"]]<-strsplit(pat,'_')[[1]][3]
  seu_raw[["primary"]]<-strsplit(pat,'_')[[1]][2]
  }else{
  seu_raw[["ID"]]<-strsplit(pat,'_')[[1]][1]
  seu_raw[["sequencing"]]<-strsplit(pat,'_')[[1]][2]
  seu_raw[["cancer"]]<-'Control'
}

# Identify doublets using scrublet
#doublet_rate_tmp<-doublet_rate[doublet_rate$sample==pat,2]
doublet_rate_tmp<- 0.0644
writeMM(seu_raw@assays$RNA@counts, paste0('data/',pat,'/matrix_',pat,'_raw.mtx'))
system(paste('python3 scrublet_code.py', pat, doublet_rate_tmp))
doublets <- read.table(paste0('data/',pat,'/doublets_',pat,'_raw.txt'),header = T)
#doublets <- read.table(paste0('data/',pat,'/matrix_RU1065C_MET_LI_raw.mtx'),header = T)
#doublets <- read.table(paste0('_raw.txt'),header = T)
seu_raw[['predicted_doublets']]<-doublets$predicted_doublets
seu_raw[['doublet_scores']]<-doublets$doublet_scores
system(paste0('rm data/',pat,'/matrix_',pat,'_raw.mtx'))
system(paste0('rm data/',pat,'/doublets_',pat,'_raw.txt'))

seu_raw <- NormalizeData(seu_raw)
seu_raw <- FindVariableFeatures(seu_raw)
seu_raw <- ScaleData(seu_raw)
seu_raw <- RunPCA(seu_raw)
seu_raw <- RunUMAP(seu_raw,dims=1:40)
seu_raw <- FindNeighbors(seu_raw, dims = 1:40)
seu_raw <- FindClusters(seu_raw,resolution=0.1)

sweep.list <- paramSweep_v3(seu_raw, PCs = 1:40, sct=FALSE)
sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK<-bcmvn %>% arrange(desc(BCmetric))
pK<-pK[1,2]
pK<-as.numeric(levels(pK[[1]]))[pK[[1]]]
nExp<-round(doublet_rate_tmp*dim(seu_raw@assays$RNA@counts)[2])
seu_raw <- doubletFinder_v3(seu_raw, PCs = 1:40,pK = pK,nExp = nExp)
seu_raw$doublet<-seu_raw@meta.data[,paste0('DF.classifications_0.25_',pK,'_',nExp)]
seu_raw$DF_score<-seu_raw@meta.data[,paste0('pANN_0.25_',pK,'_',nExp)]


### subset 
if(grepl('sclc',pat)){
  minFeature<-200
  maxFeature<- 12000
  minCount<- 800
  maxCount<- 70000
  maxMT<-15
}else{
  minFeature<-200
  maxFeature<- 5000
  minCount<- 800
  maxCount<- 15000
  maxMT<-15
}

seu <- subset(seu_raw, subset = nFeature_RNA > minFeature & 
                nFeature_RNA < maxFeature & nCount_RNA > minCount &
                nCount_RNA < maxCount & percent.mt < maxMT & 
                doublet =='Singlet')

### Workflow RNA
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu,npcs = 50)
seu <- RunUMAP(seu, dims = 1:40)
seu <- FindNeighbors(seu, dims = 1:40)
seu <- FindClusters(seu)

### cell type identification
seu_sce <- as.SingleCellExperiment(seu)

bped<-BlueprintEncodeData()
pred_bped_main <- SingleR(test = seu_sce, ref = bped, labels = bped$label.main)
pruneScores(pred_bped_main)
seu[['celltype_bped_main']]<-pred_bped_main$pruned.labels
pred_bped_fine <- SingleR(test = seu_sce, ref = bped, labels = bped$label.fine)
pruneScores(pred_bped_fine)
seu[['celltype_bped_fine']]<-pred_bped_fine$pruned.labels

iced<-DatabaseImmuneCellExpressionData()
pred_iced_main <- SingleR(test = seu_sce, ref = iced, labels = iced$label.main)
pruneScores(pred_iced_main)
seu[['celltype_iced_main']]<-pred_iced_main$pruned.labels
pred_iced_fine <- SingleR(test = seu_sce, ref = iced, labels = iced$label.fine)
pruneScores(pred_iced_fine)
seu[['celltype_iced_fine']]<-pred_iced_fine$pruned.labels

hpca<-HumanPrimaryCellAtlasData()
pred_hpca_main <- SingleR(test = seu_sce, ref = hpca, labels = hpca$label.main)
pruneScores(pred_hpca_main)
seu[['celltype_hpca_main']]<-pred_hpca_main$pruned.labels
pred_hpca_fine <- SingleR(test = seu_sce, ref = hpca, labels = hpca$label.fine)
pruneScores(pred_hpca_fine)
seu[['celltype_hpca_fine']]<-pred_hpca_fine$pruned.labels

mid<-MonacoImmuneData()
pred_mid_main <- SingleR(test = seu_sce, ref = mid, labels = mid$label.main)
pruneScores(pred_mid_main)
seu[['celltype_mid_main']]<-pred_mid_main$pruned.labels
pred_mid_fine <- SingleR(test = seu_sce, ref = mid, labels = mid$label.fine)
pruneScores(pred_mid_fine)
seu[['celltype_mid_fine']]<-pred_mid_fine$pruned.labels

### stats
stats<-as.data.frame(matrix(data=NA,nrow = 1, ncol = 12))
colnames(stats)<-c('sample','n_raw_features','n_raw_cells','n_predicted_doublets_scrublet','DF_doublets','n_features','n_cells','median_features','median_counts','cutoff_features','cutoff_counts','cutoff_mt')
rownames(stats)<-seu$patient[1]
stats$sample<-seu$patient[1]
stats$n_raw_features<-dim(seu_raw@assays$RNA@counts)[1]
stats$n_raw_cells<-dim(seu_raw@assays$RNA@counts)[2]
stats$n_predicted_doublets <-length(which(seu_raw$predicted_doublets ==T))
stats$DF_doublets <-length(which(seu_raw$doublet =='Doublet'))
stats$n_features<-dim(seu@assays$RNA@counts)[1]
stats$n_cells<-dim(seu@assays$RNA@counts)[2]
stats$median_features<-round(median(seu@meta.data$nFeature_RNA))
stats$median_counts<-round(median(seu@meta.data$nCount_RNA))
stats$cutoff_features<-paste(minFeature,maxFeature)
stats$cutoff_counts<-paste(minCount,maxCount)
stats$cutoff_mt<-paste(maxMT)

### Save objects
ifelse(!dir.exists(file.path(paste0("data/",seu$patient[1],'/'))), dir.create(file.path(paste0("data/",seu$patient[1],'/'))), FALSE)
saveRDS(seu, file = paste0("data/",seu$patient[1],'/data_',seu$patient[1],'_cb_DF.rds'))

### write pdf reports
pdf(file = paste0("data/",seu$patient[1],"/plots_", seu$patient[1],"_cb_DF.pdf"))

# stats
textplot(t(stats),cex=1.2,halign='left')

# plots raw data
ggplot(seu_raw@meta.data, aes(x=seu_raw$nCount_RNA,y = seu_raw$nFeature_RNA, col=seu_raw$percent.mt)) + 
  rasterise(geom_point(size=0.5,alpha=0.5),dpi=300)+ scale_colour_gradient(low="blue", high="green") + 
  labs(color = "Percent MT") + theme_classic() + ggtitle('Raw object')

DimPlot(seu_raw, reduction = "umap",label = T,
        group.by = paste0('DF.classifications_0.25_',pK,'_',nExp),
        raster = T,shuffle = T)
FeaturePlot(seu_raw, features = paste0('pANN_0.25_',pK,'_',nExp),order=T, raster = T)
VlnPlot(seu_raw,features = paste0('pANN_0.25_',pK,'_',nExp),
        group.by = paste0('DF.classifications_0.25_',pK,'_',nExp),pt.size = 0)
FeaturePlot(seu_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.rps',
                                  'percent.rpl',paste0('pANN_0.25_',pK,'_',nExp)), 
            order=T, raster = T)
FeatureScatter(seu_raw,feature1 = paste0('pANN_0.25_',pK,'_',nExp),feature2 = 'nFeature_RNA',
               shuffle = T,raster = T,group.by = paste0('DF.classifications_0.25_',pK,'_',nExp))

# QC plot filtered
ggplot(seu@meta.data, aes(x=seu$nCount_RNA,y = seu$nFeature_RNA, col=seu$percent.mt)) + 
  rasterise(geom_point(size=0.5,alpha=0.5),dpi=300)+ scale_colour_gradient(low="blue", high="green") + 
  labs(color = "Percent MT") + theme_classic()+ ggtitle('Filtered object')

ggplot(seu@meta.data, aes(x=seu$nCount_RNA,y = seu$nFeature_RNA, col=seu$DF_score)) + 
  rasterise(geom_point(size=0.5,alpha=0.5),dpi=300)+ scale_colour_gradient(low="lightgrey", high="darkviolet") + 
  labs(color = "DF_score") + theme_classic()+ ggtitle('Filtered object')

print(VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.rps','percent.rpl','DF_score'),
              ncol = 3,group.by = 'patient',pt.size = 0))

FeaturePlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.rps','percent.rpl','DF_score'), 
            min.cutoff = "q05", max.cutoff = 'q95',order=T, raster = T)

# PCA
print(ElbowPlot(seu))
DimPlot(seu, reduction = "pca",group.by = 'ident',raster = T,shuffle = T)

# UMAP
DimPlot(seu, reduction = "umap",label = T,group.by = 'ident',raster = T,shuffle = T)

# features
FeaturePlot(seu, features = c('ASCL1',"NEUROD1", 'POU2F3', "YAP1"), min.cutoff = "q05",
            max.cutoff = 'q95',order=T, raster = T)

FeaturePlot(seu, features = c('PLCG2',"FGFR1", 'SOX2', "REST",'ENO2','DDC','INSM1','MYC'), 
            min.cutoff = "q05",max.cutoff = 'q95',order=T, raster = T)

FeaturePlot(seu, features = c('HES6','HES1','MYCL1','DLL3','BCL2','NFIB','NCOR2'), 
            min.cutoff = "q05",max.cutoff = 'q95',order=T, raster = T)

FeaturePlot(seu, features = c('PTPRC',"CD8A", 'CD68', "MKI67"), min.cutoff = "q05",
            max.cutoff = 'q95',order=T, raster = T)

## bped
plotScoreHeatmap(pred_bped_main, clusters=seu_sce@colData$ident,fontsize = 6,main='pred_bped_main')
DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_bped_main',repel = T,label.size = 2.5,raster = T,shuffle = T) + 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=6))

## hpca
plotScoreHeatmap(pred_hpca_main, clusters=seu_sce@colData$ident,fontsize = 6,main='pred_hpca_main')
DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_hpca_main',repel = T,label.size = 2.5,raster = T,shuffle = T) + 
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=6))

dev.off()

system(paste0("aws s3 sync data/",seu$patient[1],"/ s3://sclc-seq/Seurat/",seu$patient[1],"/ --exclude '*' --include '*_cb*' --exclude '.*' --quiet"))

print(paste('End:',Sys.time()))
