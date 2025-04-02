library(decontX)
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

sce <- readRDS("scRNA/scRNA.rds")
sum(table(sce@meta.data$orig.ident))##
DefaultAssay(sce) <- "RNA"
sce@meta.data[1:4,]
#rownames(scRNA@meta.data)
# 
#counts <- sce@assays$RNA@counts
counts <- as.SingleCellExperiment(sce)
#counts1<-GetAssayData(scRNA,slot="counts")
#
decontX_results<-decontX::decontX(x=counts,
                                  #
                                  batch=counts$orig.ident,
                                  #z=dplyr::coalesce(pbmc@meta.data$seurat_annotations,"unknown"),
                                  seed=2024)
#
#zx=decontX_results$decontX_contamination
#saveRDS(zx,"01_scRNA/temp/decontX_results.rds")
#length(zx)
sce@meta.data$comtamination<-decontX_results$decontX_contamination
head(sce@meta.data)
zx=sce@meta.data
length(which(sce@meta.data$comtamination<=0.2))#4755
saveRDS(zx,"scRNA/decontX_results.rds")


