rm(list = ls())

suppressMessages(library(clusterProfiler))
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(magrittr))
suppressMessages(library(gtools))
suppressMessages(library(stringr))
suppressMessages(library(Matrix))
suppressMessages(library(tidyverse))
suppressMessages(library(patchwork))
suppressMessages(library(data.table))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggpubr))
suppressMessages(library(scales))
suppressMessages(library(ggsci))
suppressMessages(library(sctransform))
suppressMessages(library(harmony))
suppressMessages(library(tidydr))
suppressMessages(library(celldex))
suppressMessages(library(pheatmap))
suppressMessages(library(clustree))
suppressMessages(library(xlsx))
suppressMessages(library(gridExtra))
suppressMessages(library(ggthemes))
suppressMessages(library(ggnewscale))
suppressMessages(library(CellChat))
suppressMessages(library(ggpubr))
suppressMessages(library(patchwork))
suppressMessages(library(monocle))
suppressMessages(library(clusterProfiler))
suppressMessages(library(enrichplot))
suppressMessages(library(ggrepel))
suppressMessages(library(DoubletFinder))
suppressMessages(library(future))

options(stringsAsFactors = F)

#01.#########
dir_name=list.files('00_origin_datas/GEO/GSE149655_RAW/')
datalist=list()
for (i in 1:length(dir_name)){
  dir.10x = paste0("00_origin_datas/GEO/GSE149655_RAW/",dir_name[i])
  list.files(dir.10x)
  my.data <- Read10X(data.dir = dir.10x)
  datalist[[i]] = CreateSeuratObject(counts = my.data, project = dir_name[i], min.cells = 3, min.features = 200)
  Samples1=dir_name[i]
  datalist[[i]] = AddMetaData(datalist[[i]] , Samples1,col.name = "Samples")
}
names(datalist)=dir_name
rm(my.data)


####
for (i in 1:length(datalist)){
  sce <- datalist[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")# 
  #sce[["percent.Ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")# 
  datalist[[i]] <- sce
  rm(sce)
}

sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
head(sce@meta.data)
sum(table(sce@meta.data$orig.ident))##
saveRDS(sce,"scRNA/scRNA.rds")

DefaultAssay(sce) <- "RNA"
temp_de <- readRDS("scRNA/decontX_results.rds")
identical(rownames(sce@meta.data),rownames(temp_de))
sce@meta.data$comtamination<-temp_de$comtamination
head(sce@meta.data)
scRNA = subset(sce, subset= comtamination<=0.2)
sum(table(scRNA@meta.data$orig.ident))##

scRNA1 = subset(scRNA, subset= nFeature_RNA<7000 & nCount_RNA <= 1000000 & percent.mt<40)

sum(table(scRNA1@meta.data$orig.ident))##

p1= VlnPlot(scRNA1, features=c("nFeature_RNA", "nCount_RNA","percent.mt"),  pt.size=0)+NoLegend()+theme(text=element_text(family="Times"))+theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('PDFs/Vlnplot.pdf',p1,height = 6,width = 12)
ggsave('PDFs/Vlnplot.jpg',p1,height = 6,width = 12)

DefaultAssay(scRNA1)
#sc <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
#sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 1000)
#sc <- ScaleData(sc, features = rownames(sc))
scRNA_SCT = SCTransform(scRNA1, vars.to.regress="percent.mt", verbose=FALSE)
saveRDS(scRNA_SCT,"scRNA/scRNA_SCT.rds")
#scRNA_SCT <- readRDS("01_scRNA/scRNA_SCT.rds")
###############################################################################
scRNA = RunPCA(scRNA_SCT, verbose=FALSE)

#DimPlot(scRNA , reduction = "pca" )##
#VizDimLoadings(scRNA , dims = 1:20, reduction = "pca")##
#DimHeatmap(scRNA , dims = 1:20, cells = 500, balanced = TRUE)
scRNA = RunHarmony(scRNA, group.by.vars="orig.ident", max.iter.harmony=50, lambda=0.5, assay.use="SCT")


scRNA_dim <- FindNeighbors(scRNA, dims=1:30, reduction="harmony")
scRNA_dim<- FindClusters(scRNA_dim  , resolution=0.3)##这里的resolution=0.1
scRNA_umap <- RunUMAP(scRNA_dim, dims=1:30, reduction="harmony")
#scRNA_tsne <- RunTSNE(scRNA_dim, dims=1:30, reduction="harmony")

mydata<-scRNA_umap
#mydata<-scRNA_tsne

#scRNA_tsne <- RunTSNE(scRNA, dims=1:15, reduction="harmony")
#scRNA_tsne <- FindNeighbors(scRNA_tsne, dims=1:15, reduction="harmony")
#scRNA_tsne <- FindClusters(scRNA_tsne  , resolution=0.1)##
#mydata<-scRNA_tsne

######################################  
p_umap = DimPlot(scRNA_umap , reduction="umap", group.by="orig.ident", pt.size=1)+theme(legend.position="right", plot.title=element_blank())
#ggsave('PDFs/p_umap.pdf',p_umap,height = 12,width = 12)

#p_tsne = DimPlot(scRNA_tsne , reduction="tsne", group.by="Samples", pt.size=1)+theme(legend.position="right", plot.title=element_blank())
#ggsave('PDFs/p_tsne.pdf',p_tsne,height = 12,width = 12)
#ggsave('PDFs/p_tsne.jpg',p_tsne,height = 12,width = 12)

#DimPlot(mydata , reduction="umap", label=T,group.by="seurat_clusters", pt.size=1)+theme(legend.position="right", plot.title=element_blank())
#DimPlot(scRNA_tsne , reduction="tsne", label=T,group.by="seurat_clusters", pt.size=1)+theme(legend.position="right", plot.title=element_blank())
FigureS1=mg_merge_plot(p1,p_umap,ncol=2,nrow=1,labels = c('A','B'),legend ="bottom",common.legend =T,widths = c(2.5,1))
ggsave('PDFs/Figure S1.pdf',FigureS1,height = 5,width = 12)
ggsave('PDFs/Figure S1.jpg',FigureS1,height = 5,width = 12)



colors=c("#E89DA0","#B2D3A4","#88CEE6","#F5D2A8","#B383B9","#FCED82","#FE88B1","#8BE0A4","#3C77AF","#D1352B","#8FA4AE")

plot_clusters=DimPlot(mydata, reduction="umap",label = T,label.size = 10,group.by="seurat_clusters", pt.size=1)+
  theme(legend.position="right", plot.title=element_blank())+
  theme(text=element_text(family="Times"))+
  #scale_color_manual(values=colors)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed"))
plot_clusters

markers <- FindAllMarkers(mydata, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
write.table(markers, "scRNA/All_celltype_DEG.txt", col.names=T, row.names=T, quote=F, sep="\t")
#markers <- read.table("scRNA/All_celltype_DEG.txt",header = T,check.names = F,fill=T,sep = "\t")
markers["pct.diff"]=markers$pct.1-markers$pct.2

#################################################################################
### 
Top5 <- markers %>% group_by(cluster) %>% slice_max(n =50, order_by =avg_logFC )
Top51 <- markers %>% group_by(cluster) %>% slice_max(n =50, order_by =pct.diff )
Cellmarker <- read.table)
Cellmarker <- Cellmarker[,c(9,7)]
colnames(Top51)[7]="marker"
colnames(Top5)[7]="marker"
Cellmarker2 <- merge(Cellmarker,Top51,by="marker")
write.table(Cellmarker2, "scRNA/Cellmarker2.txt", col.names=T, row.names=F, quote=F, sep="\t")

VlnPlot(mydata, features=c("MZB1",	"JCHAIN"), pt.size=0)+NoLegend()+theme(axis.title.x=element_blank())

plot_clusters
#0:Epithelial cell:"CLDN3",	"EPCAM"
#1:'T cells':"CD3E","CD3D"
#2:Macrophage:"ACP5",	"AIF1"
#3:Mast cell:"TPSAB1",	"MS4A2"
#4:'Endothelial cells' :"CLDN5","FCN3"
#5:'Fibroblast cells':"DCN",	"LUM"
#6:'Dendritic cells' :"LAMP3","MFSD2A"
#7:B cell:"CD79A",	"BANK1"
#8:Epithelial cell:"CLDN3",	"EPCAM"
#9:Plasma cell:"MZB1","JCHAIN"
#10:Plasma cell:"MZB1","JCHAIN"
#11:'Fibroblast cells':"DCN",	"LUM"




genes = unique(c(
  "CLDN3",	"EPCAM"
 ,"CD3E","CD3D"
 ,"ACP5",	"AIF1"
  ,"TPSAB1",	"MS4A2"
 , "CLDN5","FCN3"
 ,"DCN",	"LUM"
 ,"LAMP3","MFSD2A"
 ,"CD79A",	"BANK1"
 ,"CLDN3",	"EPCAM"
 ,"MZB1","JCHAIN"
  ,"MZB1","JCHAIN"
  ,"DCN",	"LUM"
  
))


cell_label = c(
  "Epithelial cell",
 'T cells',
  "Macrophage",
 "Mast cell",
 "Endothelial cells",
 'Fibroblast cells',
 'Dendritic cells' ,
 "B cell",
 "Epithelial cell",
  "Plasma cell",
 "Plasma cell",
  'Fibroblast cells')
#################################################################################
## 
names(cell_label) <- levels(mydata)
mydata1 <- RenameIdents(mydata, cell_label)
mydata1[["cell_type"]] = Idents(mydata1)

plot_cell=DimPlot(mydata1, pt.size=1, label=T, label.size=4)+
  theme( plot.title=element_blank())+
  scale_color_manual(values=colors)+
  theme(text=element_text(family="Times"))+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed"))


diff.marker.dotplot1= DotPlot(object = mydata1, features = unique(genes),
                              dot.scale =6,
                              dot.min = 0,
                              scale =T)+
  RotatedAxis()+ ggtitle("Marker Genes")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_text(size = 8)) +xlab('')+ylab('')+coord_flip()+scale_color_gradientn(colors=c( "snow", "red"))
diff.marker.dotplot1
cell.prop<-as.data.frame(prop.table(table(Idents(mydata1), as.vector(mydata1$orig.ident))))
colnames(cell.prop)<-c("Cell_type","Samples","proportion")
unique(cell.prop$Samples)

bili=ggplot(cell.prop,aes(Samples,proportion,fill=Cell_type))+
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values=colors)+
  ggtitle("")+
  theme_bw()+
  theme(text=element_text(family="Times"), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(face="bold",angle=45, size=12, hjust=1), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank(), legend.position="right")


hub=c("CD200R1", "BTN2A2" , "DNASE2B" ,"STAP1",   "SAMD9"   ,"BIRC3"   ,"SEMA7A" )


umap_sig2= DotPlot(object = mydata1, features =hub ,
                   dot.scale =6,
                   dot.min = 0,
                   scale =T)+
  RotatedAxis()+ ggtitle("Hub Genes")+
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab('')+ylab('')+coord_flip()+scale_color_gradientn(colors=c( "snow", "red"))

VlnPlot(mydata1, features=hub, pt.size=0)+NoLegend()+theme(axis.title.x=element_blank())


Figure1a=mg_merge_plot(plot_clusters,plot_cell,ncol=2,nrow=1,labels = c('A','B'))
Figure1d=mg_merge_plot(diff.marker.dotplot1,bili,umap_sig2,ncol=3,nrow=1, legend = "right",labels = c("C",'D',"E"))

Figure1=mg_merge_plot(Figure1a,Figure1d,ncol=1,nrow=2, legend = "right")
ggsave('PDFs/Fig6.pdf',Figure1,height = 12,width = 16)
ggsave('PDFs/Fig6.jpg',Figure1,height = 12,width = 16)

library(scran)
View(cc.genes)
g2m_genes<-cc.genes$g2m.genes
g2m_genes<-CaseMatch(search = g2m_genes,match = rownames(mydata1))
s_genes<-cc.genes$s.genes
s_genes<-CaseMatch(search = s_genes,match = rownames(mydata1))
mydata2<-CellCycleScoring(mydata1,g2m.features = g2m_genes,s.features = s_genes)
table(mydata2$Phase)
DimPlot(mydata2,group.by = 'Phase')
cell.prop1<-as.data.frame(prop.table(table(Idents(mydata2), as.vector(mydata2$Phase))))
colnames(cell.prop1)<-c("Cell_type","Phase","proportion")
bili2=ggplot(cell.prop1,aes(Cell_type,proportion,fill=Phase))+
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values=colors)+
  ggtitle("")+
  theme_bw()+
  theme(text=element_text(family="Times"), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(face="bold",angle=45, size=12, hjust=1), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank(), legend.position="right")

