rm(list = ls())


suppressMessages(library(data.table))
suppressMessages(library(devtools))
suppressMessages(library(customLayout))
suppressMessages(library(stringr))
suppressMessages(library(ConsensusClusterPlus))
suppressMessages(library(tidydr))
suppressMessages(library(openxlsx))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tidyverse))
suppressMessages(library(clusterProfiler))
suppressMessages(library(pheatmap))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(GSVA))
suppressMessages(library(GSEABase))
suppressMessages(library(fgsea))
suppressMessages(library(corrplot))
suppressMessages(library(colorspace))
suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(maftools))
suppressMessages(library(vegan))
suppressMessages(library(forcats))
suppressMessages(library(ggpubr))
suppressMessages(library(ggplot2))
suppressMessages(library(rstatix))
suppressMessages(library(ggstatsplot))
suppressMessages(library(ggcor))
suppressMessages(library(ggstance))
suppressMessages(library(tidyverse))
suppressMessages(library(GOplot))
suppressMessages(library(caret))
suppressMessages(library(writexl))
suppressMessages(library(rcartocolor))
suppressMessages(library(ggcorrplot))
suppressMessages(library(psych))
suppressMessages(library(clusterProfiler))
suppressMessages(library(dplyr))
suppressMessages(library(cols4all))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(scales))
suppressMessages(library(oncoPredict))
suppressMessages(library(gghalves))
suppressMessages(library(cowplot))
suppressMessages(library(IOBR))
suppressMessages(library(estimate))
suppressMessages(library(UpSetR))
suppressMessages(library(ggbiplot))
suppressMessages(library(ggsci))
suppressMessages(library(WGCNA))
suppressMessages(library(circlize))
suppressMessages(library(rJava))
suppressMessages(library(xlsxjars))
suppressMessages(library(xlsx))
suppressMessages(library(glmnet))
suppressMessages(library(tidyr))
suppressMessages(library(pROC))
suppressMessages(library(ROCR))

tcga.exp=readMatrix('00_origin_datas/Preprocessed/tcga.exp1.txt')
tcga.t.exp=readMatrix('00_origin_datas/Preprocessed/tcga.t.exp1.txt')
tcga.cli=readMatrix('00_origin_datas/Preprocessed/tcga.cli.txt')

tcga.t.cli=tcga.cli[colnames(tcga.t.exp),]
identical(as.vector(tcga.cli$sampleID),as.vector(colnames(tcga.exp)))
geneSets <- read.table("00_origin_datas/efferocytosis related genes.txt",header = T,check.names = F,fill=T,sep = "\t")
geneSets <- geneSets$`efferocytosis related genes` #signature read
length(as.vector(unique(geneSets)))#167

tcga_exp=as.data.frame(tcga.exp)

ssgsea_list <- list()
ssgsea_list[['ssGSEA']] <- as.vector(unique(geneSets))
ssgsea_score <- ssGSEAScore_by_muti_group_genes(gene.exp = tcga_exp,
                                                genelist = ssgsea_list)
ssgsea_score <- as.data.frame(t(ssgsea_score))
ssgsea_score$"Group"=""
######### -01 
identical(rownames(ssgsea_score),as.vector(tcga.cli$sampleID))

inds1=which(substr(colnames(tcga.exp),14,15)=="01")
inds2=which(substr(colnames(tcga.exp),14,15)=="11")

ssgsea_score$"Group"[inds1]="Case"
ssgsea_score$"Group"[inds2]="Control"
Fig1A=mg_violin(ssgsea_score[, c("Group", "ssGSEA")]
                ,melt = T
                ,xlab = ''
                ,legend.pos = 'tl'
                ,ylab = 'ssGSEA score')
ggsave('01_WGCNA/Fig1A.pdf',Fig1A,height = 6,width = 6.5)

dev.off()
###############


#############
cox.pval=0.05
exp_m6A=as.matrix(tcga.t.exp[which(rownames(tcga.t.exp)%in%as.vector(geneSets)),])
tcga.t.cli$OS.time=as.numeric(tcga.t.cli$OS.time)
tcga.t.cli$OS=as.numeric(tcga.t.cli$OS)
identical(as.vector(tcga.t.cli$sampleID),colnames(exp_m6A))
tcga.pcd.cox=cox_batch(t(scale(t(as.matrix(exp_m6A))))
                       ,time = tcga.t.cli$OS.time/365
                       ,event = tcga.t.cli$OS)

table(tcga.pcd.cox$p.value<0.05)
table(tcga.pcd.cox$p.value<0.01)
table(tcga.pcd.cox$p.value<0.001)

tcga.pcd.cox=tcga.pcd.cox[order(tcga.pcd.cox$HR,decreasing = T),]

tcga.pcd.cox.sig=tcga.pcd.cox[which(tcga.pcd.cox$p.value<cox.pval),]
nrow(tcga.pcd.cox.sig)
pdf('PDFs/bioForest.pdf',height = 9,width =5,onefile = F)
bioForest(rt = tcga.pcd.cox.sig,col=c('#5C8980','#A5604A'))
dev.off()

write.table(tcga.pcd.cox.sig,'ConsensusClusterPlus/tcga.pcd.cox.sig.txt',sep = "\t",quote = F,row.names = T,col.names = T)

###########
clusterAlg_name=c('hc','pam','km','kmdist')[3]
distance_name=c('pearson','spearman','euclidean','binary','maximum','canberra','minkowski')[1]
tcga_consen_data=as.matrix(tcga.t.exp[which(rownames(tcga.t.exp)%in%(rownames(tcga.pcd.cox.sig))),])
tcga_consen_data=t(scale(t(tcga_consen_data),scale = F))   

#tcga_consen_data=sweep(tcga_consen_data,1,apply(tcga_consen_data, 1, median))
#tcga_consen_data=as.dist(1-cor(tcga_consen_data,method = 'spearman'))

tcga_clust_subtype <- ConsensusClusterPlus(tcga_consen_data, maxK = 10, reps = 500, pItem = 0.8, pFeature = 1, title = "ConsensusClusterPlus", clusterAlg = clusterAlg_name, distance = distance_name, plot = "pdf", writeTable = F, seed = 123456)#########
save(tcga_clust_subtype,file='ConsensusClusterPlus/tcga.subtype.RData')
load('ConsensusClusterPlus/tcga.subtype.RData')#########
k=2
colors = c("#66C5CC", "#DCB0F2", "#D3B484", "#87C55F", "#DCB0F2", "#87C55F", "#D3B484", "#F6CF71", "#9EB9F3", "#F89C74", "#66C5CC")

subtype.cols=c("#B497E7","#F89C74")
tcga.subtype <- data.frame( Samples=names(tcga_clust_subtype[[k]]$consensusClass),Subtype=tcga_clust_subtype[[k]]$consensusClass)
tcga.subtype$Subtype=paste0('C',tcga.subtype$Subtype)
write.table(tcga.subtype,file = 'ConsensusClusterPlus/Subtype.txt',sep = '\t',quote = F,row.names = F,col.names = T)


table(tcga.subtype$Subtype)
colnames(tcga.t.cli)[1]='Samples'
tcga.subtype.cli=merge(tcga.subtype,tcga.t.cli,by='Samples')




fig1f=ggplotKMCox(data.frame(time = tcga.subtype.cli$OS.time/365
                             , event = tcga.subtype.cli$OS
                             , tcga.subtype.cli$Subtype) 
                  ,add_text = '',show_confint = F,palette = subtype.cols)


fig1f
ggsave("ConsensusClusterPlus/fig1f.pdf", width = 8, height = 8)

#################PCA
tcga.subtype$Samples=as.vector(tcga.subtype$Samples)
tcga_exp_var=t(tcga.t.exp[,tcga.subtype$Samples])

tcga_exp_var=tcga_exp_var[ , which(apply(tcga_exp_var, 2, var) != 0)]##
dim(tcga_exp_var)
cluster.pca <- prcomp(tcga_exp_var, scale=T)
cluster.pca.plot <- ggbiplot(cluster.pca, scale=1, groups = tcga.subtype$Subtype,
                             ellipse = TRUE,ellipse.prob=0.3, circle = F,var.axes=F) +
  scale_color_manual(values = subtype.cols) + 
  theme_bw() +
  theme(legend.direction = 'horizontal', legend.position = 'top',
        panel.grid = element_blank(),text = element_text(family = 'Times')) +
  xlab('PCA1') + ylab('PCA2')+xlim(-3,3)+ylim(-3,3)
cluster.pca.plot
ggsave("ConsensusClusterPlus/cluster.pca.plot.pdf", cluster.pca.plot,width = 6, height = 6)

#################### 
colnames(tcga.subtype.cli)[c(4:8,13)]
tcga.subtype.cli[1:4,]
tcga.subtype.cli.cmp=list()
for(i in c(4:8,13)){
  #group.color=subtype.cols[1:4]
  
  p=plotMutiBar_tmp(table(tcga.subtype.cli[,i],tcga.subtype.cli$Subtype)
                    ,fill.color = colors
                    ,isAuto = F,showValue = F
                    ,legTitle=colnames(tcga.subtype.cli)[i])
  tcga.subtype.cli.cmp=c(tcga.subtype.cli.cmp,list(p$Bar))
}
length(tcga.subtype.cli.cmp)

fig1e=mg_merge_plot(tcga.subtype.cli.cmp,nrow = 1,ncol = length(c(4:8,13))
                    ,labels='F')
fig1e
ggsave("ConsensusClusterPlus/fig1e2.pdf", fig1e,width = 8, height = 6)

Fig1=ggpubr::ggarrange(fig1f,cluster.pca.plot,fig1e, ncol = 2, nrow = 2,labels=LETTERS[1:4])

ggsave('PDFs/Fig1.pdf',Fig1,height = 12,width = 15)


####WGCNA
######## WGCNA ssGSEA########
tcga.t.exp1=readMatrix('00_origin_datas/Preprocessed/tcga.t.exp1.txt')

tcga_wgcna=tcga.t.exp1
allowWGCNAThreads(nThreads = 36)#
enableWGCNAThreads(nThreads = 36)# 

tcga.mads=apply(tcga_wgcna, 1, mad)
tpm_T2=tcga_wgcna[which(tcga.mads>quantile(tcga.mads, probs=seq(0, 1, 0.25))[1]),]
dim(tpm_T2)
tpm_T2=t(tpm_T2)
dim(tpm_T2)
range(tpm_T2)
########
#sampleTree = hclust(dist(tpm_T2), method ="average");
#pdf('03_WGCNA/hust.pdf',width = 10,height = 5)

#plot(sampleTree, main ="Sample clustering to detect outliers",sub="", xlab="", cex = 0.8,cex.lab = 1.5, cex.axis= 1.5, cex.main = 1.5)
#abline(h = 190,col="red"); # 
#dev.off()
#clust = cutreeStatic(sampleTree, cutHeight = 275, minSize = 2)
#table(clust)	# 
#keepSamples = (clust==1)	# 
#tpm_T2 = tpm_T2[keepSamples, ]		
########
##########
pdf('01_WGCNA/FigABC.pdf',width = 8,height = 8)
tpm_T2.power=mg_wgcna_get_power(tpm_T2)
dev.off()

tpm_T2.power$cutPower
tpm_T2.module=mg_WGCNA_getModule(tpm_T2,
                                 power = tpm_T2.power$cutPower,
                                 deepSplit=2,
                                 mergeCutHeight=0.25,
                                 minModuleSize=100)

table(tpm_T2.module$Modules[,2])
length(table(tpm_T2.module$Modules[,2]))

pdf('01_WGCNA/FigD.pdf',height = 5,width = 6)
plotDendroAndColors(tpm_T2.module$Tree, tpm_T2.module$Modules,
                    c("Dynamic Module",'Merged Module'),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

zx=data.frame(tpm_T2.module$Modules)
write.table(zx,'01_WGCNA/tcga.wgcna.module.genes.txt',sep = "\t", row.names = TRUE,quote = FALSE,col.names = T)

fig2e=mg_barplot_point(labels = names(table(tpm_T2.module$Modules[,2]))
                       ,values = as.numeric(table(tpm_T2.module$Modules[,2]))
                       ,point_sizes = 2
                       ,point_cols = names(table(tpm_T2.module$Modules[,2]))
                       ,xlab = 'Number of Genes',legend.pos = NULL)
fig2e
dev.off()
ggsave('01_WGCNA/Fig2e.pdf',fig2e,height = 5,width = 3)

#### 

#### 
# Calculate eigengenes
MEs = tpm_T2.module$MEs
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
pdf('01_WGCNA/Fig2f.pdf',height = 6,width = 12,onefile = T)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
dev.off()


# ########## 
#ssgsea_score$"Class"=ifelse(ssgsea_score$Group=='Metastatic','1','0')
#ssgsea_score$"Class"=as.numeric(ssgsea_score$"Class")
#ssgsea_score1=ssgsea_score[colnames(tcga_wgcna),c(1,3)]

ssgsea_score$"Samples"=rownames(ssgsea_score)
ssgsea_score1=merge(ssgsea_score,tcga.subtype,by="Samples")
rownames(ssgsea_score1)=ssgsea_score1$Samples
ssgsea_score1$"Class"=ifelse(ssgsea_score1$Subtype=='C1','1','0')
ssgsea_score1=ssgsea_score1[rownames(tpm_T2),c(2,5)]

identical(rownames(tpm_T2),rownames(ssgsea_score1))
spms <-as.data.frame(ssgsea_score1[,c("ssGSEA","Class")])
colnames(spms)=c("ssGSEA_score","Subtype")
rownames(spms)=rownames(ssgsea_score1)
head(spms)

MEs_col<-tpm_T2.module$MEs
dim(MEs_col)

modTraitCor = cor(MEs_col[,rownames(MEDiss)[METree$order]]
                  , spms
                  ,use = 'p')
modTraitP = corPvalueStudent(modTraitCor, dim(spms)[1])

textMatrix = paste(signif(modTraitCor, 2), "\n", " (", format(modTraitP,scientific =TRUE,digits = 3), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
dim(textMatrix)


pdf('01_WGCNA/Figf-2.pdf',width =6,height = 6)
labeledHeatmap(Matrix = data.frame(modTraitCor),
               xLabels = colnames(modTraitCor),
               yLabels = rownames(modTraitCor),
               cex.lab = 1,
               ySymbols = rownames(modTraitCor), 
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = data.frame(textMatrix),
               setStdMargins = FALSE,xLabelsAngle = 0,
               cex.text = 0.8, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

rownames(modTraitCor)=gsub("ME","",rownames(modTraitCor))
rownames(textMatrix)=gsub("ME","",rownames(textMatrix))
colnames(modTraitCor)


#
geneModuleMembership <- signedKME(tpm_T2
                                  , data.frame(tpm_T2.module$MEs)
                                  , outputColumnName = "")
head(geneModuleMembership)
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership)
                                           , nrow(tpm_T2.module$MEs)))
#
geneTraitSignificance <- as.data.frame(cor(tpm_T2
                                           , spms
                                           , use = 'pairwise.complete.obs'))

GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance)
                                           , nrow(spms)))


table(tpm_T2.module$Modules[,2])
# 
module = c("magenta")
#module = c("green")

pheno = c("ssGSEA_score")
modNames<-colnames(geneModuleMembership)

# 
module_column = match(module, modNames)
pheno_column = match(pheno, colnames(spms))

# 
gene_M=as.data.frame(tpm_T2.module$Modules)
identical(rownames(geneModuleMembership),rownames(gene_M))
moduleGenes <-  which(gene_M[,2]%in%module)


zx1=as.data.frame(geneModuleMembership[moduleGenes, module_column])
rownames(zx1)=rownames(geneModuleMembership)[moduleGenes]###MM
write.table(rownames(zx1),'02_DEGs/wgcna_genes.txt',sep = "\t",quote = F,row.names = F,col.names = F)

##############DEG##########
C1_sample=as.vector(tcga.subtype.cli$Samples[which(tcga.subtype.cli$Subtype=="C1")]) 
C2_sample=as.vector(tcga.subtype.cli$Samples[which(tcga.subtype.cli$Subtype=="C2")]) 

geo.limma_C1vsC2=mg_limma_DEG(exp = tcga.t.exp[,c(C1_sample,C2_sample)],group = c(rep("C1",length(C1_sample)),rep("C2",length(C2_sample))), ulab = 'C1',dlab = 'C2')
geo.limma_C1vsC2$Summary
df.deg.sig=geo.limma_C1vsC2$DEG[which(geo.limma_C1vsC2$DEG$adj.P.Val<0.05 & abs(geo.limma_C1vsC2$DEG$logFC)>log2(1.5)),]
write.table(df.deg.sig,'02_DEGs/limma_C1vsC2.txt',sep = "\t",quote = F,row.names = T,col.names = T)
dim(df.deg.sig)
############
fig3a=my_volcano_FDR(geo.limma_C1vsC2,p_cutoff = 0.05,fc_cutoff = log2(1.5),col = c('#BC3C29FF','#20854EFF','grey'))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        text = element_text(color = "black",family = 'Times',size = 14),
        axis.text = element_text(color = "black",family = 'Times',size = 14),
        legend.position = 'top')#+xlim(-2,2)+ylim(0,20)
fig3a
ggsave("02_DEGs/fig3a.pdf", width = 6, height = 6)

tcga.deg.sig=intersect(rownames(df.deg.sig),rownames(zx1))
hub=df.deg.sig[which(rownames(df.deg.sig)%in%tcga.deg.sig),]
########
enrichment=mg_clusterProfiler(as.vector(rownames(hub)))
write.table(enrichment$Enrich_tab,file = '02_DEGs/enrichment.txt',sep = '\t',quote = F,row.names = T,col.names = T)
kegg_dot=enrichplot::dotplot(enrichment$KEGG)+scale_y_discrete(labels=function(y)str_wrap(y,width = 25))+ scale_x_continuous(limits = c(0.05, 0.15))
bp_dot=enrichplot::dotplot(enrichment$GO_BP)+scale_y_discrete(labels=function(y)str_wrap(y,width = 25))+ scale_x_continuous(limits = c(0.12, 0.2))
ggsave('02_DEGs/enrichment_kegg.pdf',kegg_dot,height = 8,width = 9)
ggsave('02_DEGs/enrichment_GO_BP.pdf',bp_dot,height = 8,width = 9)



#############lasso##############
#########03######

tcga.t.exp=readMatrix('00_origin_datas/Preprocessed/tcga.t.exp1.txt')
tcga.t.cli=readMatrix('00_origin_datas/Preprocessed/tcga.cli.txt')
tcga.t.cli=tcga.cli[colnames(tcga.t.exp),]
identical(as.vector(tcga.t.cli$sampleID),as.vector(colnames(tcga.t.exp)))

tcga.t.exp_use=as.data.frame(tcga.t.exp)
tcga.subtype.cli=tcga.t.cli
identical(colnames(tcga.t.exp_use),as.vector(tcga.subtype.cli$sampleID))
colnames(tcga.subtype.cli)[1]="Samples"
rownames(tcga.subtype.cli)=tcga.subtype.cli$Samples


tcga.cox=cox_batch(t(scale(t(tcga.t.exp_use[rownames(hub),])))
                   ,time =  tcga.subtype.cli$OS.time/365
                   ,event =tcga.subtype.cli$OS)
dim(tcga.cox)

table(tcga.cox$p.value<0.05)
table(tcga.cox$p.value<0.01)
table(tcga.cox$p.value<0.001)
writeMatrix(tcga.cox,outpath = '03_Lasso/tcga.cox.txt')



p.cutoff=0.01
tcga.cox_use=tcga.cox
tcga.cox_use$coef=log(tcga.cox_use$HR)
tcga.cox_use$Gene=rownames(tcga.cox_use)
tcga.cox_use$type=rep('None',nrow(tcga.cox_use))
tcga.cox_use$type[which(tcga.cox_use$p.value<p.cutoff & tcga.cox_use$coef>0)]='Risk'
tcga.cox_use$type[which(tcga.cox_use$p.value<p.cutoff & tcga.cox_use$coef<0)]='Protective'
table(tcga.cox_use$type)

######### lasso
tcga.gene.sig=rownames(tcga.cox)[which(tcga.cox$p.value<p.cutoff)]
length(tcga.gene.sig)

table(tcga.cox_use$type)


tcga.cox_forVis=tcga.cox_use
tcga.cox_forVis=tcga.cox_forVis[which(tcga.cox_forVis$type %in% c('Risk','Protective')),]
tcga.cox_forVis$p.value=-log10(tcga.cox_forVis$p.value)
range(tcga.cox_forVis$p.value)


#################### LASSO
table(tcga.cox$p.value<0.05)
table(tcga.cox$p.value<0.01)
table(tcga.cox$p.value<0.001)

tcga.exp.sig=tcga.t.exp_use[tcga.gene.sig,]
tcga.exp.sig=t(tcga.exp.sig)
dim(tcga.exp.sig)


dim(tcga.exp.sig)
options(ggrepel.max.hnscerlaps = Inf)
tcga.subtype.cli$Samples=as.vector(tcga.subtype.cli$Samples)
identical(rownames(tcga.exp.sig),tcga.subtype.cli$Samples)
tcga.lasso.res=mg_lasso_cox_use(tcga.exp.sig
                                , time = tcga.subtype.cli$OS.time/365
                                , event = tcga.subtype.cli$OS
                                , nfolds = 5
                                , lambda.min = T
                                , figLabels=c('B','C'))
tcga.lasso.res$Genes

tcga.lasso.res$lambda

tcga.lasso.res$plot

tcga.exp.for.cox=tcga.t.exp_use[match(tcga.lasso.res$Genes,row.names(tcga.t.exp_use)),]
dim(tcga.exp.for.cox)
identical(colnames(tcga.exp.for.cox),tcga.subtype.cli$Samples)

lst.modl=createCoxModel_use((t(tcga.exp.for.cox))
                            , time = tcga.subtype.cli$OS.time/365
                            , event = tcga.subtype.cli$OS
                            , isStep =T)
lst.modl$Cox
lst.modl$Genes
lst.modl$fmla

lst.modl.Coef=lst.modl$Coef
names(lst.modl.Coef)=lst.modl$Genes
lst.modl.Coef

tcga.risk.score=lst.modl$Score
#tcga.risk.score=scale(tcga.risk.score)[,1]
tcga.risk.score=mosaic::zscore(tcga.risk.score)

range(tcga.risk.score)

lst.modl$Coef

gene.coef=data.frame(Gene=lst.modl$Genes,Coef=lst.modl$Coef)
gene.coef$Type=ifelse(lst.modl$Coef>0,'Risk','Protective')
gene.coef$Type=factor(gene.coef$Type,levels=c('Risk','Protective'))
table(gene.coef$Type)

fig3c=gene.coef %>% 
  ggplot(aes(reorder(Gene, Coef), Coef)) +
  geom_col(aes(fill = Type)) +
  geom_text(aes(label=round(Coef,digits = 3)),color="black",hjust = "left")+
  
  coord_flip() +
  scale_fill_manual(values=pal_nejm(alpha = 0.9)(8)[c(1,2)]) +
  coord_flip() +
  labs(x = "") +
  labs(y = "Lasso Cox coefficient") +
  theme_classic()+theme(legend.position = c(0,1))
# theme(axis.text.y = element_text(angle = 0, hjust = 1),legend.position="top")


fig3AB=mg_plot_lasso_use(fit = tcga.lasso.res$Mode1
                         , cv_fit = tcga.lasso.res$Model2
                         , show_text = F
                         , figLabels = c('A', 'B'))
fig3AB
fig3abc=mg_merge_plot(fig3AB,fig3c,nrow = 1,ncol = 2,widths = c(2,1))
#savePDF('PDFs/Fig7AB.pdf',fig7A,height = 4,width = 9)
#savePDF('PDFs/fig3abc.pdf',fig3abc,height = 5,width = 15)


tcga.exp.forCox<- cbind(time=tcga.subtype.cli$OS.time/365,
                        status=tcga.subtype.cli$OS,
                        t(tcga.t.exp_use)[rownames(tcga.subtype.cli), lst.modl$Genes])



dim(tcga.exp.forCox)

fmla <- as.formula(paste0("Surv(time, status) ~",paste0(lst.modl$Genes,collapse = '+')))
cox <- coxph(fmla, data =as.data.frame(tcga.exp.forCox))
fig3d=survminer::ggforest(cox,data=tcga.exp.forCox,noDigits = 3)
fig3abd=mg_merge_plot(fig3AB,fig3d,nrow = 1,ncol =2,widths = c(2,1))

#savePDF('PDFs/fig3abc.pdf',fig3abcd,height = 5,width = 16)


############### TCGA
#tcga.cutoff <- survminer::surv_cutpoint(data.frame(time=tcga.subtype.cli$OS.time/365, event=tcga.subtype.cli$OS, risk=tcga.risk.score),time = "time", event = "event",variables = c("risk"))
#tcga.cutoff=tcga.cutoff$cutpoint$cutpoint
#tcga.cutoff=median(tcga.risk.score)
tcga.cutoff=0
identical(colnames(tcga.exp.for.cox),tcga.subtype.cli$Samples)
risk.group.color=c("#EE0000FF","#3B4992FF")
names(risk.group.color)=c('High','Low')
tcga.roc=plotCoxModel_Batch_use(riskScore = tcga.risk.score
                                ,dat = t(tcga.exp.for.cox[match(lst.modl$Genes, row.names(tcga.exp.for.cox)), ])
                                , time = tcga.subtype.cli$OS.time/365
                                , event = tcga.subtype.cli$OS
                                , cutoff = tcga.cutoff
                                , labs = c('High','Low')
                                , title = 'RiskType'
                                , hetColor = c('#3B4992FF', 'white', '#EE0000FF')
                                
                                , pal = risk.group.color
                                , mks = c(1,2,3,4,5))
tcga.roc1=tcga.roc[[1]]
#pdf('PDFs/fig3e2.pdf',height = 6,width = 6)
tcga.roc1

dev.off()

tcga.group=ifelse(tcga.risk.score>tcga.cutoff,'High','Low')
tcga.group=data.frame(tcga.group)
colnames(tcga.group)='group'
table(tcga.group$group)
write.table(cbind(tcga.risk.score,tcga.group),file = '03_Lasso/tcga.group.txt',sep='\t',quote = F)


###########################
gse.t.exp <- read.table("00_origin_datas/GEO_expression_GSE31210.txt",header = T,check.names = F,fill=T,sep = "\t")
gse.t.cli <- read.table("00_origin_datas/GEO_cli_GSE31210.txt",header = T,check.names = F,fill=T,sep = "\t")
identical(colnames(gse.t.exp),as.vector(gse.t.cli$Samples))
gse.GSE31210.t.exp=gse.t.exp
match(lst.modl$Genes,row.names(gse.GSE31210.t.exp))
length(lst.modl$Genes)
gse.GSE31210.t.cli.os=gse.t.cli
gse.GSE31210.t.cli.os$Samples =as.vector(gse.GSE31210.t.cli.os$Samples )
identical(gse.GSE31210.t.cli.os$Samples , colnames(gse.GSE31210.t.exp)) ####
gse.GSE31210.model.dat=gse.GSE31210.t.exp[match(lst.modl$Genes,row.names(gse.GSE31210.t.exp)),]

gse.GSE31210.risk.score=predictScoreByCoxModel(coxModel = lst.modl
                                              ,(t(gse.GSE31210.model.dat)))
#gse.GSE31210.risk.score=scale(gse.GSE31210.risk.score)
gse.GSE31210.risk.score=mosaic::zscore(gse.GSE31210.risk.score)

lst.modl$fmla
#lst.vd.mod1$fmla

#gse.GSE31210.cutoff <- survminer::surv_cutpoint(data.frame(time=as.numeric(gse.GSE31210.t.cli.os$OS1),
# event=gse.GSE31210.t.cli.os$os_status1,
# risk=gse.GSE31210.risk.score), time =                                                   "time", event = "event", variables = c("risk"))
#gse.GSE31210.cutoff=gse.GSE31210.cutoff$cutpoint$cutpoint
gse.GSE31210.cutoff=0


test.roc=plotCoxModel_Batch_use(riskScore = gse.GSE31210.risk.score
                               , dat = t(gse.GSE31210.t.exp[intersect(lst.modl$Genes, row.names(gse.GSE31210.t.exp)),])
                               , time = as.numeric(gse.GSE31210.t.cli.os$OS.time/365) 
                               , event = as.numeric(gse.GSE31210.t.cli.os$OS)
                               , cutoff = gse.GSE31210.cutoff
                               , labs = c('High','Low')
                               , title = 'RiskType'
                               , hetColor = c('#3B4992FF', 'white', '#EE0000FF')
                               , pal = risk.group.color
                               , mks = c(1:5))
test.roc1=test.roc[[1]]
#pdf('PDFs/fig3g2.pdf',height = 6,width = 6)

test.roc1
dev.off()
gse.GSE31210.group=ifelse(gse.GSE31210.risk.score>gse.GSE31210.cutoff,'High','Low')
gse.GSE31210.group=data.frame(gse.GSE31210.group)
colnames(gse.GSE31210.group)='group'
table(gse.GSE31210.group)

write.table(cbind(gse.GSE31210.risk.score,gse.GSE31210.group),file = '03_Lasso//geo.group.txt',sep='\t',quote = F)


Fig4_ROC=mg_merge_plot(tcga.roc1,test.roc1,ncol=2,nrow=1)
Fig4=ggpubr::ggarrange(fig3abc,Fig4_ROC, ncol = 1, nrow = 2,heights = c(1,2))

ggsave('PDFs/Fig4.pdf',Fig4,height = 12,width = 15)
ggsave('PDFs/Fig4.jpg',Fig4,height = 12,width = 15)

identical(rownames(as.data.frame(tcga.risk.score)),rownames(tcga.subtype.cli))
tcga.risktype.cli=data.frame(tcga.subtype.cli,Riskscore=tcga.risk.score)
tcga.risktype.cli$Risktype=ifelse(tcga.risktype.cli$Riskscore>tcga.cutoff,'High','Low')
write.table(tcga.risktype.cli,file = '03_Lasso/tcga.risktype.cli.txt',sep='\t',quote = F)
#
#######################
#####################06.###################
#  #############
tcga.risktype.cli=read.table('03_Lasso/tcga.risktype.cli.txt',header = T,check.names = F,fill=T,sep = "\t")
tcga.t.exp=readMatrix('00_origin_datas/Preprocessed/tcga.t.exp1.txt')
tcga.t.cli=readMatrix('00_origin_datas/Preprocessed/tcga.cli.txt')
tcga.t.cli=tcga.cli[colnames(tcga.t.exp),]
identical(as.vector(tcga.t.cli$sampleID),as.vector(colnames(tcga.t.exp)))
identical(as.vector(tcga.risktype.cli$Samples),colnames(tcga.t.exp))

risk.group.color1=c("#EE0000FF","#3B4992FF")
names(risk.group.color1)=c('High','Low')

library(estimate)
#### ESTIMATE
#tcga.exp.estimate<-deconvo_estimate(eset=tcga.t.exp)
#save(tcga.exp.estimate,file='06_imm/tcga.exp.estimate.RData')
load('06_imm/tcga.exp.estimate.RData')
tcga.exp.estimate=get.IOBR.immu.format(tcga.exp.estimate)


############ TIMER 
#tcga.exp.timer<-deconvo_timer(eset=as.matrix(tcga.t.exp),indications=rep('LUAD',ncol(tcga.t.exp)))
#save(tcga.exp.timer,file='06_imm/tcga.exp.timer.RData')
load('06_imm/tcga.exp.timer.RData')
tcga.exp.timer=get.IOBR.immu.format(tcga.exp.timer)

library(EPIC)
############ EPIC 
#tcga.exp.epic<-deconvo_epic(eset=as.matrix(tcga.t.exp),tumor = TRUE)
#save(tcga.exp.epic,file='06_imm/tcga.exp.epic.RData')
load('06_imm/tcga.exp.epic.RData')
tcga.exp.epic=get.IOBR.immu.format(tcga.exp.epic)

############ MCP-counter 
#tcga.exp.mcp<-deconvo_mcpcounter(eset=as.matrix(tcga.t.exp))
#save(tcga.exp.mcp,file='06_imm/tcga.exp.mcp.RData')
load('06_imm/tcga.exp.mcp.RData')
tcga.exp.mcp=get.IOBR.immu.format(tcga.exp.mcp)

### CIBERSORT
#tcga.exp.cibersort<-deconvo_cibersort(eset=tcga.t.exp,arrays=F)
#save(tcga.exp.cibersort,file='06_imm/tcga.exp.cibersort.RData')
load('06_imm/tcga.exp.cibersort.RData')
tcga.exp.cibersort=get.IOBR.immu.format(tcga.exp.cibersort)

#######sssGSEA#######
#geo.immu.ssgsea=immu_ssgsea(exp = tcga.t.exp)
#save(geo.immu.ssgsea,file='06_imm/geo.immu.ssgsea.RData')
#load('06_imm/geo.immu.ssgsea.RData')



#
tcga.t.estimate=tcga.exp.estimate[rownames(tcga.risktype.cli),1:3]
tcga.t.timer=tcga.exp.timer[rownames(tcga.risktype.cli),]
tcga.t.epic=tcga.exp.epic[rownames(tcga.risktype.cli),]
tcga.t.mcp=tcga.exp.mcp[rownames(tcga.risktype.cli),]
tcga.t.cibersort=tcga.exp.cibersort[rownames(tcga.risktype.cli),1:22]
#tcga.t.ssGSEA28=as.data.frame(geo.immu.ssgsea[tcga.risktype.cli$Samples,])




fig5a=get_PlotMutiBoxplot(tcga.t.estimate,tcga.risktype.cli
                          ,group_cols = risk.group.color1
                          ,legend.pos = NULL
                          ,ylab = 'Score'
                          ,group.val = 'Risktype',xangle=45)+labs(color='Risktype')

fig5a=groupViolin(tcga.t.estimate,
                  tcga.risktype.cli$Risktype,
                  ylab = 'ssgsea Immune Score',
                  group_col=risk.group.color1)

fig5b=get_PlotMutiBoxplot(tcga.t.timer,tcga.risktype.cli
                          ,group_cols = risk.group.color1
                          ,legend.pos = NULL
                          ,ylab = 'Score'
                          ,group.val = 'Risktype',xangle=45)+labs(color='Risktype')
fig5b=groupViolin(tcga.t.timer,
                       tcga.risktype.cli$Risktype,
                       ylab = 'Score',
                       group_col=risk.group.color1)


fig5c=get_PlotMutiBoxplot(tcga.t.epic,tcga.risktype.cli
                          ,group_cols = risk.group.color
                          ,legend.pos = NULL
                          ,ylab = 'Score'
                          ,group.val = 'Risktype',xangle=45)+labs(color='Risktype')
fig5c

fig5d=get_PlotMutiBoxplot(tcga.exp.mcp,tcga.risktype.cli
                          ,group_cols = risk.group.color1
                          #,legend.pos = NULL
                          ,ylab = 'Score'
                          ,group.val = 'Risktype',xangle=45)+labs(color='Risktype')
fig5d=groupViolin(tcga.exp.mcp,
                       tcga.risktype.cli$Risktype,
                       ylab = 'Score',
                       group_col=risk.group.color1)







fig5e=get_PlotMutiBoxplot(tcga.t.cibersort,tcga.risktype.cli
                          ,group_cols = risk.group.color
                          ,legend.pos = NULL
                          ,ylab = 'Score'
                          ,group.val = 'Risktype',xangle=45)+labs(color='Risktype')
fig5e



################imm############

# 05######
library('oncoPredict')
tcga.t.exp=readMatrix('00_origin_datas/Preprocessed/tcga.t.exp1.txt')
tcga.t.cli=readMatrix('00_origin_datas/Preprocessed/tcga.cli.txt')
tcga.t.cli=tcga.cli[colnames(tcga.t.exp),]
identical(as.vector(tcga.t.cli$sampleID),as.vector(colnames(tcga.t.exp)))
identical(tcga.t.cli$sampleID,colnames(tcga.t.exp))
#drug_exp=as.matrix(tcga.t.exp)

#GDSC2_Expr = readRDS(file=file.path(dir,'Training Data/GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
#GDSC2_Res = readRDS(file = file.path(dir,"Training Data/GDSC2_Res.rds"))
#GDSC2_Res <- exp(GDSC2_Res)
#calcPhenotype(trainingExprData = as.matrix(GDSC2_Expr),
#              trainingPtype = as.matrix(GDSC2_Res),
#              testExprData = as.matrix(drug_exp),
#              batchCorrect = 'eb',  #   "eb" for ComBat
#              powerTransformPhenotype = TRUE,
#              removeLowVaryingGenes = 0.2,
#              minNumSamples = 10,
#              printOutput = TRUE,
#             removeLowVaringGenesFrom = 'rawData' )

tcga.drug.ic50=read.csv('06_imm/calcPhenotype_Output/DrugPredictions.csv',row.names = 1)
dim(tcga.drug.ic50)
tcga.exp.forCox=cbind(tcga.risktype.cli$Riskscore
                      ,t(tcga_exp)[as.vector(tcga.risktype.cli$Samples),lst.modl$Genes])
colnames(tcga.exp.forCox)[1]='Riskcsore'

all(rownames(tcga.drug.ic50)==rownames(tcga.exp.forCox))
all(rownames(as.data.frame(tcga.risk.score))==rownames(tcga.exp.forCox))
tcga.drug.ic501=tcga.drug.ic50[match(rownames(tcga.exp.forCox),rownames(tcga.drug.ic50)),]
all(rownames(tcga.drug.ic501)==rownames(tcga.exp.forCox))
######## 

mat1=tcga.exp.forCox
mat2=tcga.drug.ic501
all(rownames(mat2)==rownames(mat1))

#df=cor(risk.score,tcga.drug.ic501)

# 
correlation <- matrix(NA, nrow = ncol(mat1), ncol = ncol(mat2))
p_values <- matrix(NA, nrow = ncol(mat1), ncol = ncol(mat2))
rownames(correlation)=colnames(mat1)
colnames(correlation)=colnames(mat2)
rownames(p_values)=colnames(mat1)
colnames(p_values)=colnames(mat2)
# 
for (i in 1:ncol(mat1)) {
  for (j in 1:ncol(mat2)) {
    result <- cor.test(mat1[, i], mat2[, j])
    correlation[i, j] <- result$estimate
    p_values[i, j] <- result$p.value
  }
}


p_values1=as.data.frame(t(as.data.frame(p_values)))
p_values1$"to"=rownames(p_values1)
p_values2 <- melt(p_values1, id.vars = "to", variable.name = "from", value.name = "p.adj")
correlation1=as.data.frame(t(as.data.frame(correlation)))
correlation1$"to"=rownames(correlation1)
correlation2 <- melt(correlation1, id.vars = "to", variable.name = "from", value.name = "cor")
identical(rownames(p_values1),rownames(correlation1))

df=merge(p_values2,correlation2,by=c("to","from"))
df=df[,c(2,1,3,4)]

df=data.frame(df)
df <- df %>%
  mutate(col = cut(cor, breaks = c(-1, 0, 1),
                   labels = c("negative", "positive")),
         p.signif = cut(p.adj, breaks = c(0,0.001, 0.01, 0.05,1),
                        labels = c("***", "**","*",""),
                        right = FALSE, include.lowest = TRUE))
df=data.frame(df)
head(df)



inds=which(df$from=='Riskcsore')
drug.selected=df$to[inds][which(df$p.adj[inds]<0.05 & abs(df$cor[inds])>0.4)]
# View(df[inds,])
length(drug.selected)
table(df$col[inds])
unique(drug.selected)
#### 
# df=df[df$to %in% (names(which(table(df$to,df$p.signif)[,4]<2))),]
#########
df=df[df$to %in% drug.selected,]
head(df)
unique(df$to)
########
df1=reshape2::dcast(df[,1:3],from ~ to)
head(df1)
df1[1:3,1:4]
rownames(df1)=df1$from
df1=df1[,-1]
head(df1)

hc.r=hclust(dist(df1[lst.modl$Genes,]))
hc.r$order
plot(hc.r)

pw.order=c(lst.modl$Genes[hc.r$order],'Riskcsore')
df$from=factor(df$from,pw.order)

##########################
drug_plot=ggplot(df,aes(y = from,x = to))+
  geom_point(aes(size=(-log10(p.adj)),color=cor))+
  geom_text(aes(y = from,x= to,label=p.signif))+
  theme_bw()+
  # scale_color_continuous(type = circlize::colorRamp2(c(-1, 0, 1), c(risk.group.color[2], 'white', risk.group.color[1])))+
  # scale_colour_gradient2(low = risk.group.color[2],high =risk.group.color[1] )+
  scale_colour_gradient2(low = '#3B4992FF',high ='#EE0000FF' )+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=1))+
  ylab('')+xlab('')+labs(color="Corr",size="-log10(FDR)")

##########################
unique(drug.selected)
drug.selected_use=c("Doramapimod_1042")
setdiff(drug.selected_use,df$to)
all(rownames(tcga.drug.ic501)==rownames(tcga.risktype.cli))

p.all=list()
#for(drug in drug.selected_use){
drug <- drug.selected_use
df=data.frame(RiskType=tcga.risktype.cli$Risktype,RiskScore=tcga.risktype.cli$Riskscore,tcga.drug.ic501[,drug])
colnames(df)[3]=drug
head(df)
p2=plot_cor_point(x = df$RiskScore,y = df$Doramapimod_1042,
                  xlab = 'RiskScore',ylab = 'Doramapimod_1042 IC50',
                  right_col ="#9998FF",top_col = "#C99BFF",
                  marginal.type='density',method = 'spearman')


library(gghalves)
library(ggplot2)
library(cowplot)
p1=ggplot(data = df, aes(x = RiskType, y = Doramapimod_1042, fill = RiskType)) +
  geom_half_violin(side='R',position = position_nudge(x = 0.2, y = 0), alpha = 1) +
  geom_point(aes(y = Doramapimod_1042, color = RiskType), position = position_jitter(width = 0.15), size = 1, alpha = 1) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 1) +
  labs(y = 'Doramapimod_1042 IC50', x = 'RiskType') +
  guides(fill = FALSE, color = FALSE) +
  scale_fill_manual(values = risk.group.color) +
  scale_colour_manual(values = risk.group.color) +theme_classic2()+ggpubr::stat_compare_means(comparisons = list(c(1,2)),method = 'wilcox.test',label= "p.signif")
p1
fig6a=mg_merge_plot(fig5a,fig5b,fig5d,nrow = 1,ncol = 3,labels = LETTERS[1:3],common.legend = T ,widths = c(1,1.2,1.5))
fig6b=mg_merge_plot(drug_plot,p1,p2,nrow = 1,ncol = 3,labels = LETTERS[4:6],widths = c(1.5,1,1))

fig6=mg_merge_plot(fig6a,fig6b,nrow = 2,ncol =1)

ggsave("PDFs/Fig5.pdf", fig6, width = 14, height = 10)

