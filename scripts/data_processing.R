rm(list = ls())
library(magrittr)
library(dplyr)
ann=readMatrix('gencode.v22.annotation.gene.txt')
ann.pcg=ann[which(ann$gene_type=='protein_coding'),]###
dim(ann.pcg)# [1] 
tcga.cli = read.delim("00_origin_datas/TCGA/TCGA-LUAD.clinical.tsv")
dim(tcga.cli)# [1]    721  89
tcga.cli1=readMatrix('00_origin_datas/TCGA/TCGA-LUAD.survival.tsv',row = F)########### TCGA
dim(tcga.cli1) #645   4
tcga.cli1 <- tcga.cli1[,c(1,4,3,2)]
head(tcga.cli1)
colnames(tcga.cli1)[1:2]=c("sampleID","PATIENT")
surv_data = dplyr::select(tcga.cli, c(sample,
                                      submitter_id,
                                      gender.demographic,
                                      #tumor_grade.diagnoses,
                                      ajcc_pathologic_stage.diagnoses,
                                      ajcc_pathologic_t.diagnoses,
                                      ajcc_pathologic_n.diagnoses,
                                      ajcc_pathologic_m.diagnoses,
                                      age_at_index.demographic
                                      #vital_status.demographic
                                      #days_to_diagnosis.diagnoses,
                                      #days_to_last_follow_up.diagnoses
)) 
head(surv_data)
colnames(surv_data)[1]=c("sampleID")
meta_data=merge(surv_data,tcga.cli1,by="sampleID")
colnames(meta_data)[c(3:8)]=c("Gender" ,"Stage","pathologic_T","pathologic_N","pathologic_M","age_at")
head(meta_data)

meta_data[which(meta_data=="not reported",arr.ind = T)]=NA
meta_data[which(meta_data=="NA",arr.ind = T)]=NA
meta_data[which(meta_data=="",arr.ind = T)]=NA
#meta_data <- lapply(meta_data, as.vector)





table(meta_data$pathologic_N)
meta_data$pathologic_N=substr(meta_data$pathologic_N,0,2)
meta_data$pathologic_N[which(meta_data$pathologic_N=='NX')]=NA
table(meta_data$pathologic_N)

table(meta_data$pathologic_M)
meta_data$pathologic_M[which(meta_data$pathologic_M=="MX")]=NA
#meta_data$pathologic_M[which(meta_data$pathologic_M=="cM0 (i+)")]=NA
meta_data$pathologic_M=gsub("[abcd]","",meta_data$pathologic_M)
meta_data$pathologic_M=as.vector(meta_data$pathologic_M)
table(meta_data$pathologic_M)

table(meta_data$pathologic_T)
meta_data$pathologic_T=gsub("[abcd]","",meta_data$pathologic_T)
meta_data$pathologic_T[which(meta_data$pathologic_T=="TX")]=NA
table(meta_data$pathologic_T)

table(meta_data$Stage)
meta_data$Stage=gsub("[ABC]","",meta_data$Stage)
#meta_data$Stage[which(meta_data$Stage=="stge x")]=NA
table(meta_data$Stage)

table(meta_data$age_at)
meta_data$age_at=as.numeric(meta_data$age_at)
meta_data$Age=ifelse(meta_data$age_at>60,'>60','<=60')
meta_data=as.data.frame(meta_data)
table(meta_data$Age)

meta_data=as.data.frame(meta_data)
head(meta_data)
meta_data=meta_data[meta_data$OS.time>30,]#######


###TCGA-LUAD FPKM
tcga.exp.fpkm=read.table("00_origin_datas/TCGA/TCGA-LUAD.star_fpkm.tsv"
                         ,row.names = 1,sep="\t",header=T,as.is=T,quote="\""
                         ,fill=T,check.names = F,stringsAsFactors = F)
range(tcga.exp.fpkm)
dim(tcga.exp.fpkm)# # [1] 60483   585
tcga.exp.fpkm=2^(tcga.exp.fpkm)-1
tcga.exp.tpm=mg_FPKM2TPMs(tcga.exp.fpkm)
tcga.exp.tpm=log2(tcga.exp.tpm+1)
rownames(tcga.exp.tpm)=unlist(lapply(strsplit(as.vector(rownames(tcga.exp.tpm)),"\\."),function(x) x[1]))###

tcga.exp=exp_ensg2symbol(tcga.exp.tpm,method = "max")
#######保存

range(tcga.exp)

tcga.exp=tcga.exp[rownames(tcga.exp) %in% ann.pcg$gene_name,]
dim(tcga.exp)
#  18177   448
table(substr(colnames(tcga.exp),14,15))
#01  02  11 
#528   2  59 
#colnames(tcga.exp)=substr(colnames(tcga.exp),0,15)

#tcga.exp=tcga.exp[,-which(substr(colnames(tcga.exp),14,15)=="06")]
######### -01 
inds1=which(substr(colnames(tcga.exp),14,15)=="01")
inds2=which(substr(colnames(tcga.exp),14,15)=="11")
#table(substr(colnames(tcga.exp),14,16))

sample_com=intersect(colnames(tcga.exp)[c(inds1,inds2)], as.vector(meta_data$sampleID))
table(substr(sample_com,14,15))
sample_com=sample_com[which(substr(sample_com,16,16)=="A")]
sample_com=sample_com[which(substr(sample_com,15,15)=="1")]
table(substr(sample_com,14,16))
rownames(meta_data)=meta_data$sampleID
tcga.exp=tcga.exp[,sample_com]
meta_data=meta_data[sample_com,]
#colnames(tcga.exp)=substr(colnames(tcga.exp),0,15)
#rownames(tcga.cli)=substr(rownames(tcga.cli),0,15)

identical(as.vector(meta_data$sampleID),as.vector(colnames(tcga.exp)))
table(substr(colnames(tcga.exp),14,15))
#01  11 
#366  32
writeMatrix(tcga.exp,outpath = '00_origin_datas/Preprocessed/tcga.exp.txt')
writeMatrix(meta_data,outpath = '00_origin_datas/Preprocessed/tcga.cli.txt')
dim(tcga.exp)


######### -01 
inds1=which(substr(colnames(tcga.exp),14,15)=="01")
inds2=which(substr(colnames(tcga.exp),14,15)=="11")

tcga.t.exp=tcga.exp[,inds1]
dim(tcga.t.exp)
# [1] 18177   490
tcga.n.exp=tcga.exp[,inds2]
dim(tcga.n.exp)
# [1] 18177    58

writeMatrix(tcga.t.exp,outpath = '00_origin_datas/Preprocessed/tcga.t.exp.txt')
writeMatrix(tcga.n.exp,outpath = '00_origin_datas/Preprocessed/tcga.n.exp.txt')





#GSE31210 ############
GSE31210 <- getGEOExpData('GSE31210')
GSE31210_cli <- GSE31210$Sample
table(GSE31210_cli$tissue)
GSE31210_cli <- GSE31210_cli[!startsWith(GSE31210_cli$tissue, 'normal'), ]

GSE31210_cli <- GSE31210_cli[, c("Acc", "days before death/censor", "death")]
colnames(GSE31210_cli) <- c("Samples", "OS.time", "OS")
rownames(GSE31210_cli) <- GSE31210_cli$Samples

# GSE31210_cli <- GSE31210_cli[GSE31210_cli$OS.time != 'Not available', ]
table(GSE31210_cli$OS)
GSE31210_cli$OS <- ifelse(GSE31210_cli$OS == 'alive', 0, 1)
GSE31210_cli <- crbind2DataFrame(GSE31210_cli)


# GSE31210_cli$OS.time <- GSE31210_cli$OS.time * 30
table(GSE31210_cli$OS)
# GSE31210_cli$OS <- ifelse(GSE31210_cli$OS == 'no', 0, 1)
GSE31210_cli$DataSet <- 'GSE31210'

GSE31210_exp <- GSE31210$Exp$GPL570_54675_Data_col1

GSE31210_exp <- exp_probe2symbol_v2(GSE31210_exp[, GSE31210_cli$Samples],
                                    anno = GPL570_anno[, c(5,3)])
GSE31210_exp <- log2(GSE31210_exp + 1)
boxplot(GSE31210_exp[, 1:5])
# intersect(genes_lncRNA, rownames(GSE19188_exp))
write.table(GSE31210_exp,file="00_origin_datas/GEO_expression_GSE31210.txt",sep="\t",quote=F,col.names = TRUE)
write.table(GSE31210_cli,file="00_origin_datas/GEO_cli_GSE31210.txt",sep="\t",quote=F,col.names = TRUE)

tcga.exp=readMatrix('00_origin_datas/Preprocessed/tcga.exp.txt')
tcga.t.exp=readMatrix('00_origin_datas/Preprocessed/tcga.t.exp.txt')
gse.t.exp <- read.table("00_origin_datas/GEO_expression_GSE31210.txt",header = T,check.names = F,fill=T,sep = "\t")
comm_gene=intersect(rownames(tcga.exp),rownames(gse.t.exp))
tcga.exp=tcga.exp[comm_gene,]
tcga.t.exp=tcga.t.exp[comm_gene,]

writeMatrix(tcga.t.exp,outpath = '00_origin_datas/Preprocessed/tcga.t.exp1.txt')
writeMatrix(tcga.exp,outpath = '00_origin_datas/Preprocessed/tcga.exp1.txt')
