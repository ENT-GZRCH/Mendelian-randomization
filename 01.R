if(! require("GEOquery")) BiocManager::install("GEOquery",update=F,ask=F)

library(GEOquery)
#GSE36830--------------------------
setwd("F:/stroke")
gset=getGEO('GSE36830',destdir = '.',getGPL = F)
#
exprSet=exprs(gset[[1]])
pdata=pData(gset[[1]])

group_list=c(rep('CRSwNP',12),rep('CT',6))
group_list=factor(group_list,levels = c('CT','CRSwNP'))

##————————————标准化和log2转化——————————————————#
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
max(exprSet)
exprSet=log2(exprSet+1)
exprSet=as.data.frame(exprSet)


library(hgu133plus2.db)
probeset <- rownames(exprSet)
SYMBOL <-  annotate::lookUp(probeset,"hgu133plus2.db", "SYMBOL")
symbol = as.vector(unlist(SYMBOL))
probe2symbol = data.frame(probeset,symbol,stringsAsFactors = F)
library(dplyr)
library(tibble)
exprSet=as.data.frame(exprSet)
exprSet <- exprSet %>% 
  rownames_to_column(var="probeset") %>% 
  #合并探针的信息
  inner_join(probe2symbol,by="probeset") %>% 
  #去掉多余信息
  dplyr::select(-probeset) %>% 
  #重新排列
  dplyr::select(symbol,everything()) %>% 
  #求出平均数(这边的点号代表上一步产出的数据)
  mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% 
  #去除symbol中的NA
  filter(symbol != "NA") %>% 
  #把表达量的平均值按从大到小排序
  arrange(desc(rowMean)) %>% 
  # symbol留下第一个，也就是最大的
  distinct(symbol,.keep_all = T) %>% 
  #反向选择去除rowMean这一列
  dplyr::select(-rowMean) %>% 
  # 列名变成行名
  column_to_rownames(var = "symbol")

dim(exprSet)

exprSet_bioc=exprSet

save(exprSet,pdata,file ='GSE36830.Rdata')
write.table(exprSet,file ='GSE36830.txt',sep = '\t',col.names = NA,quote = F)

design=model.matrix(~ group_list)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
allDiff=topTable(fit,adjust='fdr',coef=2,number=Inf,p.value=1) 


load("deg.Rdata")
gene=df$mRNAs

alldiff_mito=allDiff[rownames(allDiff) %in% gene,]

rt=exprSet[rownames(alldiff_mito),]
pheatmap::pheatmap(rt,scale = 'row',cluster_cols = F,border_color = NA,
                   breaks = seq(-2,2,length.out = 100))

#GSE23552--------------------------
gset=getGEO('GSE23552',destdir = '.',getGPL = F)

exprSet=exprs(gset[[1]])
pdata=pData(gset[[1]])

pdata1=pdata[pdata$characteristics_ch1.1=='tissue: polyp',]
pdata2=pdata[pdata$characteristics_ch1=='disease state: Normal',]

pdata <- rbind(pdata1, pdata2)

exprSet=as.data.frame(exprSet)
exprSet1=exprSet[,rownames(pdata1)]
exprSet2=exprSet[,rownames(pdata2)]
exprSet1=as.data.frame(exprSet1)
exprSet2=as.data.frame(exprSet2)

group_list2=c(rep('CRSwNP',10),rep('CT',13))
group_list2=factor(group_list2,levels = c('CT','CRSwNP'))


#按列合并exprSet1,2
exprSet=cbind(exprSet1,exprSet2)

boxplot(exprSet,outline=FALSE, notch=T,col=group_list2, las=2)
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list2, las=2)
max(exprSet)
exprSet=log2(exprSet+1)
max(exprSet)
exprSet=as.data.frame(exprSet)

library(data.table)
anno = fread("GPL5175-3188.txt",data.table = F,sep = '\t') 
probe2symbol <- anno[,c("ID","gene_assignment")]  



library(data.table)
symbol=tstrsplit(anno$gene_assignment,"//",fixed=TRUE)[[2]]  #拆分
symbol <- trimws(symbol,which = c("both","left","right"),whitespace = "[ \t\r\n]") #去除空格
probe2symbol["symbol"] <- symbol
probe2symbol <-probe2symbol[,c("ID","symbol")] #去掉gene_assigment列
anno <- probe2symbol
#看一下列名
colnames(anno) 
#!!!!
anno <- anno[!is.na(as.character(anno$symbol)), ]
#无脑跑代码即可
tmp = rownames(exprSet)%in% anno[,1] 
exprSet = exprSet[tmp,] 
dim(exprSet) 
match(rownames(exprSet),anno$ID) 
anno = anno[match(rownames(exprSet),anno$ID),] 
match(rownames(exprSet),anno$ID) 
dim(exprSet) 
dim(anno) 
tail(sort(table(anno[,2])), n = 12L) 
#注释到相同基因的探针，保留表达量最大的那个探针 
{ 
  MAX = by(exprSet, anno[,2],  
           function(x) rownames(x)[ which.max(rowMeans(x))]) 
  MAX = as.character(MAX) 
  exprSet = exprSet[rownames(exprSet) %in% MAX,] 
  rownames(exprSet) = anno[ match(rownames(exprSet), anno[,1] ),2] 
} 
dim(exprSet) 


save(exprSet,pdata,file ='GSE23552.Rdata')
write.table(exprSet,file ='GSE23552.txt',sep = '\t',col.names = NA,quote = F)
load("GSE23552.Rdata")
#######去批次
library(limma)
library(sva)
# 运行以下代码去批次
mergeFile="merge.preNorm.txt"            #合并后的文件名称
normalizeFile="merge.normalzie.txt"      #矫正后的文件名称
# 一定是在昨天输出这两个txt的文件夹下
files=c('GSE36830.txt','GSE23552.txt')       #输入文件名称

length(files)
library(data.table)

geneList=list()
for(i in 1:length(files)){
  fileName=files[i]
  rt=fread(fileName,data.table = F)
  header=unlist(strsplit(fileName, "\\.|\\-"))
  geneList[[header[1]]]=as.vector(rt[,1])
}
# 共同基因
intersectGenes=Reduce(intersect, geneList)
library(limma)
#数据合并
allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
  fileName=files[i]
  header=unlist(strsplit(fileName, "\\.|\\-"))
  #读取输入文件，并对输入文件进行整理
  rt=fread(fileName,data.table = F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
  rt=avereps(data)
  colnames(rt)=paste0(header[1], "_", colnames(rt))
  #对数值大的数据取log2
  #qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  #LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
  #if(LogC){
  #rt[rt<0]=0
  #rt=log2(rt+1)}
  #rt=normalizeBetweenArrays(rt)
  #数据合并
  if(i==1){
    allTab=rt[intersectGenes,]
  }else{
    allTab=cbind(allTab, rt[intersectGenes,])
  }
  batchType=c(batchType, rep(header[1],ncol(rt)))
}

allTabOut=rbind(geneNames=colnames(allTab), allTab)
write.table(allTabOut, file=mergeFile, sep="\t", quote=F, col.names=F)
#对数据进行批次矫正，输出矫正后的结果
library(sva)
normalizeTab=ComBat(allTab, batchType, par.prior=TRUE)
normalizeTab=rbind(geneNames=colnames(normalizeTab), normalizeTab)
write.table(normalizeTab, file=normalizeFile, sep="\t", quote=F, col.names=F)

#合并前的PCA###################
#install.packages('ggplot2')
library(ggplot2)        #引用包
#读取输入文件,提取数据
rt=read.table('merge.preNorm.txt', header=T, sep="\t", check.names=F, row.names=1)
data=t(rt)
Project=gsub("(.*?)\\_.*", "\\1", rownames(data))
rownames(data)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data))

#PCA分析
data.pca=prcomp(data)
pcaPredict=predict(data.pca)
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Type=Project)

#定义颜色
bioCol=c("#bc2c29",'#0072b5',"#20854e","#ef7c1c","#EE4C97","#FF9900","#20854E","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121")
bioCol=bioCol[1:length(levels(factor(Project)))]

#绘制PCA图

p=ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = Type)) +
  scale_colour_manual(name="",  values=bioCol)+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)

#合并后PCA###################
#读取输入文件,提取数据
rt=read.table('merge.normalzie.txt', header=T, sep="\t", check.names=F, row.names=1)
data=t(rt)
Project=gsub("(.*?)\\_.*", "\\1", rownames(data))
rownames(data)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data))

#PCA分析
data.pca=prcomp(data)
pcaPredict=predict(data.pca)
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Type=Project)

#定义颜色
bioCol=c("#bc2c29",'#0072b5',"#20854e","#ef7c1c","#EE4C97","#FF9900","#20854E","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121")
bioCol=bioCol[1:length(levels(factor(Project)))]

#绘制PCA图

p=ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = Type)) +
  scale_colour_manual(name="",  values=bioCol)+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)

### 富集（运气做）
rt=read.table('merge.normalzie.txt', header=T, sep="\t", check.names=F, row.names=1)
# 牢记样本顺序
load('GSE36830.Rdata')
pdata1=pdata
group_list1=c(rep('CT',6),rep('CRSwNP',12))
load('GSE23552.Rdata')
pdata2=pdata
group_list2=c(rep('CRSwNP',10),rep('CT',13))

group_list=c(group_list1,group_list2)
group_list=factor(group_list,levels=c('CT','CRSwNP'))
#3的数据集，提前2个
a=grep('GSE36830',colnames(rt))
b=grep('GSE23552',colnames(rt))

exprSet=rt[,c(a,b)]

load('deg.Rdata')
df <- df[df$p_val_adj<0.05,]
df <- df[df$avg_log2FC>0.5,]
#!df <- subset(df, p_val_adj < 0.05, avg_log2FC > 0.5)
gene=df$mRNAs

mito=list(gene)
names(mito)='ScRNA'

library(GSVA)
exprSet=as.matrix(exprSet)
score=gsva(exprSet,mito,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
score
max(score)
min(score)


score=as.data.frame(t(score))

score$group=group_list
library(ggpubr)

ggboxplot(score,x='group',y='ScRNA',
          color = "black",
          fill = "group",
          palette=c("#20854e","#E18727"))+stat_compare_means()











