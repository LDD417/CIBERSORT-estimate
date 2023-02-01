#####1.CIBERSORT 算法####
rm(list = ls())

library("limma")            #加载voom函数
load("data/TCGA_tumor.Rdata")

v <-voom(data, plot = F, save.plot = F)
out=v$E
out=rbind(ID=colnames(out),out)
write.table(out,file="data/symbol.txt",sep="\t",quote=F,col.names=F)       

source("CIBERSORT.R")
results=CIBERSORT("data/ref.txt", "data/symbol.txt", perm=100, QN=TRUE)

#####1.estimate 算法####
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos = rforge, dependencies = TRUE)
library(estimate)
rm(list = ls())

load("data/TCGA_tumor.Rdata")
out=rbind(ID=colnames(data),data)
write.table(out,file="data/tumorsymbol.txt",sep="\t",quote=F,col.names=F)

filterCommonGenes(input.f="data/tumorsymbol.txt", 
                  output.f="data/commonGenes.gct", 
                  id="GeneSymbol")

estimateScore(input.ds = "data/commonGenes.gct",
              output.ds="data/estimateScore.gct", 
              platform="illumina")

#输出每个样品的打分
scores=read.table("data/estimateScore.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.","\\-",rownames(scores))
out=rbind(ID=colnames(scores),scores)
write.table(out,file="data/tcgascores.txt",sep="\t",quote=F,col.names=F)
