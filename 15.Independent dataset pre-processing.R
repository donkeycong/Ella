#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)               #Citation package
expFile="geneMatrix.txt"     #Expression data files
conFile="s1.txt"             #Sample file for the control group
treatFile="s2.txt"           #Sample file of the experimental group
setwd("15.normalize")      #Setting up a working directory

#Read input files and collate input files
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#Read the control sample file and extract the expression data from the control group
s1=read.table(conFile, header=F, sep="\t", check.names=F)
sampleName1=as.vector(s1[,1])
conData=data[,sampleName1]

#Read the sample file of the experimental group and extract the expression data of the experimental group
s2=read.table(treatFile, header=F, sep="\t", check.names=F)
sampleName2=as.vector(s2[,1])
treatData=data[,sampleName2]

#Data rectification
data=cbind(conData, treatData)
data=normalizeBetweenArrays(data)

#Correct the data and output the corrected expression
conNum=ncol(conData)
treatNum=ncol(treatData)
Type=c(rep("con",conNum),rep("treat",treatNum))
outData=rbind(id=paste0(colnames(data),"_",Type),data)
write.table(outData, file="diffGeneExp.txt", sep="\t", quote=F, col.names=F)
