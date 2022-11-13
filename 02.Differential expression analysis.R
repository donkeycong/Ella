#Citation package
library(limma)
library(pheatmap)

inputFile="geneMatrix.txt"       #input

logFCfilter=log2(1.5)                   #logFC
P.Val.Filter=0.05                    #P.value
setwd("02.diff")      #Setting up a working directory

#Read input files and collate input files
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data=normalizeBetweenArrays(data)


#Read all files in the directory ending with "s1.txt"

sampleName1=c()
files=dir()
files=grep("s1.txt$", files, value=T)
for(file in files){
    rt=read.table(file, header=F, sep="\t", check.names=F)      #Reading input files
    geneNames=as.vector(rt[,1])      #Name of extracted gene
    uniqGene=unique(geneNames)       #Genes taken unique
    sampleName1=c(sampleName1, uniqGene)
}

#Read all files in the directory ending with "s2.txt"
sampleName2=c()
files=dir()
files=grep("s2.txt$", files, value=T)
for(file in files){
    rt=read.table(file, header=F, sep="\t", check.names=F)      #Reading input files
    geneNames=as.vector(rt[,1])      #Name of extracted gene
    uniqGene=unique(geneNames)       #Genes taken unique
    sampleName2=c(sampleName2, uniqGene)
}

#Extraction of data from the experimental and control groups
conData=data[,sampleName1]
treatData=data[,sampleName2]
data=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

#Analysis of variances
Type=c(rep("con",conNum),rep("treat",treatNum))
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("con","treat")
fit <- lmFit(data,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)
allDiffOut=rbind(id=colnames(allDiff),allDiff)
write.table(allDiffOut, file="all.txt", sep="\t", quote=F, col.names=F)

#Output corrected expressions
outData=rbind(id=paste0(colnames(data),"_",Type),data)
write.table(outData, file="normalize.txt", sep="\t", quote=F, col.names=F)



#Exporting the results of differential expression analysis
diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & P.Value < P.Val.Filter )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut, file="diff.txt", sep="\t", quote=F, col.names=F)

#Exporting differential gene expression
diffGeneExp=data[row.names(diffSig),]
diffGeneExpOut=rbind(id=paste0(colnames(diffGeneExp),"_",Type), diffGeneExp)
write.table(diffGeneExpOut, file="diffGeneExp.txt", sep="\t", quote=F, col.names=F)


#Heat mapping of differential genes
geneNum=50
diffSig=diffSig[order(as.numeric(as.vector(diffSig$logFC))),]
diffGeneName=as.vector(rownames(diffSig))
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(2*geneNum)){
    hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
    hmGene=diffGeneName
}
hmExp=data[hmGene,]
Type=c(rep("Con",conNum),rep("Treat",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf", width=10, height=8)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         cluster_rows = T,
         show_colnames = T,
         show_rownames = T,
         scale="row",
         fontsize = 8,
         fontsize_row=7,
         fontsize_col=8)
dev.off()

