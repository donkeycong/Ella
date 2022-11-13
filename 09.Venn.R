##install.packages("venn")


library(venn)                   #Citation package
outFile="interGenes.txt"        #Output file name
setwd("09.venn")    #Setting up a working directory
geneList=list()

#Read the results file of the lasso regression
rt=read.table("LASSO.gene.txt", header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])       #Name of extracted gene
uniqGene=unique(geneNames)        #Genes taken unique
geneList[["LASSO"]]=uniqGene

#Read the SVM results file
rt=read.table("SVM-RFE.gene.txt", header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])       #Name of extracted gene
uniqGene=unique(geneNames)        #Genes taken unique
geneList[["SVM-RFE"]]=uniqGene

#Plotting venn diagram
mycol=c("blue2","red2")
pdf(file="venn.pdf", width=5, height=5)
venn(geneList,col=mycol[1:length(geneList)],zcolor=mycol[1:length(geneList)],box=F,ilabels=F)
dev.off()

#Preservation of intersecting genes
intersectGenes=Reduce(intersect,geneList)
write.table(file=outFile, intersectGenes, sep="\t", quote=F, col.names=F, row.names=F)

