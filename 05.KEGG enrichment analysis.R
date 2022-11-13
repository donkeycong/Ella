
#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#Citation package
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05       #p-value filter condition
qvalueFilter=0.05       #Corrected p-value filter condition

#Defining colours
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}
	
setwd("05.KEGG")          #Setting up a working directory
rt=read.table("diff.txt", header=T, sep="\t", check.names=F)     #Reading input files

#Gene name to gene id conversion
genes=as.vector(rt[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #Removal of the gene with gene id NA

#kegg enrichment analysis
kk=enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$id[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
#Preservation of enrichment results
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#Define the number of display pathways
showNum=30
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}

#Bar chart
pdf(file="barplot.pdf", width=8, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, color=colorSel)
dev.off()

#Air bubble chart
pdf(file="bubble.pdf", width = 8, height = 7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", color=colorSel)
dev.off()

