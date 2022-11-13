#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#Citation package
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

inputFile="all.txt"         #Input files
gmtFile="c2.cp.kegg.v7.4.symbols.gmt"      #Gene set files
setwd("06.GSEA")      #Setting up a working directory

#Reading files and sorting input files
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
logFC=as.vector(rt[,2])
names(logFC)=as.vector(rt[,1])
logFC=sort(logFC, decreasing=T)

#Read in gene set files
gmt=read.gmt(gmtFile)

#Enrichment analysis
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result.txt",sep="\t",quote=F,row.names = F)

#Output the graph of the experimental group enrichment
termNum=5      #Number of display pathways
kkUp=kkTab[kkTab$NES>0,]
if(nrow(kkUp)>=termNum){
	showTerm=row.names(kkUp)[1:termNum]
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in Treat")
	pdf(file="GSEA.treat.pdf", width=7, height=5.5)
	print(gseaplot)
	dev.off()
}

#Output of normal group enrichment graphs
termNum=5      #Number of display pathways
kkDown=kkTab[kkTab$NES<0,]
if(nrow(kkDown)>=termNum){
	showTerm=row.names(kkDown)[1:termNum]
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in Control")
	pdf(file="GSEA.con.pdf", width=7, height=5.5)
	print(gseaplot)
	dev.off()
}


