#install.packages('e1071')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore")


inputFile="normalize.txt"       #Input files
setwd("10.CIBERSORT")       #Setting up a working directory
source("10.source-CIBERSORT.R")        #Citation package

#Immune cell infiltration analysis
outTab=CIBERSORT("ref.txt", inputFile, perm=1000, QN=TRUE)

#Filtering of immune infiltration results and preservation of immune cell infiltration results
outTab=outTab[outTab[,"P-value"]<0.05,]
outTab=as.matrix(outTab[,1:(ncol(outTab)-3)])
outTab=rbind(id=colnames(outTab),outTab)
write.table(outTab, file="CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)


