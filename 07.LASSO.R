
#install.packages("glmnet")
  
  
  set.seed(123)
library(glmnet)                   #Citation package
inputFile="diffGeneExp.txt"       #Input files
setwd("07.lasso")       #Setting up a working directory

#Reading input files
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
rt=t(rt)

#Building the model
x=as.matrix(rt)
y=gsub("(.*)\\_(.*)", "\\2", row.names(rt))
fit=glmnet(x, y, family = "binomial", alpha=1)
cvfit=cv.glmnet(x, y, family="binomial", alpha=1,type.measure='deviance',nfolds = 10)
pdf(file="cvfit.pdf",width=6,height=5.5)
plot(cvfit)
dev.off()

#Export of screened hub genes
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
lassoGene=row.names(coef)[index]
lassoGene=lassoGene[-1]
write.table(lassoGene, file="LASSO.gene.txt", sep="\t", quote=F, row.names=F, col.names=F)

