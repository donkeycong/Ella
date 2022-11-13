#install.packages("glmnet")
#install.packages("pROC")


#Citation package
library(glmnet)
library(pROC)

expFile="diffGeneExp.txt"      #Expression data files
geneFile="interGenes.txt"      #List file of intersecting genes
setwd("16.ROC")    #Setting up a working directory

#Reading input files
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)

#Obtain grouping information for samples
y=gsub("(.*)\\_(.*)", "\\2", colnames(rt))
y=ifelse(y=="con", 0, 1)

#Read the list file of genes
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)

#Loop on intersecting genes and plot ROC curves
bioCol=rainbow(nrow(geneRT), s=0.9, v=0.9)    #Defining the colour of a graphic
aucText=c()
k=0
for(x in as.vector(geneRT[,1])){
  k=k+1
  #Plotting the ROC curve
  roc1=roc(y, as.numeric(rt[x,]))     #Obtain the parameters of the ROC curve
  if(k==1){
    pdf(file="ROC.genes.pdf", width=8, height=8)
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="")
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }else{
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", add=TRUE)
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }
}
#Plot the legend to obtain the area under the ROC curve
legend("bottomright", aucText, lwd=2, bty="n", col=bioCol[1:(ncol(rt)-1)])
dev.off()

#Building a logic model
rt=rt[as.vector(geneRT[,1]),]
rt=as.data.frame(t(rt))
logit=glm(y ~ ., family=binomial(link='logit'), data=rt)
pred=predict(logit, newx=rt)     #Get the scoring of the model

#Plotting the model's ROC curve
roc1=roc(y, as.numeric(pred))      #Obtain the parameters of the model ROC curve
ci1=ci.auc(roc1, method="bootstrap")     #To obtain the fluctuation range of the area under the ROC curve
ciVec=as.numeric(ci1)
pdf(file="ROC.model.pdf", width=8, height=8)
plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main="Model")
text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
dev.off()