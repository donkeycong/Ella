#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")
#install.packages("ggExtra")


#Citation package
library(limma)
library(reshape2)
library(ggpubr)
library(ggExtra)

gene="STK3"                       #Gene Name
expFile="normalize.txt"             #Expression data files
immFile="CIBERSORT-Results.txt"     #Immune cell infiltration results document
setwd("13.immuneCor")     #Setting up a working directory

#Read gene expression files and process the data
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#Obtain target gene expression
data=t(data[gene,,drop=F])
data=as.data.frame(data)

#Read immune cell result files and collate data
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)

#Data merge
sameSample=intersect(row.names(immune), row.names(data))
rt=cbind(immune[sameSample,,drop=F], data[sameSample,,drop=F])

#Circulating immune cells and plotting correlation scatter plots
outTab=data.frame()
for(i in colnames(rt)[1:(ncol(rt)-1)]){
	x=as.numeric(rt[,gene])
	y=as.numeric(rt[,i])
	if(sd(y)==0){y[1]=0.00001}
	cor=cor.test(x, y, method="spearma")
	
	outVector=cbind(Gene=gene, Cell=i, cor=cor$estimate, pvalue=cor$p.value)
	outTab=rbind(outTab,outVector)
	
	if(cor$p.value<0.05){
		outFile=paste0("cor.", i, ".pdf")
		df1=as.data.frame(cbind(x,y))
		p1=ggplot(df1, aes(x, y)) + 
				  xlab(paste0(gene, " expression")) + ylab(i)+
				  geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
				  stat_cor(method = 'spearman', aes(x =x, y =y))
		p2=ggMarginal(p1, type="density", xparams=list(fill = "orange"), yparams=list(fill = "blue"))

#Correlation graphics
		pdf(file=outFile, width=5.2, height=5)
		print(p2)
		dev.off()
	}
}
#Export of immunisation function and p-value table files
write.table(outTab,file="cor.result.txt",sep="\t",row.names=F,quote=F)

