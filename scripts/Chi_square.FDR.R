args = commandArgs(trailingOnly=TRUE)
library(rstatix)

input <- args[1]
output <- args[2]
cutoff <- args[3]

MSdata <- read.table(input, header=TRUE)
SubData <- MSdata[,4:8]
pval <- c()
for(i in 1:nrow(SubData))
{
	if(apply(SubData[i,],1,sum) > 0)
	{
		Sub <- as.vector(t(SubData[i,]))
		if(sd(Sub) <= cutoff)
		{
			pval[i] <- "NA"
		}
		else
		{
			a <- chisq_test(SubData[i,],simulate.p.value=TRUE)
			pval[i] <- a$p
		}
	}
	else
	{
		pval[i] <- "NA"
	}
	 
}
pval <- as.numeric(pval)
padj.fdr <- p.adjust(pval, method="fdr")	###	fdr == BH

MSdata$pval <- pval
MSdata$FDR <- padj.fdr
sig <- MSdata %>% filter(FDR != "NA" & FDR < 0.05)
write.table(sig,output,sep="\t",quote=FALSE,na="NA",row.names=FALSE)
