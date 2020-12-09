args = commandArgs(trailingOnly=TRUE)
library(ggplot2)

input <- args[1]
output <- args[2]

MSdata <- read.table(input, header=T, col.names=c("Genomes","Genes","Class"))
MSdata$Genomes=as.factor(MSdata$Genomes)
p <- ggplot(MSdata, aes(x=Genomes, y=Genes, fill=Class)) + geom_boxplot(position="identity",alpha=.5) + theme_bw()
svg(output)
print(p)
dev.off()
