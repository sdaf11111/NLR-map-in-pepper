library(dplyr)
library(ggtree)
library(ggplot2)
library(svglite)
library(scales)

setwd("scripts/Figure2.Data")
info <- read.csv("Annuum.Intact.NBARC.group.trimal92.csv")
tree <- read.tree("Annuum.Intact.NBARC.tree.trimal92.nwk")

###	Assign group information for each gene	###
IDgroup <- select(info, Group)
row.names(IDgroup) <- info$ID
groupInfo <- split(row.names(IDgroup), IDgroup$Group)
tree<- groupOTU(tree, groupInfo)

###	find node for MRCA	###
a <- c()
for (i in groupInfo)
{
	a=c(a,MRCA(tree,i))
}
b <- c()
for( j in a)
{
	for(k in a)
	{
		b=c(b, MRCA(tree, j, k))
	}
}
d<-as.numeric(unique(c(a,b)))

###	Color code	###
heatmap.colours <- c("#be9fe1","#8ac6d1","#e1ccec","#fddb3a",
                    "#C0C0C0","#c9b6e4","#d5c455","#ffb6b9","#fae3d9",
                    "#9aceff","#d7cde6","#bbded6","#ede59a","#4f98ca",
                    "#4a69bb","#f5cdaa","NA",
	        "#c3d14a","#63b637","#FF0000","#008000","#FF1493","#FF4500")
names(heatmap.colours) <- c("G1","G2","G3","G4",
		"G5","G6","G7","G8","G9",
		"G10","G11","G12","G13","GT",
		"GR","G14","NG",
		"CHIL","Know","CANN","CECW","CZUN","CASF")
###	Plot	###
p <- ggtree(tree, layout='circular', size=0.2) %<+% info +	### ggtree(tree, layout='circular', aes(color=group))
  geom_aline(linetype="solid", size=0.5, aes(color=group), alpha=0.5) +	
  scale_colour_manual(values=heatmap.colours, breaks=c("G1","G2","G3","G4","G5","G6","G7","G8","G9","G10","G11","G12","G13","GT","GR","G14","CHIL","Know","CANN","CECW","CZUN","CASF"), name="Group") +
  geom_tippoint(aes(color=Species), size=0.2)+ theme(legend.position="bottom") +
  geom_nodepoint(aes(subset = node %in% node[d] & as.numeric(label) >= 90), fill="red", shape=23, size=1, na.rm = TRUE)
q <- flip(p, 3001, 4090) %>% rotate(3507)

ggsave(file="NLR_tree.92.svg", plot=q)

