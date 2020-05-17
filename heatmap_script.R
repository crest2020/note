
getwd()
setwd("D://Rdocuments")
heat=read.csv("heat_tgf.csv",header = T)
heat
#install.packages("pheatmap")
library(pheatmap)
#?pheatmap
rownames(heat)=heat[,1]
annotation_col=data.frame(group=factor(c(rep("WT",3),rep("KO",3))))
annotation_col
rownames(annotation_col)=colnames(heat[,-1])
annotation_col
pheatmap(heat[,2:7],scale = "row",cluster_row = F,cluster_cols = F,show_colnames = F,annotation_col = annotation_col)


