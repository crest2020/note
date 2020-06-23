
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

t <- t(scale(t(t)))
p1 <- pheatmap(t,
 annotation_col = col,
 annotation_legend = T,
 border_color = 'black',#设定每个格子边框的颜色，border=F则无边框
 cluster_rows = T, #对行聚类
 cluster_cols = T, #队列聚类
 show_colnames = F, #是否显示列名
 show_rownames = F #是否显示行名
)

table(abs(t)>2)
t[t>=2]=2
t[t<=-2]=-2
table(is.na(t))
p2 <- pheatmap(t,
 annotation_col = col,
 annotation_legend = T,
 border_color = 'black',#设定每个格子边框的颜色，border=F则无边框
 cluster_rows = T, #对行聚类
 cluster_cols = T, #队列聚类
 show_colnames = F, #是否显示列名
 show_rownames = F #是否显示行名
)




