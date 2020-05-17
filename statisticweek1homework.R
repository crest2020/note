library(gplots)
library(devtools)
library(Biobase)
library(RSkittleBrewer)
library(org.Hs.eg.db)
library(AnnotationDbi)

install.packages("RSkittleBrewer")
library(BiocInstaller)
biocLite("alyssafrazee/RSkittleBrewer")

detach("package:AnnotationHub",unload=TRUE)

library(RSkittleBrewer)
# Make the colors pretty
trop = RSkittleBrewer("tropical")
palette(trop)
par(pch=19)
trop

con = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata=pData(bm)
edata=exprs(bm)
fdata = fData(bm)
ls()
con
bm
gene_na=rowSums(is.na(edata))
table(gene_na)

par(mfrow=c(1,2))
hist(log2(edata[,1]+1),col=2)
hist(log2(edata[,2]+1),col=2)
plot(density(log2(edata[,1]+1),col=2))
lines(density(log2(edata[,2]+1)),col=3)
lines(0,1,col=9)
?line
abline(0,1)

mm=log2(edata[,1]+1)-log2(edata[,2]+1)
aa=log2(edata[,1]+1)+log2(edata[,2]+1)
plot(aa,mm,col=2)
edata1=as.data.frame(edata)
filt_data=filter(edata1,rowMeans(edata1)>1)
dim(filt_data)
boxplot(log2(filt_data+1),col=2)
dim(edata)
boxplot(log2(edata+1),col=2)
head(filt_data)

install.packages("dendextend")
library(dendextend)

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata=pData(bm)
edata=as.data.frame(exprs(bm))
fdata = fData(bm)
ls()

edata=edata[rowMeans(edata)>5000,]
edata=log2(edata+1)
dim(edata)
head(edata,n=20)
dist1=dist(t(edata))
colramp=colorRampPalette(c(3,"white",2))(9)
heatmap(as.matrix(dist1),col=colramp,Colv = NA, Rowv = NA)

kmeans2=kmeans(edata,centers = 3)
kmeans3=kmeans(edata,centers=4)
table(kmeans2$cluster,kmeans3$cluster)


matplot(t(kmeans3$centers),col=1:3,lwd=3)
 kmeans3$centers
hclust1=hclust(dist1)
plot(hclust1)
plot(hclust1,hang=-1)
dend=as.dendrogram(hclust1)
dend=color_labels(hclust1,4,col=1:4)
plot(dend)
labels_colors(dend)=c(rep(1,10),rep(2,9))
plot(dend)

kmeans1=kmeans(edata,centers=3)
names(kmeans1)
?kmeans
matplot(t(kmeans1$centers),col=1:3,type="l",lwd = 3)
?matplot
 
table(kmeans1$cluster,kmeans2$cluster)
heatmap(as.matrix(edata)[order(kmeans1$cluster),],col=colramp,Colv=NA,Rowv = NA)


a=matrix(0,6,6)
a
rowSums(a)
