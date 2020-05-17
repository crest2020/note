library(Biobase)
library(genefu)
data(nkis)
nkiES = ExpressionSet(t(data.nkis), phenoData=AnnotatedDataFrame(demo.nkis),
                      featureData=AnnotatedDataFrame(annot.nkis))
nkiES

head(data.nkis)
sum(duplicated(featureNames(nkiES)))
names(fData(nkiES))
nrow(fData(nkiES))
length(fData(nkiES)$EntrzGene.ID)
head(fData(nkiES))
head(annot.nkis)

length(annot.nkis$EntrezGene.ID)

sum(is.na(annot.nkis$EntrezGene.ID))    
x<-annot.nkis$EntrezGene.ID
x<-x[!is.na(x)]
sum(duplicated(x))
head(x)
nrow(annot.nkis)
y<-c(1,2,3,3,3,3)


x[which(duplicated(x)==TRUE)]
duplicated(x)
data.nkis[1,1]
?mad
datanew<-t(data.nkis)

yy<-apply(datanew,1,mad)

min(yy)

which.min(yy)

which.min(apply(exprs(nkiES),1,mad))

?which.min
install.packages("impute")
library(BiocInstaller)
biocLite("impute")

library(impute)
ied = impute.knn(exprs(nkiES))
mean(is.na(ied$data))

which.min(apply(ied$data,1,mad))

exprs(nkiES) = ied
?impute.knn

library(BiocInstaller)
biocLite("IRanges")
library("IRanges")
ir<-IRanges(5,10)
ir

browseVignettes("IRanges")
seq(10, 0.001, length = 5)
rep(0.001, 4)
?rpois
?rle
start(ir)
end(ir)
width(ir)
IRanges(5,width=6)
IRanges(start=c(3,5,17),end=c(10,8,20))
?"intra-range-methods"
shift(ir,-2)
narrow(ir,start=2)
narrow(ir,end=5)
ir
flank(ir,width=3,start=TRUE,both=TRUE)
ir
plotir(ir,1)
irr<-IRanges(start=c(3,5,17),end=c(10,8,20))
range(irr)
irr
reduce(irr)
?"inter-range-methods"
disjoin(ir)
disjoin(irr)
library(BiocInstaller)
biocVersion()

library(devtools)
install_github("genomicsclass/ERBS")

library(ERBS)

data(HepG2)
class(HepG2)
length(HepG2)
HepG2
median(HepG2$signalValue)
?mcols
median(mcols(HepG2)$signalValue)
max(mcols(HepG2)$signalValue)
which(mcols(HepG2)$signalValue==max(mcols(HepG2)$signalValue))
mcols(HepG2)$seqnames[120]
HepG2$seqnames[2]
values(HepG2)
HepG2[120]
chr<-seqnames(HepG2)
table(chr)

median(ranges(HepG2))
median(width(HepG2))
hist(width(HepG2),nclass=25)
?hist
hist(width(HepG2))
m<-IRanges(101,200)
m*2
m
narrow(m,start=20)
m+25
n<-IRanges(c(1,11,21),c(3,15,27))
n
width(n)
sum(width(n))
x<-IRanges(c(101,106,201,211,221,301,306,311,351,361,401,411,501),c(150,160,210,270,225,310,310,330,390,380,415,470,510))
gaps(x)
sum(width(gaps(x)))
disjoin(x)
par(mfrow=c(2,1))
plotRanges(x)
??plorRanges
?resize
resize(x,1)
Plotranges(x)
library(GenomicRanges)
library(ph525x)
?plotranges

gr<-GRanges("chrZ",IRanges(start=c(5,10),end=c(35,45)),strand="+",seqlengths = c(chrZ=100L))
gr
shift(gr,10)
shift(gr,80)
mcols(gr)
mcols(gr)$value<-c(1,7)
gr
gr2<-GRanges("chrZ",IRanges(11:13,51:53))
mcols(gr)$value<-NULL
gr1<-GRangesList(gr,gr2)
gr1

gr3<-GRanges("chrZ",IRanges(c(1,11,21,31,41),width=5))
gr4<-GRanges("chrZ",IRanges(c(19,33),c(38,35)))

gr3
gr4

ir
fo<-findOverlaps(gr3,gr4)
fo
gr3%over%gr4
gr3[gr3%over%gr4]
?rle
r<-Rle(c(1,1,1,0,0,-2,-2,-2,rep(-1,20)))
r
rep(-1,20)
Views(r,start = c(4,2),end=c(7,6))
?Views

library(GenomicRanges)
x = GRanges("chr1", IRanges(c(101,201,401,501),c(150,250,450,550)), strand="+")

y = GRanges("chr1", IRanges(c(101,221,301,401,541),c(150,250,350,470,550)), strand="+")
library(ph525x)
par(mfrow=c(2,1))
plotGRanges = function(x) plotRanges(ranges(x))
GRangesList(x,y)
plotGRanges(x)
plotGRanges(y)

c(x,y)
x %over% y
x
y
disjoin(c(x,y))
plotGRanges(c(x,y))
z = GRanges("chr1", IRanges(c(101,201,401,501),c(150,250,450,550)), strand="-")

x%over%z
irr
findOverlaps(x,y)
library(ERBS)
data(HepG2)
data(GM12878)
HepG2[17,]
?distanceToNearest()
index<-which.min(distanceToNearest(GM12878,HepG2[17,])$distance)
mm<-distanceToNearest(GM12878,HepG2[17,])
mm
mm$value
mm$distance
index<-which.min(mcols(mm)$distance)
index
GM12878[945,]
queryHits(mm)

d = distanceToNearest(HepG2[17],GM12878)
dd<-distanceToNearest(HepG2,GM12878)
head(dd)
dis<-mcols(dd)$distance
mean(dis<2000)
library(Homo.sapiens)
source("https://bioconductor.org/biocLite.R")
biocLite("Rattus.norvegicus")
library(Rattus.norvegicus)
ghs<-genes(Rattus.norvegicus)
ghs
?precede
library(rtracklayer)
ghshuman<-genes(Homo.sapiens)
ghshuman
ghs
?genes
library(BiocInstaller)
biocLite("rtrackerlayer")

source("https://bioconductor.org/biocLite.R")
biocLite("rtracklayer")
library(rtracklayer)

library(ph525x)
?genes
library(Rattus.norvegicus)
rat<-genes(Rattus.norvegicus)
rat
library(Homo.sapiens)
ghshuman<-genes(Homo.sapiens)
ghshuman

which.max(mcols(ghshuman)$GENEID)
hist(log10(width(ghshuman)))
qqnorm(width(ghshuman))
abline(0,1)
GRanges(ghshuman)
library(ERBS)
erbs<-HepG2

median(width(ghshuman))
