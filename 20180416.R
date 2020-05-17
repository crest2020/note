library(limma)
library(CLL)
library(BiocInstaller)
biocLite("CLL")
data(package="CLL")
data("sCLLex")
exprset=exprs(sCLLex)

head(exprset)
samples=sampleNames(sCLLex)
pdata=pData(sCLLex)
head(pdata)
gruop_list=as.character(pdata[,2])
head(gruop_list)
table(gruop_list)
dim(exprset)
exprset[1:5,1:5]
par(cex=.7)
n.sample=ncol(exprset)
n.sample
cols=rainbow(n.sample*1.2)
boxplot(exprset,col=cols,mian="expression value",las=2)
design=model.matrix(~0+factor(gruop_list))
head(design)
colnames(design)=levels(factor(gruop_list))
rownames(design)=colnames(exprset)
design
contrast.matrix=makeContrasts(paste0(unique(gruop_list),collapse = "-"),levels=design)
contrast.matrix
fit=lmFit(exprset,design)
fit
fit2=contrasts.fit(fit,contrast.matrix)
fit2=eBayes(fit2)
tempoutput=topTable(fit2,coef = 1,n=Inf)
nrdeg=na.omit(tempoutput)
head(nrdeg)
bb=c(rep(0,5),rep(1,5))
bb
dd=bb<1
dd
mean(dd)
mean(dd[1:7])
dat=read.csv("mice_pheno.csv")
dat
head(dat)
hf=dat[dat$Diet=="hf"&dat$Sex=="F",]
hf
tail(dat)
table(dat$Diet)
chow=dat[dat$Diet=="chow"&dat$Sex=="F",]
chow
head(hf)
head(chow)
chow=na.omit(chow)
hf=na.omit(hf)
mean(chow$Bodyweight)
mean(hf$Bodyweight)
head(na.omit(chow))
mean(hf$Bodyweight,na.rm = T)
head(is.na(hf))
table(is.na(hf))
head(na.omit(hf)
head(hf)
hf
hf
chow=na.omit(chow)
chow
hf=na.omit(hf)
t.test(sample(chow$Bodyweight,20),sample(hf$Bodyweight,20))$p.value

chow$Bodyweight
na.omit(chow)
?na.omit
reject=function(N,alpha=0.05){
  pval=t.test(sample(chow$Bodyweight,N),sample(hf$Bodyweight,N))$p.value
  pval<alpha
}
reject(10)
reject(5)
replicate(10,reject(5))
Ns=seq(5,50,5)
N
rejections=replicate(10,reject(N))
rejections
power=sapply(Ns,function(N) {
  rejections=replicate(1000,reject(N))
  mean(rejections)
})
power
power
plot(Ns,power,type = "b")
Ns=30
alpha=c(0.1,0.05,0.01,0.001,0.0001)
power=sapply(alpha,function(alpha) {
  rejections=replicate(1000,reject(Ns,alpha))
  mean(rejections)
})
plot(alpha,power,type="b",log = "x")
t.test(sample(chow$Bodyweight,N),sample(hf$Bodyweight,N))$conf.int/mean(chow$Bodyweight)*100
?t.test
pp=t.test(sample(chow$Bodyweight,N),sample(hf$Bodyweight,N))
confint(pp)
?confint
library(limma)
?voom
?rlm
?lm
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb=TxDb.Hsapiens.UCSC.hg19.knownGene
txdb
gr=GRanges(seqnames = "chr1",strand = "+",ranges = IRanges(start=11874,end=14409))
gr
micro=microRNAs(txdb)
subsetByOverlaps(microRNAs(txdb),gr)
subsetByOverlaps(genes(txdb),gr,ignore.strand=T)
subsetByOverlaps(transcripts(txdb),gr)
subsetByOverlaps(exons(txdb),gr)
subsetByOverlaps(exonsBy(txdb,by="tx"),gr)

subsetByOverlaps(cds(txdb),gr)
subsetByOverlaps(cdsBy(txdb,by="tx"),gr)
subset(transcriptLengths(txdb,with.cds_len = T),gene_id=="100287102")
mean(width(micro))
head(width(micro),n=100)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
head(select(org.Hs.eg.db,keys = "ORMDL3",keytype = "SYMBOL",columns = c("PMID","GENENAME")))
library(GO.db)
GO.db
help('select')
k5=keys(GO.db)[1:5]
k5
cgo=columns(GO.db)
head(cgo)
select(GO.db,keys=k5,columns = cgo[1:3])
con=GO_dbconn()
dbListTables(con)
dbGetQuery(con,'select _id,go_id, term from go_term limit 5')
dbGetQuery(con,"select * from go_bp_parents where _id=27")
dbGetQuery(con,"select _id, go_id, term from go_term where _id=5938")
dbGetQuery(con,"select * from go_bp_parents where _id=5938")
library(BiocInstaller)
biocLite("Mus.musculus")
library(KEGGREST)
brca2k=keggGet("hsa:675")
names(brca2k[[1]])
brpat=keggGet("path:hsa05212")
names(brpat[[1]])
brpat[[1]]$GENE[seq(1,132,2)]
library(png)
library(grid)
brpng=keggGet("hsa05212","image")
grid.raster(brpng)
library(Homo.sapiens)
library()
class(Homo.sapiens)
Homo.sapiens
tx=transcripts(Homo.sapiens)
head(tx)
tx
keytypes(Homo.sapiens)
columns(Homo.sapiens)
head(select(Homo.sapiens,keys="TP53",keytype = "SYMBOL",columns = "OMIM"))

p=10^-7
N=5*10^6
winners=rbinom(1000,N,p)
head(winners)
sum(winners==1)
tab=table(winners)
tab
length(winners)
N*p
plot(tab)
prop.table(tab)
N=100000
lambdas=2^seq(1,16,len=N)
head(lambdas)
y=rpois(N,lambdas)
x=rpois(N,lambdas)
ind=which(y>0&x>0)
head(ind)
library(rafalib)
head(y)
y=y[ind]
x=x[ind]
  plot(log2(lambdas[ind]),log2(y/x))
?splot
library(parathyroidSE)
data("parathyroidGenesSE")
se=parathyroidGenesSE
head(se)
exprs(se)
x=assay(se)[,23]
y=assay(se)[,24]
ind=which(y>0&x>0)
splot(log2(x)+log2(y),log2(x/y),subset=ind)
library(matrixStats)
vars=rowVars(assay(se)[,c(2,8,16,21)])
means=rowMeans(assay(se)[,c(2,8,16,21)])
splot(means,vars,log="xy",subset=which(means>0&vars>0))
abline(0,1,col=2,lwd=2)
datadir="http://www.biostat.jhsph.edu/bstcourse/bio751/data"
x=read.csv(file.path(datadir,"hcmv.csv"))[,2]
head(x)
breaks=seq(0,4000*round(max(x)/4000),4000)
head(breaks)
length(breaks)
?seq
head(4000*round(max(x)/4000))
seq(0,10,2)
tmp=cut(x,breaks)
tmp
length(breaks)
head(x)
length(x)
?cut
head(breaks)
length(tmp)
length(breaks)
cut(rep(1,5), 4)
tmp

counts=table(tmp)
hist(counts)
l=function(lambda) sum(dpois(counts,lambda,log=T))
lamdas=seq(3,7,len=100)
ls=exp(sapply(lamdas,l))
head(ls)
?dpois
plot(lamdas,ls,type="l")
mle=optimize(l,c(0,10),maximum = T)
abline(v=mle$maximum)
print(c(mle$maximum,mean(counts)))
theoretical=qpois((seq(0,99)+0.5)/100,mean(counts))
qqplot(theoretical,counts)
abline(0,1)
counts
plot(counts)
?dpois
x=seq(0,10)
dpois(x,0.01)
ppis(1,0.01)
pp
ppois(1,0.01)
qpois(0.9,1)
library(Biobase)
library(maPooling)
data("maPooling")
pd=pData(maPooling)
strain=factor(as.numeric(grepl("b",rownames(pd))))
table(strain)
grep("b",rownames(pd))
pooled=which(rowSums(pd)==12&strain==1)
pooled
techreps=exprs(maPooling[,pooled])
individual=which(rowSums(pd)==1&strain==1)
techreps
individual
individual=individual[-grep("tr",names(individual))]
bioreps=exprs(maPooling)[,individual]
library(matrixStats)
techsds=rowSds(techreps)
biosds=rowSds(bioreps)
library(rafalib)
shist(biosds,unit=0.1,col=1,xlim=c(0,1.5))
shist(techsds,unit=0.1,col=2,add=T)
legend("topright",c("Biological","Technical"),col=c(1,2),lty=c(1,1))
qqnorm(biosds)
qqline(biosds)
library(rafalib)
mypar(3,3)
sds=seq(0,2,len=100)
hist(sds)
for(d in c(1,5,10)){
  for(s0 in c(0.1,0.2,0.3)){
    tmp=hist(biosds,main=paste("s_0=",s0,"d=",d),xlab = "sd",ylab="density",freq = F,nc=100,xlim = c(0,1))
    dd=df(sds^2/s0^2,11,d)
    k=sum(tmp$density)/sum(dd)
    lines(sds,dd*k,type="l",col=2,lwd=2)
  }
}
dd

?df
k
library(limma)
estimates=fitFDist(biosds^2,11)
theoretical=sqrt(qf((seq(0,999)+0.5)/1000,11,estimates$df2)*estimates$scale)
observed=biosds
mypar(1,2)
qqplot(theoretical,observed)
abline(0,1)
seq(1:10)
seq(along=1:10)
seq(1,10,by=2)
?seq
seq_len(1:10)
x=1:10
z=NULL
for(i in seq(1:10)){
  if(x[i]<5){
    z=c(z,x[i]-1)
    
  }
  else {
    z=c(z,x[i]/x[i])
  }
}
z
z=NULL
z
for(i in seq(x)){
  if(x[i]<5){
    z=c(z,x[i]-1)
  }else{
    stop("value need to be <5")
  }
}
z
z=0
while(z<5){
  z=z+2
  print(z)
}
apply(iris[,1:3],1,mean)
x
test=function(x){
  if(x<5){
    x-1
  }
  else{
    x/x
  }
}
apply(as.matrix(x),1,test)
as.matrix(x)
apply(as.matrix(x),1,function(x){if(x<5){x-1}else{x/x}})
table(factor(iris[,5]))
tapply((iris[,4]),factor(iris[,5]),mean)
aggregate(iris[,1:4],list(iris$Species),mean)
mylist=as.list(iris[1:3,1:3])
mylist
install.packages("TeachBayes")
library(TeachBayes)
data(package="TeachBayes")
??TeachBayes
?table
?table
?table
d=c(1,1,2,4,5)
table(d)
prob.table(d)
prop.table(d)
?

plot(x<-rnorm(40,2e+07,sd=1e+07),y=rep(1,times=40),type="h",col="blue")
rnorm(40,2e+07,sd=1e+07)
y=rep(1,times=40)
y
abline(h=0.78,col="green",lwd=12)
lines(a<-rnorm(5,2e+07,sd=1e+07),b<-rep(1,times=5),type="h",col="red",lwd=2)
rnorm(5,2e+07,sd=1e+07)
plot(1:10)
lines(1:10,lty=1)
plot(x=1:10,y=1:10)
abline(-1,1,col="green")
abline(1,1,col="lightgray")
abline(v=5,col="brown")
library(grid)
library(ggplot2)
dsmall=diamonds[sample(nrow(diamonds),1000),]
sample(90,10)
?sample
a=ggplot(dsmall,aes(color,price/carat,color=color))+geom_jitter(size=4,alpha=I(1/1.5))
a
?geom_jitter
b=ggplot(dsmall,aes(color,price/carat,color=color))+geom_boxplot()
b
c=ggplot(dsmall,aes(color,price/carat,fill=color))+geom_boxplot()+theme(legend.position = "none")
      ?labs
?theme()
c

y <- as.data.frame(matrix(runif(30), ncol=3, dimnames=list(letters[1:10], LETTERS[1:3])))
y
cat("palette():\n")
palette()
palette()
palette(rainbow(20,start=0.1,end=0.2))
par(mfrow=c(1,1,2,2))
barplot(as.matrix(y[1:8,]),beside=T,col=1:8)
y
plot(1:10,col=1:10)
?mfrow
layout(matrix(c(1,1,2,3),2,2,byrow = T))
library(rafalib)
library(tissuesGeneExpression)
data(tissuesGeneExpression)
colind=tissue%in%c("kidney","colon","liver")
data(package="tissuesGeneExpression")
mat=e[,colind]
head(mat)
dim(mat)
group=factor(tissue[colind])
head(colind)
head(group)
table(group)
s=svd(mat-rowMeans(mat))
PC1=s$d[1]*s$v[,1]
PC2=s$d[2]*s$v[,2]
mypar(1,1)
plot(PC1,PC2,pch=21,bg=as.numeric(group))
?plot
legend("bottomright",levels(group),col=seq(along=levels(group)),pch=15,cex=1.5)
library(tissuesGeneExpression)
data(tissuesGeneExpression)
d=dist(t(e))
library(rafalib)
mypar()
hc=hclust(d)
hc
plot(hc,labels=tissue,cex=0.5)
table(tissue)
head(e)
myplclust(hc,labels = tissue,lab.col = as.fumeric(tissue),cex=0.5)
abline(h=120)
hclusters=cutree(hc,h=120)
table(true=tissue,cluster=hclusters)
hckusters=cutree(hc,k=8)
table(true=tissue,cluster=hckusters)
set.seed(1)
km=kmeans(t(e[1:2,]),centers = 7)
names(km)
mypar(1,2)
plot(e[1,],e[2,],col=as.fumeric(tissue),pch=16)
plot(e[1,],e[2,],col=km$cluster,pch=16)
table(true=tissue,cluster=km$cluster)
km=kmeans(t(e),centers = 7)
km
head(d)
mds=cmdscale(d)
mypar(1,2)
plot(mds[,1],mds[,2],col=km$cluster,pch=16)
table(true=tissue,cluster=km$cluster)
library(RColorBrewer)
hmcol=colorRampPalette(brewer.pal(9,"GnBu"))(100)
head(hmcol)
library(genefilter)
rv=rowVars(e)
idx=order(-rv)[1:40]
library(gplots)
cols=palette(brewer.pal(8,"Dark2"))[as.fumeric(tissue)]
head(cbind(colnames(e),cols))
heatmap.2(e[idx,],labCol = tissue,trace = "none",col = hmcol)
heatmap.2(e[idx,],labCol = tissue,trace="none")
ee=e[,1:30]
heatmap.2(ee[idx,],trace="none")
head(tissue)
head(e,n=2)
table(tissue)
library(tissuesGeneExpression)
data(package="tissuesGeneExpression")
head(tab)
table(colnames(e)==tab$filename)
table(tissue==tab$Tissue)
library(Biobase)
library(SpikeIn)
install.packages("SpikeIn")
library(BiocInstaller)
biocLite("hgu95acdf")
n
library(SpikeIn)
library(hgu95acdf)
data("SpikeIn95")
mypar()
plot(X,Y,data=SpikeIn95)
head(SpikeIn95)
SpikeIn95
head(assayData(SpikeIn95))
?maplot
n <- 10000
signal <- runif(n,4,15)
bias <- (signal/5 - 2)^2
x <- signal + rnorm(n)
y <- signal + bias + rnorm(n)
maplot(x,y)
head(signal)
library(AnnotationDbi)
library(org.Hs.eg.db)
columns(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
kk=head(keys(org.Hs.eg.db,keytype = "ENTREZID"))
res=mapIds(org.Hs.eg.db,keys=kk,column = "ENSEMBL",keytype = "ENTREZID")
res
k
dim(res)
length(res)
res=select(org.Hs.eg.db,keys=kk,columns=c("ENTREZID","ENSEMBL","SYMBOL"),keytype = "ENTREZID")
res
idx=match(kk,res$ENTREZID)
idx
res[idx,]
library(biomaRt)
m=useMart("ensembl",dataset = "hsapiens_gene_ensembl")
map=getBM(mart=m,attributes = c("ensembl_gene_id","entrezgene"),filters = "ensembl_gene_id",values = k)
map
k
library(GenomicRanges)
z=GRanges("chr1",IRanges(1000001,1001000),strand = "+")
library(BSgenome.Hsapiens.UCSC.hg19)
granges
dnastringset=getSeq(Hsapiens,z)
z
head(dnastringset)
library(Biostrings)
dnastringset=readDNAStringSet("transcripts.fa")
library(Rsamtools)
which=GRanges("chr1",IRanges(1000001,1001000))
which
what=c("rname","strand","pos","qwidth","seq")
param=ScanBamParam(which=which,what=what)
biocVersion()
library(BiocInstaller)
biocVali
)
biocValid()
?mean
?mad
example(mad)
help(package="genefilter",help_type = "html")
??genefilter
class(6)
library(Biobase)
?ExpressionSet
?"ExpressionSet-class"
methods(class="ExpressionSet")
methods(class="lm")
library(rafalib)
data("sample.ExpressionSet")
whatMethods(sample.ExpressionSet)
whatMethods(sample.ExpressionSet)
whatm
library(BSgenome.Hsapiens.UCSC.hg19)
Hsapiens$chr17

library(GSE5859Subset)
library(Biobase)
.oldls=ls()
data(GSE5859Subset)
.newls=ls()
head(.oldls)
newstuff=setdiff(.newls,.oldls)
newstuff
cls=lapply(newstuff,function(x) class(get(x)))

?get
get(newstuff)
class(get(newstuff[2]))
newstuff
head(get(newstuff))
head(geneAnnotation)
?apply
  ?lapply
names(cls)=newstuff
cls
newstuff
boxplot(date~factor(ethnicity),data=sampleInfo,ylab="Chip date")
all.equal(colnames(geneExpression),sampleInfo$filename)
all.equal(rownames(geneExpression),geneAnnotation$PROBEID)
ind=which(geneAnnotation$SYMBOL=="BRCA2")
head(ind)
boxplot(geneExpression[ind,]~sampleInfo$ethnicity,ylab="BRCA2 by hgfocus")
es1=exprs(geneExpression)
es1=ExpressionSet(geneExpression)
head(es1)
pData(es1)=sampleInfo
fData(es1)=geneAnnotation
es1
boxplot(es1$date~es1$ethnicity)
es1
head(geneExpression)
library(GSE5859)
data(GSE5859)
e
library(annotate)
p=pmid2MIAME("15269782")
p
getClass("ExpressionSet")
getClass("MIAME")
library(S4Vectors)
sa2=DataFrame(sampleInfo)
sa2
getClass("DataFrame")
library(RNAseqData.HNRNPC.bam.chr14)
bfp=RNAseqData.HNRNPC.bam.chr14_BAMFILES
length(bfp)
bfp[1]
library(Rsamtools)
bfl=BamFileList(file=bfp)
bfl
seqinfo(bfl)
hnrnpcloc=GRanges("chr14",IRanges(21677296,21737638))
hnrnpcloc
library(GenomicAlignments)
library(BiocParallel)
register(SerialParam())
hnse=summarizeOverlaps(hnrnpcloc,bfl)
hnse
head(assay(hnse))
dim(assay(hnse))
?Views
?sample
cc=c("A","T","C","G")
cc
aa=replicate(100,sample(cc,1))
aa
paste0(aa)
length(aa)
paste(c("a","b","c"),sep="")
letters(c(1:26))
letter(c(1:26))
letters[1:26]
LETTERS
paste(aa,collapse = "")
?toString
library(RNAseqData.HNRNPC.bam.chr14)
bfp=RNAseqData.HNRNPC.bam.chr14_BAMFILES
library(Rsamtools)
bfl=BamFileList(file=bfp)
hnrnpcloc=GRanges("chr14",IRanges(21677296,21737638))
library(GenomicAlignments)
library(BiocParallel)
register(SerialParam())
hnse=summarizeOverlaps(hnrnpcloc,bfl)
hnse
assay(hnse)
rowRanges(hnse)
seqinfo(hnse)
metadata(hnse)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb=TxDb.Hsapiens.UCSC.hg19.knownGene
gr14=genes(txdb,vals=list(tx_chrom="chr14"))
gr14
gg=genes(txdb)
keepSeqlevels(gg,"chr14")
library(org.Hs.eg.db)
gr14$symbol=mapIds(org.Hs.eg.db,keys=gr14$gene_id,keytype = "ENTREZID",column = "SYMBOL")
Homo.sapiens
library(BSgenome.Hsapiens.UCSC.hg19)
gr14
char=rep(c("hela_wt","hela_hkd"),each=4)
char
library(GenomicFiles)
bff=GenomicFiles(files=path(bfl))
?GenomicFiles
bff
colData(bff)$condition=char
sid=c(1,1,1,1,2,2,3,3)
sid
bff$sample=sid
bff
hnse=summarizeOverlaps(gr14[c(1:4,305)],files(bff))
gr14[c(1:4,305)]
bfp
length(bfp)
dim(bfp)
class(bfp)
colData(hnse)=cbind(colData(hnse),colData(bff))
hnse
assay(hnse)
par(mfrow=c(2,2))
for(i in 2:5){
  boxplot(assay(hnse)[i,]~hnse$condition,ylab=rowRanges(hnse)$symbol[i])
}\
library(ALL)
data(ALL)
allse=makeSummarizedExperimentFromExpressionSet(ALL)
allse
library(erma)
ef=dir(system.file("bed_tabix",package = "erma"),pattern  = "bed.gz$")
length(ef)
ef
mm=makeErmaSet()
mm
head(colData(mm))
names(files(mm))=colData(mm)$Epigenome.Mnemonic
colData(mm)$Epigenome.Mnemonic
library(GenomicFiles)
head(files(mm))
library(png)
library(grid)
im=readPNG(system.file("pngs/emparms.png",package="erma"))
grid.raster(im)
grid.r
stateProfile(mm,"BRCA2")
stateProfile(mm,"EOMES")
gm=promoters(range(genemodel("BRCA2")))
gm
library(BiocParallel)
register(SnowParam(workers=2))
library(GenomicFiles)
ans=reduceByFile(gm,files(mm),MAP=function(range,file){table(import(file,genome="hg19",which=range)$name)})
ans=unlist(ans,recursive = F)
names(ans)=colData(mm)$Epigenome.Mnemonic
ans[1:4]
library(rtracklayer)
ans = reduceByFile( gm, files(mm), MAP=function(range,file) {
  table( import(file, genome="hg19", which=range)$name ) } )
colData(mm)$Epigenome.Mnemonic
colData(mm)$Epigenome.Mnemonic
files(mm)
ans
?import
library(IRanges)
IRanges(start = c(3,5,17),end=c(10,8,20))
ir=IRanges(5,10)
ir
start(ir)
end(ir)
width(ir)
shift(ir,-2)
shift(ir,2)
ir
narrow(ir,start=2)
narrow(ir,end=4)
flank(ir,width=3,start=T,both=F)
flank(ir,width=3,start=F,both=F)
flank(ir,width=3,start = T,both=T)
flank(ir,width=3,start=T,both=T)
ir*2
ir*-2
ir

library(Biostrings)
library(biomaRt)
library(BSgenome)
ag=available.genomes()
?available.genomes
length(ag)
head(ag)
?grep
grep("TAIR",ag)
ag
library(BSgenome.Hsapiens.UCSC.hg19)
Hsapiens
genome(Hsapiens)
Hsapiens$chr17
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb=TxDb.Hsapiens.UCSC.hg19.knownGene
txdb
ghs=genes(txdb)
ghs
f1=dir(system.file("extdata",package = "ERBS"),full=T)[1]
readLines(f1,4)
library(rtracklayer)
imp=import(f1,format = "bedGraph")
imp
genome(imp)
export(imp,"denoex.bed")
cat(readLines("denoex.bed",n=5),sep="\n")
library(AnnotationHub)
ah=AnnotationHub()
ah
query(ah,"HepG2")
mean(sample(1:10,10))
mean(sample(1:10,10))
mean(sample(1:100,10))
x=rnorm(10)
ks.test(x,"pnorm")
set.seed(3)
x=runif(n=20,min=0,max=20)
y=runif(n=20,min=0,max=20)
x
y
plot(ecdf(x),do.points=F,verticals = T,xlim=c(0,20))
lines(ecdf(y),lty=3,do.points=F,verticals=T)
ks.test(x,y)
boxplot(mtcars$mpg~mtcars$am,ylab="mpg",names=c("automatic","manual"))
wilcox.test(mpg~am,data=mtcars)
ftable=table(mtcars$am,mtcars$gear)
ftable
mosaicplot(ftable,main="Number of Forward Gears within automatic and manual cars",color=T)
chisq.test(ftable)
boxplot(mtcars$mpg~factor(mtcars$gear),xlab="gear",ylab="mpg")
library(stats)
oneway.test(mtcars$mpg~factor(mtcars$gear))
mtcars.aov=aov(mtcars$mpg~as.factor(mtcars$gear))
summary(mtcars.aov)
model.tables(mtcars.aov,"means")
mtcars_posthoc=TukeyHSD(mtcars.aov)
mtcars_posthoc
plot(mtcars_posthoc)
par(mfrow=c(1,2))
library(stats)
boxplot(mtcars$mpg~mtcars$gear,subset=(mtcars$am==0),xlab="gear",ylab="mpg",main="automatic")
boxplot(mtcars$mpg~mtcars$gear,subset=(mtcars$am==1),xlab="gear",ylab="mpg",main="manual")
boxplot(mtcars$mpg~factor(mtcars$gear)*factor(mtcars$am),xlab="gear * transmission",ylab="mpg",main="boxplot of mpg by gear * transmission")
interaction.plot(mtcars$gear,mtcars$am,mtcars$mpg,type="b",col=c(1:3),leg.bty ="o",leg.bg = "beige",lwd=2,pch=c(18,24,22),xlab="number of gears",ylab="mean miles per gallon",main="interacion plot")
mpg_anova2=aov(mtcars$mpg~factor(mtcars$gear)*factor(mtcars$am))
summary(mpg_anova2)
TukeyHSD(mpg_anova2)
par(mfrow=c(1,2))
plot(TukeyHSD(mpg_anova2)
)
library(car)                                                                        
data("Quartet")
str(Quartet)
plot(Quartet$x,Quartet$y1)
lmfit=lm(y1~x,data=Quartet)
abline(lmfit,col="red")
lmfit
summary(lmfit)
coefficients(lmfit)
confint(lmfit,level = 0.95)
fitted(lmfit)
residuals(lmfit)
anova(lmfit)
vcov(lmfit)
influence(lmfit)
newdata=data.frame(x=c(3,6,12))
newdata
predict(lmfit,newdata,interval = "confidence",level = 0.95)
predict(lmfit,newdata,interval = "predict")
par(mfrow=c(2,2))
plot(lmfit)
plot(cooks.distance(lmfit))
par(mfrow=c(1,1))
plot(Quartet$x,Quartet$y2)
plot(Quartet$x,Quartet$y2)
lmfit=lm(Quartet$y2~I(Quartet$x)+I(Quartet$x^2))
abline(lmfit,col="red")
lines(sort(Quartet$x),lmfit$fit[order(Quartet$x)],col="red")

plot(Quartet$x,Quartet$y3)
library(MASS)
lmfit1=rlm(Quartet$y3~Quartet$x)
abline(lmfit1,col="red")
plot(Quartet$x,Quartet$y2)
lines(sort(Quartet$x),lmfit$fit[order(Quartet$x)],col="red")
lmfit$fit
lmfit$fit[order(Quartet$x)]
sort(Quartet$x)
order(Quartet$x)
Quartet$x
str(SLID)
par(mfrow=c(2,2))
plot(SLID$wages~SLID$language)
plot(SLID$wages~SLID$age)
plot(SLID$wages~SLID$education)
plot(SLID$wages~SLID$sex)
lmfit=lm(wages~.,data = SLID)
summary(lmfit)
lmfit=lm(wages~age+sex+education,data = SLID)
summary(lmfit)
par(mfrow=c(2,2))
plot(lmfit)
lmfit=lm(log(wages)~age+sex+education,data=SLID)
plot(lmfit)
vif(lmfit)
sqrt(vif(lmfit))>2
library(lmtest)
bptest(lmfit)
?bptest
library(mgcv)
library(MASS)
attach(Boston)
str(Boston)
fit=gam(dis~s(nox))
summary(fit)
plot(nox,dis)
x=seq(0,1,length=500)
y=predict(fit,data.frame(nox=x))
lines(x,y,col="red",lwd=2)
plot(fit)
head(seq)
head(x)
plot(x,y)
gam.check(fit)
library(C50)
install.packages("C50")

data(churn)
str(churnTrain)
churnTrain=churnTrain[,!names(churnTrain)%in%c("state","area_code","account_length")]
set.seed(2)
sample(10)
ind=sample(2,nrow(churnTrain),replace = T,prob=c(0.7,0.3))
sample(6,6,replace = T)
sample(6,6,replace = F)
trainset=churnTrain[ind==1,]
testset=churnTrain[ind==2,]
dim(trainset)
dim(testset)
head(names(trainset))
head                                                                (trainset)

library(rpart)
churn.rp=rpart(churn~.,data=trainset)
churn.rp
printcp(churn.rp)
plotcp(churn.rp)
summary(churn.rp)
plot(churn.rp,margin = 0.1)
text(churn.rp,all=T,use.n = T)
plot(churn.rp,uniform = T,branch = 0.6,margin = 0.1)
text(churn.rp,all=T,use.n = T)
predictions=predict(churn.rp,testset,type="class")
table(testset$churn,predictions)
library(caret)
install.packages("caret")
confusionMatrix(table(predictions,testset$churn))
?confusionMatrix
table(predictions)
table(testset$churn)
confusionMatrix(table(testset$churn,predictions))
table(testset$churn,predictions)
min(churn.rp$cptable[,"xerror"])
which.min(churn.rp$cptable[,"xerror"])
churn.cp=churn.rp$cptable[7,"CP"]
churn.rp
?prune
prune.tree=prune(churn.rp,cp=churn.cp)
plot(prune.tree,margin = 0.1)
text(prune.tree,all=T,use.n = T)
predictions=predict(prune.tree,testset,type="class")
table(testset$churn,predictions)
confusionMatrix(table(predictions,testset$churn))
library(class)
?knn
fit=glm(churn~.,data=trainset,family = binomial)
summary(fit)
fit=glm(churn~international_plan+voice_mail_plan+total_intl_calls+number_customer_service_calls,data=trainset,family = binomial)
summary(fit)
pred=predict(fit,testset,type="response")
class=pred>.5
summary(class)
tb=table(testset$churn,class)
tb
churn.mod=ifelse(testset$churn=="yes",1,0)
head(churn.mod)
head(testset$churn,n=100)
head(pred,n=100)
pred_class=churn.mod
length(pred_class)
head(churn.mod,n=100)
pred_class[pred<.5]=1-pred_class[pred<=.5]


library(e1071)
model=svm(churn~.,data=trainset,kernel="radial",cost=1,gamma=1/ncol(trainset))
summary(model)
iris.subset=subset(iris,select=c("Sepal.Length","Sepal.Width","Species"),Species%in%c("setosa","virginica"))
?subset
head(iris.subset)
plot(x=iris.subset$Sepal.Length,y=iris.subset$Sepal.Width,col=iris.subset$Species,pch=18)
svm.model=svm(Species~.,data=iris.subset,kernel="linear",cost=1,scale = F)
points(iris.subset[svm.model$index,c(1,2)],col="red",cex=2)
w=t(svm.model$coefs)%*%svm.model$SV
w
b=-svm.model$rho
abline(a=-b/w[1,2],b=-w[1,1]/w[1,2],col="red",lty=5)
svm.model=svm(Species~.,data=iris.subset,type="C-classification",kernel="linear",cost=10000,scale=F)
data(iris)
model.iris=svm(Species~.,iris)
plot(model.iris,iris,Petal.Width~Petal.Length,slice=list(Sepal.Width=3,Sepal.Length=4))
plot(model,trainset,total_day_minutes~total_intl_charge)
svm.pred=predict(model,testset[,!names(testset)%in%c("churn")])
head(testset)
names(testset)
rownames(testset)
table(duplicated(rownames(testset)))
sample(2,10,replace = T)
?sample
svm.table=table(svm.pred,testset$churn)
svm.table
classAgreement(svm.table)
library(caret)
confusionMatrix(svm.table)


library(car)
data("Quartet")
model.regression=svm(Quartet$y1~Quartet$x,type="eps-regression")
predict.y=predict(model.regression,Quartet$x)
predict.y
plot(Quartet$x,Quartet$y1,pch=19)
points(Quartet$x,predict.y,pch=15,col="red")
tuned=tune.svm(churn~.,data=trainset,gamma = 10^(-6:-1),cost=10^(1:2))
summary(tuned)
model.tuned=svm(churn~.,data=trainset,gamma=tuned$best.parameters$gamma,cost=tuned$best.parameters$cost)
summary(model.tuned)
svm.tuned.pred=predict(model.tuned,testset[,!names(testset)%in%c("churn")])
svm.tuned.table=table(svm.tuned.pred,testset$churn)
svm.tuned.table
classAgreement(svm.tuned.table)
confusionMatrix(svm.tuned.table)
data(iris)
ind=sample(2,nrow(iris),replace = T,prob=c(0.7,0.3))
trainset=iris[ind==1,]
testset=iris[ind==2,]
install.packages("neuralnet")
library(neuralnet)
head(trainset)
head(testset)
trainset$setosa=trainset$Species=="setosa"
trainset$virginica=trainset$Species=="virginica"
trainset$versicolor=trainset$Species=="versicolor"
network=neuralnet(versicolor+virginica+setosa~Sepal.Length+Sepal.Width+Petal.Length+Petal.Width,trainset,hidden=3)
network$result.matrix
head(network$generalized.weights[[1]])
plot(network)
install.packages("nnet")
par(mfrow=c(2,2))
gwplot(network,selected.covariate = "Sepal.Length")


net.predict=compute(network,testset[-5])$net.result
tail(testset)
?head
library(nnet)
data(iris)
net.prediction=c("versicolor","virginica","setosa")[apply(net.predict,1,which.max)]
predict.table=table(testset$Species,net.prediction)
predict.table
classAgreement(predict.table)
classAgreement(predict.table)
?classAgreement
library(e1071)
library(caret)
confusionMatrix(predict.table)
predict.table2=table(net.prediction,testset$Species)
confusionMatrix(predict.table2)
iris.nn=nnet(Species~.,data=trainset,size=2,rang=0.1,decay=5e-4,maxit=200)
summary(iris.nn)
iris.predict=predict(iris.nn,testset,type="class")
library(FSelector)
library(randomForest)
install.packages("FSelector")
weigts=random.forest.importance(churn~.,trainset,importance.type = 1)
library(C50)
data(churn)
set.seed(2)
ind=sample(2,nrow(churnTrain),replace = T,prob=c(0.7,0.3))
churnTrain=churnTrain[,!names(churnTrain)%in%c("state","area_code","account_length")]
trainset=churnTrain[ind==1,]
testset=churnTrain[ind==2,]
dim(trainset)
dim(testset)
weigts
subset=cutoff.k(weigts,5)
f=as.simple.formula(subset,"Class")
f
data("swiss")
head(Swiss)
swiss
swiss=swiss[,-1]
swiss.pca=prcomp(swiss,center=T,scale=T)
swiss.pca
summary(swiss.pca)
predict(swiss.pca,newdata=head(swiss,1))
screeplot(swiss.pca,type="lines")
plot(swiss.pca$x[,1],swiss.pca$x[,2],xlim=c(-4,4))
biplot(swiss.pca)

library(C50)
library(caret)
library(e1071)
data(churn)
str(churnTrain)
churnTrain=churnTrain[,!names(churnTrain)%in%c("state","area_code","account_length")]
set.seed(2)
ind=sample(2,nrow(churnTrain),replace = T,prob=c(0.7,0.3))
trainset=churnTrain[ind==1,]
testset=churnTrain[ind==2,]
dim(trainset)
dim(testset)
tuned=tune.svm(churn~.,data=trainset,gamma=10^-2,cost=10^2,tunecontrol=tune.control(cross = 10))
summary(tuned)
tuned$performances
svmfit=tuned$best.model
svmfit
table(trainset[,c("churn")],predict(svmfit))
control=trainControl(method = "repeatedcv",number=10,repeats = 3)
control
model=train(churn~.,data=trainset,method="rpart",preProcess="scale",trControl=control)
model
summary(model)
importance=varImp(model,scale=F)
importance
plot(importance)
library(rpart)
model.rp=rpart(churn~.,data=trainset)
model.rp$variable.importance
install.packages("xgboost")
library(devtools)
install_github("xgboost")
??devtools
library(rminer)
install.packages("xgboost")
install.packages("xgboost", repos="http://dmlc.ml/drat/", type = "source")

new_train=trainset[,!names(churnTrain)%in%c("churn","international_plan","voice_mail_plan")]
cor_mat=cor(new_train)
cor_mat
highlycorrelated=findCorrelation(cor_mat,cutoff=0.75)
highlycorrelated
names(new_train)[highlycorrelated]
library(e1071)
tuned=tune.svm(churn~.,data=trainset,gamma = 10^-2,cost=10^2,tunecontrol=tune.control(cross = 10))
summary(tuned)
tuned$performances
svmfit=tuned$best.model
summary(svmfit)
table(trainset[,c("churn")],predict(svmfit))
library(caret)
control=trainControl(method="repeatedcv",number=10,repeats=3)
control
model=train(churn~.,data=trainset,method="rpart",preProcess="scale",trControl=control)
warnings()
model

library(BSgenome)
available.genomes()
library(BiocInstaller)
biocLite("TxDb.Rnorvegicus.UCSC.rn6.refGene")
biocUpdatePackages("BiocInstaller")
library(TxDb.Rnorvegicus.UCSC.rn5.refGene)
source("https://bioconductor.org/biocLite.R")
biocValid()
source("https://bioconductor.org/biocLite.R")
biocLite("TxDb.Rnorvegicus.UCSC.rn6.refGene")
n
remove.packages("BiocInstaller")
source("https://bioconductor.org/biocLite.R")
source("https://bioconductor.org/biocLite.R")
library(BiocInstaller)
library(BiocInstaller)
biocLite("S4Vectors")
n
library(S4Vectors)
pwd()
wd()
getwd()
findOverlaps()
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(TxDb.Rnorvegicus.UCSC.rn5.refGene)
?import
gr=import("M2_rn5.bed",format="BED")
head(gr)
gr
txdb=TxDb.Rnorvegicus.UCSC.rn5.refGene
txdb
??promoters()
?promoters
pro20k=promoters(txdb,upstream = 10000,downstream = 10000)
pro2k=promoters(txdb,upstream = 2000,downstream = 2000)
help(promoters())
args(promoters)
pro2k=trim(pro2k)
pro20k
ov=findOverlaps(pro2k,gr)
ov
ovranges=pro2k[queryHits(ov)]
head(ovranges$tx_name)
name=select(org.Rn.eg.db,keys=ovranges$tx_name,columns=c("SYMBOL","ACCNUM"),keytype="ACCNUM")
head(name,n=100)
library(org.Rn.eg.db)
keytypes(org.Rn.eg.db)
dim(name)
name_unique=unique(name$SYMBOL)
length(name_unique)
name_u=name[!duplicated(name$SYMBOL),]
dim(name_u)
write.csv(name_u,"m2chiprn5_2kp0001.csv")
?write.csv
library(org.Hs.eg.db)
browseVignettes("org.Hs.eg.db")
x=org.Hs.egENSEMBLTRANS
head(x)
length(x)
mapped_genes=mappedkeys(x)
head(mapped_genes)
xx=as.list(x[mapped_genes])
length(xx)
x=org.Hs.egENSEMBL
mapped_genes=mappedkeys(x)
xx=as.list(x[mapped_genes])
length(x)
length(xx)
table(unlist(lapply(xx,length)))
head(xx)
class(xx)
head(unlist(lapply(xx,length)))
head(xx,n=100)
table(xx)
x=org.Hs.egCHR
mapped_genes=mappedkeys(x)
xx=as.list(x[mapped_genes])
length(x)
length(xx)
table(unlist(xx))
barplot(table(unlist(xx)))
library(GenomicRanges)
z=GRanges("chr1",IRanges(1000001,1001000),strand="+")
z
start(z)
end(z)
width(z)
strand(z)
mcols(z)
range(z)
ranges(z)
seqnames(z)
seqlevels(z)
seqlengths(z)
promoters
library(BSgenome.Hsapiens.UCSC.hg19)
dnastringset=getSeq(Hsapiens,z)
dnastringset
library(Biostrings)
substr(dnastringset,1,10)
subseq(dnastringset,1,10)
library(GenomicFeatures)
Views(dnastringset,1,10)
Views
library(IRanges)
complement(dnastringset)
reverseComplement(dnastringset)
vmatchPattern("GGGT",dnastringset)
class(dnastringset)
gr=IRanges(1,10)
Views(dnastringset,gr)
vmatchPattern
?vmatchPattern
letterFrequency(dnastringset,"CG")
alphabetFrequency(dnastringset,as.prob = T)
library(rtracklayer)
gr=import("fo.bed",format="BED")
gr
library(TxDb.Hsapiens.UCSC.hg18.knownGene)
txdb=TxDb.Hsapiens.UCSC.hg18.knownGene
pro=promoters(txdb,upstream = 2000,downstream = 2000)
pro=trim(pro)
pro
ov=findOverlaps(pro,gr)
ov
ovrange=pro[queryHits(ov)]
ovrange
library(org.Hs.eg.db)
name=select(org.Rn.eg.db,keys=ovranges$tx_id,columns=c("SYMBOL","ENTREZID"),keytype="ENTREZID")
keytypes(org.Hs.eg.db)
select(org.Hs.eg.db,keys="uc001ypz.1",columns=c("GENENAME","UCSCKG"),keytype="UCSCKG")
as.character(ovrange$tx_id)
ovrange$tx_id
columns(org.Hs.eg.db)
library(org.Rn.eg.db)
columns(org.Rn.eg.db)
library(TxDb.Rnorvegicus.UCSC.rn5.refGene)
txsbrat=TxDb.Rnorvegicus.UCSC.rn5.refGene
promoters(txsbrat)
gr
seqlevelsStyle(gr)="UCSC"
gr
pro
seqnames()
\
ovrange
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
promoters(TxDb.Hsapiens.UCSC.hg38.knownGene)

aa=available.packages(contriburl = "https://cran.rstudio.com/bin/windows/contrib/3.2/")
head(aa)
library(BiocInstaller)
biocLite("cummeRbund")
library(cummeRbund)
cuff_data=readCufflinks("diff_out")
csDensity(genes(cuff_data))
csScatter(genes(cuff_data),'C1','C2')
csVolcano(genes(cuff_data), 'C1', 'C2')
mygene=getGene(cuff_data,'regucalcin')
mygene
expressionBarplot(isoforms(mygene))
cuff_data
genes_diff=diffData(genes(cuff_data))
head(genes_diff)
sig_genes_diff=subset(genes_diff,significant=="yes")

head(sig_genes_diff)
dim(sig_genes_diff)
nrow(sig_genes_diff)
head(sig_genes_diff)
isoform_diff_data=diffData(isoforms(cuff_data))

head(isoform_diff_data)
sig_isoform_diff_data=subset(isoform_diff_data,significant=="yes")
nrow(sig_isoform_diff_data)
geneID=sig_genes_diff$gene_id
tss_diff_data=diffData(TSS(cuff_data),'C1','C2')
nrow(tss_diff_data)
sig_tss_diff_data=subset(tss_diff_data,significant=="yes")
nrow(sig_tss_diff_data)
splicing_diff_data=distValues(splicing(cuff_data))
sig_splicing_diff_data=subset(splicing_diff_data,significant=="yes")
nrow(sig_splicing_diff_data)
head(sig_isoform_diff_data)
head(sig_splicing_diff_data)
runInfo(cuff_data)
replicates(cuff_data)
aa=(annotation(genes(cuff_data)))
subset(aa,gene_short_name=="regucalcin")
nrow(aa$gene_short_name)
class(aa$gene_short_name)
?grep
head(repFpkmMatrix(genes(cuff_data)))
MAplot(genes(cuff_data),"C1","C2")
?read.table
bb=read.table("loctogenename.txt",sep=";")

nrow(aa)
rownames(aa)=aa$gene_id
subset(aa,aa$gene_id==c("XLOC_000001","XLOC_000002","XLOC_000006"))
aa[aa$gene_id==c("XLOC_000005","XLOC_000002","XLOC_000006"),]
head(aa)
head(aa$gene_id)
cuff_name
head(aa)
head(cuff_name)
length(geneID)
class(geneID)
head(geneID)
head(cuff_name)
head(aa)
colnames(aa)=c("gene_id","short_name")
head(aa)
aa[2,2]
cc=c("XLOC_000001","XLOC_000002","XLOC_000003")
subset(aa,aa$gene_id==c("XLOC_000001","XLOC_000002","XLOC_000004"))
head(aa)
aa[which(aa$gene_id==c("XLOC_000001","XLOC_000002","XLOC_000006","XLOC_000007")),]
getwd()
setwd("D:/Rdocuments/htcount")

aa["XLOC_000001",]



ls()
list.files()
list.dirs()
library("DESeq2")
??deseq2
biocLite("pasilla")
library(pasilla)
directory <- system.file("extdata", package="pasilla", mustWork=TRUE)
n
sampleFiles <- grep("treated",list.files(directory),value=TRUE)
directory
sampleFiles
getwd()
samplef=list.files()
samplef
sampleCondition <- sub("(.*treated).*","\\1",sampleFiles)
sampleCondition
samolecondition=c(rep("C1",3),rep("C2",3))
samolecondition
shortname=substr(samplef,1,5)
shortname
sampletable=data.frame(samplename=shortname,filename=samplef,conditon=samolecondition)
sampletable
ddshtseq=DESeqDataSetFromHTSeqCount(sampleTable = sampletable,design = ~conditon)
ddshtseq
cds=estimateSizeFactors(ddshtseq)
cds
sizeFactors(cds)
cds=estimateDispersions(cds)
vsd=varianceStabilizingTransformation(cdsB)
p=plotPCA(vsd,intgroup="conditon")
p
res=nbinomWaldTest(cds)
plotMA(res)
head(res)

ddshtseq=ddshtseq[rowSums(counts(ddshtseq))>1,]
ddshtseq$conditon
dds=DESeq(ddshtseq)
res=results(dds)
head(res)
res_sig=subset(res,pvalue<0.05)
head(res_sig)
nrow(res_sig)
summary(res)






head(rownames(res_sig))
head(sig_genes_diff)
?merge
head(merge(sig_genes_diff,aa,by="gene_id"),n=50)
head(aa)
head(sig_genes_diff)
fbgtogenename=read.table("fbgtogenename.txt",sep=";")
head(fbgtogenename)
colnames(fbgtogenename)=c("shortname","fbg")
res1=cbind(res_sig,rownames(res_sig))
mcols(res_sig)=rownames(res_sig)
head(res_sig)
?cbind
class(res_sig)
res3=as.data.frame(res_sig)
head(res2)
res3=cbind(rownames(res3),res3)
head(res3)
colnames(res3)[1]=c("fbg")
nrow(res)
ncol(res2)
head(res2)
head(fbgtogenename)
head(merge(res3,fbn,by="fbg"),n=100)
fbntogenename=read.table("fbntogenename.txt",sep=";")
head(fbntogenename)
colnames(fbntogenename)=c("fbg","shortname")
head(fbntogenename)
class(fbntogenename)
res2=data.frame(res2)
head(res2)
res2[which(res2$fbg==c("FBgn0031213")),]
head(res3)
head(fbntogenename)
rownames(res3)=NULL
head(res3)
fbntogenename[which(fbntogenename$fbg==c(" FBgn0031213")),]
class(fbntogenename)
fbn=trim(fbntogenename)
install.packages("raster")
library(raster)
head(fbn)
ee=merge(res3,fbn,by="fbg")$shortname
ff=merge(sig_genes_diff,aa,by="gene_id")$short_name
ee=trim(ee)

class(ee)
class(ff)
ee=as.vector(ee)
ff=trim(ff)
intersect(ee,ff)
head(ee)


shortname=c("c77","c78","t79","t80")
condition=c(rep("control",2),rep("treated",2))
condition
setwd("SRRcount")
list.files()
samplename=list.files()
samplename
library(DESeq2)
sampletable=data.frame(samplename=shortname,filename=samplename,condition=condition)
sampletable
dds=DESeqDataSetFromHTSeqCount(sampleTable = sampletable,design = ~condition)
dds=dds[rowSums(counts(dds))>1,]
result=DESeq(dds)
results=results(result)
head(results[results$pvalue<0.01,])

setwd("D://Rdocuments")
list.files(".//SRRcount")

browseVignettes("DESeq2")
library(BiocInstaller)
biocLite("edgeR")
library(edgeR)
x=read.table("tableCounts.txt")
dim(x)
rownames(x)=x[,1]
head(x)
xx=x[,c(2,3,4,5)]
head(xx)
head(x)
group=factor(c(1,1,2,2))
y=DGEList(counts = xx,group = group)
y=calcNormFactors(y)
class(y)
y
design=model.matrix(~group)
design
y=estimateDisp(y,design)
y
fit=glmFit(y,design)
lrt=glmLRT(fit,coef=2)
topTags(lrt)
lrt
yy=DGEList(counts=xx,group=group)
fityy=glmQLFit(yy,design)
yy=estimateDisp(yy,design)
qlf=glmQLFTest(fityy,coef=2)
topTags(qlf)
topTags(lrt)
qlft=qlf$table
sum(qlft$PValue<0.05)
lrtt=lrt$table
sum(lrtt$PValue<0.05)




library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)
pheno_data=read.csv("geuvadis_phenodata.csv")
pheno_data
bg_chrX=ballgown(dataDir = "ballgown",samplePattern = "ERR",pData = pheno_data)
bg_chrX
bg_chrX_filt=subset(bg_chrX,"rowVars(texpr(bg_chrX))>1",genomesubset=TRUE)
head(texpr(bg_chrX))
results_transcripts=stattest(bg_chrX_filt,feature = "transcript",covariate = "sex",adjustvars = c("population"),getFC = T,meas="FPKM")
head(results_transcripts)
results_genes=stattest(bg_chrX_filt,feature = "gene",covariate = "sex",adjustvars = c("population"),getFC = TRUE,meas="FPKM")
results_genes=data.frame(geneNames=ballgown::geneNames(bg_chrX_filt),geneIDs=ballgown::geneIDs(bg_chrX_filt),results_genes)
head(results_transcripts,n=100)
head(results_genes)
results_transcripts=arrange(results_transcripts,pval)
results_genes=arrange(results_genes,pval)
head(results_transcripts,n=20)
all(pheno_data$ids==list.files("ballgown"))
head(results_genes)
browseVignettes("ballgown")

library(UsingR)
x=father.son$fheight
sample(x,10)
samallest=floor(min(x))
largetest=ceiling(max(x))
values=seq(samallest,largetest,len=300)
heitghtecdf=ecdf(x)
heitghtecdf
head(x)
length(x)
head(values)
plot(values,heitghtecdf(values),type="l",xlab="height",ylab="P<a")
heitghtecdf(60)
?ecdf

library(biomaRt)
m=useMart("ensembl",dataset = "dmelanogaster_gene_ensembl")
m
lf=listFilters(m)
lf
head(lf)
lf[grep("flybase",lf$description,ignore.case = T),]
map=getBM(mart=m,attributes=c("ensembl_gene_id","flybasename_gene"),filters="flybasename_gene",values="lgs")
att=(listAttributes(m))
att[grep("g",att$name,ignore.case = T),]

?norm
?qnorm
qnorm(0.369)
pnorm(2)
qnorm(2)
qnorm(0.631)
library(GSE5859Subset)
data(GSE5859Subset)
dim(geneExpression)
dim(sampleInfo)
head(sampleInfo)
table(sampleInfo$group)
head(geneExpression)
match(sampleInfo$filename,colnames(geneExpression))
dim(geneAnnotation)
head(geneAnnotation)
head(rownames(geneExpression))
library(parathyroidSE)
data("parathyroidGenesSE")
se=parathyroidGenesSE
dim(se)
head(Se)
head(se)
se
head(assay(se))
x=assay(se)[,23]
y=assay(se)[,24]
head(y)
ind=which(x>0&y>0)
head(ind)
??splot
library(rafalib)
?splot
splot(log(x)+log(y),log(x)/log(y),subset=ind)
library(matrixStats)
vars=rowVars(assay(se)[,c(2,8,16,21)])
means=rowMeans(assay(se)[,c(2,8,16,21)])
splot(means,vars,log="xy",subset=which(means>0&vars>0))
abline(0,1,col=2,lwd=2)
rbinom(1500,1,0.05)
?rbinom
?rnorm
mean(rnorm(10000,mean=0,sd=1))
library(tissuesGeneExpression)
data(tissuesGeneExpression)
d=dist(t(e))
head(d)
library(rafalib)
mypar()
hc=hclust(d)
hc
myplclust(hc,labels = tissue,lab.col=as.fumeric(tissue),cex=0.5)
tissue
data(package="tissuesGeneExpression")
set.seed(1)
km=kmeans(t(e[1:2,]),centers = 7)
names(km)
km
mypar(1,2)
plot(e[1,],e[2,],col=as.fumeric(tissue),pch=16)
plot(e[1,],e[2,],col=km$cluster,pch=16)
head(tissue)
table(tissue)
?cmdscale
km=kmeans(t(e),centers = 7)
mds=cmdscale(d)
mypar(1,2)
plot(mds[,1],mds[,2])
plot(mds[,1],mds[,2],col=km$cluster,pch=16)
table(km$cluster)
km$cluster
idx=1:20
idx
library(gplots)
heatmap.2(e[idx,],labCol = tissue,trace="none")

shortname=c("59","60","61","62")
condition=c("control","treated","control","treated")
setwd("count180205")
list.files()
samplename=list.files()
samplename
sampletable=data.frame(samplename=shortname,filename=samplename,condition=condition)
sampletable
library(DESeq2)
dds=DESeqDataSetFromHTSeqCount(sampleTable = sampletable,design = ~condition)
dds=dds[rowSums(counts(dds))>1,]
ddsresult=DESeq(dds)
result=results(ddsresult)
result_pvalue=result[result$pvalue<0.05,]
nrow(result_pvalue)
names=rownames(result_pvalue)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
library(org.Mm.eg.db)
keytypes(org.Mm.eg.db)
symbol=select(org.Mm.eg.db,keys=names,keytype = "ACCNUM",columns = "SYMBOL")
tail(symbol,n=200)
library(BiocInstaller)
biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore"))
install.packages("WGCNA")        
library("WGCNA")
?sapply
?stripchart
2**(-0.3)
2**(-1.4)
