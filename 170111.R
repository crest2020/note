library(GSE5859Subset)
data(GSE5859Subset)
month=format()
data(package="GSE5859Subset")
head(sampleInfo)
head(geneAnnotation)
head(geneExpression)
dim(geneExpression)
month=format(sampleInfo$date,"%m")
head(month)
head(sampleInfo)
table(sampleInfo$group,month)
library(qvalue)
res=rowttests(geneExpression,as.factor(sampleInfo$group))
??rowttests
library(genefilter)
?rowttests
mypar(1,2)
library(rafalib)
mypar(1,2)
hist(res$p.value[which(!geneAnnotation$CHR%in%c("chrX","chrY"))],main="",ylim=c(0,1300))
head(res$p.value)
plot(res$dm,-log10(res$p.value))
head(res$statistic)
head(res$dm)
chr=geneAnnotation$CHR
?rowttests
browseVignettes("genefilter")
points(res$dm[which(chr=="chrX")],-log(res$p.value[which(chr=="chrX")]),col=1,pch=16)
points(res$dm[which(chr=="chrY")],-log(res$p.value[which(chr=="chrY")]),col=2,pch=16,xlab="Effect size",ylab="-log10(p-value)")
legend("bottomright",c("chrX","chrY"),col=1:2,pch=16)
qvals=qvalue(res$p.value)$qvalue
index=which(qvals<0.1)
abline(h=-log10(max(res$p.value[index])))
cat("Total genes with q-value < 0.1: ",length(index),"\n",
    "Number of selected genes on chrY: ",sum(chr[index]=="chrY",na.rm=T),"\n",
    "Number of selected genes on chrX: ",sum(chr[index]=="chrX",na.rm=T),sep=" ")
library(psych)
install.packages("psych")
fa.parallel(USJudgeRatings[,-1],fa="pc",n.iter=100,show.legend=FALSE,main="Scree plot with parallel analysis")
hist(res[,1])
fa.parallel(Harman23.cor$cov,n.obs=302,fa="both",n.iter = 100,show.legend = T,main="Scree plot with parallel analysis")
?fa.parallel
pa=principal(USJudgeRatings[,-1],nfactors = 1)
pa
pc=principal(Harman23.cor$cov,nfactors = 2,rotate="none")
pc
rc=principal(Harman23.cor$cov,nfactors = 2,rotate="varimax")
rc
pc=principal(USJudgeRatings[,-1],nfactors = 1,score=T)
pc
head(pc$scores)
cor(USJudgeRatings$CONT,pc$scores)
round(unclass(rc$weights),2)
options(digits = 2)
covariances=ability.cov$cov
correlation=cov2cor(covariances)
correlation
fa.parallel(correlation,n.obs = 12,fa="both",n.iter = 100,main="Scree plots with parallel analysis")
fa=fa(correlation,nfactors = 2,rotate="none",fm="pa")
fa
fa.varimax=fa(correlation,nfactors=2,rotate="varimax",fm="pa")
fa.varimax
fa.promax=fa(correlation,nfactors=2,rotate="promax",fm="pa")
fa.promax
factor.plot(fa.varimax,labels=rownames(fa.varimax$loadings))
fa.varimax$weights
fa.promax$weights
library(lattice)
histogram(~height|voice.part,data=singer,main="Distribution of Heights by Voice Pitch",xlab="Height(inches)")
attach(mtcars)
gear=factor(gear,level=c(3,4,5),labels = c("3 gear","4 gear","5 gear"))
cyl=factor(cyl,level=c(4,6,8),labels=c("4 cylinders","6 cylinders","8 cylinders"))
densityplot(~mpg,main="Density Plot",xlab = "Miles per Gallon")
densityplot(~mpg|cyl,main="Density plot by number of cylinders",xlab="Miles per Gallon")
bwplot(cyl~mpg|gear,main="Box Plots by Cylinders and gears",xlab="Miles per gallon")
xyplot(mpg~wt|cyl*gear,main="scatter plot by cylinders and gears",xlab="Car weight",ylab="Miles per gallon")
cloud(mpg~wt*qsec|cyl,main="3D scatter plot by cylinders")
dotplot(cyl~mpg|gear,main="dot plots by number of gears and cylinders",xlab="Miles per Gallon")
splom(mtcars[c(1,3,4,5,6)],main="scatter plot matrix for mtcars data")
displacement=equal.count(mtcars$disp,number=3,overlap=0)
xyplot(mpg~wt|displacement,data=mtcars,main="Miles per Gallon by Engine displacement",xlab="Weight",ylab = "Miles per Gallon",layout=c(3,1),aspect = 1.5)
mypanel=function(x,y){
  panel.xyplot(x,y,pch=19)
  panel.rug(x,y)
  panel.grid(h=-1,v=-1)
  panel.lmline(x,y,col="red",lwd=1,lty=2)
}
mypanel
xyplot(mpg~wt|displacement,data=mtcars,layout=c(3,1),aspect=1.5,main="Miles per Gallon",xlab="Weight",ylab="Miles per gallon",panel=mypanel)
x=rep(1:10,2)
x
y=c(2,4,6)
x%in%y
y%in%x
duplicated(x)
x[!duplicated(x)]
unique(x)
z=c(1,2,3,3,4,5)
unique(z)
z
duplicated(z)
z[!duplicated(z)]
duplicated(z)
index=z%in%z[duplicated(z)]
z[!index]
y=1:200
interval=cut(y,right = F,breaks = c(1,2,6,11,21,51,101,length(y)+1),
             labels = c("1","2-5","6-10","11-20","21-50","51-100",">=101"))
table(interval)
?cut
plot(interval,ylim=c(0,110),xlab="intervals",ylab="counts",col="green")
text(labels=as.character(table(interval)),x=seq(0.7,8,by=1.2),y=as.vector(table(interval))+2)
dataf=matrix(c(1:24),4,6)
dataf
?data.matrix
rownames(dataf)
row.names(dataf)
rownames(dataf)=1:4
rownames(dataf)
row.names(dataf)
a=c(1,2,3,4,6)
b=c(6,7,8,9,2)
match(a,b)
a%in%b
summary(dataf)
cor(dataf[,1:2])
cor(dataf[1:2,])
myDF=as.data.frame(matrix(rnorm(100000),10000,10))
myCol=c(1,1,1,2,2,2,3,3,4,4)
myDFmean=t(aggregate(t(myDF),by=list(myCol),FUN=mean,na.rm=T))
list(myCol)
head(myDFmean)
head(myDF)
myList=tapply(colnames(myDF),myCol,list)
myList
names(myList)=sapply(myList,paste,collapse="_")
names(myList)
myDFmean=sapply(myList, function(x) mean(as.data.frame(t(myDF[,x]))))
myList
myDFmean[1:4,]
x=paste(rep("A",times=12),1:12,sep="")
x
y=paste(rep("B",12),1:12,sep="")
append(x,y)
Z <- array(1:12, dim=c(12,8))
z
Z
month.name(1:12)
?gsub
x <- data.frame(matrix(rep(c("P","A","M"),20),10,5))
x
count("p",x)
sum(x=="P")
table(iris$Sepal.Length,exclude = NULL)[iris$Sepal.Length[1]]
length(iris$Sepal.Length)
iris$Sepal.Length
x=c(1,1,2,2,3,5)
x
table(x)[x]
x=c("a","a","b","e")
table(x)
table(x)[x]
my_counts <- table(iris$Sepal.Length, exclude=NULL)[iris$Sepal.Length]
my_counts
myvec <- c("a", "a", "b", "c", NA, NA)
table(factor(myvec, levels=c(unique(myvec), "z"), exclude=NULL))
labels=paste("sample",1:5,sep="")
labels
combn(labels,m=2,FUN=paste,collapse="-")
seq(along=labels)
?seq
y=as.data.frame(matrix(runif(30),ncol=3,dimnames = list(letters[1:10],LETTERS[1:3])))
y
plot(y[,1],y[,2])
plot(y[,1],y[,2],type="n",main="plot of lablels")
text(y[,1],y[,2],rownames(y))
plot(y[,1],y[,2],pch=20,col="red",main="plot of symbols and labels")
text(y[,1]+0.03,y[,2],rownames(y))
grid(5,5,lwd=2)
op=par(mar=c(2,2,2,2),bg="lightblue")plot(y[,1],y[,2],type='p',col="red",cex.lab=1.2,cex.axis=1.2,cex.main=1.2,cex.sub=1,lwd=4,pch=20,xlab = "x label",ylab="y label",main="my main",sub="mu sub")
par(op)
plot(y[,1],y[,2])
myline=lm(y[,2]~y[,1],data=y[,1:2])
abline(myline,lwd=2)
summary(myline)
text(y[,1],y[,2],expression(sum(frac(1,srqt(x^2*pi)))),cex=1.3)
plot(y)
pairs(y)
library(affy)
library(BiocInstaller)
biocLite("affy")
myseq <- c("ATGCAGACATAGTG", "ATGAACATAGATCC", "GTACAGATCAC") # Creates a sample sequence data set.
grep("ATG",myseq)
myseq[grep("ATG",myseq)]
pos1=regexpr("AT",myseq)
pos1
pos2=gregexpr("AT",myseq)
as.numeric(pos2[[1]])
gsub("^ATG","atg",myseq)
labels=paste("Sample",1:5,sep="")
labels
combn(labels,m=2,FUN=paste,collapse="-")
aggregate(iris[,1:4],by=list(iris$Species),FUN=mean,na.rm=T)
myMA=matrix(rnorm(100000),10000,10,dimnames=list(1:10000,paste("C",1:10,sep="")))
myMA
apply(myMA,1,mean)[1:4]
my_frame=data.frame(Month=month.name,N=12)
my_frame
month.name
month.abb
month.name
my_query=c("May","August")
my_frame[my_frame$Month%in%my_query,]
my_frame
my_frame$Month
frame1=iris[sample(1:length(iris[,1]),30),]
frame1
sample(1:10,5)
dim(iris)
dim(frame1)
my_result=merge(frame1,iris,by.x=0,by.y=0,all=T)
dim(my_result)
head(my_result)
head(frame1)
?merge
b=data.frame(c(1:10),c(11:20))
b
colnames(b)=c("A","B")
colnames(cc)=c("B","C")
cc=data.frame(c(11:20),c(21:30))
cc
merge(b,cc,by="B")
merge(b,cc)
dd=merge(b,cc,by.x=0,by.y = 0)
sort(merge(b,cc,by.x=0,by.y = 0))
colnames(dd)
colnames(dd)=c(1:5)
dd
order(dd$`1`)
?sort
dd[order(dd[,1],decreasing = F),]
library(help=lattice)
library(lattice)
xyplot(1:10~1:10)
xyplot(1:10~1:100|rep(LETTERS[1:5],each=2),as.table=T)
myplot=xyplot(Petal.Width~Sepal.Width|Species,data=iris,layout=c(3,1,1))
print(myplot)
xyplot(Petal.Width~Sepal.Width|Species,data=iris,layout=c(3,1,1))
show.settings()
default=trellis.par.get()
mytheme=default
names(mytheme)
mytheme["background"][[1]][[2]]="grey"
trellis.par.set(mytheme)
show.settings()
xyplot(1:10~1:10|rep(LETTERS[1:5],each=2),as.table=T,layout=c(1,5,1),col=c("red","blue"))
trellis.par.set(default)
library(ggplot2)
ggplot(iris,aes(Sepal.Length,Sepal.Width))+geom_point()
ggplot(iris,aes(Sepal.Length,Sepal.Width))+geom_point(aes(color=Species),size=4)
ggplot(iris,aes(Sepal.Length,Sepal.Width))+geom_point(aes(color=Species),size=4)+ylim(2,4)+xlim(4,8)+scale_color_manual(values = rainbow(10))
ggplot(iris,aes(Sepal.Length,Sepal.Width,label=1:150))+geom_text()
??opts
??ggplot2
?ggplot
ggplot(iris,aes(Sepal.Length,Sepal.Width,label=1:150))+geom_point()+geom_text(hjust=0.5,vjust=-0.5)
ggplot(iris,aes(Sepal.Length,Sepal.Width))+geom_point()+ggtitle(panel.background=theme_rect(fill = "white", colour = "black"))
ggplot(iris,aes(Sepal.Length,Sepal.Width))+geom_point()+stat_smooth(method = "lm",se=F)
ggplot(iris,aes(Sepal.Length,Sepal.Width))+geom_point()+coord_trans(x="log2",y="log2")
library(biomaRt)
m=useMart("ensembl",dataset = "dmelanogaster_gene_ensembl")
lf=listFilters(m)
head(lf)
lf[grep("name",lf$description,ignore.case = T),]
map=getBM(mart=m,attributes = c("ensembl_gene_id","flybasename_gene"),filters = "flybasename_gene",values = "lgs")
map
library(GenomicFeatures)
grl=exonsBy(TxDb.Dmelanogaster.UCSC.dm3.ensGene, by="gene")
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
grl
gene=grl[[map$ensembl_gene_id[1]]]
gene
rg=range(gene)
plot(c(start(rg),end(rg)),c(0,0),type="n",xlab=seqnames(gene)[1],ylab="")
arrows(start(gene),rep(0,length(gene)),end(gene),rep(0,length(gene)),lwd=3,length=.1)
arrows(start(gene),rep(0,length(gene)),end(gene),rep(0,length(gene)),lwd=3,length=0.1,code=ifelse(as.character(strand(gene)[1])=="+",2,1))
strand(gene)[1]
library(BiocInstaller)
library(Gviz)
gtrack=GenomeAxisTrack()
atrack=AnnotationTrack(gene,name="Gene Model")
plotTracks(list(gtrack,atrack))
xgr=as(xcov,"GRanges")
biocLite("pasillaBamSubset")
library(pasillaBamSubset)
fl1=untreated1_chr4()
fl2=untreated3_chr4()
fl1
library(Rsamtools)
x=readGAlignments(fl1)
library(GenomicAlignments)
xcov=coverage(x)
z=GRanges("chr4",IRanges(456500,466000))
z
xcov[z]
xcov$chr4[ranges(z)]
xnum=as.numeric(xcov$chr4[ranges(z)])
plot(xnum)
y=readGAlignmentPairs(fl2)
ycov=coverage(y)
ynum=as.numeric(ycov$chr4[ranges(z)])
plot(xnum,type="l",col="blue",lwd=2)
lines(ynum,col="red",lwd=2)
plot(xnum,type="l",col="blue",lwd=2,xlim=c(6200,6600))
lines(ynum,col="red",lwd=2)
xgr=as(xcov,"GRanges")
ygr=as(ycov,"GRanges")
dtrack1=DataTrack(xgr[xgr%over%z],name="sample1")
help(%over%)
?over
xgr%over%z
xgr[xgr%in%z]
z
dtrack2=DataTrack(ygr[ygr%over%z],name="sample2")
plotTracks(list(gtrack,atrack,dtrack1,dtrack2))
plotTracks(list(gtrack,atrack,dtrack1,dtrack2),type="polygon")
library(ggbio)
autoplot(gene)
autoplot(fl1,which=z)
autoplot(fl2,which=z)
library(ph525x)
dfHclust(mtcars)
system.time(lapply(1:8, function(x) Sys.sleep(1)))
library(parallel)
detectCores()
options(mc.cores=4)
system.time(mclapply(1:8,function(x)Sys.sleep(1)))
?options
library(BiocInstaller)
biocLite("RNAseqData.HNRNPC.bam.chr14")
library(RNAseqData.HNRNPC.bam.chr14)
file2=dir(system.file("extdata",package="RNAseqData.HNRNPC.bam.chr14"))
file1=file2[grep("bam$",file2)]
file1
fns=RNAseqData.HNRNPC.bam.chr14_BAMFILES
fns
system.file(package = "RNAseqData.HNRNPC.bam.chr14")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb=TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb,force=T)="chr14"
ebg=exonsBy(txdb,by="gene")
ebg
library(GenomicAlignments)
s1=system.time(i1<-summarizeOverlaps(ebg,fns[3]))
s1
BiocParallel::bpparam()
library(BiocParallel)
register(SnowParam(workers=3))
system.time(mclapply(1:8,function(x)Sys.sleep(1)))
library(GO.db)
GO.db$conn@dbname
dbListTables(GO.db$conn)
dbGetQuery(GO.db$conn,"select * from go_term limit 5")
sri=read.csv("SraRunInfo.csv",stringsAsFactors = F)
sri
head(sri)
keep=grep("CG8144|Untreated-",sri$LibraryName)
keep
sri=sri[keep,]
sri$LibraryName
sri$LibraryName=gsub("S2_DRSC_","",sri$LibraryName)
sri$LibraryName
samples=unique(sri[,c("LibraryName","LibraryLayout")])
samples
for(i in seq_len(nrow(samples))){
  rw=(sri$LibraryName==samples$LibraryName[i])
  if(samples$LibraryLayout[i]=="PAIRED"){
    samples$fastq1[i]=paste0(sri$Run[rw],"_1.fastq",collapse = ",")
    samples$fastq2[i]=paste0(sri$Run[rw],"_2.fastq",collapse = ",")
  }else{
    samples$fastq1[i]=paste0(sri$Run[rw],".fastq",collapse = ",")
    samples$fastq2[i]=""
  }
}
samples$fastq1
samples[,c("fastq1","LibraryName")]

sri$LibraryName
rw
sri$Run[rw]
sri$LibraryName==samples$LibraryName[1]
samples$LibraryName[1]
sri$Run
sri$Run[TRUE,FALSE,]
samples$condition="CTRL"
samples$condition[grep("RNAi",samples$LibraryName)]="KD"
samples$condition
samples$shortname=paste(substr(samples$condition,1,2),substr(samples$LibraryLayout,1,2),seq_len(nrow(samples)),sep=".")
samples$shortname
seq_len(nrow(samples))
substr(samples$condition,1,2)
substr(samples$LibraryLayout,1,2)
samples$LibraryLayout
?substr
samples
?paste0
nth <- paste0(1:12, c("st", "nd", "rd", rep("th", 9)))
nth
paste(1:12, c("st", "nd", "rd", rep("th", 9)))
library(BiocInstaller)
biocLite("h5vc")
library(h5vc)
library(rhdf5)
tallyfile=system.file("extdata","example.tally.hfs5",package="h5vcData")
h5ls(tallyfile)
sampleData=getSampleData(tallyfile,"/ExampleStudy/16")
sampleData
dim(sampleData)
H5close()
data=h5readBlock(filename = tallyfile,
                 group="/ExampleStudy/16",
                 names=c("Coverages","Counts"),
                 range=c(29000000,29001000))
str(data)
position=29979628
windows=30
samples=sampleData$Sample[sampleData$Patient=="Patient8"]
samples
data1=h5readBlock(filename=tallyfile,
                  group="/ExampleStudy/16",
                  names=c("Coverages","Counts","Deletions","Reference"),
                  range=c(position-windows,position+windows))

p=mismatchPlot(data=data1,sampledata = sampleData,samples=samples,windowsize = windows,position = position)
print(p)
benchOOM
library(ph525x)
benchOOM()
benchOOM
ggshot()
biocLite("geuvStore")
library(geuvStore)
m=makeGeuvStore()
class(m)
m
library(gQTLBase)
utl=system.time(l1<-storeApply(m,length))
utl
library(doParallel)
registerDoParallel(cores=2)
ut2=system.time(l2<-storeApply(m,length))
ut2
print(sum(unlist(l2)))
all.equal(unlist(l1),unlist(l2))
library(harbChIP)
library(BiocInstaller)
biocLite("harbChIP")
data(harbChIP)
harbChIP
abstract(harbChIP)
mind=which(sampleNames(harbChIP)=="MBP1")
head(mind)
qqnorm(exprs(harbChIP)[,mind],main="MBP1 binding")
qqline(0,1)
abline(0,1)
exp=exprs(harbChIP)
dim(exp)
head(exp,2)
topb=featureNames(harbChIP)[order(exprs(harbChIP)[,mind],decreasing = TRUE)[1:5]]
topb
library(org.Sc.sgd.db)
biocLite("org.Sc.sgd.db")
select(org.Sc.sgd.db,key=topb,keytype = "ORF",columns = "COMMON")
exp
fdata=fData(harbChIP)
head(fdata)
pdata=pData(harbChIP)
head(pdata)
length(pdata)
dim(pdata)
biocLite("yeastCC")
library(yeastCC)
data(package="yeastCC")
