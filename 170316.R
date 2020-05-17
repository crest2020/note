library(VIM)
install.packages("VIM")
library(mice)
install.packages("mice")
data(sleep,package ="VIM")
da=(sleep[!complete.cases(sleep),])
cotableunts("NA",da)
table(da[1:10,])
dim(da)
duplicated(as.list(da))

da=da[1:5,1:5]
da
class(as.matrix(da))
da1=as.matrix(da)
da2=as.vector(da1)
da2
table(da2,useNA="ifany")
table
?table
da2
sum(is.na(sleep$Dream))
mean(is.na(sleep$Dream))
is.na(sleep$Dream)
table(is.na(sleep$Dream))
library(mice)
md.pattern(sleep)
aggr(sleep,prop=F,numbers=T)
matrixplot(sleep)
marginplot(sleep[c("Gest","Dream")],pch=c(20),col=c("darkgray","red","blue"))
x=as.data.frame(abs(is.na(sleep)))
x
abs(TRUE)
head(sleep,n=5)
head(x,n=5)
y=x[which(sd(x)>0),]
sd(x)
y
newdata=sleep[complete.cases(sleep),]
cor(newdata)
fit=lm(Dream~Span+Gest,data=newdata)
summary(fit)
imp=mice(sleep,seed=1234)
fit1=with(imp,lm(Dream~Span+Gest))
pooled=pool(fit1)
summary(pooled)
imp
summary(fit1)
cor(sleep,use="pairwise.complete.obs")

library(coin)
states=as.data.frame(state.x77)
set.seed(1234)
spearman_test(Illiteracy~Murder,data=states,distribution=approximate(B=9999))
library(MASS)
wilcox_test(U1~factor(U2),data=UScrime,distribution="exact")
head(UScrime)
?wilcoxsign_test()
U2
head(factor(UScrime$U2))
library(vcd)
arthritis=transform(Arthritis,improved=as.factor(as.numeric(Improved)))
set.seed(1234)
chise_test=chisq_test(Treatment~Improved,data=Arthritis,distribution=approximate(B=9999))
chise_test
table(Arthritis$Treatment,Arthritis$Improved)
library(lmPerm)
set.seed(1234)
fit=lmp(weight~height,data=women,perm="prob")
summary(fit)
fit2=lmp(weight~height+I(height^2),data=women,perm = "prob")
summary(fit2)
fit3=lmp(Murder~Population+Illiteracy+Income+Frost,data=states,perm="Prob")
head(states)
summary(fit3)
library(multcomp)
set.seed(1234)
fit44=aovp(response~trt,data=cholesterol,perm = "prob")
summary(fit44)
fit5=aovp(weight~gesttime+dose,data=litter,perm = "Prob")
summary(fit5)
?aovp
fit6=aovp(len~supp*dose,data = ToothGrowth,perm="Prob")
summary(fit6)
library(devtools)
install_github("coloncancermeth/genomicsclass")
library(coloncancermeth)
library(SpikeInSubset)
data(rma95)
fac=factor(rep(1:2,each=3))
fac
library(genefilter)
rtt=rowttests(exprs(rma95),fac)

mask=with(rtt,abs(dm)<.2&p.value<0.01)
spike=rownames(rma95)%in%colnames(pData(rma95))
cols=ifelse(mask,"red",ifelse(spike,"dodgerblue","black"))
with(rtt,plot(-dm,-log10(p.value),cex=.8,pch=16,xlim=c(-1,1),ylim=c(0,5),xlab="difference in means",col=cols))
abline(h=2,v=c(-.2,.2),lty=2)
rtt$s=apply(exprs(rma95),1,function(row) sqrt(.5*(var(row[1:3])+var(row[4:6]))))
with(rtt,plot(s,-log(p.value),cex=.8,pch=16,log="x",xlab="estimate of standard deviation",col=cols))
library(limma)
fit=lmFit(rma95,design=model.matrix(~fac))
colnames(coef(fit))
fit=eBayes(fit)
fit
fit=topTable(fit,coef = 2)
dim(fit)
?topTable
library(BiocInstaller)
biocLite("package.skeleton")
library(OrganismDbi)
library(org.Sc.sgd.db)
biocLite("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")
gd=list(join1=c(GO.db="GOID",org.Sc.sgd.db="GO"),join2=c(org.Sc.sgd.db="ENTREZID",TxDb.Scerevisiae.UCSC.sacCer3.sgdGene="GENEID"))
gd
if(!file.exists("Sac.Cer3")){
  makeOrganismPackage(pkgname = "Sac.Cer3",
  graphData=gd,organism="Saccharomyces cerevisiae",version="1.0.0",maintainer = "xuhongzhan@pku.edu.cn",author="xuhongzhan",destDir = ".",license = "Artistic-2.0")}
install.packages("Sac.cer3", repos=NULL, type="source")
getwd()
if (!file.exists("Sac.cer3")) # don't do twice...
  makeOrganismPackage(pkgname="Sac.cer3",  # simplify typing!
                      graphData=gd, organism="Saccharomyces cerevisiae",
                      version="1.0.0", maintainer="Student <ph525x@harvardx.edu>",
                      author="Student <ph525x@harvardx.edu>",
                      destDir=".",
                      license="Artistic-2.0")


library(Sac.cer3)
?makeOrganismPackage
myseq <- c("ATGCAGACATAGTG", "ATGAACATAGATCC", "GTACAGATCAC")
myseq[grep("ATG",myseq)]
pos1=regexpr("ATG",myseq)
pos1
as.numeric(pos1)
attributes(pos1)$match.length
substring(myseq[1],c(1,3))
myseq
myseq[1]
?substring
substring(myseq, c(1,4,7), c(2,6,10))
myseq_comp <- chartr("ATGC", "TACG", myseq)
myseq_comp
myseq
reverse(myseq_comp)
library(Biostrings)
reverseComplement(myseq1)
?reverseComplement
myseq1=DNAString(myseq[1])
AAdf <- read.table(file="http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/AA.txt", header=T,"\t")
AAdf
aav=AAdf[,2]
names(aav)=AAdf[,1]
y <- c("ATGCATTGGACGTTAG")
y
y=gsub("(...)","\\1_",y)
y
?gsub
y=unlist(strsplit(y,"_"))
y
y=y[grep("^...$",y)]
y
aav
aav[y]
sapply(1:100,function(x) paste(sample(c("A","T","C","G"),20,replace = T),collapse = ""))
sapply(1:100, function(x) paste(sample(1:40,20,replace=T),collapse = " "))
myseq=sapply(1:12,function(x) paste(sample(c("A","T","C","G"),40,replace = T),collapse = ""))
myseq
writeLines(myseq)
writeLines(strtrim(myseq,10))
writeLines(strwrap(myseq,indent = 20))
myname=paste(">",month.name,sep="")
myname
writeLines(as.vector(t(cbind(myname,myseq))),"myseq.fasta")
?DNAString
library(rafalib)
library(tissuesGeneExpression)
data(tissuesGeneExpression)
colind=tissue%in%c("kidney","colon","liver")
mat=e[,colind]
head(mat)
group=factor(tissue[colind])
group
head(tissue)
head(colind)
dim(mat)
s=svd(mat-rowMeans(mat))
PC1=s$d[1]*s$v[,1]
PC2=s$d[2]*s$v[,2]
plot(PC1,PC2,pch=21,bg=as.numeric(group))
as.numeric(group)
legend("bottomright",levels(group),col=seq(along=levels(group)),pch=15,cex=1)
plot(s$d^2/sum(s$d^2))
d=dist(t(mat))
mds=cmdscale(d)
plot(mds[,1],mds[,2],bg=as.numeric(group),pch=21,xlab="first dimension",ylab="second dimension")
legend("bottomleft",levels(group),col=seq(along=levels(group)),pch=15,cex=1)
dir()
install.packages("xlsx")
library(xlsx)
library(rJava)
dir()
aa=read.xlsx("CufflinksFPKM.xlsx",sheetName = "WT-TG-FPKM",header = T)
?read.xlsx
?gsub
a=c("Chr1","chr2","chr3")
b=c("10.1-17246","2.3-57295857","5.7-48919","MT.1-481189","X.6-28318")
cc=gsub("(.*)[T]*\\..*","chr\\1",b)
print(a)
gsub("chrMT","chrM",cc)

?dnorm
?rbinom
dbinom(1,10,0.05)
rnorm(2,0,1)
pnorm(0,0,1)
dnorm(0,0,1)
i <- array(c(1:5,5:1),dim=c(3,2))
i
array(c(1:5,5:1))
array(1:8, dim=c(12,8))
my_frame <- data.frame(y1=rnorm(12), y2=rnorm(12), y3=rnorm(12), y4=rnorm(12))
my_frame
names(my_frame)
?names
summary(my_frame)
myDF <- as.data.frame(matrix(rnorm(200), 20, 10))
myDF
myList <- tapply(colnames(myDF), c(1,1,1,2,2,2,3,3,4,4), list)
myList
colnames(myDF)
 sapply(myList, paste, collapse="_")
 x <- 1:10; x <- x[1:12]; z <- data.frame(x,y=12:1)
 x
 y
 z
 is.na(x)
 is.na(z)
 z[is.na(z)]=0
 z
 as.numeric(1:10<=5)
 ?regexp
 month.name[grep("A",month.name)]
 gsub('(i.*a)','\\1_xxx',iris$Species,perl=T)
 table(iris$Species)
 ?gsub
 gsub("x","hz","xhaoihti4ha")
 x=rep(1:10,2)
 x
 sample(1:10,5,replace = T)
 sample(1:10,10,replace=T)
 unique(x)
 as.integer(runif(100,min=1,max=5))
 sample(1:5,100,replace = T)
 runif(10)
 runif(10)*5
 ?runif
 sprintf("ID:%09d",1:10)
 1:30
 letters
 LETTERS
 month.name
 month.abb
 seq(1,30,by=0.5)
 rep(LETTERS[1:8],times=5)
 rep(1,5)
 paste(LETTERS[1:8],1:12,sep=":")
 x=paste(rep("A",times=12),1:12,sep="")
 x
 y=paste(rep("B",times=12),1:12,sep="")
 y
 append(x,y)
 data.frame(x,y)
 x=1:100
 x
 x[2:23]
 x[-(2:23)]
 x[5]=99
 x
which(rep(letters,2)=="c")
match(c("c","g"),rep(letters,2)
)
x=rep(1:10,2)
u=c(2,4,6)
x
y
u
x[x%in%u]
y=1:200
interval=cut(y,right=F,breaks=c(1,2,6,11,21,51,101,length(y)+1),labels=c("1","2-5","6-10","11-20","21-50","51-100",">=101"))
interval
table(interval)
plot(interval,ylim=c(0,110),xlab="Intervals",ylab="Counts",col="green")
text(labels=as.character(table(interval)),x=seq(1,8,by=1),y=as.vector(table(interval))+2)
seq(0.7,8,by=1.2)
?text
?plot
x=matrix(1:30,3,10,byrow=T)
x
dim(x)
dim(x)=c(3,5,2)
x
dim(x)
dim(x)=c(5,3,2)
x
y=matrix(1:30,3,10)
y
x=array(1:25,dim=c(5,5))
x
y=c(x)
y
x[c(1:5),3]
mean(x[c(1:5),3])
i=array(c(1:5,5:1),dim=c(3,2))
i
i=array(c(1:5,5:1))
i
dim(i)
library(reshape2)
iris_mean=aggregate(iris[,1:4],by=list(Species=iris$Species),FUN=mean)
iris_mean
head(iris)
iris_mean
df_mean=melt(iris_mean,id.vars = c("Species"),variable.name = "Samples")
df_mean
dcast(df_mean,formula = Species~...)
my_list=list(name="Fred",wife="Mary",no.children=3,child.ages=c(4,7,9))
attributes(my_list)
names(my_list)
my_list[2]
my_list[[2]]
my_list$wife
my_list[[4]][2]
my_list$wife=1:12
my_list
my_list=c(my_list,list(my_title2=month.name[1:12]))
unlist(my_list)
data.frame(unlist(my_list))
x=1:10
x=x[1:12]
x
z=data.frame(x,y=12:1)
z
is.na(x)
x=letters[1:10]
x
print9x
print(x)
x[!is.na(x)]
apply(z,1,mean,na.rm=T)
z[is.na(z)]=0
z
unique(iris$Sepal.Length)
length(unique(iris$Sepal.Length))
myvec=c("a","a","b","c",NA,NA)
table(factor(myvec,levels=c(unique(myvec),"z"),exclude = NULL))
factor(c("good","bad"))
factor(c("good","bad"),levels=c("good","bad"))
labels=paste("Sample",1:5,sep="")
combn(labels,m=2,FUN=paste,collapse="-")
seq(along=labels)
lables
labels
month.name[1]
x=gsub("(a)","\\1_",month.name[1],perl=T)
x
strsplit(x,"_")
paste(rev(unlist(strsplit(x,NULL))),collapse = "")
?rev
cat(month.name,file="zzz.txt",sep="\n")
x=readLines("zzz.txt")
x
x=x[c(grep("^J",as.character(x),perl = T))]
x
c(grep("^J",as.character(x),perl = T))
strsplit(x,"u")
x=shell("dir",intern = T)
x
?shell.exec
?read.table
seq(1:10,along=1)
seq(1:10)
seq(1:10,along.with = 5)
?seq
seq(along.with=5)
seq(from=1,to=10,along.with = 3)
?seq_along
seq(0, 1, length=11)
myseq <-DNAStringSet(c("CCCATGCAGACATAGTG", "CCCATGAACATAGATCC", "CCCGTACAGATCACGTG"))
myseq
names(myseq)=letters[1:3]
myseq
trimLRPatterns(Lpattern = "CCC",Rpattern = "AGTG",subject = myseq,max.Lmismatch = 0.33,max.Rmismatch = 1)
?trimLRPatterns
trimLRPatterns(Lpattern = "CCC",Rpattern = "AGTG",subject = myseq,max.Lmismatch = c(0,0,0),max.Rmismatch = c(1,0,0),ranges = T)
library(Biostrings)
myseq1 <- readDNAStringSet("ftp://ftp.ncbi.nih.gov/genbank/genomes/Bacteria/Halobacterium_sp_uid217/AE004437.ffn", "fasta")
library(help=Biostrings)
library(Biostrings)
??biostrings
BString
readDNAStringSet()
GENETIC_CODE
IUPAC_CODE_MAP
data("BLOSUM80")
?substitution.matrices
data("phiX174Phage")
mymm=which(rowSums(t(consensusMatrix(phiX174Phage))[,1:4]==6)!=1)
mymm
consensusMatrix(phiX174Phage)
head(t(consensusMatrix(phiX174Phage)))
which(rowSums(t(consensusMatrix(phiX174Phage)))==0)
sqrt(4)
sqrt(0.443**2+0.25**2)
20*0.0254
library(Biostrings)
dna=DNAString("CCAGGAACATCCGCATCTT")
reverseComplement(dna)
library(reshape2)
library(ggplot2)
mat=matrix(rnorm(20),5)
mat
m=melt(mat)
m
g=ggplot(m,aes(x=Var1,y=Var2,fill=value))+xlab("X-labels")+ylab("Y-labels")+opts(title="Heatmap Example")
??ggplot
