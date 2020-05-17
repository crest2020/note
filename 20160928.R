library(gplots)
library(devtools)
library(Biobase)
library(RSkittleBrewer)
library(org.Hs.eg.db)
library(AnnotationDbi)
gdURL <- "http://www.stat.ubc.ca/~jenny/notOcto/STAT545A/examples/gapminder/data/gapminderDataFiveYear.txt"
gDat <- read.delim(file = gdURL)
head(gDat)
str(gDat)
opar=par(pch=19)
length(gDat)
rowsum(gDat)
nrow(gDat)
ncol(gDat)
dim(gDat)[2]
a=sample(1:nrow(gDat),8)
a
jdat=gDat[a,]
jdat
jdat=jdat[order(jdat$gdpPercap),]
jdat
jxlim=c(460,60000)
jylim=c(47,82)
plot(lifeExp~gdpPercap,jdat,log='x',col=jcolor,pch=19,main="start my engines")
??plot
?plot
jdat
?plot
?pch
head(colors())
with(jdat,text(x=gdpPercap,y=lifeExp,labels=jcolor,pos=1))
rep(c(1, 3, 1), c(5, 1, 2)) 

rep(c(5,1,2))
?rep
library(RColorBrewer)
display.brewer.pal(n=8,name="Dark2")
jcolor=brewer.pal(n=8,name="Dark2")
library(RSkittleBrewer)
trop=RSkittleBrewer("tropical")
palette(trop)
par(pch=19)
con = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm=bodymap.eset
bm
pdata=pData(bm)
edata=exprs(bm)
fdata=fData(bm)
ls()
pdata
table(pdata$gender)
args(table)
example(table)
?table
showMethods("table")
table(pdata$gender,pdata$race)
summary(edata)
table(pdata$age,"any")
table(is.na(pdata$age))
sum(pdata$age==" ")
sum(is.na(edata))
sum(is.na(pdata$age))
gene_na=rowSums(is.na(edata))
table(gene_na)
head(is.na(edata))
dim(fdata)
dim(pdata)
dim(edata)
edata
head(edata)
rownames(pdata)
colnames(edata)==rownames(pdata)
fdata
head(fdata)
head(edata)
table(rownames(edata)==rownames(fdata))
boxplot(log2(edata+1),col=2)
args(boxplot)
?boxplot
par(mfrow=c(1,2))
plot(density(log2(edata[,1]+1)),col=2)
pallete()
palette()
qqplot(log2(edata[,1]+1),log2(edata[,2]+1),col=3)
mm=log2(edata[,1]+1)-log2(edata[,2]+1)
aa=log2(edata[,1]+1)+log2(edata[,2]+1)
plot(aa,mm,col=2)
edata=as.data.frame(edata)
filt_edata=filter(edata,rowMeans(edata)>1)
plot(density(log2(filt_data[,2]+1)),col=2)
edata
dim(edata)
head(edata)
dim(filt_data)
boxplot(as.matrix(log2(filt_edata+1)),col=2)
dim(edata)
dim(filt_edata)
head(filt_edata)
head(filt_data)
head(edata)
example(filter)
rep(1,3)
head(rowMeans(edata)>1)
?filter
filt_edata=edata[which(rowMeans[edata]>1)]
x <- 1:100
stats::filter(x, rep(1, 3))
x
pexp(6,1/5,lower.tail = F)
args(pexp)
qexp(0.25,1/5)
??axes
attach(mtcars)
opar=par(no.readonly = TRUE)
par(fig=c(0,0.8,0,0.8))
plot(mtcars$wt,mtcars$mpg,xlab="Miles Per Gallon",ylab="Car Weight")
par(fig=c(0,0.8,0.55,1),new=T)
boxplot(mtcars$wt,horizontal=)
table(is.na(pdata$race))
table(pdata$gender==" ")
sum(is.na(edata))
sum(is.na(pdata))
sum(is.na(fdata))
head(pdata)
tail(pdata)
dim(pdata)
tail(pdata,n=10)
pdatanona=pdata[!is.na(pdata)]
class(pdata)
dim(na.omit(pdata))
an=read.table("GTEx.txt",fill = T)
dim(an)
head(an)
an[1,1]
an[1,2]
b=head(an)
write.table(b,"b.csv")
dim(an)
colnames(an)=an[1,]
an2=an[-1,]
dim(an2)
colnames(an2)
an[1,]
an[2,]
smtsd=an$`1397`
write.csv(smtsd,"smtsd.csv")
bind1k=read.csv("binding_tss_genes2.csv",header=T)
dim(bind1k)
head(bind1k)
whole=read.csv("whole name.csv",header=T)
whole1=whole[,1:4]
head(whole1)
bind1kMergeWhole=merge(bind1k,whole1,by="SYMBOL")
head(bind1kMergeWhole)
dim(bind1kMergeWhole)
write.csv(bind1kMergeWhole,file="bind1kMergeWhole.csv")
library(downloader)
library(devtools)
url="https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleMiceWeights.csv"
filename="femaleMiceWeights.csv"
download(url,destfile=filename)
args(download)
?download
dat=read.csv(filename)
install_github("genomicsclass/dagdata")
dir=system.file(package = "dagdata")
list.files(dir)
list.files(file.path(dir,"data"))
filename=file.path(dir,"extdata/femaleMiceweights.csv")
dat=read.csv(filename)
head(dat)
dim(dat)
boxplot(dat)
table(dat)
boxplot(dat[dat$Diet=="chow",2],dat[dat$Diet=="hf",2])
points()
dat[which(dat$Diet=="chow"),]$Bodyweight
?sample
sample(13:24,1)
n=100000
null=vector("numeric",n)
mean(null>=3.02)
mean(dat[dat$Diet=="chow",]$Bodyweight)-mean(dat[dat$Diet=="hf",]$Bodyweight)
library(UsingR)
dir=system.file(package = "UsingR")
list.files(dir)
list.files(file.path(dir,"errata"))
head(father.son)
sample(x,10)
round(sample(x,10),1)
class(x)
?father.son
?round
head(dat)
boxplot(dat$Bodyweight~dat$Diet)
points(dat$Bodyweight~jitter(as.numeric(dat$Diet)),col=as.numeric(dat$Diet))
dat1=dat$Bodyweight
library(rafalib)
B <- 250
N=6
n=0
mypar()
plot(mean(dat1)+c(-7,7),c(1,1),type="n",
     xlab="weight",ylab="interval",ylim=c(1,B))
abline(v=mean(dat1))
for (i in 1:B) {
  chow <- sample(dat1,N,replace = T)
  se <- sd(chow)/sqrt(N)
  interval <- c(mean(chow)-2.57*se, mean(chow)+2.57*se)
  covered <- 
    mean(dat1) <= interval[2] & mean(dat1) >= interval[1]
  color <- ifelse(covered,1,2)
  if(color==1) {n=n+1}
  lines(interval, c(i,i),col=color)
}
table(color)
color(1)
color()
head(colors())
?col
??color
head(palette())
if(color==1) {n=n+1}
n
n=0
n/B
qt(1-0.05/2,df=5)
qnorm(1-0.01/2)
qnorm(0.95)
head(dat)
View(dat)
control=subset(dat[,2],dat$Diet=="chow")
head(control)
table(dat)
treatment=subset(dat[,2],dat$Diet=="hf")
table(dat$Diet)
obsdiff=mean(treatment)-mean(control)
obsdiff/mean(control)
library(downloader)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename <- "femaleControlsPopulation.csv"
if (!file.exists(filename)) download(url,destfile=filename)
population <- read.csv(filename)
population <- unlist(population) # turn it into a numeric
dim(population)
head(population)
n=10000
null=vector("numeric",n)
for(i in 1:n){
  control=sample(population,12)
  treatment=sample(population,12)
  null[i]=mean(treatment)-mean(control)
}
ttestgenerator=function(n){
  treatment1=rnorm(n,mean(treatment),sd(treatment))
  control1=rnorm(n,mean(control),sd(control))
  tstat=(mean(treatment1)-mean(control1))/sqrt(var(treatment1)/n+var(control1)/n)
  return(tstat)
}
ttest=replicate(10000,ttestgenerator(5))
mean(ttest<qt(1-0.05/2,df=2*5-2))
ps=(seq(0,9999)+0.5)/10000
head(ps)

qqplot(ttest,qt(ps,df=8),xlim=c(-6,6),ylim=c(-6,6))
abline(0,1)

library(UsingR)
library(rafalib)
x=father.son$fheight
ps=(seq(0,99)+0.5)/100
qs=quantile(x,ps)
head(cbind(x,qs))
normalqs=qnorm(ps,mean(x),popsd(x))
plot(normalqs,qs)
abline(0,1)
qqnorm(x)
qqline(x)
df=c(3,6,12,30)
mypar(2,2)
for(df1 in df){
  x=rt(1000,df1)
  qqnorm(x);qqline(x)
}

hist(x)
mypar(1,2)
hist(exec.pay)
qqnorm(exec.pay)
qqline(exec.pay)
boxplot(exec.pay)
x=father.son$fheight
y=father.son$sheight
cor(x,y)
?signif
groups=split(y,round(x))
boxplot(groups)
class(groups)
groups
groups$`59`
groups[[1]]
boxplot(control,treatment)
stripchart(control,treatment,method="jitter")
?stripchart
dat=list(control,treatment)
boxplot(dat,outline=F)
stripchart(dat,method="jitter",vertical = T,pch=16,add=T,col=c(1,2))
mypar(1,1)
?boxplot
barplot(dat)
datd=InsectSprays
head(datd)
class(datd)
table(datd$spray)
boxplot(datd)
table(datd)
?boxplot
datdd=split(datd,datd$spray)
datdd$A
boxplot(datdd$A)
head(datd)
boxplot(datd$count~datd$spray)
library(GSE5859Subset)
data(package="GSE5859Subset")
vignette("BSgenomeForge")
help(package="genefilter", help_type="html")
methods(class="ExpressionSet")
methods(class="lm")
whatMethods(pdata)
library(rafalib)
help(package="rafalib")
library("methods")
showMethods(pdata)
?cmdscale
b=c(1,0.5,6,7,9,10,3,5,4,2,6,8,3,7,0.9)
b[order(-b)]
idx
length(b)

library(BiocInstaller)
biocLite("erma")
?system.file
library(erma)
dir=system.file(package = "erma")
list.files(dir)
ef=list.files(file.path(dir,"bed_tabix"),pattern = "bed.gz$")
length(ef)
?dir
head(ef)
?makeErmaSet()
mm
head(colData(mm))
?install.packages
library(RCurl)
browseVignettes("RCurl")
??RCurl
