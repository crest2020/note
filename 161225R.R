con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
mp
pdata=pData(mp)
edata=exprs(mp)
fdata=fData(mp)
head(pdata)
head(edata)
head(fdata)
edatano=edata-rowMeans(edata)
svd1=svd(edatano)
names(svd1)
head(svd1$d)
(svd1$d[1])^2/sum(svd1$d^2)
svd0=svd(edata)
(svd0$d[1])^2/sum(svd0$d^2)
edatalog=log2(edata+1)
svd2=svd(edatalog)
(svd3$d[1])^2/sum(svd3$d^2)
edatalogc=edatalog-rowMeans(edatalog)
svd3=svd(edatalogc)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
bm
edata=exprs(bm)
head(rowSums(edata==0))
head(edata)
?rowSums
head(rowSums(edata))
hist(rowSums(edata==0))
set.seed(3)
prev=1/20
m=50
n=30
d=rbinom(n*m,1,prev)
?rbinom
head(d)
dbinom(n*m,1,prev)
set.seed(1)
qbinom(5,1,0.5)
pnorm(2,0,1)
0.05399097/2
accuracy=0.9
test=rep(NA,n*m)
test[d==1]=rbinom(sum(d==1),1,p=accuracy)
test[d==0]=rbinom(sum(d==0),1,p=1-accuracy)
head(test)
head(d)
table(test==d)
159/(159+1341)
x=rnorm(10,0,1)
rank(x)
x
set.seed(779)

x[1] <- 5
x[2] <- 7

cat("t-test pval:",t.test(x,y)$p.value)
cat("Wilcox test pval:",wilcox.test(x,y)$p.value)
library(rafalib)
mypar(1,2)
stripchart(list(x,y),vertical=TRUE,ylim=c(-7,7),ylab="Observations",pch=21,bg=1)
abline(h=0)
x
y
list(x,y)
xrank=rank(c(x,y))[seq(along=x)]
seq(along=x)
seq(25)
rank(c(x,y))
xrank=rank(c(x,y))[seq(along=x)]
xrank
yrank=rank(c(x,y))[-seq(along=y)]
yrank
stripchart(list(xrank,yrank),vertical = T,ylab="Rank",pch=21,bg=1,cex = 1.25)
ws=sapply(x, function(z) rank(c(z,y))[1]-1)
text(rep(1.05,length(ws)),xrank,ws,cex=0.8)
rank(c(x,y))[1]
ws
W <-sum(ws) 
W
n1<-length(x);n2<-length(y)
Z <- (mean(ws)-n2/2)/ sqrt(n2*(n1+n2+1)/12/n1)
print(Z)
wilcox.test(x,y)
library(AnnotationHub)
ah=AnnotationHub()
query(ah,c("nrf1","narrowpeak"))[["AH25951"]]
da=da[[1]]
da
da[da$qVaue<0.05,]
query(ah,"atf")
prox=read.table("GSM1477597_PROX1vsInput_hg18.bed.peak.txt")
head(prox)
library(GenomicRanges)
gr=GRanges(seqnames = prox$V1,ranges = IRanges(start = prox$V2,end=prox$V3))
mcols(gr)$score=prox$V5
gr
library(BiocInstaller)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Rnorvegicus.UCSC.rn5.refGene)
biocLite("TxDb.Hsapiens.UCSC.hg18.knownGene")
library(TxDb.Hsapiens.UCSC.hg18.knownGene)
txdb18=TxDb.Hsapiens.UCSC.hg18.knownGene
prom=promoters(genes(txdb18))
prom
ov=subsetByOverlaps(prom,gr)
ov
rownames(prom)=mcols(prom)$gene_id
??MethylSet
class(gr)
entrzid=mcols(ov)$gene_id

library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
se=select(org.Hs.eg.db,keys=entrzid,keytype="ENTREZID",columns = c("ENTREZID","GENENAME","SYMBOL"))
write.csv(se,file="prox1chipseq.csv")
?write.csv
library(Biostrings)
dnastring=DNAString("CCACCGCCGCGTTTGCAGAGAAAAGGCGC")
reverseComplement(dnastring)
dn=DNAString("GCGTTTGCAGAGAAAAGGCGCACCACC")
length(dn)
library(ggplot2)
p=ggplot(diamonds,aes(carat,price,colour=cut))
p
summary(p)
library(ggplot2)
qplot(displ,hwy,data=mpg,colour=factor(cyl))
p=ggplot(diamonds,aes(carat,price,colour=cut))
p=p+layer(geom="point",position = "identity",stat="identity")
p
p=ggplot(diamonds,aes(x=carat))
p=p+layer(geom="bar",params=list(fill="steelblue",binwidth=0.5),stat="bin",position="identity")
?layer
p
ggplot(msleep,aes(sleep_rem/sleep_total,awake))+geom_point()
qplot(sleep_rem/sleep_total,awake,data=msleep)+geom_smooth()
qplot(sleep_rem/sleep_total,awake,data=msleep,geom = c("point","smooth"))
p=ggplot(msleep,aes(sleep_rem/sleep_total,awake))
summary(p)
library(scales)
bestfit=geom_smooth(method="lm",se=F,colour=alpha("steelblue",0.5),size=2)
qplot(sleep_rem,sleep_total,data=msleep)+bestfit
qplot(awake,brainwt,data=msleep,log="y")+bestfit
qplot(bodywt,brainwt,data=msleep,log="xy")+bestfit
p=ggplot(mtcars,aes(mpg,wt,colour=cyl))+geom_point()
p
mtcars=transform(mtcars,mpg=mpg^2)
p %+% mtcars
aes(x=weight,y=heigt,colour=age)
p=ggplot(mtcars)
summary(p)
p=p+aes(wt,hp)
p
summary(p)
p+geom_point()
p+geom_point(aes(colour=factor(cyl)))
p+geom_point(aes(y=disp))
p
p=ggplot(mtcars,aes(mpg,wt))
p+geom_point(aes(colour="darkblue"))
p=ggplot(Oxboys,aes(age,height,group=Subject))+geom_line()
library(nlme)
p
p+geom_smooth(aes(group=Subject),method = "lm",se=F)
p+geom_smooth(aes(group=1),method="lm",se=F)
boysbox=ggplot(Oxboys,aes(Occasion,height))+geom_boxplot()
boysbox+geom_line(aes(group=Subject),colour="#3366FF")
xgrid=with(df,seq(min(x),max(x),length=50))
head(df)
df=data.frame(x=c(3,1,5),y=c(2,4,6),label=c("a","b","c"))
df
p=ggplot(df,aes(x,y))+xlab(NULL)+ylab(NULL)
p+geom_point()+labs(title="geom_point")
p+geom_bar(stat="identity")+labs(title="geom_bar(stat=\"identity\")")
p+geom_line()+labs(title="geom_line")
p+geom_area()+labs(title="geom_area")
p+geom_path()+labs(title="geom_path")
p+geom_text(aes(label=label))+labs(title="geom_text")
p+geom_tile()+labs(title="geom_tile")
p+geom_polygon()+labs(title="geom_polygon")
depth_dist=ggplot(diamonds,aes(depth))+xlim(58,68)
depth_dist+geom_histogram(aes(y=..density..),binwidth = 0.1)+facet_grid(.~cut)
depth_dist+geom_histogram(aes(fill=cut),binwidth = 0.1,position="fill")
depth_dist+geom_freqpoly(aes(y=..density..,colour=cut),binwidth=0.1)
library(plyr)
qplot(cut,depth,data=diamonds,geom="boxplot")
qplot(carat,depth,data=diamonds,geom="boxplot",group=round_any(carat,0.1,floor),xlim=c(0,3))
qplot(class,cty,data=mpg,geom="jitter")
qplot(class,drv,data=mpg,geom="jitter")
qplot(depth,data=diamonds,geom="density",xlim=c(54,70))
qplot(depth,data=diamonds,geom="density",xlim=c(54,70),fill=cut,alpha=I(0.2))
df=data.frame(x=rnorm(2000),y=rnorm(2000))
norm=ggplot(df,aes(x,y))
norm+geom_point()
norm+geom_point(shape=1)
norm+geom_point(shape=".")
norm+geom_point(colour="black",alpha=1/3)
norm+geom_point(colour="black",alpha=1/10)
td=ggplot(diamonds,aes(table,depth))+xlim(50,70)+ylim(50,70)
td+geom_point()
td+geom_jitter()
jit=position_jitter(width = 0.5)
td+geom_jitter(position = jit,colour="black",alpha=1/50
               
library(Biostrings)
dna=DNAString("CGAGCCAGGAGGAGCGTGTTGAAAAGGCGC")
reverseComplement(dna)



library(ggplot2)
td=ggplot(diamonds,aes(table,depth))+xlim(50,70)+ylim(50,70)
td+geom_point()
td+geom_jitter()
jit=position_jitter(width=0.5)
td+geom_jitter(position = jit,colour="black",alpha=1/200)
d=ggplot(diamonds,aes(carat,price))+xlim(1,3)+theme(legend.position="none")
d+stat_bin2d()
d+stat_bin2d(bins=100)
d+stat_bin2d(binwidth=c(0.02,200))
d+stat_binhex()
install.packages("hexbin")
library(hexbin)
d+stat_binhex(bins=10)
d+stat_binhex(binwidth = c(0.02,200))
d=ggplot(diamonds,aes(carat,price))+xlim(1,3)+theme(legend.position="none")
d+geom_point()+geom_density2d()
d+stat_density2d(geom = "tile",aes(size=..density..),contour=F)+scale_size_area()
last_plot()+scale_fill_gradient(limits=c(1e-5,8e-4))
 library(maps)
install.packages("maps")
data(package="maps")
data(us.cities)
head(us.cities)
big_cities=subset(us.cities,pop>500000)
qplot(long,lat,data=big_cities)+borders("world",size=0.5)
?borders
library(maps)
library(ggplot2)
data("us.cities")
big_cities=subset(us.cities,pop>500000)
big_cities
qplot(long,lat,data=big_cities)+borders("state",size=0.5)
tx_cities=subset(us.cities,country.etc=="TX")
ggplot(tx_cities,aes(long,lat))+borders("county","texas",colour = "grey70")+geom_point(colour="black",alpha=0.5)
library(maps)
states=map_data("state")
arrests=USArrests
names(arrests)=tolower(names(arrests))
arrests$region=tolower(rownames(USArrests))
choro=merge(states,arrests,by="region")
head(choro)
qplot(long,lat,data=choro,group=group,fill=assault,geom="polygon")
qplot(long,lat,data=choro,group=group,fill=assault/murder,geom="polygon")
library(plyr)
ia=map_data("county","iowa")
mid_range=function(x) mean(range(x,na.rm=T))
centres=ddply(ia,.(subregion),colwise(mid_range,.(lat,long)))
?ddply
ggplot(ia,aes(long,lat))+geom_polygon(aes(group=group),fill=NA,colour="grey60")+
geom_text(aes(label=subregion),data=centres,size=2,angle=45)
library(BiocInstaller)
biocLite("gmapR")
n
library()
install.packages("")
source("https://bioconductor.org/biocLite.R")
biocLite("gmapR",type="source")
?biocLite
biocLite("Rtools")
n
library(devtools)
install_git("Rtools")
install_github("Rtools")
library(gmapR)
library(BiocInstaller)
install.packages('gmapR_1.16.0.tar.gz',repos = NULL,type="source")
y
biocLite("DESeq2")
library(DESeq2)
library(BiocInstaller)
biocLite("gmapR")
biocLite("VariantAnnotation")
library(VariantAnnotation)
library(Biostrings)
dna=DNAString("ATCGATCGATCG")
method(dna)
methods(dna)
?DNAString
translate(dna)
library(VariantAnnotation)
library(org.Hs.eg.db)
library(BiocInstaller)
biocLite("cgdv17")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
biocLite('PolyPhen.Hsapiens.dbSNP131')
file=system.file('vcf','NA06985_17.vcf.gz',package = 'cgdv17')
hdr=scanVcfHeader(file)
head(hdr)
hdr
info(hdr)
geno(hdr)
colnames(meta(hdr))
meta(hdr)$META
genesym=c("TRPV1","TRPV2","TRPV3")
geneid=select(org.Hs.eg.db,keys = genesym,keytype="SYMBOL",columns = "ENTREZID")
geneid
txdb=TxDb.Hsapiens.UCSC.hg19.knownGene
txdb
txdb=renameSeqlevels(txdb,gsub("chr","",seqlevels(txdb)))
txdb
txdb=keepSeqlevels(txdb,"17")
txbgene=transcriptsBy(txdb,"gene")
txbgene
gnrng=unlist(range(txbgene[geneid$ENTREZID]),use.names = F)
gnrng
names(gnrng)=geneid$SYMBOL
gnrng
param=ScanVcfParam(which=gnrng,info="DP",geno = c("GT","cPd"))
param
vcf=readVcf(file,"hg19",param)
vcf
head(fixed(vcf))
geno(vcf)
cds=locateVariants(vcf,txdb,CodingVariants())
five=locateVariants(vcf,txdb,FiveUTRVariants())
splice=locateVariants(vcf,txdb,SpliceSiteVariants())
intron=locateVariants(vcf,txdb,IntronVariants())
all=locateVariants(vcf,txdb,AllVariants())
table(sapply(split(mcols(all)$GENEID,mcols(all)$QUERYID),function(x) length(unique(x))>1))
head(all)
head(split(mcols(all)$GENEID,mcols(all)$QUERYID))
a=c(1,2,3)
b=c(4,5,6)
split(a,b)
table(sapply(split(a,b),function(x) length(unique(x))>1))
idx=sapply(split(mcols(all)$QUERYID,mcols(all)$GENEID),unique)
idx
sapply(names(idx),function(nm){
  d=all[mcols(all)$GENEID %in% nm,c("QUERYID","LOCATION")]
  table(mcols(d)$LOCATION[duplicated(d)==FALSE])
})
names(idx)
seqlevelsStyle(vcf)="UCSC"
seqlevelsStyle(txdb)='UCSC'
aa=predictCoding(vcf,txdb,Hsapiens)
table(sapply(split(mcols(aa)$GENEID,mcols(aa)$QUERYID),function(x) length(unique(x))>1))
idx=sapply(split(mcols(aa)$QUERYID,mcols(aa)$GENEID,drop=T),unique)
idx
sapply(idx,length)
sapply(names(idx),function(nm){
  d=aa[mcols(aa)$GENEID %in% nm,c("QUERYID","CONSEQUENCE")]
  table(mcols(d)$CONSEQUENCE[duplicated(d)==F])
})

biocLite("ensemblVEP")
library(ensemblVEP)
dest=tempfile()
writeVcf(vcf,dest)
gr=ensemblVEP(file=dest)
pnorm(-2)
pnorm(2)
pnorm(-2)+1-pnorm(2)
pnorm(-1.96)
library(dplyr)
dat=read.csv('mice_pheno.csv')
controlPopulation=filter(dat,Sex=="F"&Diet=='chow') %>% select(Bodyweight) %>% unlist
hfPopulation=filter(dat,Sex=="F"&Diet=="hf") %>%select(Bodyweight) %>% unlist
head(controlPopulation)
head(hfPopulation)
library(rafalib)
mypar(1,2)
hist(hfPopulation,xlim = c(10,50),breaks = 6)
hist(controlPopulation,xlim=c(10,50),breaks = 6)
?hist
qqnorm(hfPopulation)
qqline(hfPopulation)
qqnorm(controlPopulation)
qqline(controlPopulation)
library(downloader)
dat=na.omit(dat)
length(dat)
head(dat)
dim(dat)
type(dat)
class(dat)
1-2*pnorm(-1)
1-2*pnorm(-2)
1-2*pnorm(-3)
1-2*pnorm(-3*popsd(controlPopulation))
sum(controlPopulation-mean(controlPopulation))/popsd(controlPopulation)
mean(controlPopulation)
popsd(controlPopulation)
head(dat)
y=filter(dat,Sex=="M"&Diet=="chow")%>%select(Bodyweight)%>%unlist
head(y)
popsd(y)
qqplot(y)
qqnorm(y)
qqline(y)
length(y)
avgs=replicate(10000,mean(sample(y,25)))
hist(avgs)
qqnorm(avgs)
qqline(avgs)
sd(avgs)
popsd(y)
4.42/25
4.42/5
sdd=replicate(10000,sd(sample(y,25)))
sum(sdd<3.5)
mean(sdd<3.5)
x=seq(0.0001,0.9999,len=300)
qqnorm(x)
mypar(2,2)
plot(qnorm(x),qt(x,df=100))
abline(0,4)
ablien
xx=seq(-0.5,0.5,len=300)
plot(xx,qt(x,df=100))
abline(0,1)
qqline(xx)
qqnorm(xx)
qqline(xx)
plot(qnorm(xx),qt(xx,df=20))
head(qnorm(xx))
library(biomaRt)
listMarts()
ensembl=useMart("ensembl",dataset="hsapiens_gene_ensembl")
ensembl
attr=listAttributes(ensembl)
head(attr)
?listAttributes
dim(attr)
head(attr,n=30)
filters=listFilters(ensembl)
head(filters)
entrez=c('673','837')
getBM(attributes = c('entrezgene','go_id'),filters = 'entrezgene',values=entrez,mart=ensembl)
goids


write.csv(filters,"filters.csv")
go=c("GO:0000114")
chrom=c(17,20)
hh=getBM(attributes = c("hgnc_id"),filters = c('chromosome_name','go_id'),values=list(chrom,go),mart = ensembl)
head(hh)
head(hh)
getBM(attributes = )
chrom
go
head(filters)
find(listAttributes(ensembl),'')
?find
?grep
grep('hgnc_symbol',listAttributes(ensembl),value = T)

head(listAttributes(ensembl),n=50)
write.csv(listAttributes(ensembl),'attributes.csv')
affyids=c("202763_at","209310_s_at","207500_at")
getBM(attributes = c("affy_hg_u133_plus_2","entrezgene"),filters='affy_hg_u133_plus_2',values=affyids,mart=ensembl)
getBM(attributes = c("affy_hg_u133_plus_2","hgnc_symbol","chromosome_name","start_position","end_position","band"),filters='affy_hg_u133_plus_2',values=affyids,mart=ensembl)
getBM(c("entrezgene","hgnc_symbol"),filters='go_id',values='GO:0004707',mart=ensembl)
refsedids=c("NM_005359","NM_000546")
ipro=getBM(attributes = c('refseq_mrna','interpro','interpro_description'),filters = 'refseq_mrna',values=refsedids,mart=ensembl)
ipro
getBM(c('affy_hg_u133_plus_2','ensembl_gene_id'),filters = c('chromosome_name','start','end'),values=list(16,1100000,1250000),mart=ensembl)
entrez=c('673','7157','837')
?getSequence
getSequence(id=entrez,type='entrezgene',seqType = "peptide",upstream = 100,mart=ensembl)
listMarts(archive = TRUE)
ensemblold=useMart('ensembl_mart_46',dataset = 'hsapiens_gene_ensembl',archive = TRUE)
library(dnastring)
library(Biostrings)
dna=DNAString("CACACCAGCGCGGCTCGGTGaaaAAGCGC")
reverseComplement(dna)
library(pwr)
install.packages("pwr")
library(pwr)
pwr.t.test(d=0.8,sig.level = 0.05,power=0.9,type="two.sample")
pwr.t.test(n=20,d=0.5,sig.level = 0.01,type="two.sample")
pwr.t2n.test(n1=20,n2=39,sig.level = 0.05,power=0.9,alternative ="two.sided" )
?pwr.t2n.test
pwr.anova.test(k=5,f=0.25,sig.level = 0.05,power=0.8)
?pwr.anova.test
pwr.r.test(r=0.25,sig.level = 0.05,power=0.9,alternative = "greater")
pwr.f2.test(u=3,f2=0.0769,sig.level = 0.05,power=0.9)                                                                                                                            
library(vcd)
install.packages("vcd")
mosaic(Titanic,shade = TRUE,legend=TRUE)
score=c(40,57,45,55,58,57,64,55,62,65)
treatment=seq("A",5)
?seq
seq("A",5)
treatment=c(rep("A",5),rep("B",5))
data=data.frame(score,treatment)
plot(data~treatment)
plot(data,pch=19)
plot(score~treatment)
treatment=as.factor(treatment)
hist(data)
library(coin)
cholesterol
install.packages("coin")
library(lmPerm)
install.packages("lmPerm")
library(coin)
score
t.test(score~treatment,data=data,var.equal=T)
oneway_test(score~treatment,data=data,distribution=approximate(B=1000))
?oneway_test
library(MASS)
UScrime=transform(UScrime,So=factor(So))
UScrime
wilcox_test(Prob~So,data=UScrime,distribution="exact")
library(multcomp)
set.seed(1234)
a=data("cholesterol")
a
oneway_test(response~trt,data=cholesterol,distribution=approximate(B=9999))
anova(response~trt)
library(UsingR)
library(rafalib)
x=father.son$fheight
ps=(seq(0,99)+0.5)/100
ps
qs=quantile(x,ps)
?quantile
qs
normalqs=qnorm(ps,mean(x),popsd(x))
ps
normalqs
plot(normalqs,qs,xlab="Normal percemtiles",ylab="Height quantiles")
abline(0,1)
qqnorm(x)
abline(0,1)
qqline(x)
n=1000
x=rnorm(n)
qqnorm(x)
qqline(x)
dfs=c(3,6,12,30)
mypar(2,2)
for(df in dfs){
  x=rt(1000,df)
  qqnorm(x,xlab="t quantiles",main=paste0("d.f=",df),ylim=c(-6,6))
  
}
myoar(1,2)
hist(exec.pay)
qqnorm(exec.pay)
qqline(exec.pay)
boxplot(exec.pay,ylab="10,000s of dollars",ylim=c(0,400))
library(multcomp)
attach(cholesterol)
table(trt)
head(cholesterol)
aggregate(response,by=list(trt),FUN=sd)
fit=aov(response~trt)
summary(fit)
library(gplots)
plotmeans(response~trt,xlab="Treatment",ylab="Response",main="Mean plot\nwith 95% CI")
detach(cholesterol)
TukeyHSD(fit)
?TukeyHSD
library(multcomp)
par(mar=c(5,4,6,2))
tuk=glht(fit,linfct=mcp(trt="Tukey"))
?glht
?mcp
plot(cld(tuk,level=0.05),col="lightgrey")
?cld
library(car)
library("car")
install.packages("car")
library(car)
qqPlot(lm(response~trt,data=cholesterol),simulate=TRUE,main="Q-Q Plot",labels=FALSE)
bartlett.test(response~trt,data=cholesterol)
contrast=rbind("no drug vs. drug"=c(3,-1,-1,-1))
contrast
data(litter,package="multcomp")
attach(litter)
table(dose)
aggregate(weight,by=list(dose),FUN=mean)
head(litter)
fit=aov(weight~gesttime+dose)
summary(fit)
library(effects)
install.packages("effects")
effect("dose",fit)
TukeyHSD(fit)
summary(glht(fit,linfct=mcp(dose=contrast)))
fit2=aov(weight~gesttime*dose,data=litter)
summary(fit2)
library(HH)
install.packages("HH")
library(HH)
ancova(weight~gesttime+dose,data=litter)
detach("package:HH")
detach("HH")
?detach
install.packages("titanic")
install.packages("faraway")
library(titanic)
library(faraway)
data("worldcup")
library(ggplot2)
ggplot(data=titanic_test,aes(x=Fare))+geom_histogram()
head(Titanic)
ggplot(data=titanic_test)+geom_histogram(aes(x=Fare))
ggplot()+geom_histogram(data=titanic_test,aes(x))
ggplot(worldcup,aes(x=Time,y=Passes))+geom_point()
ggplot(worldcup,aes(x=Time,y=Passes,color=Position,size=Shots))+geom_point()
head(worldcup,n=20)
args(ggplot)
?ggplot
?mutate()
library(dplyr)
noteworthy_players=worldcup%>%filter(Shots==max(Shots)|Passes==max(Passes))%>%mutate(point_label=paste(Team,Position,sep=", "))
head(noteworthy_players)
ggplot(worldcup,aes(x=Passes,y=Shots))+geom_point()+geom_text(data=noteworthy_players,aes(label=point_label),vjust="inward",hjust="inward")
?geom_text
ggplot(worldcup,aes(x=Time))+geom_histogram(binwidth=10)+geom_vline(xintercept = 90*0:6,color="blue",alpha=0.5)
ggplot(worldcup,aes(x=Time,y=Passes))+geom_point(color="darkgreen")
group=factor(c(1,1,2,2))
model.matrix(~group)
model.matrix(formula(~group))
g=c(1,1,2,2)
model.matrix(~g)
group=factor(c("control","control","highfat","highfat"))
model.matrix(~group)
group=factor(c(1,1,2,2,3,3))
model.matrix(~group)
model.matrix(~group+0)
diet=factor(c(1,1,1,1,2,2,2,2))
sex=factor(c("f","f","m","m","f","f","m","m"))
table(diet,sex)
model.matrix(~diet+sex+0+0)
model.matrix(~diet+sex+diet:sex)
model.matrix(~diet*sex)
group=factor(c(1,1,2,2))
group
group=relevel(group,"2")
model.matrix(~group)
group=factor(c("control","control","hf","hf","hf"))
model.matrix(~group)
group=relevel(group,"hf")
group
model.matrix(~group)
group=factor(group,levels = c("control","hf"))
model.matrix(~group)
group=
group
model.matrix(~group,data=as.factor((group=5:8)))
tt=seq(0,3.4,len=4)             
tt
model.matrix(~tt+I(tt^2))
dat=read.csv("femaleMiceWeights.csv")
stripchart(dat$Bodyweight~dat$Diet,vertical=TRUE,method="jitter",main="Body weight over diet")
levels(dat$Diet)
levels(dat$Bodyweight)
x=model.matrix(~Diet,data=dat)
x
y=dat$Bodyweight
solve(crossprod(x))%*%crossprod(x,y)
s=split(dat$Bodyweight,dat$Diet)
s
mean(s[["chow"]])
mean(s[["hf"]])
fit=lm(Bodyweight~Diet,data=dat)
fit
summary(fit)
summary(fit)$coefficients
ttest=t.test(s[["hf"]],s[["chow"]],var.equal=TRUE)
summary(fit)$coefficients[2,3]
ttest$statistic
set.seed(1)
B=1000
h0=56.67
v0=0
g=9.8
n=25
tt=seq(0,3.4,len=n)
x=cbind(1,tt,tt^2)
x
tt
A=solve(crossprod(x))%*%t(x)
A
y=h0+v0*tt-0.5*g*tt^2+rnorm(n,sd=1)
y
betahats=A%*%y
betahats
betahat=replicate(B,{
 y=h0+v0*tt-0.5*g*tt^2 +rnorm(n,sd=1)
 betahats=A%*%y
 return(betahats[3])
})

head(betahat)
library(rafalib)
mypar(1,2)
hist(betahat)
qqnorm(betahat)
qqline(betahat)
mean(betahat)
sd(betahat)
1/sqrt(25)



library(UsingR)
x=father.son$fheight
y=father.son$sheight
n=length(y)
n
N=50
B=1000
betahat=replicate(B,{
  index=sample(n,N)
  sampledat=father.son[index,]
  x=sampledat$fheight
  y=sampledat$sheight
  lm(y~x)$coef
})
head(betahat)
dim(betahat)
betahat=t(betahat)
head(betahat)
qqnorm(betahat[,1])
qqline(betahat[,1])
qqnorm(betahat[,2])
qqline(betahat[,2])
cor(betahat[,1],betahat[,2])
mean((betahat[,1]-mean(betahat[,1]))*(betahat[,2]-mean(betahat[,2])))


n=nrow(father.son)
N=50
index=sample(n,N)
sampledat=father.son[index,]
x=sampledat$fheight
y=sampledat$sheight
X=model.matrix(~x)
head(X)
N=nrow(X)
p=ncol(X)
XtXinv=solve(crossprod(X))
resid=y-X%*%XtXinv%*%crossprod(X,y)
s=sqrt(sum(resid^2)/(N-p))
ses=sqrt(diag(XtXinv))*s
summary(lm(y~x))$coef[,2]
ses
head(s)
x=father.son$fheight
beta=c(34,0.5)
var(beta[1]+beta[2]*x)
n=length(tt)
n
y=h0+v0*tt-0.5*g*tt^2+rnorm(n,sd=1)
var(y)
spider=read.csv("spider_wolff_gorb_2013.csv",skip=1)
boxplot(spider$friction~spider$type*spider$leg,col=c("grey90","grey40"),las=2,main="Comparison of friction coefficients of different leg pairs")
table(spider$leg,spider$type)
boxplot(spider$friction~spider$type*spider$leg,col=c("grey90","grey40"),las=2,main="Comparison of friction coefficients of different leg pairs")
spider.sub=spider[spider$leg=="L1",]
head(spider.sub)
fit=lm(friction~type,data=spider.sub)
summary(fit)$coefficients
s=split(spider.sub$friction,spider.sub$type)
s
mean(s[["pull"]])
mean(s[["push"]])
X=model.matrix(~type,data=spider.sub)
head(X)
length(X)
tail(X)
library(rafalib)
imagemat(X,main="Model matrix for linear model")
X=model.matrix(~type+leg,data=spider)
colnames(X)
head(X)
tail(X)
imagemat(X,main="Model matrix for liner model with two factors")
fitTl=lm(friction~type+leg,data=spider)
summary(fitTl)$coef
Y=spider$friction
X=model.matrix(~type+leg,data=spider)
beta.hat=solve(crossprod(X))%*%crossprod(X,Y)
t(beta.hat)
library(contrast)
L3vsL2=contrast(fitTl,list(leg="L3",type="pull"),list(leg="L2",type="pull"))
L3vsL2
sigma.hat=sum(fitTl$residuals^2)/(nrow(X)-ncol(X))*solve(crossprod(X))
signif(sigma.hat,2)
fitX=lm(friction~type+leg+type:leg,data=spider)
summary(fitX)
coef(fitX)
L2push.vs.pull=contrast(fitX,list(leg="L2",type="push"),list(leg="L2",type="pull"))
L2push.vs.pull
coefs=coef(fitX)
coefs
coefs[2]+coefs[6]
library(multcomp)
c=matrix(c(0,0,0,0,0,-1,1,0),1)
c
B<-10000
minpval <- replicate(B, min(runif(10000,0,1))<0.01)
mean(minpval>=1)
table(min(runif(10000,0,1))<0.1)
p=10^-7
p
N=5*10^6
winners=rbinom(1000,N,p)
winners
tab=table(winners)
tab
plot(tab)
prop.table(tab)
N=10000
lamdas=2^seq(1,16,len=N)
lamdas
y=rpois(N,lamdas)
x=rpois(N,lamdas)
sum(y>0)
sum(x>0)
ind=which(y>0 & x>0)
library(rafalib)
splot(log2(lamdas),log2(y/x),subset=ind)
library(parathyroidSE)
data("parathyroidGenesSE")
se=parathyroidGenesSE
dim(se)
class(se)
se
x=assay(se)[,23]
y=assay(se)[,24]
ind=which(y>0 & x>0)
splot((log(x)+log(y))/2,log(x/y),subset=ind)
library(matrixStats)
vars=rowVars(assay(se)[,c(2,8,16,21)])
means=rowMeans(assay(se)[,c(2,8,16,21)])
splot(means,vars,log="xy",subset=which(means>0&vars>0))
abline(0,1,col=2,lwd=2)
datadir="http://www.biostat.jhsph.edu/bstcourse/bio751/data"
x=read.csv(file.path(datadir,"hcmv.csv"))[,2]
breaks=seq(0,4000*round(max(x)/4000),4000)
tmp=cut(x,breaks)
counts=table(tmp)
counts
hist(counts)
l=function(lambda) sum(dpois(counts,lambda,log=TRUE))
lambdas=seq(3,7,len=100)
ls=exp(sapply(lambdas,l))
plot(lambdas,ls,type="l")
mle=optimize(l,c(0,10),maximum = T)
abline(v=mle$maximum)
?dpois
head(ls)
theoretical=qpois((seq(0,99)+0.5)/100,mean(counts))
qqplot(theoretical,counts)
abline(0,1)
hist(counts)
mean(counts)
sd(counts)
library(Biobase)
library(maPooling)
data("maPooling")
pd=pData(maPooling)
head(pd)
strain=factor(as.numeric(grepl("b",rownames(pd))))
head(strain)
a=c("chaoeeia","gvharg","ghr","j","k")
grep("vh",a)
grepl("vh",a)
dim(pd)
pooled=which(rowSums(pd)==12&strain==1)
pd
pooled
thchreps=exprs(maPooling[,pooled])
individuals=which(rowSums(pd)==1&strain==1)
individuals=individuals[-grep("tr",names(individuals))]
bioreps=exprs(maPooling)[,individuals]
library(matrixStats)
techsds=rowSds(thchreps)
biosds=rowSds(bioreps)
library(rafalib)
shist(biosds,unit=0.1,col=1,xlim=c(0,1.5))
shist(techsds,unit=0.1,col=2,add=T)
legend("topright",c("Biological","Technical"),col=c(1,2),lty=c(1,1))
?shist
qqnorm(biosds)
qqline(biosds)
mypar(3,3)
sds=seq(0,2,len=100)
for( d in c(1,5,10)){
  for(s0 in c(0.1,0.2,0.3)){
    tmp=hist(biosds,main=paste("s_0=",s0,"d= ",d),xlab="sd",ylab="density",freq=F,nc=100,xlim=c(0,1))
    dd=df(sds^2/s0^2,11,d)
    k=sum(tmp$density)/sum(dd)
    lines(sds,dd*k,type="l",col=2,lwd=2)
  }
}
library(limma)
estimates=fitFDist(biosds^2,11)
theoretical=sqrt(qf((seq(0,999)+0.5)/1000,11,estimates$df2)*estimates$scale)
qf(1,11,12)
qf(0.99,11,12)
observed=biosds
mypar(1,2)
qqplot(theoretical,observed)
abline(0,1)
tmp=hist(biosds,main=paste("s_0=",signif(estimates[[1]],2),"d =",signif(estimates[[2]],2)),
         xlab="sd",ylab="density",freq=F,nc=100,xlim=c(0,1),ylim=c(0,9))
dd=df(sds^2/estimates$scale,11,estimates$df2)
k=sum(tmp$density)/sum(dd)
lines(sds,dd*k,type="l",col=2,lwd=2)
head(dd)
head(sds)
length(sds)
sds
sd(sds)
sum(dd)
sum(tmp$density)
lines(sds,dd,col=3,lwd=3)
