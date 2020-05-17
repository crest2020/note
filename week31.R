m<-matrix(1:10,2,5)
m[-1,]
m[0,-1]
m[,-1]

data(litter,package = "multcomp")
attach(litter)
table(dose)

aggregate(weight,by=list(dose),FUN=mean)
fit<-aov(weight~gesttime+dose)
summary(fit)

library(multcomp)
contrast<-rbind("no drug vs.drug"=c(3,-1,-1,-1))

contrast
summary(glht(fit,linfct=mcp(dose=contrast)))

contrast1<-rbind("no drug vs.drug"=c(1,-1,-1,-1))
summary(glht(fit,linfct=mcp(dose=contrast1)))
head(litter)
fit
man(mtcars)
library(BiocInstaller)



set.seed(100)
B = 1000
plessthan = replicate(B,{
  pvals = replicate(20,{
    cases = rnorm(10,30,2)
    controls = rnorm(10,30,2)
    t.test(cases,controls)$p.value
  })
  sum(pvals<=0.05)
})
table(plessthan) ##just for illustration
mean(plessthan)



library(Biobase)
library(genefu)
data(nkis)
dim(demo.nkis)
head(demo.nkis)[,1:8]
head(annot.nkis)[,1:8]

nkes = ExpressionSet(t(data.nkis), phenoData=AnnotatedDataFrame(demo.nkis),
                     featureData=AnnotatedDataFrame(annot.nkis))
head(nkes)
install.packages("devtools")
library(devtools)
install_github("genomicsclass/ph525x")
install.packages("Homo.sapiens")


library(GSE5859Subset)
library(Biobase)
.oldls = ls()
data(GSE5859Subset)
.newls = ls()
head(.oldls)
install.packages("BiocInstaller")
library(BiocInstaller)

biocLite("BSgenome.Hsapiens.UCSC.hg19")
biocLite("SNPlocs.Hsapiens.dbSNP.20120608")
biocLite("Homo.sapiens")
library(devtools)
install_github("genomicsclass/ph525x")

library("BSgenome.Hsapiens.UCSC.hg19")
library("Homo.sapiens")
Hsapiens$chr16
nchar(Hsapiens$chr16)
library(BiocInstaller)
biocLite("genefu")

library(genefu)
data(sig.gene70)
dim(sig.gene70)
head(sig.gene70)[,1:6]

ncbi<-sig.gene70$NCBI.gene.symbol
ncbi

write.table(ncbi,"a.txt")
write.csv(ncbi,"b.csv")
write.csv(sig.gene70,"siggene70.csv")
sig.gene70
colnames(sig.gene70)
grep("kinase",sig.gene70$Description,fixed=TRUE)
library(BiocInstaller)
biocLite("COPDSexualDimorphism.data")

library(COPDSexualDimorphism.data)
data(lgrc.expr.meta)
head(expr.meta)
nrow(expr.meta)
sum(expr.meta$gender=="2-Female")
table(expr.meta$gender)
median(expr.meta$pkyrs)
qqplot(expr.meta$pkyrs)
qqnorm(expr.meta$pkyrs)
abline(0,1)
boxplot(pkyrs~gender, data=expr.meta)

expr.meta$pyp1 = expr.meta$pkyrs+1
library(MASS)
lm1 = lm(pyp1~gender, data=expr.meta)
boxcox(lm1)


library(BSgenome.Hsapiens.UCSC.hg19)
chr11seq<-BSgenome.Hsapiens.UCSC.hg19[["chr11"]]
subseq(chr11seq,start=10^6,width=25)
?countPattern
pattern<-c("ATG","TGA","TAA","TAG")

for (i in pattern){print(i) }
?print
countPattern("TAG",chr11seq)

?alphabetFrequency
chr7seq<-BSgenome.Hsapiens.UCSC.hg19[["chr7"]]
alphabetFrequency(chr7seq)
example("alphabetFrequency")
alphabetFrequency(chr7seq, as.prob=TRUE)
library(BiocInstaller)
biocLite("SNPlocs.Hsapiens.dbSNP144.GRCh37")


library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
snps144 <- SNPlocs.Hsapiens.dbSNP144.GRCh37
s17 <- snplocs(snps144, "ch17")
head(s17)
location<-s17$loc[s17$RefSNP_id==73971683]
location

biocLite("gwascat")
library(gwascat)
data(ebicat37)
ebicat37
write.csv(sort(table(ebicat37$DISEASE.TRAIT),decreasing=TRUE),"c.csv")
library(tissuesGeneExpression)

library(devtools)
install_github("genomicsclass/tissuesGeneExpression")

library(tissuesGeneExpression)
data(tissuesGeneExpression)
head(e[,1:5])
table(tissue)
qq<-e[rownames(e)=="209169_at",]

boxplot(e["209169_at",]~tissue,las=2)
IDs = c("201884_at", "209169_at", "206269_at", "207437_at", "219832_s_at", "212827_at")
data(tissuesGeneExpression)
boxplot(e["212827_at",]~tissue,las=2)

library(BiocInstaller)
biocLite("hgu133aprobe")
library(hgu133aprobe)

head(hgu133aprobe)
hgu133aprobe[hgu133aprobe$Probe.Set.Name=="206269_at",]

library(GSE5859Subset)
library(devtools)
install_github("genomicsclass/GSE5859Subset")
library(GSE5859Subset)
data(GSE5859Subset)

?ls()
newstuff = setdiff(.newls, .oldls)
newstuff


library(GSE5859Subset)
library(Biobase)
.oldls = ls()
data(GSE5859Subset)
.newls = ls()
newstuff = setdiff(.newls, .oldls)
newstuff

cls = lapply(newstuff, function(x) class(get(x)))
names(cls) = newstuff
cls
.oldls
lapply(.oldls, function(x) class(get(x)))

j<-rbind(.oldls,lapply(.oldls, function(x) class(get(x))))
j
i<-cbind(.oldls,lapply(.oldls, function(x) class(get(x))))
class(i)
is.data.frame(i)
is.character(i)
is.matrix(i)
i<-as.data.frame(i)
i
is.matrix(i)
is.data.frame(i)
boxplot(date~factor(ethnicity), data=sampleInfo, ylab="chip date")

library(GSE5859Subset)
data()
all.equal(rownames(geneExpression), geneAnnotation$PROBEID)

ind = which(geneAnnotation$SYMBOL=="BRCA2")
boxplot(geneExpression[ind,]~sampleInfo$ethnicity, ylab="BRCA2 by hgfocus")

?ExpressionSet
es1 = ExpressionSet(geneExpression)
pData(es1) = sampleInfo
fData(es1) = geneAnnotation
es1
boxplot(es1$date~es1$ethnicity)
library(BiocInstaller)
biocLite("GSE5859")
library(devtools)
install_github("genomicsclass/GSE5859")
library(GSE5859)
data(GSE5859)
e

library(annotate)
install.packages("annotate")
library(BiocInstaller)
biocLite("annotate")

library(annotate)
p = pmid2MIAME("15269782")
p
abstract(p)
e
experimentData(e)<-p
e

experimentData(e)
