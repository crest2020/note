library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb=TxDb.Hsapiens.UCSC.hg19.knownGene
txdb
genes(txdb)
transcripts(txdb)
cds(txdb)
exons(txdb)
microRNAs(txdb)
tRNAs(txdb)
library(BiocInstaller)
biocLite("FDb.UCSC.tRNAs")
promoters(txdb)
transcriptsBy(txdb,by=c("gene","exon","cds"))
??transcriptsBy
transcriptsBy(txdb,)
cdsBy(txdb)
exonsBy(txdb)
intronsByTranscript(txdb)
fiveUTRsByTranscript(txdb)
threeUTRsByTranscript(txdb)
gr=GRanges(seqnames = "chr1",strand="+",ranges=IRanges(start=11874,end=14409))
gr

subsetByOverlaps(genes(txdb),gr,ignore.strand=T)
subsetByOverlaps(transcripts(txdb),gr)
subsetByOverlaps(exons(txdb),gr)
subsetByOverlaps(exonsBy(txdb,by="tx"),gr)

subsetByOverlaps(cds(txdb),gr)
subsetByOverlaps(cdsBy(txdb,by="tx"),gr)
subset(transcriptLengths(txdb,with.cds_len = T),gene_id=="100287102")
extractTranscriptSeqs(transcripts(txdb),gr)

library(BSgenome.Hsapiens.UCSC.hg19)

genome=BSgenome.Hsapiens.UCSC.hg19
cds=cdsBy(txdb,by="tx",use.names=T)
cds_seqs=extractTranscriptSeqs(genome,cds)
cds_seqs
cds_seqs["uc009vjk.2",]
cds_seqs$
colnames(cds_seqs)
gr=GRanges(seqnames = "chr1",strand=c("+","-","+"),ranges = IRanges(start = c(1,3,5),width = 3))
gr
gr
flank(gr,2,start=F)
seqinfo(gr)
gr
seqlengths(gr)=c("chr1"=10)
seqinfo(gr)
gr
gr
seqlevels(gr)=c("chr1","chr2")
seqnames(gr)=c("chr1","chr2","chr1")
gr
sort(gr)
seqlevels(gr)=c("chr2","chr1")
gr
sort(gr)
genome(gr)="hg19"
gr
gr2=gr
genome(gr2)="hg18"
gr2
findOverlaps(gr,gr2)
values(gr)=DataFrame(score=c(0.1,0.5,0.3))
gr
gr$score
gr$score2=c("1","2","3")
gr
gr2=GRanges(seqnames = c("Chr1","chr2","chr1"),strand="*",ranges = IRanges(start=c(1,3,5),width=3))
gr2
gr2
gr2
grindex=findOverlaps(gr,gr2)
queryHits(grindex)
gr[queryHits(grindex),]
subsetByOverlaps(gr,gr2)
df=data.frame(chr="chr1",start=1:3,end=4:6,score=7:9)
df
makeGRangesFromDataFrame(df,keep.extra.columns = T)
resize(gr,width=1,fix="center")
gr
resize(gr2,width=2,fix="center")
gr2
gr=GRanges(seqnames=c("chr1","chr2"),ranges=IRanges(start=1:2,end=4:5))
gr
gr$score=c(5,7)
gr
genome(gr)="hg19"
gr
seqlengths(gr)=c(100,400)
gr
seqlengths(gr)
seqlevels(gr)
keepSeqlevels(gr,"chr2")
newStyle=mapSeqlevels(seqlevels(gr),"NCBI")
newStyle
gr=renameSeqlevels(gr,newStyle)
gr
newStyle=mapSeqlevels(seqlevels(gr),"UCSC")
newStyle
gr=renameSeqlevels(gr,newStyle)
gr
seqinfo(gr)
library(AnnotationHub)
ah=AnnotationHub()


shortname=c("C1_R1","C1_R2","C1_R3","C2_R1","C2_R2","C2_R3")
condition=c(rep("control",3),rep("treated",3))
condition
setwd("htcount180424")
samplename=list.files()
library(DESeq2)
sampletable=data.frame(samplename=shortname,filename=samplename,condition=condition)
sampletable
dds=DESeqDataSetFromHTSeqCount(sampleTable = sampletable,design = ~condition)
dds
dds=dds[rowSums(counts(dds))>1,]
dds
ddsresult=DESeq(dds)
ddsresult
result=results(ddsresult)
result
table(result$pvalue<0.05)
res=result[result$pvalue<0.05,]
sum(is.na(result$pvalue))
nrow(result)
table(result$pvalue<0.05)
is.na(result$pvalue)
index=is.na(result$pvalue)
result_omitna=result[!index,]
nrow(result_omitna)
nrow(result)
result_pvalue=result_omitna[result_omitna$pvalue<0.05,]
nrow(result_pvalue)
result_pvalue
rownames(result_pvalue)
write.csv(row.names(result_pvalue),"result_genename.csv")
packages <- c("HGNChelper", "rjson")
install.packages(packages, dependencies = T)
library(rjson)
library(HGNChelper)
library(HGNChelper)
dir.create("TCGA")
setwd("D://TCGA-Assembler")
list.files()
source("Module_A.R")
source("Module_B.R")
sCancer="BRCA"
vPatientID <- c("TCGA-A7-A13F", "TCGA-AO-A12B", "TCGA-AR-A1AP", "TCGAAR-A1AQ", "TCGA-AR-A1AS", "TCGA-AR-A1AV", "TCGA-AR-A1AW", "TCGA-BHA0BZ", "TCGA-BH-A0DD", "TCGA-BH-A0DG")
vPatientID
sPath1 <- "./QuickStartExample/Part1_DownloadedData"
sPath2 <- "./QuickStartExample/Part2_BasicDataProcessingResult"
sPath3 <- "./QuickStartExample/Part3_AdvancedDataProcessingResult"
path_copynumber=DownloadCNAData(cancerType = sCancer,assayPlatform = "cna_cnv.hg19",inputPatientIDs = vPatientID,saveFolderName = sPath1)
path_copynumber
??DownloadCNAData
list.files()
path_copynumber
path_methylation_450=DownloadMethylationData(cancerType = sCancer,assayPlatform = "methylation_450",inputPatientIDs = vPatientID,saveFolderName = sPath1)
path_miRNAExp=DownloadmiRNASeqData(cancerType = sCancer,assa="mir_HiSeq.hg19.mirbase20",inputPatientIDs = vPatientID,saveFolderName = sPath1)
path_geneExp=DownloadRNASeqData(cancerType = sCancer,assayPlatform = "gene.normalized_RNAseq",inputPatientIDs = vPatientID,saveFolderName = sPath1)
path_protein_iTRAQ=DownloadCPTACData(cancerType = sCancer,assayPlatform = "proteome_iTRAQ",inputPatientIDs = vPatientID,saveFolderName = sPat
                                     
path_somaticMutation=DownloadSomaticMutationData(cancerType = sCancer,assayPlatform = "somaticMutation_DNAseq",inputPatientIDs = vPatientID,saveFolderName = sPath1)                                    
list_copyNumber=ProcessCNAData(inputFilePath = )


library(genefilter)
library(GSE5859Subset)
data(GSE5859Subset)
g=factor(sampleInfo$group)
g
sampleInfo
results=rowttests(geneExpression,g)
head(results)
GSE5859Subset
data(package="GSE5859Subset")
head(geneAnnotation)
head(geneExpression)
pvals=results$p.value
head(pvals)
head(results)

m=nrow(geneExpression)
m
n=ncol(geneExpression)
n
randomData=matrix(rnorm(m*n),m,n)
nullpvals=rowttests(randomData,g)$p.value
head(randomData)
head(geneExpression)
plot(results$dm,-log10(results$p.value))
head(results)
plot(results$dm,-log10(results$p.value),xlab="Effect size",ylab="-log(base 10) p-values")
library(rafalib)
mypar(1,2)
hist(nullpvals,ylim=c(0,1400))
hist(pvals,ylim=c(0,1400))
replicate(10,sum(rnorm(100)))
permg=sample(g)
permg
perresults=rowttests(geneExpression,permg)
hist(perresults$p.value)
randomdata1=matrix(rnorm(m*n),m,n)
g
library(genefilter)
nullpvalue=rowttests(randomdata1,g)$p.value
hist(nullpvalue,ylim=c(0,1400))
m
n
library(Biobase)
library(GSE5859)
data(GSE5859)
ge=exprs(e)
ge
dim(ge)
library(rafalib)
boxplot(ge,range=0,names=1:ncol(e))
boxplot(ge[,1],ge[,2])
boxplot(ge[,c(1,2,3,4,5,6)])
x=ge[,1]

y=ge[,2]
mypar(1,2)
plot(x,y)
plot((x+y)/2,x-y)
sd(y-x)
library(tissuesGeneExpression)
data(package="tissuesGeneExpression")
dim(e)
head(e)
head(tab)
head(tissue)
table((colnames(e)==tab$filename))
dim(e)
table(tab$Tissue)
d=dist(t(e))
d
class(d)
as.matrix(d)[1,2]
as.matrix(d)[1,87]



library(rafali
library(MASS)
n=100
y=t(mvrnorm(n,c(0,0),matrix(c(1,0.95,0.95,1),2,2)))
y
s=svd(y)
s
sqrt(2)*s$u
PC1=s$d[1]*s$v[,1]
PC2=s$d[2]*s$v[,2]
plot(PC1,PC2,xlim=c(-3,3),ylim=c(-3,3))
set.seed(1)
ind=sample(nrow(e),500)
Y=t(apply(e[ind,],1,scale))
head(Y)
head(e)
dim(Y)
yy=t(Y)
dim(yy)
d=apply(e[ind,],1,scale)
head(d)
dim(d)
d=apply(e[ind,],1,mean)
head(d)
dim(D)
dim(d)
class(d)
s=svd(Y)
U=s$u
V=s$v
D=diag(s$d)
yhat=U%*%D%*%t(V)
head(yhat)
resid=Y-yhat
max(abs(resid))
plot(s$d)
??mvrnorm

library(BiocInstaller)
biocLite("RTCGA")
biocLite("RTCGA.clinical")
biocLite("RTCGA.mutations")
library(RTCGA)

??viginette
?vignette
browseVignettes("RTCGA")
library(RTCGA.clinical)
head(checkTCGA("DataSets",cancerType="ACC"))

library(RTCGA.mRNA)
biocLite("cgdsr")
library(cgdsr)
browseVignettes("cgdsr")
head(getCancerStudies("BRCA"))
mycgds = CGDS("http://www.cbioportal.org/")
test(mycgds)
getCancerStudies(mycgds)[,c(1,2)]
getGeneticProfiles(mycgds,'brca_tcga')[,c(1:2)]
??getGeneticProfiles
getCaseLists(mycgds,'gbm_tcga')[,c(1:2)]
getProfileData(mycgds, c("p53","actb"), "brca_tcga_mrna", "brca_tcga_all")[c(1:5),]
getClinicalData(mycgds, "ova_all")[c(1:5),]
df = getProfileData(mycgds, "NF1", c("gbm_tcga_gistic","gbm_tcga_mrna"), "gbm_tcga_all")

head(df)
boxplot(df[,2]~df[,1])
stripchart(df[,2] ~ df[,1], vertical=T, add=T, method="jitter",pch=1,col='red')
df
boxplot(df[,2]~df[,1])
stripchart(df[,2]~df[,1],vertical=T,add=T,method="jitter",pch=1,col="red")
df=getProfileData(mycgds,c("MDM2","MDM4"),"gbm_tcga_mrna","gbm_tcga_all")
head(df)
plot(df)
abline(0,1)
plot(df,main="MDM2 and MDM4 gene expression",xlab="MDM2",ylab="MDM4")
plot(mycgds, "gbm_tcga", c("MDM2","MDM4"), "gbm_tcga_mrna" ,"gbm_tcga_all")
df.met=getProfileData(mycgds,"PTEN","prad_mskcc_mrna_median_Zscores","prad_mskcc_mets")
df.pri=getProfileData(mycgds,"PTEN","prad_mskcc_mrna_median_Zscores","prad_mskcc_primary")
boxplot(list(t(df.pri),t(df.met)),main="PTEN expression in primary and metastatic tumors",xlab="Tumor type",ylab="Gene expression level",names=c("primary","metastatic"),outpch=NA)
head(df.met)s
df.met
stripchart(list(t(df.pri),t(df.met)),add=T,method = "jitter",vertical = T,pch=1,col = "red")
??ttest
t.test(df.pri,df.met)
?boxplot
boxplot(list(t(df.pri),t(df.met)),main="PETN expression in primary and metastatic tumours",xlab="Tumour type",ylab="Gene expression level",names=c("primary","metastatic"))
dim(df.pri)
library(BiocInstaller)
biocLite("RTCGA.mRNA")
library(RTCGA.mRNA)
library(RTCGA)
library(RTCGA.clinical)
library(RTCGA.mutations)
library(RTCGA.rnaseq)
library()
biocLite('RTCGA.mRNA.20160128')
expr=expressionTCGA()

source("https://bioconductor.org/biocLite.R")
biocLite("RTCGA.mRNA")
remove.packages("BiocInstaller", lib=.libPaths())
source("https://bioconductor.org/biocLite.R")
biocValid()
remove.packages("BiocInstaller", lib=.libPaths())
source("https://bioconductor.org/biocLite.R")
library(BiocInstaller)
biocLite("preprocessCore")
library("affy")
library(gcrma)
library(limma)
data(package="gcrma")

?gcrma
?rma
targets=readTargets("targets.txt")
targets <- readTargets()

diabets=c("Type1","Type2","Type1","Type1")
diabets
as.numeric(factor(diabets,ordered = T))
?factor

library(vcd)
dose=c(20,30,40,45,60)
drugA=c(16,20,27,40,60)
drugB=c(15,18,25,31,40)
opar=par(no.readonly = T)
opar
par(lwd=2,cex=1.5,font.lab=2)
plot(dose,drugA,type="b",pch=15,lty=1,col="red",ylim=c(0,60),main="Drug A vs Drug B",xlab="Drug Dosage",ylab="Drug Response")
lines(dose,drugB,type="b",pch=17,lty=2,col="blue")
abline(h=c(30),lwd=1.5,lty=2,col="gray")
?abline
library(Hmisc)
minor.tick(nx=3,ny=3,tick.ratio = 0.5)
legend("topleft",inset=.01,title="Drug Type",c("A","B"),lty=c(1,2),pch=c(15,17),col=c("red","blue"))
par(opar)
help(legend)
library(vcd)
counts=table(Arthritis$Improved,Arthritis$Treatment)
counts
barplot(counts,main="Stacked Bar plot",xlab="Treatment",ylab = "Frequecncy",col=c("red","yellow","green"))
barplot(counts,main="Grouped bar plot",xlab="treatment",ylab="frequency",col=c("red","yellow","green"),beside = T)
legend(locator(1),inset=.005,fill=c("red","yellow","green"),title="tretment",c("None","Some","Marked"),cex=0.7,text.width=.5)
