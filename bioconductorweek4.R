library(BiocInstaller)
biocLite("ShortRead")
library(ArrayExpress)
library(GEOquery)
library(ShortRead)
library(yeastRNASeq)
fastqFilePath <- system.file("reads", "wt_1_f.fastq.gz", package = "yeastRNASeq")
reads=readFastq(fastqFilePath)
reads
FastqFile(fastqFilePath)
dn=DNAStringSet(sread(reads)[1:2]) 
quality(reads)[1:2]
id(reads)[1:2]
dn[5]
substr(dn,5,5)
cha=substr(sread(reads),5,5)
table(cha,as.prob=TRUE)
363841/(363841+160678+268298+7088+200095)
quality=as(quality(reads),"matrix")
mean(quality[,5])
library(leeBamViews)
bamFilePath <- system.file("bam", "isowt5_13e.bam", package="leeBamViews")
bamFilePath
bamfile=BamFile(bamFilePath)
bamfile
seqinfo(bamfile)
seqlevels(bamfile)
aln=scanBam(bamfile)
length(aln)
class(aln)
aln=aln[[1]]
aln
names(aln)
list.files(system.file("bam", package="leeBamViews"
lapply(aln, function(xx) xx[1])        
yieldSize(bamfile)=NA
open(bamfile)
scanBam(bamfile)[[1]]$seq
close(bamfile)
gr=GRanges(seqnames = "Scchr13",ranges=IRanges(start=800000,end=801000))
param=ScanBamParam(which=gr,what=scanBamWhat())
aln=scanBam(bamfile,param = param)
names(aln)
table(duplicated(aln[[1]]$pos))
browseVignettes("Rsamtools")
x=IRanges(start=c(100,200,300),width=50)
coverage(x,shift=-100)
library(Rsamtools)
scanBamWhat()
library(BiocInstaller)
biocLite("TxDb.Dmelanogaster.UCSC.dm3.refGene")
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb=TxDb.Dmelanogaster.UCSC.dm3.ensGene
txdb
exonby=exonsBy(txdb,by="tx",use.names=T)
exonby
head(seqnames(exonby))
exon_strand=unlist(runValue(strand(exonby)),use.names = F)
head(exon_strand)
tx=transcripts(txdb,columns=c("tx_name","gene_id"))
head(tx)
df=mcols(tx)
head(df)
exvytx2gene=as.character(df$gene_id)
head(exvytx2gene)
exvytx2gene=factor(exvytx2gene,levels = unique(exvytx2gene))
head(exvytx2gene)
exvytx2gene=exvytx2gene[names(exvytx2)]
head(exvytx2gene)
library(GenomicAlignments)
flag0=(isDuplicate=FALSE, isNotPassingQualityControls=FALSE)
head(flag0)
?scanBamFlag
library(minfi)
library(GEOquery)
getGEOSuppFiles("GSE68777")
head(list.files("GSE68777/"))
untar("GSE68777//GSE68777_RAW.tar",exdir = "GSE68777/idat")
untar("GSE68777/GSE68777_RAW.tar", exdir = "GSE68777/idat")
?untar
library(utils)
untar("GSE68777/GSE68777_RAW.tar", exdir = "GSE68777/idat")
getGEOSuppFiles("GSE68777")
library(leeBamViews)
library(Rsamtools)
bamfilpath=system.file("bam","isowt5_13e.bam", package="leeBamViews")
bamfile=BamFile(bamfilpath)
bamfile
seqinfo(bamfile)
gr=GRanges(seqnames = "Scchr13",ranges = IRanges(start = 807762,end=808068))
gr
params=ScanBamParam(which=gr,what=scanBamWhat())
aln=scanBam(bamfile,param = params)
names(aln)
sum(table(aln[[1]]$pos)==1)
length(unique(aln[[1]]$pos))
table(duplicated(aln[[1]]$pos))
327-198
bpaths <- list.files(system.file("bam", package="leeBamViews"), pattern = "bam$", full=TRUE)
bpaths
bamview=BamViews(bpaths)
aln=scanBam(bamview)
names(aln)
bamRanges(bamview)=gr
aln=scanBam(bamview)
aln
names(aln)
(aln[[2]]$pos)
aln[[1]]
aln
a2=aln[3]
a=a2[[1]]
length(a$pos)
a$pos
length(a$`Scchr13:807762-808068`$pos)
len=0
for (i in seq(1:8)){
  a=aln[i]
  len=len+length(a[[1]]$`Scchr13:807762-808068`$pos)
}
len/8

length(a$pos)
library(oligo)
library(GEOquery)
getGEOSuppFiles("GSE38792")
list.files("GSE38792")
head(list.files("GSE68777/idat",pattern="idat"))
idatFiles <- list.files("GSE68777/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip,overwrite=TRUE)
?gunzip
library(GEOquery)
library(minfi)
rgset=read.450k.exp("GSE68777/idat")
rgset
pData(rgset)
head(sampleNames(rgset))
geoMat <- getGEO("GSE68777")
geoMat
pdall=pData(geoMat[[1]])
head(pdall)
pd=pdall[,c("title","geo_accession","characteristics_ch1.1","characteristics_ch1.2")]
head(pd)
names(pd)[c(3,4)]=c("group","sex")
head(pd)
pd$group=sub("^diagnosis: ","",pd$group)
pd$sex=sub("^Sex: ","",pd$sex)
head(pd)
table(sampleNames(rgset)==pd$title)
sampleNames(rgset)=sub(".*_5","5",sampleNames(rgset))
rownames(pd)=pd$title
pd=pd[sampleNames(rgset),]
pData(rgset)=pd
rgset
grset=preprocessQuantile(rgset)
grset
granges(grset)
getBeta(grset)[1:3,1:3]
head(getIslandStatus(grset))
library(minfiData)
data(package="minfiData")
data(RGsetEx)
rg=RGsetEx
rg
grset=preprocessFunnorm(rg)
pdata=pData(rg)
nornames=rownames(pData(rg))[pdata$status=="normal"]
cancernames=rownames(pdata)[pdata$status=="cancer"]
cancernames
nornames
grset
head(rownames(beta)==rownames(status))
beta=getBeta(grset)
getBeta(grset)[1:3,1:3]
status=getIslandStatus(grset)
head(status)
opensea=beta[status=="OpenSea",]
dim(beta)
length(status)
dim(opensea)
meanc=apply(opensea,2,mean)
meanc
mean(meanc[nornames])-mean(meanc[cancernames])
library(AnnotationHub)     
ah=AnnotationHub()
query(ah,c("narrow","dnase","Caco2"))
caco=ah[["AH22442"]]
caco
gr=granges(grset)
gr
head(status)
grcpg=gr[status=="Island",]
grcpg
ov=findOverlaps(caco,gr)
length(unique(queryHits(ov)))
aa=c(1,2,2,3,5,5,6)
unique(aa)
123048\
library(BiocInstaller)
biocLite("zebrafishRNASeq")
