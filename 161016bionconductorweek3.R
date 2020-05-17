library(GenomicRanges)
library(BiocInstaller)
biocLite("minfiData")
library(ALL)
data(package="ALL")
data(ALL)
ALL
expall=exprs(ALL)
head(expall)
library(hgu95av2.db)
mean(expall[,5])
library(biomaRt)
head(listMarts())
mart <- useMart(host='feb2014.archive.ensembl.org', biomart = "ENSEMBL_MART_ENSEMBL")
head(listDatasets(mart),n=50)
ensemblhsapiens=useDataset("hsapiens_gene_ensembl",mart)
fdata=fData(ALL)
fdata
ALL
pdata=pData(ALL)
head(pdata)
head(featureNames(ALL))
values <- c("202763_at","209310_s_at","207500_at")
chrom_some=seq(1:22)
ensemblid=getBM(attributes = c("ensembl_gene_id", "affy_hg_u95av2"),
      filters = c("affy_hg_u95av2","chromosome_name"), values = list(featureNames(ALL),chrom_some), mart = ensemblhsapiens)

grep("affy*",aaa[,1])
aaa=listAttributes(ensemblhsapiens)
class(aaa)
?grep
grep
aaa[seq(92,109,1),]
seq(92,109,1)
seq(2:7)
?seq
ALL
head(ensemblid)
head(featureNames(ALL))
length(unique(ensemblid[,1]))
b=c(1,1,1,2,3)
sum(table(ensemblid[,2])>=1)
library(minfiData)
data(package="minfiData")
data("MsetEx")
MsetEx
data("RGsetEx")
RGsetEx
ass=assayData(MsetEx)
ass
ass[1]
methods(ass)
a=c(1,2,3)
showMethods(a)
methods(a)
methods(class = "MethylSet")
MsetEx
getClass("MethylSet")
meth=getMeth(MsetEx)

head(meth)
colnames(meth[,2])
rownames(meth[1,])
class(meth)
getClass("meth")
methods("meth")
?matrix
dimnames(meth)[[2]]
methdataframe=as.data.frame(meth)
aa=methdataframe[,2]
dim(meth)
class(aa)
mean(aa)
?mean
library(GEOquery)
elist=getGEO("GSE788")
class(elist)
length(elist)
names(elist)
edata=elist[[1]]
head(edata)
exp=exprs(edata)
head(exp)
mean(exp[,2])
library(airway)
data(airway)
airway
library(GenomicRanges)
colData(airway)
air=assay(airway,"counts")
head(air)
mean(air)
mean(apply(air,2,mean))
mean(airway$avgLength)
sum(air[,"SRR1039512"]>=1)
gr=rowRanges(airway)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb=TxDb.Hsapiens.UCSC.hg19.knownGene
subsetByOverlaps(cds(txdb), gr)
findOverlaps(transcripts(txdb),unlist(gr))
gr
gr
unlist(gr)
intronsByTranscript(txdb)
browseVignettes("GenomicRanges")
?GRangesList
trans=transcripts(txdb)
trans
subsetByOverlaps(airway,trans)

gr <- GRanges(seqnames = "1", ranges = IRanges(start = 1, end = 10^7))
subsetByOverlaps(airway, reduce(trans))


reduce(trans)
head(IRanges(gr))
gr
rowRanges(airway)[6]
gr <- GRanges(seqnames = "chr1", strand = c("+", "+", "+"),
              ranges = IRanges(start = c(1,3,5), width = 3))
gr
gr2 <- GRanges(seqnames = c("chr1", "chr2", "chr1"), strand = "*",
               ranges = IRanges(start = c(1, 30, 4), width = 1))
gr2
subsetByOverlaps(gr,gr2)
intersect(gr2,gr,ignore.strand=T)
newstyle=mapSeqlevels(seqlevels(gr),"dbSNP")
newstyle
?mapSeqlevels()
air1=rowRanges(airway)
air1
newstyle=mapSeqlevels(seqlevels(trans),"NCBI")
trans=renameSeqlevels(trans,newstyle)
newstyle
trans
seqle=paste0("chr",seq(1:22))
seqle
air22=keepSeqlevels(air1,as.character(seq(1:22)))
air22
newstyle=mapSeqlevels(seqlevels(air22),"UCSC")
air22=renameSeqlevels(air22,newstyle)
air22
trans
subsetByOverlaps(air22,cds)
intron=intronsByTranscript(txdb)
aaa=subsetByOverlaps(air22,exons(txdb))
30391-24221
cds=cds(txdb)
subsetByOverlaps(intron,trans)
trans
intron
trans=keepSeqlevels(trans,seqle)
trans
air22
intron=keepSeqlevels(intron,seqle)
intron
subsetByOverlaps(trans,intron)
air22[[1]]$exon_id
mcols(air22[1])

as.list(air22$exon_id)==0
class(mcols(airway)$exon_id)
unlist(air22)
table((unlist(air22)$exon_id)==0)
air22
median(width(airway))
air
head(air)
srr=air[,"SRR1039508"]
head(srr)
srr[1]
dim(srr1)
rownames(srr1)=rownames(air)
head(srr1)
head(names(ov))
ov
?GRangesList
nameov=names(ov)
nameov=as.data.frame(nameov)
dim(nameov)
dim(srr1)
air22
a






ense=data.frame(srr1,row=rownames(srr1))
head(ense)
?merge
mer=merge(ense,nameov,by="row")
names(nameov)
names(nameov)="row"
head(mer)
dim(mer)
sum(mer$srr)
head(ense)
dim(ense)
sum(ense$srr)
18582828/20637971
mer

library(AnnotationHub)
ah=AnnotationHub()
query(ah,c("H3K4me3","e096","narrow"))
hek4=ah[["AH30596"]]
hek4

prom=promoters(txdb)
promhek4=subsetByOverlaps(prom,hek4)
promhek4
transhek4prom=subsetByOverlaps(exons(txdb),promhek4)
transhek4prom
newair=subsetByOverlaps(air22,transhek4prom)
newairname=names(newair)
newairname=as.data.frame(newairname)
dim(newairname)
sum(duplicated(newairname))
head(newairname)
names(newairname)="row"
mernew=merge(srr12,newairname,by="row")
head(mernew)
dim(mernew)
median(mernew$srr)
?promoters
head(srr)
srr1=as.data.frame(srr)
srr12=data.frame(srr1,row=rownames(air))
which(mernew$srr==226)
srr12
head(srr12)
dim(srr12)
length(unique(srr12$row))
query(ah,c("homo","cpg"))
cpg=ah[["AH5086"]]
cpg
cpg22=keepSeqlevels(cpg,"chr22")
cpg22
library(BSgenome.Hsapiens.UCSC.hg19)
cpg22seq=Views(Hsapiens,cpg22)
cpg22seq
di=dinucleotideFrequency(cpg22seq)[,"CG"]
cc=(letterFrequency(cpg22seq,"C"))
gg=(letterFrequency(cpg22seq,"G"))
len=width(cpg22seq)
head(len)
mean(di*len/cc/gg)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb=TxDb.Hsapiens.UCSC.hg19.knownGene
promm=promoters(txdb,upstream = 900,downstream = 100)
promm
prom22=keepSeqlevels(promm,"chr22")
prom22cds=subsetByOverlaps(prom22,cds(txdb))
prom22seq=Views(Hsapiens,prom22cds)
prom22seqplus=pro22seq[strand(pro22seq)=="+",]
prom22seqminus=pro22seq[strand(pro22seq)=="-",]
table(vcountPattern("TATAAA",DNAStringSet(prom22seqplus)))
trans=transcripts(txdb)
cds=cds(txdb)
transov=subsetByOverlaps(trans,cds)
pro=promoters(transov,upstream=900,downstream=100)
pro
pro22=keepSeqlevels(pro,"chr22")
pro22
pro22seq=Views(Hsapiens,pro22)
table(vcountPattern("TTTATA",DNAStringSet(pro22seq)))
pro22ov=subsetByOverlaps(pro22,prom22)
pro22ov
pro22
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb=TxDb.Hsapiens.UCSC.hg19.knownGene
pp=promoters(genes(txdb))
pp
pp1=promoters(transcripts(txdb))
table(pp1==pp)
trans=transcripts(txdb)
cds=cds(txdb)
transov=subsetByOverlaps(trans,cds)
transov
pr=promoters(transov,upstream = 900,downstream = 100)
pr
pr22=keepSeqlevels(pr,"chr22")
pr22
int=intersect(pr22,pr22,ignore.strand=T)
sum(width(int))

geneov=subsetByOverlaps(genes(txdb),cds(txdb))
geneov
pro=promoters(geneov,upstream = 900,downstream = 100)
pro222=keepSeqlevels(pro,"chr22")
pro222
in1=intersect(pro222,pro222,ignore.strand=T)
sum(width(in1))
trans
