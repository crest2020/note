library(ERBS)
data("GM12878")
data("HepG2")
class(HepG2)
HepG2
granges(HepG2)
seqinfo(HepG2)="hg3819"
?seqinfo
values(HepG2)
HepG2[1:10,]
seqnames(HepG2)
chr=seqnames(HepG2)
as.character(chr)
table(chr)
table(chr)[1:24]
HepG2[chr=="chr20"]
x<-HepG2[order(HepG2),]
seqnames(x)
browseVignettes("IRanges")
library(GenomicRanges)
gr=GRanges(seqnames=(c("chr1","chr2")),ranges=IRanges(start=1:2,end=4:5))
gr
seqlevels(gr,force=TRUE)=c("chr1")
gr
gr2=GRanges(seqnames=(c("chr1","chr2")),ranges=IRanges(start=1:2,end=4:5))
gr2
dropSeqlevels(gr2,"chr2")
library(AnnotationHub)
ah=AnnotationHub()
ah
ah[1]
ah1=subset(ah,species=="Homo sapiens")
ah1
ah2=display(ah)
qhs=query(ah,"CpG Islands")
qhs
qhsh=query(ah,c("CpG Islands","Homo sapiens"))
qhsh
qhsh[[1]]
cpg=qhsh[[1]]
cpg
cpg[seqnames(cpg)=="chr4",]
cpg[seqnames(cpg)=="chrX",]
cpg[seqnames(cpg)=="chrY",]
28691-896-181
table(seqnames(cpg))
h3k4me3=query(ah,c("H3K4me3","H1"))
h3k4me3
write.dcf(h3k4me3,"h3k4me3")
mcols(h3k4me3)$description
less(h3k4me3)
head(h3k4me3,20)
library(rtracklayer)
h19=h3k4me3[[8]]
h19
width(ranges(h19))
args(seqnames)
?seqnames
library(BSgenome.Hsapiens.UCSC.hg19)
chr22 = Hsapiens[["chr22"]]
s = subseq(chr22,start=23456789,width=1000)
print( as.character(s) )


ah19 = subset(ah,ah$genome=="hg19")
query(ah19,"genes")
query(ah19,"CpG Islands")
cgi = ah[["AH5086"]]
class(cgi)
library(BSgenome.Hsapiens.UCSC.hg19)
cgiseq= getSeq(Hsapiens,cgi)
genome(cgi)[1:24]
genome(Hsapiens)[1:24]
cgiseq

median(letterFrequency(cgiseq,"G",as.prob = TRUE))
median(vcountPattern("GC",cgiseq)/width(cgiseq)/letterFrequency(cgiseq,"C",as.prob = TRUE)/letterFrequency(cgiseq,"G",as.prob = TRUE))
width(cgiseq)
bb=vcountPattern("CG",cgiseq)/width(cgiseq)
bb
cgiseq
head(width(cgiseq))
median(bb/(0.3418675*0.342723))

chr2use = seqlevels(cgi)[1:24]
index = which( seqnames(cgi) %in% chr2use)
noncgi = shift(cgi[index],20000)
cgi
cgi[index]
noncgiseq= getSeq(Hsapiens,noncgi)
nullres = alphabetFrequency(noncgiseq)
head(nullres)
keepIndex=nullres[,"G"]>0 &  nullres[,"C"]>0 & nullres[,"N"]==0
nullres = nullres[keepIndex,]
nullres
head(nullres)
noncgiseq=noncgiseq[keepIndex,]
median(vcountPattern("CG",noncgiseq)/width(noncgiseq)/letterFrequency(noncgiseq,"C",as.prob = TRUE)/letterFrequency(noncgiseq,"G",as.prob = TRUE))
chr22 = Hsapiens[["chr22"]]


lengths = c(100,200,300,100,100)
mat = cbind(c(1,1,0,1,0),c(1,1,1,1,1),c(0,1,1,0,1))
mat
lengths%*%mat
counts = c(125,350,300,125,100)
theta.hat = c(1, 2, 3) / 10000
w=1000
mat %*% theta.hat * lengths * w
LHS = counts/(lengths * w)
LHS
lm.fit(mat, LHS)$coefficients
counts1<-c(60,320,420,60,140)
LHS1<-counts1/(lengths*w)
LHS1
lm.fit(mat,LHS1)$coefficients
library(BiocInstaller)
biocValid()
