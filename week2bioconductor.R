library(GenomicRanges)
library(rtracklayer)
library(AnnotationHub)
ah=AnnotationHub()
table(ah$dataprovider)
ah[1]
table(table(ah$species))
unique(ah$dataprovider)
unique(ah$genome)
ahhsapiens=subset(ah,species=="Homo sapiens")
ahhsapiens
query(ah,"h3k4me3")
qhs=query(ah,"h3k4me3")
mcols(ah[2])
ahs=query(ah,c("nrf1"))
qhs=query(ahhsapiens,c('h3k4me3',"gm12878"))
?query
help(package="AnnotationHub")
ahs$title
ahs[11]
qhs
qhs$title
gr1=qhs[4]
gr1=gr1[[1]]
gr2=qhs[7]
gr2=gr2[[1]]
gr2
table(width(gr1))
refseqhg19=query(ahhsapiens,"refseq")
refseqhg19$genome
genome(gr1)
refseqhg19=refseqhg19[[1]]
refseqhg19
table(table(refseqhg19$name))
promoterhg19=promoters(refseqhg19)
table(width(promoterhg19))
args(promoters)
ov=findOverlaps(promoterhg19,gr1)
ov
length(unique(subjectHits(ov)))/length(promoterhg19)
overlap=unique(subjectHits(ov))
overlap
refseqhg19
refseqhg19[1,]$name
index=unique(subjectHits(ov))
overlapgene=refseqhg19[index,]$name
max(index)
library(Biostrings)
library(rtracklayer)
dna1=DNAString("ACGT-N")
dna1
dna2=DNAStringSet(c("ACGT","GTCA","GCTA"))
dna2
IUPAC_CODE_MAP
dna1[2:4]
dna1
dna2[2:3]
dna2
dna2[[2]]
names(dna2)=paste0("seq",1:3)
dna2
width(dna2)
sort(dna2)
rev(dna2)
rev(dna1)
translate(dna2)
reverseComplement(dna1)
dna1
alphabetFrequency(dna1)
alphabetFrequency(dna2)
letterFrequency(dna2,"GC")
consensusMatrix(dna2,as.prob = T)
library(BSgenome)
library(BiocInstaller)
biocLite("BSgenome.Scerevisiae.UCSC.sacCer2")
available.genomes()
library(BSgenome.Scerevisiae.UCSC.sacCer2)
Scerevisiae
seqlengths(Scerevisiae)
seqnames(Scerevisiae)
Scerevisiae$chrI
letterFrequency(Scerevisiae$chrI,"CG",as.prob = T)
param=new("BSParams",X=Scerevisiae,FUN=letterFrequency)
?new
head(bsapply(param,letters="GC"))
param=new("BSParams",X=Scerevisiae,FUN=letterFrequency,simplify=T)
bsapply(param,letters="GC")
sum(bsapply(param,letters="GC"))/sum(seqlengths(Scerevisiae))
dnaseq=DNAString("ACGTACGT")
matchPattern(dnaseq,Scerevisiae$chrI)
countPattern(dnaseq,Scerevisiae$chrI)
vmatchPattern(dnaseq,Scerevisiae)
dnaseq
dnaseq==reverseComplement(dnaseq)
library(BSgenome)
library(rtracklayer)
library(BSgenome.Scerevisiae.UCSC.sacCer2)
available.genomes()
library(AnnotationHub)
dnaseq=DNAString("ACGTACGT")
vi=matchPattern(dnaseq,Scerevisiae$chrI)
vi
ranges(vi)
Scerevisiae$chrI[start(vi):end(vi)]
alphabetFrequency(vi)
letterFrequency(vi,"GC")
dinucleotideFrequency(vi)
args(dinucleotideFrequency)
shift(vi,10)
gr=vmatchPattern(dnaseq,Scerevisiae)
vi2=Views(Scerevisiae,gr)
vi2
gr
qh=query(ah,c("saccer2","genes"))
qh
qhrat=query(ah,c("rattus","gene"))
qhrat$title
genes=qh[[]]
which(qh$title=="SGD Genes")
genens=qh[[1]]
genens
qh
qhrat
genesrat=qhrat[[2]]
qhrat$title
genesratrgd=qhrat[[which(qhrat$title=="RGD Genes")]]
genesratrgd
prom=promoters(genens)
head(prom,n=3)
prom=trim(prom)
promviews=Views(Scerevisiae,prom)
gcprom=letterFrequency(promviews,"GC",as.prob = T)
head(gcprom)
params=new("BSParams",X=Scerevisiae,FUN=letterFrequency,simplify=T)
gccontent=bsapply(params,letters="GC")
gcpercentage=sum(gccontent)/sum(seqlengths(Scerevisiae))
gcpercentage
plot(density(gcprom))
abline(v=gcpercentage,col="red")
sessionInfo()
library(GenomicRanges)
library(AnnotationHub)
library(rtracklayer)
r1=Rle(c(1,1,1,1,2,2,3,3,3,2,2))
r1
runLength(r1)
as.numeric(r1)
runValue(r1)
ir=IRanges(start = c(2,6),width=2)
aggregate(r1,ir,FUN=mean)
ir=IRanges(start=1:10,width=3)
r11=coverage(ir)
r11
ir
?coverage
slice(r11,2)
vi=Views(r11,start = c(3,7),width=3)
vi
mean(vi)
?aggregate
r1
as.numeric(r11)
ir
r1
r11
library(AnnotationHub)
mean(vi)
vi
gr=GRanges(seqnames = "chr1",ranges = IRanges(start = 1:10,width=3))
gr
rl=coverage(gr)
slice(r1,2)
as.numeric(rl)
grview=GRanges("chr1",ranges = IRanges(start = 2,end=7))
vi=Views(rl,as(grview,"RangesList"))
vi
vi[[1]]
sum(width(vi))
r1
grr=GRanges("chr1",ranges = IRanges(start=2,end=100))
vi=Views(rl,as(grr,"RangesList"))
vi
vi[[1]]
rl
sum(vi)
widths(vi)
width(vi)
length(vi)
args(width)
?width
library(rtracklayer)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb=TxDb.Hsapiens.UCSC.hg19.knownGene
txdb
genes(txdb)
transcripts(txdb)
cds(txdb)
transcriptsBy(txdb,by=c("gene","exon","cds"))
head(transcriptLengths(txdb))
exons(txdb)
library(BiocInstaller)
biocLite("mirbase.db")
library(mirbase.db)
microRNAs(txdb)
promoters(txdb)[1]
 args(promoters)
12073-9874
intronsByTranscript(txdb)
five=fiveUTRsByTranscript(txdb)
extractTranscriptSeqs(txdb,five[1])
?extractTranscriptSeqs
five[1]
gr=GRanges(seqnames="chr1",strand="+",ranges=IRanges(start = 11874,end=14409))
subsetByOverlaps(exonsBy(txdb,by="tx"),gr)
subsetByOverlaps(cdsBy(txdb,by="tx"),gr)
subset(transcriptLengths(txdb,with.cds_len = T),gene_id=="100287102")
library(Rsamtools)
library(AnnotationHub)
table(ah$rdataclass)
library(AnnotationHub)
library(rtracklayer)
ahub.gr=subset(ah,rdataclass=="GRanges"&species=="Homo sapiens")
gr=ahub.gr[[1]]
ahub.gr[1]
seqinfo(gr)
ahub.bw=subset(ah,rdataclass=="BigWigFile"&species=="Homo sapiens")
head(mcols(ahub.bw))
bw=ahub.bw[[1]]
library(Biostrings)
dna=DNAString("CAGCGCCTGCGCCTCATACTCCACCCTGCGCATGCGTAGC")
reverseComplement(dna)
dnaf=DNAString("AAACTCATACTCCACCCTGAAAATGCGTAGC")
dnar=DNAString("TTTCAGGGTGGAGTATGAGTTTCAGGCGCTG")
length(dnaf)
length(dnar)
matchPattern(dna,dnaf)
diff(dnaf,dna)
match(dnaf,dna)
pairwiseAlignment(dnaf,dna)
pairwiseAlignment(dnar,reverseComplement(dna))
library(ERBS)
library(devtools)
?install_github
data(package="ERBS")
hepg2=data(HepG2)
hepg2
HepG2
data(GM12878)
GM12878
seqnames(GM12878)
library(GenomicRanges)
seqnames(GM12878)
table((seqlengths(GM12878))==(seqlengths(HepG2)))
head(seqlengths(HepG2))
library(qqmath)
library(BiocInstaller)
biocLite("qqmath")
library(devtools)
install_github("qqmath")
install.packages("qqmath")
library(BSgenome.Hsapiens.UCSC.hg19)
letterFrequency(Hsapiens$chr22,"GC")/(length(Hsapiens$chr22)-letterFrequency(Hsapiens$chr22,"N"))
length(Hsapiens$chr22)
ah
ahhsapiens
query(ahhsapiens,c("h3k27me3","H1","wig","e003"))
h3k=h3k27me3[[1]]
h3k27me3
h3k
h3kChr22=h3k[seqnames(h3k)=="chr22",]
h3kChr22GC=Views(Hsapiens,h3kChr22)
h3kChr22GCcontent=letterFrequency(h3kChr22GC,"GC",as.prob = T)
letterFrequency(h3kChr22GC[1:3],"GC",as.prob = T)
mean(h3kChr22GCcontent)
seqlengths(Hsapiens)
seqnames(Hsapiens)
head(h3kChr22$signalValue)
cor(h3kChr22GCcontent,h3kChr22$signalValue)
head(mcols(h3k27me3)$tags)
h3k
mcols(query(ahhsapiens,c("h3k27me3","H1","wig","e003")))
foldchange
foldchange[1]
library(GenomicRanges)
library(rtracklayer)
foldchange$resource
fold=import(foldchange,format = "BigWig",as="GRanges")
?import
 import.wig(foldchange)
import.bw(query(ahhsapiens,c("h3k27me3","H1","wig","e003"))[[1]])
import
traceback()
import(BigWigFile("https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/foldChange/E003-H3K27me3.fc.signal.bigwig"),as="GRanges")
library(GenomicRanges)
library(rtracklayer)
import.bw(foldchange)
out.gr=import(foldchange,which=h3k[1:3])
h3k[1:3]
gr.chr22 <- GRanges(seqnames = "chr22",
ranges = IRanges(start = 1, end = seqlengths(gr)["chr22"]))
out.chr222 <- import(foldchange, which = gr.chr22)
out.chr222
cor(length(out.chr222)$score,length(h3k$signalValue)
reduce(out.chr222)


gr22=h3k[seqnames(h3k)=="chr22",]
gr22
out.chr222
chr22=import(foldchange,which=gr22)
chr22
gr22
aggregate()
rlechr22=out.chr22[["chr22"]]
rlechr222
head(runValue(rlechr22))
runLength(rlechr22[5])
head(which(runValue(rlechr22)>=1))
runLength(rlechr22[1])
runValue(rlechr22[2])
runLength(rlechr22)[]
runValue(rlechr22)

sum(runLength(rlechr22)[which(runValue(rlechr22)>=1)])
h3k$chr22
h3k[seqnames(h3k)=="chr22",]
h3kview=Views(h3k,foldchange)
?Views
Scerevisiae
library(BSgenome.Scerevisiae.UCSC.sacCer2)
Scerevisiae
Views(Scerevisiae,"ATGCGCGc")
bb=import(foldchange,which=h3k)
bb
aggregate(bb,h3k,FUN=mean)
h3k
seqnames(h3k)
h3k
table(seqnames(h3k))
bbb=Views(bb,h3k)
?Views
col(runValue(bb))
bb
h3k
Views(foldchange,h3k)
Scerevisiae$chrI
out.chr222
browseVignettes("rtracklayer")
chr22rle=import(foldchange,which=gr22,as="Rle")
chr22rle$chr22
chr22
?Views
Views(Rle(foldchange),h3k)
coverage(chr22)
chr22$chr22
gr22
chr22rle$chr22
findOverlaps(chr22,gr22)
library(Biostrings)
a=import(foldchange,as="Rle")
a=Views()
?reduce
findOverlaps(bb,h3k)
ov=findOverlaps(h3k,bb)
bb[2155670]
h3k[1]
ov
queryHits(ov)
queryLength(ov)
ov[queryHits(ov)==1]
h3k
mean(h3k[1:2181]$signalValue)
bb[1]
table(queryHits(ov)==1)
ov
duplicated(queryHits(ov))
table(queryHits(ov))
rle=Rle(queryHits(ov))
rle
(runLength(rle))[1]
index=runLength(rle)
index[2]
head(index)

a=vector(mode = "numeric",length=0)
n=0
for (i in 1:123497) {
  if (i==1){
    a[i]=mean(bb[1:index[i]]$score)
  }
  else {
    n=index[i-1]+n
    a[i]=mean(bb[(n+1):(n+index[i])]$score)}
  }
a
index[123497]
a
n

h3k[1]
h3k[2]
a=vector(mode = "numeric",length=0)
i=123497
bb[(index[i-1]+1):(index[i-1]+index[i])]
ov  
   
cor(h3k[1:12300]$signalValue,a)   
mean(bb[1:])
a
rle
bb
h3k
ov
h3k[1]
bb[2155670]
mean(bb[1:2181]$score)
bb[10790607]
h3k[1]
h3k[2]
h3k[3]
bb[10790607]
h3k[123497]
ov
bb[2155670]
subjectHits(ov)[1]

a=vector(mode = "numeric",length=0)
n=0
for (i in 1) {
  
    a[i]=mean(bb[subjectHits(ov)[i]:(subjectHits(ov)[i]+index[i])]$score)
  
}
mean(c(1,2,4))
for (i in 1:200) {
   print(mean(bb[2155669+i]$score)) }
print a
ov
bb[2155670]
cor(a,h3k[1:1230]$signalValue)
index[1]
ov[2181]
ov[2182]
ov[1]
2157850-2155670
subjectHits(ov)[1]
subjectHits(ov)[1]+index[1]
index[1]
vv=c(1,2,3)
print(vv)
dataa=c(1,1,1,2,3,3,3,3,3)
datab=c(7,8,9,10,23,45,67,88,90)
data1=data.frame(dataa,datab)
data1

datas=split(data1,data1$dataa)
aa=datas[[1]]$datab

class(aa)
bb=c(7,89,0)
class(bb)
mean(datas[1]$datab)
ov
dataindex=data.frame(queryHits(ov),subjectHits(ov))
head(dataindex)
library(GenomicRanges)
library(rtracklayer)
colnames(dataindex)=c("query","subject")
datasplit=split(dataindex,dataindex$query)
bb[c(2155670,2155671)]
dataindex

gr.chr22 <- GRanges(seqnames = "chr22",
                    ranges = IRanges(start = 1, end = seqlengths(gr)["chr22"]))

query(ah,c("E003","narrowpeak","h3k27me3"))
foldchange=ah[["AH32033"]]
bb=import(foldchange,which=h3k22)

bbrle=import(foldchange,which=grchr22,as="Rle")
sum(runLength(bbrle$chr22)[runValue(bbrle$chr22)>=1])

ahub.gr <- subset(ah, rdataclass == "GRanges" & species == "Homo sapiens")
gr <- ahub.gr[[1]]
grchr22=GRanges(seqnames="chr22",ranges=IRanges(start=1,end=seqlengths(gr)["chr22"]))
h3k=ah[["AH29892"]]
h3k22=h3k[seqnames(h3k)=="chr22",]
h3k22
bb
ov=findOverlaps(h3k22,bb)
ov
data1=(as.data.frame(ov))
datas=split(data1,data1$queryHits)
length(datas)
table(datas)
str(datas)
index=vector(mode="numeric")
a=vector(mode="numeric")
for (i in 1:4068){
  index=(datas[[i]]$subjectHits)
  a[i]=mean(bb[index]$score)}
head(a)
cor(h3k22$signalValue,a)
ov
h3k9=h3k[seqnames(h3k)=="chr9",]
bb9=import(foldchange,which=h3k9)
ov9=findOverlaps(h3k9,bb9)
data2=as.data.frame(ov9)
datas9=split(data2,data2$queryHits)
ov9
index=vector(mode="numeric")
a=vector(mode="numeric")
for (i in 1:5497){
  index=(datas9[[i]]$subjectHits)
  a[i]=mean(bb9[index]$score)}

head(a)
cor(h3k9$signalValue,a)
query(ah,c("E055","h3k27me3","fc"))
e055=ah[["AH32470"]]
library(Rsamtools)
library(GenomicRanges)
library(rtracklayer)
e055chr22=import(e055,which=grchr22)
e055chr22
e055chr22G2=e055chr22[e055chr22$score>=2,]
e055chr22G2
e003chr22=import(foldchange,which=grchr22)
e003chr22L05=e003chr22[e003chr22$score<=0.5,]
e003chr22L05
ovv=findOverlaps(e003chr22L05,e055chr22G2)
sum(width(e003chr22L05[queryHits(ovv),]))
sum(width(e055chr22G2[subjectHits(ovv),]))
sum(width(intersect(e003chr22L05,e055chr22G2)))
letterFrequency(Scerevisiae,"CG",as.prob = T)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)
library(BSgenome.Scerevisiae.UCSC.sacCer2)
params <- new("BSParams", X = Hsapiens, FUN = letterFrequency, simplify = TRUE)
gccontent <- bsapply(params, letters = "CG")
gccontent
dinucleotideFrequency(Hsapiens$chr22,step=3,as.prob = T)["CG"]
letterFrequency(Hsapiens$chr22,"C")
seqnames(Hsapiens)
?dinucleotideFrequency
0.01656698/0.16326/0.1631285
16745219/51304566
length(Hsapiens$chr22)

8369235*8375984/51304566
578097/1366361*2
ta=DNAString("TATAAA")
countPattern(at,Hsapiens$chr22)
at=reverseComplement(ta)
at
13636+13627
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb=TxDb.Hsapiens.UCSC.hg19.knownGene
txdb
gr=promoters(genes(txdb),upstream = 900,downstream = 100)
gr22=gr[seqnames(gr)=="chr22",]
vi22=Views(Hsapiens,gr22)
DNAStringSet(vi22)[1:3]
library(GenomicRanges)
library(rtracklayer)
library(Biostrings)
table(vcountPattern(at,DNAStringSet(vi22)))
genes(txdb)
seqlengths(Hsapiens)["chr22"]
dinucleotideFrequency(Hsapiens$chr22)
(1985208+578097)*51304566
8375984*8369235/51304566

letterFrequency(Hsapiens$chr22,"G")
ratvier=import("rat6chr1.wig")
ratvier
grr6=GRanges(seqnames="chr1",IRanges(start = 5000,end=5050))
Views(ratvier,grr6)
Views(ratvier)
export.bed(ratvier)
?Views
rat=as(ratvier,"XIntegerViews")
rat
?intersect
class(Hsapiens)
Hsapiens
library(AnnotationHub)
ah=AnnotationHub()
query(ah,c("homo sapiens","cpg"))$genome
cpg=ah[["AH5086"]]
cpg
cpgchr22=cpg[seqnames(cpg)=="chr22",]
cpgchr22
cpgchr22seq=Views(Hsapiens,cpgchr22)
cpgchr22seq
dnachr22cpg=DNAStringSet(cpgchr22seq)
sum(dinucleotideFrequency(dnachr22cpg)[,"CG"])
sum(seqlengths(dnachr22cpg))
562123*58265/198100/197345
cpgchr22seq
mcols(cpgchr22seq)$GC=
dinucleotideFrequency(cpgchr22seq)[,"CG"]/100*width(cpgchr22seq)/letterFrequency(cpgchr22seq,"C")/letterFrequency(cpgchr22seq,"G")
mean(mcols(cpgchr22seq)$GC)
cpgchr22seq
width(dnachr22cpg)
cpgchr22seq
sum(letterFrequency(dnachr22cpg,"G"))
txdb=TxDb.Hsapiens.UCSC.hg19.knownGene
transcribe=transcripts(txdb)
cds=cds(txdb)
cdschr22=cds[seqnames(cds)=="chr22",]
transcribechr22=transcribe[seqnames(transcribe)=="chr22",]
ov=findOverlaps(transcribechr22,cdschr22)
transcribechr22clean=transcribechr22[queryHits(ov)]
prom=promoters(transcribechr22clean,upstream = 900,downstream = 100)
prom
dnaprom=Views(Hsapiens,prom)
dnapromset=DNAStringSet(dnaprom)
table(vcountPattern("TATAAA",dnapromset))

tx_id=mcols(dnaprom)$tx_id
tx_id=unique(tx_id)
dnapromclean=dnaprom[datafclean$index]
dnapromclean
a=c(1,2,3,3,3,5,6,7)
seq(1:3)
unique(a)
a[!duplicated(a)]
length(mcols(dnaprom)$tx_id)
dataf=data.frame(index=seq(1:18356),tx=mcols(dnaprom)$tx_id)
head(dataf)
datafclean=dataf[!duplicated(dataf$tx),]
dim(datafclean)
dnapromcleanset=DNAStringSet(dnapromclean)
sum(vcountPattern("TTTATA",dnapromcleanset))
reverseComplement(DNAString("TATAAA"))

gene=genes(txdb)
gene
table(strand(gene)=="+")
dnapromcleanplus=dnapromclean[strand(dnapromclean)=="+",]
table(vcountPattern("TATAAA",DNAStringSet(dnapromcleanplus)))
dnapromcleanminus=dnapromclean[strand(dnapromclean)=="-",]
table(vcountPattern("AAATAT",DNAStringSet(dnapromcleanminus)))
139+82
inter=subsetByOverlaps(transcribechr22,cdschr22)
inter
indexa=seq(1,1521,by=2)
indexb=seq(2,1520,by=2)

intera=inter[indexa]
interb=inter[indexb]
sum(width(intersect(intera,interb,ignore.strand=TRUE)))

interignore=subsetByOverlaps(transcribechr22,cdschr22,ig)
transcribechr22
promnew=promoters(inter,upstream = 900,downstream = 100)
promnew=promnew[seqnames(promnew)=="chr22",]
promnewnew=promoters(genes(txdb),upstream=900,downstream=100)
promnew
dnanew=Views(Hsapiens,promnew)
dnanew
dnaplus=dnanew[strand(dnanew)=="+",]
dnaplus=dnanew[strand(dnanew)=="-",]
table(vcountPattern("TATAAA",DNAStringSet(dnaplus)))
table(vcountPattern("AAATAT",complement(DNAStringSet(dnaplus))))
table(vcountPattern("TATAAA",DNAStringSet(dnanew)))
dna1=DNAString("AGGCTTTGGCAAAA")
complement(dna1)
head(DNAStringSet(dnaplus))
prom
promnewnew
