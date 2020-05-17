library(Biostrings)
library(BSgenome)
available.genomes()

grep("mask", grep("Drerio", available.genomes(), value=TRUE), invert=TRUE, value=TRUE)
grep("Drerio", available.genomes())
library(BiocInstaller)
biocLite("BSgenome.Hsapiens.UCSC.hg19.masked")

library(BSgenome.Hsapiens.UCSC.hg19.masked)
c17m = BSgenome.Hsapiens.UCSC.hg19.masked$chr17
class(c17m)
c17m
c22m<-BSgenome.Hsapiens.UCSC.hg19.masked$chr22
c22m
?masks()
22m = BSgenome.Hsapiens.UCSC.hg19.masked$chr22
round(100*sum(width(masks(c22m)$AGAPS))/length(c22m),0)
?AGAPS
masks(c22m)
library(rtracklayer)
liftchain
ch=import.chain("D:/Rdocuments/hg38ToHg19.over.chain")
ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz


library(rtracklayer)
data(targets)
targets
class(targets)
genome(targets)
head(targets)
table(targets$chrom)

ch
download.file("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
              "hg19ToHg38.over.chain.gz")
install.packages("R.utils")
library(R.utils)
gunzip("hg19ToHg38.over.chain.gz")

library(ERBS)
data(HepG2)
library(rtracklayer)
ch = import.chain("hg19ToHg38.over.chain") 
nHepG2 = liftOver(HepG2, ch)
nHepG2
HepG2
unlist(nHepG2)
start(HepG2[1])
library(ph525x)
library(Homo.sapiens)
length(transcriptsBy(Homo.sapiens, by="gene")$"2099")
modPlot("ESR1", useGeneSym=FALSE, collapse=FALSE)


library(GenomicRanges)
mtar = with(targets,
            GRanges(chrom, IRanges(start,end), strand=strand,
                    targets=target, mirname=name))

cat(export(mtar[1:5], format="bed"), sep="\n")
cat("\n")
cat(export(mtar[1:5], format="gff3"), sep="\n")


library(AnnotationHub)
ah = AnnotationHub()
ah
mah = mcols(ah)
names(mah)
sort(table(mah$species), decreasing=TRUE)[1:10]

query(ah, "HepG2")


library(Homo.sapiens)
g = genes(Homo.sapiens)
library(ERBS)
data(HepG2)
kp = g[resize(g,1) %over% HepG2]
nn = names(kp)
m = select(Homo.sapiens, keys=nn, keytype="ENTREZID",
           columns=c("SYMBOL", "GENENAME", "TERM", "GO"))
m
grep("apoptosis",m,fixed = TRUE)
install.packages("DT")
library(DT)
library(BiocInstaller)
biocLite("DT")
datatable(m)

library(devtools)
install_github("genomicsclass/maPooling")
library(maPooling)
data()
library(Biobase)
data("maPooling")


?qvalue
head(pData(maPooling))
head(maPooling)
?pData
data(maPooling)
pd=pData(maPooling)
pooled=which(rowSums(pd)==12)
pooled

individuals=which(rowSums(pd)==1)
individuals=individuals[-grep("tr",names(individuals))]

pool = exprs(maPooling)[,pooled];indiv = exprs(maPooling)[,individuals]
strain= ifelse(grepl("a",rownames(pData(maPooling))),0,1)
g_pool = strain[pooled]
g_indiv = strain[individuals]

g_pool
g_indiv
sdd<-sd(exprs(maPooling),g_pool)
library(qvalue)
install.packages("qvalue")
biocLite("qvalue")
library(qvalue)
?qvalue


library(BiocInstaller)
biocLite("BiocParallel")
library(BiocParallel)
register(SnowParam(workers=3))
system.time( bplapply(1:3, function(x) Sys.sleep(1) ) )

biocLite("GO.db")
library(GO.db)
con = GO.db$conn
head(con[1:10,])
con
dbGetQuery(con, "select count(*) from go_term")

biocLite("microbenchmark")
library(microbenchmark)
m1 = microbenchmark(
  dbGetQuery(GO.db$conn, "select term from go_term"), times=10L, unit="ms")
m2 = microbenchmark(
  keys(GO.db, keytype="TERM"), times=10L, unit="ms")

m1
m2
dbListTables(con)
library(GenomicAlignments)
?ScanBamParam()
