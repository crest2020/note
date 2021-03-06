peaks=read.table("M2_14pvaluee_4_summits.bed",header = F)
head(peaks)
library(BiocInstaller)
biocLite("TxDb.Rnorvegicus.UCSC.rn6.refGene")
library(BSgenome.Rnorvegicus.UCSC.rn6)
nrow(peaks)
library(GenomicRanges)
library(TxDb.Rnorvegicus.UCSC.rn6.refGene)
txdb=TxDb.Rnorvegicus.UCSC.rn6.refGene
head(peaks)
gr=GRanges(seqnames = peaks$V1,ranges = IRanges(start = peaks$V2,end = peaks$V3))
gr
promoter=promoters(txdb,upstream = 500,downstream = 500)
promoter=trim(promoter)
promoter
overlap=findOverlaps(gr,promoter)
overlap
nrow(peaks)
head(subjectHits(overlap))
pro_ov=promoter[subjectHits(overlap),]
library(org.Rn.eg.db)
symbol=select(org.Rn.eg.db,keys = pro_ov$tx_name,keytype = "ACCNUM",columns =c("GENENAME","SYMBOL","ENSEMBL"))

columns(org.Rn.eg.db)
write.csv(symbol,"500.csv")
