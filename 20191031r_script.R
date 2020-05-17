library(org.Rn.eg.db)
select(org.Rn.eg.db,keys="Mknk2",keytype = "SYMBOL",columns = c("ENSEMBL","GENENAME"))
fpkm=read.csv(file="fpkm.csv",header = T)
head(fpkm)

symbol=select(org.Rn.eg.db,keys=as.character(fpkm$geneID),keytype = "ENSEMBL",columns = c("GENENAME","SYMBOL"))

columns(org.Rn.eg.db)
head(fpkm)
head(symbol)
fpkm_merge=merge(fpkm,symbol,by.x="geneID",by.y="ENSEMBL",all.x=T)
head(fpkm_merge)
write.csv(fpkm_merge,"fpkm_merge.csv")









library(BiocInstaller)
biocLite("clusterProfiler")
biocLite("ChIPseeker")
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(clusterProfiler)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
files=getSampleFiles()
print(files)
peak=readPeakFile(files[[4]])
peak
covplot(peak,weightCol = "V5")
promoter=getPromoters(TxDb=txdb,upstream = 3000,downstream = 3000)
tagmatrix=getTagMatrix(peak,windows = promoter)
?plotAvgProf
plotAvg
plotAvgProf(tagmatrix, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
peak_distance=read.csv("peakdis.csv",header=T)
head(peak_distance)
d=density(peak_distance$distance)
plot(d)
head(d)
write.csv(d$x,"x.csv")
write.csv(d$y,"y.csv")
tail(d)
tail(peak_distance)
hist(d)
hist(peak_distance$distance,freq = F,ylim = c(0,0.0010),breaks = 100)
lines(density(peak_distance$distance),col="blue",lwd=2)
abline(v=0)
peak=read.table("M2_14pvalue0.001_peaks.bed",header = F)
head(peak)
library(TxDb.Rnorvegicus.UCSC.rn6.refGene)
rn6=TxDb.Rnorvegicus.UCSC.rn6.refGene
library(ChIPseeker)
library(clusterProfiler)
covplot("M2_14pvalue0.001_peaks.bed")
?covplot
promoter <- getPromoters(TxDb=rn6, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix("M2_14pvalue0.001_peaks.bed", windows=promoter)
?getTagMatrix
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
peakAnno <- annotatePeak("M2_14pvalue0.001_peaks.bed", tssRegion=c(-3000, 3000),
                         TxDb=rn6, annoDb="org.Rn.eg.db")
library(org.Rn.eg.db)
head(peakAnno)
peakAnno
plotAnnoPie(peakAnno)
peakHeatmap("M2_14pvalue0.001_peaks.bed",  TxDb=rn6, 
            upstream=3000, downstream=3000, 
            color="blue")
plotDistToTSS(peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")
library(Reac)
library(BiocInstaller)
biocLite("ReactomePA")
library(ReactomePA)
pathway1 <- enrichPathway(as.data.frame(peakAnno)$geneId,organism = "rat")
head(as.data.frame(peakAnno)$geneId)
?enrichPathway
head(pathway1, 10)
gene <- seq2gene(peak, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
?seq2gene
dotplot(pathway1)
peak1=peak=read.table("M2_14pvalue0.001_peaks.bed",header = F)
head(peak1)
gr=GRanges(seqnames = peak1$V1,ranges = IRanges(start=peak1$V2,end=peak1$V3),score=peak1$V5)
?GRanges
gr
gene <- seq2gene(gr, tssRegion = c(-100,100), flankDistance = 500, TxDb=rn6)
pathway2 <- enrichPathway(gene,organism = "rat",pvalueCutoff=0.1)
pathway2
dotplot(pathway2)
?seq2gene
?enrichPathway
fpkm=read.csv("fpkm.csv",header = T)
head(fpkm)
fpkm_merge=read.csv("fpkm_merge.csv",header = T)
head(fpkm_merge)
table(is.na(fpkm_merge$SYMBO))
fpkm_naomit=fpkm_merge[complete.cases(fpkm_merge$SYMBOL),]                                    
table(duplicated(fpkm_naomit$SYMBOL))
d=c(1,2,3,1,4,2)
duplicated(d)
d[!duplicated(d)]
fpkm_naomit_noduplicated=fpkm_naomit[!duplicated(fpkm_naomit$SYMBOL),]
table(duplicated(fpkm_naomit_noduplicated$SYMBOL))

head(fpkm_naomit_noduplicated)
write.csv(fpkm_naomit_noduplicated,"fpkm_siRNA.csv")



bc_bn=read.csv(file="over_new/bc_bn.csv",header = T)
head(bc_bn)
bc_bn_ensembl=bc_bn$ensembl
head(bc_bn_symbol)

vc_bg=read.csv(file="over_new/vc_bg.csv",header = T)
head(vc_bg)
vc_bg_ensembl=vc_bg$ensembl
head(vc_bg_symbol)
intersect_symbol=intersect(bc_bn_symbol,vc_bg_symbol)
length(intersect_symbol)
write.csv(intersect_symbol,"intersect_symbol.csv")
mg53_diff_nrf1=setdiff(vc_bg_symbol,bc_bn_symbol)
length(mg53_diff_nrf1)
write.csv(mg53_diff_nrf1,"mg53_diff_nrf1.csv")
nrf1_diff_mg53=setdiff(bc_bn_symbol,vc_bg_symbol)
length(nrf1_diff_mg53)
write.csv(nrf1_diff_mg53,"nrf1_diff_mg53.csv")


vc_bg=read.csv(file="over_new/vc_bg.csv",header = T)
head(vc_bg)
vc_bg_ensembl=vc_bg$ensembl
length(vc_bg_ensembl)
intersect_ensembl=intersect(bc_bn_symbol,vc_bg_symbol)
length(intersect_symbol)
write.csv(intersect_symbol,"intersect_symbol.csv")
mg53_diff_nrf1=setdiff(vc_bg_symbol,bc_bn_symbol)
length(mg53_diff_nrf1)
write.csv(mg53_diff_nrf1,"mg53_diff_nrf1.csv")
nrf1_diff_mg53=setdiff(bc_bn_symbol,vc_bg_symbol)
length(nrf1_diff_mg53)
write.csv(nrf1_diff_mg53,"nrf1_diff_mg53.csv")
length(bc_bn_ensembl)
length(vc_bg_ensembl)
intersect_ensembl=intersect(vc_bg_ensembl,bc_bn_ensembl)
length(intersect_ensembl)
dim(vc_bg)
dim(bc_bn)
write.csv(intersect_ensembl,"intersect.csv")


mg53_diff_nrf1_ensembl=setdiff(vc_bg_ensembl,bc_bn_ensembl)
length(mg53_diff_nrf1_ensembl)
vc_bg_ensembl[duplicated(vc_bg_ensembl)==T]

write.csv(mg53_diff_nrf1_ensembl,"mg53_diff_nrf1_ensembl.csv")

nrf1_diff_mg53_ensembl=setdiff(bc_bn_ensembl,vc_bg_ensembl)
length(nrf1_diff_mg53_ensembl)
write.csv(nrf1_diff_mg53_ensembl,"nrf1_diff_mg53_ensembl.csv")
library(org.Rn.eg.db)
select(org.Rn.eg.db,keys="Nrf1",keytype = "SYMBOL",columns = c("SYMBOL","ENSEMBL"))
columns(org.Rn.eg.db)
tbx_en

library(DESeq2)
counts=read.csv("sirna.csv",header = T)
head(counts)
shortname=c("C1","C2","C3","C6","M2","M4","M5","M6")
data.frame(shortname)
condition=factor(c(rep("Control",4),rep("MG53",4)))
condition
count_mg53=counts[,2:9]
rownames(count_mg53)=counts$geneID
head(count_mg53)
colData=data.frame(row.names=colnames(count_mg53),condition)
colData
condition
?DESeqDataSetFromMatrix
dds=DESeqDataSetFromMatrix(countData = count_mg53,colData = colData,design = ~condition)
DataFrame(condition)
dds
dds=dds[rowSums(counts(dds))>1,]
ddsresult=DESeq(dds)
result=results(ddsresult)

head(result)



result_frame=data.frame(result)
result_omitNA=result_frame[complete.cases(result_frame$padj),]
head(result_omitNA)
head(result_frame)
result_omitNA$logpvalue=-log10(result_omitNA$padj)
write.csv(result_omitNA,"MG53vsControl_siRNA.csv")
browseVignettes("DESeq2")
-log10(0.05)

###NRF1
library(DESeq2)
counts=read.csv("sirna.csv",header = T)
head(counts)

data.frame(shortname)
condition=factor(c(rep("Control",4),rep("NRF1",4)))
condition
head(counts)
count_nrf1=counts[,c(2:5,10:13)]
head(count_nrf1)
rownames(count_nrf1)=counts$geneID
head(count_nrf1)
colData=data.frame(row.names=colnames(count_nrf1),condition)
colData
condition
?DESeqDataSetFromMatrix
dds=DESeqDataSetFromMatrix(countData = count_nrf1,colData = colData,design = ~condition)
DataFrame(condition)
dds
dds=dds[rowSums(counts(dds))>1,]
ddsresult=DESeq(dds)
result=results(ddsresult)

head(result)



result_frame=data.frame(result)
result_omitNA=result_frame[complete.cases(result_frame$padj),]
head(result_omitNA)

result_omitNA$logpvalue=-log10(result_omitNA$padj)
write.csv(result_omitNA,"NRF1vsControl_siRNA.csv")
browseVignettes("DESeq2")
-log10(0.05)


###MG53 overexpression
library(DESeq2)
counts=read.table("vc_bg.count",header = F)
head(counts)


condition=factor(c(rep("Control",4),rep("MG53",4)))
condition
head(counts)
count_mg53=counts[,c(2:9)]
head(count_mg53)
rownames(count_mg53)=counts$V1
head(count_mg53)
tail(count_mg53)
colData=data.frame(row.names=colnames(count_mg53),condition)
colData
condition
?DESeqDataSetFromMatrix
dds=DESeqDataSetFromMatrix(countData = count_mg53,colData = colData,design = ~condition)

dds
dds=dds[rowSums(counts(dds))>1,]
ddsresult=DESeq(dds)
result=results(ddsresult)

head(result)



result_frame=data.frame(result)
result_omitNA=result_frame[complete.cases(result_frame$padj),]
head(result_omitNA)

result_omitNA$logpvalue=-log10(result_omitNA$padj)
head(result_omitNA)
write.csv(result_omitNA,"MG53vsControl_overexpression.csv")


####NRF1

library(DESeq2)
counts=read.table("bc_bn.count",header = F)
head(counts)


condition=factor(c(rep("Control",4),rep("NRF1",4)))
condition
head(counts)
count_nrf1=counts[,c(2:9)]
head(count_nrf1)
rownames(count_nrf1)=counts$V1
head(count_nrf1)
tail(count_nrf1)
colData=data.frame(row.names=colnames(count_nrf1),condition)
colData
condition
?DESeqDataSetFromMatrix
dds=DESeqDataSetFromMatrix(countData = count_nrf1,colData = colData,design = ~condition)

dds
dds=dds[rowSums(counts(dds))>1,]
ddsresult=DESeq(dds)
result=results(ddsresult)

head(result)

library(txdb)
library(BiocInstaller)


result_frame=data.frame(result)
result_omitNA=result_frame[complete.cases(result_frame$padj),]
head(result_omitNA)

result_omitNA$logpvalue=-log10(result_omitNA$padj)
head(result_omitNA)
write.csv(result_omitNA,"NRF1vsControl_overexpression.csv")


library(pheatmap)
heat=read.csv("heat_sirna_gly.csv",header = T)
heat
heat_sirna=heat[,3:10]
heat_sirna
rownames(heat_sirna)=heat$X.1
heat_sirna
pheatmap(heat_sirna,cluster_cols = F,cluster_rows = F,scale = "row",show_colnames =F)
?pheatmap

heat=read.csv("heat_ov_gly.csv",header = F)
heat
heat_ov=heat[,3:10]
heat_ov
rownames(heat_ov)=heat$V2
heat_ov
pheatmap(heat_ov,cluster_cols = F,cluster_rows = F,scale = "row",show_colnames =F,legend_breaks = c(-2:2))
?pheatmap
library(TxDb.Rnorvegicus.UCSC.rn6.refGene)
library(org.Rn.eg.db)
txtb=TxDb.Rnorvegicus.UCSC.rn6.refGene
org=org.Rn.eg.db
select(org,keys="ENSRNOG00000016995",keytype ="ENSEMBL",columns = c("SYMBOL","GENENAME","ENSEMBL") )



####length of exon
library(BiocInstaller)
biocLite("TxDb.Mmusculus.UCSC.mm10.knownGene")
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

library(GenomicRanges)
?GRanges
txdb=TxDb.Mmusculus.UCSC.mm10.knownGene
exon_txdb=exons(txdb)
genes_txdb=genes(txdb)
exon_txdb
genes_txdb
overlap=findOverlaps(exon_txdb,genes_txdb)
overlap
t1_raw=exon_txdb[queryHits(overlap)]
t1
t2=genes_txdb[subjectHits(overlap)]
t2

t1=as.data.frame(t1)
t1
head(t1)
t1$geneid=mcols(t2)[,1]
head(t1)
t2
length(t2)
length(t1)
length(t1_raw)
?split
head(split(t1,t1$geneid))
unlist(c(2:9))
c(2:9)
unlist(c(2:9),c(7:2))
c(7:2)
cds=cds(TxDb.Mmusculus.UCSC.mm10.knownGene)
cds
library(BiocManager)
library(BiocInstaller)
biocLite("BiocManager")
biocLite("tximportData")
biocLite("tximport")
biocLite("rhdf5")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
setwd("./counts")
list.files()

getwd()
setwd("./RNA_seq_over_001")
list.files()
bc_bn=read.csv("bc_bn.csv",header = T)
head(bc_bn)
bc_bn_fpkm=read.table("bc_bn.fpkm",header = F)
tail(bc_bn_fpkm)

colnames(bc_bn_fpkm)=c("ensembl","C1","C2","C3","C4","N1","N2","N3","N4")
head(bc_bn_fpkm)
sum(bc_bn_fpkm$N4)
?merge
colnames(bc_bn)
colnames(bc_bn_fpkm)
bc_bn_merge_fpkm=merge(bc_bn,bc_bn_fpkm,by.x = "ensembl",by.y = "ensembl",all.x = T)
head(bc_bn_merge_fpkm)
dim(bc_bn_merge_fpkm)
dim(bc_bn)
write.csv(bc_bn_merge_fpkm,"bc_bn_merge_fpkm.csv")
bc_bn_fpkm_partial=bc_bn_merge_fpkm[,9:18]
write.csv(bc_bn_fpkm_partial,"bc_bn_fpkm_partial.csv")
write.csv(bc_bn_fpkm,"bc_bn_fpkm.csv")
#####vc_bg

vc_bg=read.csv("vc_bg.csv",header = T)

head(vc_bg)
vc_bg=vc_bg[,1:10]
head(vc_bg)
vc_bg_fpkm=read.table("vc_bg.fpkm",header = F)
tail(vc_bg_fpkm)

colnames(vc_bg_fpkm)=c("ensembl","C1","C2","C3","C4","M1","M2","M3","M4")
head(vc_bg_fpkm)
sum(vc_bg_fpkm$M3)
?merge
colnames(vc_bg)
colnames(vc_bg_fpkm)
vc_bg_merge_fpkm=merge(vc_bg,vc_bg_fpkm,by.x = "ensembl",by.y = "ensembl",all.x = T)
head(vc_bg_merge_fpkm)
dim(vc_bg_merge_fpkm)
dim(vc_bg)
write.csv(vc_bg_merge_fpkm,"vc_bg_merge_fpkm.csv")
vc_bg_fpkm_partial=vc_bg_merge_fpkm[,9:18]
write.csv(vc_bg_fpkm_partial,"vc_bg_fpkm_partial.csv")
write.csv(vc_bg_fpkm,"vc_bg_fpkm.csv")


####yao

library(org.Mm.eg.db)
keytypes(org.Mm.eg.db)
keys
select(org.Mm.eg.db,keys=c("100009600"),keytype = "ENTREZID",columns = c("SYMBOL","ACCNUM"))
setwd("./yao")
list.files()
IR_24H_2.genes=read.table("IR_24H_2.genes.results.counts",header = T)
head(IR_24H_2.genes)
IR_24H_3.genes=read.table("IR_24H_3.genes.results.counts",header = T)
IR_24H_4.genes=read.table("IR_24H_4.genes.results.counts",header = T)
IR_24H_5.genes=read.table("IR_24H_5.genes.results.counts",header = T)
IR_24H=cbind(IR_24H_2.genes,IR_24H_3.genes,IR_24H_4.genes,IR_24H_5.genes)
?cbind
head(IR_24H)
table(IR_24H[,1]==IR_24H[,7])

IR_4H_3.genes=read.table("IR_4H_3.genes.results.counts",header = T)
IR_4H_4.genes=read.table("IR_4H_4.genes.results.counts",header = T)
IR_4H_5.genes=read.table("IR_4H_5.genes.results.counts",header = T)
IR_4H_7.genes=read.table("IR_4H_7.genes.results.counts",header = T)
IR_4H=cbind(IR_4H_3.genes,IR_4H_4.genes,IR_4H_5.genes,IR_4H_7.genes)
head(IR_4H)
table(IR_4H[,1]==IR_4H[,5])

SHAM1.genes=read.table("SHAM1.genes.results.counts",header = T)
SHAM3.genes=read.table("SHAM3.genes.results.counts",header = T)
SHAM4.genes=read.table("SHAM4.genes.results.counts",header = T)
SHAM5.genes=read.table("SHAM5.genes.results.counts",header = T)
SHAM=cbind(SHAM1.genes,SHAM3.genes,SHAM4.genes,SHAM5.genes)
head(SHAM)
table(SHAM[,1]==IR_24H[,7])
counts=cbind(SHAM[,c(2,4,6,8)],IR_4H[,c(2,4,6,8)],IR_24H[,c(2,4,6,8)])
head(counts)
head(IR_4H)
###count_shamvsir4h
library(DESeq2)


head(counts)
condition=factor(c(rep("SHAM",4),rep("IR_4H",4)))
condition

count_shamvsir4h=counts[,c(1:8)]

rownames(count_shamvsir4h)=SHAM[,1]
head(count_shamvsir4h)
tail(count_shamvsir4h)
colData=data.frame(row.names=colnames(count_shamvsir4h),condition)
colData
condition
?DESeqDataSetFromMatrix
dds=DESeqDataSetFromMatrix(countData =round(count_shamvsir4h) ,colData = colData,design = ~condition)

dds
dds=dds[rowSums(counts(dds))>1,]
ddsresult=DESeq(dds)
result=results(ddsresult)

head(result)
result


result_frame=data.frame(result)
result_omitNA=result_frame[complete.cases(result_frame$padj),]

write.csv(result_omitNA,"count_shamvsir4h.csv")
###count_shamvsir24h
head(counts)
condition=factor(c(rep("SHAM",4),rep("IR_24H",4)))
condition

count_shamvsir24h=counts[,c(1:4,9:12)]

rownames(count_shamvsir24h)=SHAM[,1]
head(count_shamvsir24h)
tail(count_shamvsir24h)
colData=data.frame(row.names=colnames(count_shamvsir24h),condition)
colData
condition
?DESeqDataSetFromMatrix
dds=DESeqDataSetFromMatrix(countData =round(count_shamvsir24h) ,colData = colData,design = ~condition)

dds
dds=dds[rowSums(counts(dds))>1,]
ddsresult=DESeq(dds)
result=results(ddsresult)

head(result)
result


result_frame=data.frame(result)
result_omitNA=result_frame[complete.cases(result_frame$padj),]

write.csv(result_omitNA,"count_shamvsir24h.csv")
####4-24h


head(counts)
condition=factor(c(rep("IR_4H",4),rep("IR_24H",4)))
condition

count_ir4hvsir24h=counts[,c(5:8,9:12)]

rownames(count_ir4hvsir24h)=SHAM[,1]
head(count_ir4hvsir24h)
tail(count_ir4hvsir24h)
colData=data.frame(row.names=colnames(count_ir4hvsir24h),condition)
colData
condition
?DESeqDataSetFromMatrix
dds=DESeqDataSetFromMatrix(countData =round(count_ir4hvsir24h) ,colData = colData,design = ~condition)

dds
dds=dds[rowSums(counts(dds))>1,]
ddsresult=DESeq(dds)
result=results(ddsresult)

head(result)
result


result_frame=data.frame(result)
result_omitNA=result_frame[complete.cases(result_frame$padj),]

write.csv(result_omitNA,"count_ir4hsir24h.csv")
head(counts)
table(counts[,4]<1&counts[,4]>0)


library(pheatmap)
??pheatmap
sirna_tca=read.csv("heat_new_siRNA.csv",header=F)
pheatmap(sirna_tca[,3:10],scale = "row",cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F)
ov_tca=read.csv("heat_new_ov.csv",header = F)
ov_tca
pheatmap(ov_tca[1:40,3:10],scale = "row",cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,legend_breaks=c(-2,-1,0,1,2))


setwd("./RNA_seq_over_001")
list.files()
fpkm=read.csv("vc_bg_fpkm.csv",header = T)
head(fpkm)
library(org.Rn.eg.db)
symbol=select(org.Rn.eg.db,keys=as.character(fpkm$ensembl),keytype = "ENSEMBL",columns = c("ENSEMBL","SYMBOL","GENENAME"))
head(symbol)
fpkm_merge=merge(fpkm,symbol,by.x="ensembl",by.y="ENSEMBL",all.x=T)
?merge
head(fpkm)
head(fpkm_merge)
dim(fpkm_merge)
fpkm_merge_omitna=fpkm_merge[complete.cases(fpkm_merge$SYMBOL),]
head(fpkm_merge_omitna)
dim(fpkm_merge_omitna)
table(unique(fpkm_merge_omitna$ensembl))
table(duplicated(fpkm_merge_omitna$ensembl))
fpkm_merge_omitna
?order
fpkm_merge_omitna_order=fpkm_merge_omitna[order(-fpkm_merge_omitna$M1),]
head(fpkm_merge_omitna_order)
fpkm_merge_omitna_order_unique=fpkm_merge_omitna_order[!duplicated(fpkm_merge_omitna_order$ensembl),]
table(duplicated(fpkm_merge_omitna_order_unique$ensembl))
?sample
z=c(1,2,3,2,4,5,3,5)
table(duplicated(z))
y=matrix(1:20,5,4)
y
sample(y,length(y))
x=1:10
sample(x,10)
dim(fpkm_merge_omitna_order_unique)
x=sample(1:20162,20162)
head(x,n=100)
length(x)
table(duplicated(x))
fpkm_merge_omitna_order_unique_random=fpkm_merge_omitna_order_unique[x,]
write.csv(fpkm_merge_omitna_order_unique_random,"fpkm_over_gsea.csv")

getwd()
heat=read.csv("heat_tgf.csv",header = T)
heat
install.packages("pheatmap")
library(pheatmap)
?pheatmap
rownames(heat)=heat[,1]
pheatmap(heat[,2:7],scale = "row",cluster_row = F,cluster_cols = F,show_colnames = F)

library(biomaRt)
listMarts()
ensembl=useMart("ensembl")
datasets=listDatasets(ensembl)
head(datasets)


install.packages('Seurat')
library(Seurat)
?update.packages
remove.packages("Seurat")
install.packages("cowplot")
detach("Matrix")
library(BiocInstaller)
biocLite("cowplot")
