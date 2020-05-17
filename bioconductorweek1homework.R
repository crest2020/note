library(AnnotationHub)
ah=AnnotationHub()
ah1=subset(ah,species=="Homo sapiens")
qhsnarr=query(ah1,c("H3K4me3","EpigenomeRoadMap","H1","narrowPeak"))
qhsnarr1=qhsnarr[[1]]
qhsnarr1trim=keepSeqlevels(qhsnarr1,seqlevels(qhsnarr1)[1:22])
qhsnarr1trim
sum(width(reduce(qhsnarr1trim)))

qhsnarr27=query(ah1,c("H3K27me3","EpigenomeRoadMap","H1","narrowPeak"))
qhsnarr27
qhsnarr271=qhsnarr27[[1]]
qhsnarr271trim=keepSeqlevels(qhsnarr271,seqlevels(qhsnarr271)[1:22])
qhsnarr271trim
mean(mcols(qhsnarr271trim)$signalValue)



qhs[2]
sort(table(qhs$description), decreasing=TRUE)
library(rtracklayer)
qhs[[1]]
qhs[1]
qhs[2]
qhs[3]
ov=findOverlaps(qhsnarr1trim,qhsnarr271trim)
ov
sum(width(intersect(qhsnarr1trim,qhsnarr271trim,ignore.strand=TRUE)))
length(unique(queryHits(ov)))
length(unique(subjectHits(ov)))
qhsnarr1trim[queryHits(ov)]

overlap=intersect(qhsnarr1trim,qhsnarr271trim)
intersect(overlap,cpg1trim)

ocp=findOverlaps(overlap,cpg1)
cpg1trim=keepSeqlevels(cpg1,seqlevels(cpg1)[1:22])

cpg1trim

length(unique(queryHits(ocp)))
a=GRanges(seqnames=c("chr1"),ranges=IRanges(start=1:3,end=4:6))
a
b=GRanges(seqnames=c("chr1"),ranges=IRanges(start=2:4,end=6:8))
b
c=findOverlaps(a,b,ignore.strand=TRUE)
c

sum(width(intersect(overlap,cpg1trimre)))

cpg1trimre=flank(cpg1trim,10000,TRUE,TRUE)
cpg1trimre

cpg1trimre=union(flank(cpg1trim,10000,TRUE,FALSE),flank(cpg1trim,10000,FALSE,FALSE))
cpg1trimre=union(cpg1trimre,cpg1trim)
cpg=query(ah1,c("CpG Islands"))
cpg1trimre

cpg[1]
cpg$genome
cpg1=cpg[[1]]
cpg1
sum(width(intersect(overlap,cpg1)))
