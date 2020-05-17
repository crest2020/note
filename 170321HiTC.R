library(harbChIP)
data("harbChIP")
harbChIP
es=exprs(harbChIP)
head(es)
dim(es)
abstract(harbChIP)
sampleNames(harbChIP)
mind=which(sampleNames(harbChIP)=="MBP1")
qqnorm(exprs(harbChIP)[,mind],main = "MBP1 binding")
qqline(0,1)
abline(0,1)
topb=featureNames(harbChIP)[order(exprs(harbChIP)[,mind],decreasing = T)[1:10]]
topb
library(org.Sc.sgd.db)
select(org.Sc.sgd.db,keys=topb,keytype = "ORF",columns = "COMMON")
harbChIP
library(yeastCC)
data("spYCCES")

alp=spYCCES[,spYCCES$syncmeth=="alpha"]
plot(exprs(alp)[topb[1],]~alp$time,lty=1,type="l",ylim=c(-1.5,1.5),lwd=2,ylab="Expression",xlab="Minutes elapsed")
for(i in 2:5) lines(exprs(alp)[topb[i],]~alp$time,lty=i,lwd=2)
legend(75,-0.5,lty=1:10,legend=topb[1:5],lwd=2,cex=.6,seg.len = 4)
library(gwascat)
data("gwrngs19")
gwrngs19[100]
mcols(gwrngs19[100])[,c(2,7,8,9,10,11)]
library(ERBS)
data(GM12878)
fo=findOverlaps(GM12878,gwrngs19)
fo
sort(table(gwrngs19$Disease.Trait[subjectHits(fo)]),decreasing = T)[1:5]
fo=findOverlaps(GM12878,reduce(gwrngs19))
fo
ovrng=reduce(gwrngs19)[subjectHits(fo)]
phset=lapply(ovrng, function(x) unique(gwrngs19[which(gwrngs19 %over% x)]$Disease.Trait))
ovrng
sort(table(unlist(phset)),decreasing = T)[1:5]
?match
a=c(1,2,3,4,5,3)
b=c(2,3,4,5,6)
a %in% b
match(a,b)
library(BiocInstaller)
biocLite("CAGEr")
library(CAGEr)
data("FANTOM5humanSamples")
head(FANTOM5humanSamples)
dim(FANTOM5humanSamples)
length(grep("heart",FANTOM5humanSamples$description))
FANTOM5humanSamples[890,]
astrocyteSamples=FANTOM5humanSamples[grep("Astrocyte",FANTOM5humanSamples$description),]
astrocyteSamples
data("FANTOM5mouseSamples")
head(FANTOM5mouseSamples)
nrow(FANTOM5mouseSamples)
astrocyteSamples$sample
astrocyteCAGEset=importPublicData(source="FANTOM5",dataset="human",sample=astrocyteSamples[1:3,"sample"])
astrocyteCAGEset
grep("heart",FANTOM5mouseSamples$description)
FANTOM5mouseSamples[180,]
getCTSS(astrocyteCAGEset)
corr.m=plotCorrelation(astrocyteCAGEset,samples = "all",method="pearson")
librarySizes(astrocyteCAGEset)
plotReverseCumulatives(astrocyteCAGEset,fitInRange = c(5,1000),onePlot = T)
?

library(BiocInstaller)
biocLite("HiTC")
biocLite("HiCDataHumanIMR90")
library(HiCDataHumanIMR90)
data(package="HiCDataHumanIMR90")
data(Dixon2012_IMR90)
show(hic_imr90_40)
class(intdata(hic_imr90_40$chr1chr1))
object.size(hic_imr90_40)
isComplete(hic_imr90_40)
isPairwise(hic_imr90_40)
seqlevels(hic_imr90_40)
detail(hic_imr90_40$chr6chr6)
head(summary(hic_imr90_40))
sset=reduce(hic_imr90_40,chr=c("chr5","chr6","chr7"))
sset
imr90_500=HTClist(mclapply(sset,binningC,binsize=500000,bin.adjust=F,method="sum",step=1))
mapC(imr90_500)
?HTClist
mapC(forcePairwise(imr90_500),maxrange=200)
biocLite("BSgenome.Hsapiens.UCSC.hg18")
resfrag=getRestrictionFragmentsPerChromosome(resSite="AAGCTT",overhangs5 = 1,chromosomes = "chr6",genomePack = "BSgenome.Hsapiens.UCSC.hg18")
resfrag
require(rtracklayer)
map_hg18 <- import("wgEncodeCrgMapabilityAlign100mer_chr6.bw",format="BigWig")
library(AnnotationHub)
ah=AnnotationHub()
?ah
?AnnotationHub
query(ah,c("h3k4","homo","wgEncode","wig"))$title
map_hg18 <- NULL
cutsites=getAnnotatedRestrictionSites(resSite = "AAGCTT",overhangs5 = 1,chromosomes = "chr6",genomePack = "BSgenome.Hsapiens.UCSC.hg18",wingc=200,mappability = map_hg18,winmap = 500)
resfrag
imr90_500_chr6annot=setGenomicFeatures(imr90_500$chr6chr6,cutsites)
y_intervals(imr90_500_chr6annot)
imr90_500=normLGF(imr90_500_chr6annot)
?normLGF
?wget
?curl
library(curl)
bw=curl(url="http://hgdownload.cse.ucsc.edu/goldenPath/hg18/encodeDCC/wgEncodeMapability/")
bw
imr90_500_ICE=normICE(imr90_500,max_iter=10)
mapC(HTClist(imr90_500_ICE$chr6chr6),trim.range=0.95,col.pos=c("white","orange","red","black"))
hox=extractRegion(hic_imr90_40$chr6chr6,chr=c("chr6"),from=50e6,to=58e6)
plot(hox,maxrange=50,col.pos=c("white","orange","red","black"))
di<-directionalityIndex(hox)
barplot(di, col=ifelse(di>0,"darkred","darkgreen"))
