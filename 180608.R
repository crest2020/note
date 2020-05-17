library(RTCGA)
library(RTCGA.clinical)
library(RTCGA.mRNA)
library(RTCGA.mutations)
all_TCGA_cancers=infoTCGA()
dim(all_TCGA_cancers)
rownames(all_TCGA_cancers)

expr=expressionsTCGA(BRCA.mRNA,OV.mRNA,LUSC.mRNA,extract.cols=c("GATA3","PTEN","XBP1","ESR1","MUC1"))
??expressionsTCGA
expressionsTCGA

browseVignettes("RTCGA")
remove.packages("RTCGA")
library(BiocInstaller)
biocLite("RTCGA")
library(RTCGA)
library(Rcpp)
?biocLite
install.packages("Rcpp",type = "binary")
remove.packages("Rcpp")
biocLite("RTCGA.clinical")
library(IRanges)
library(GenomicRanges)
ranges=IRanges(start=1,end=10)
ranges
?IRanges
?GRanges
grange=GRanges(seqnames="chr1",ranges=ranges,strand="+")
grange

library(RTCGA)
library(BiocInstaller)
biocLite("ggpubr")
library(RTCGA.clinical)
library(RTCGA.mRNA)
library(RTCGA.mutations)
expr=expressionsTCGA(BRCA.mRNA,OV.mRNA,LUSC.mRNA,extract.cols=c("GATA3","PTEN","XBP1","ESR1","MUC1"))
expr[1,]
table(substr(expr$bcr_patient_barcode,14,16))
unlist(strsplit(expr[1,]$bcr_patient_barcode,"-",fixed = T))[4]
boxplot(expr$GATA3)
?boxplot
table(expr$dataset)
stripchart(expr$GATA3,add=T,method="jitter",vertical = T,col="red")
library(ggpubr)
ggboxplot(expr,x="dataset",y="XBP1",title="XBP1",ylab="Expression",color="dataset",palette="jco")
expr
library(GEOquery)
gset=getGEO("GSE12345")
?getGEO
expr
nb_samples=table(expr$dataset)
nb_samples
expr$dataset=gsub(pattern = ".mRNA",replacement = "",expr$dataset)
table(expr$dataset)
expr$bcr_patient_barcode=paste0(expr$dataset,c(1:590,1:561,1:154))
expr
library(ggpubr)
ggboxplot(expr,x="dataset",y="ESR1",title="ESR1",ylab="Expression",color="dataset",palette = "simpsons")
my_conparisons=list(c("BRCA","OV"),c("OV","LUSC"))
ggboxplot(expr,x="dataset",y="GATA3",title="GATA3",ylab="Expression",color="dataset",palette="jco")+stat_compare_means(comparisons =my_conparisons)
compare_means(c(GATA3,PTEN,XBP1)~dataset,data=expr)
label.select.criteria=list(criteria="`y`>3.1 & `x`%in%c('BRCA','OV')")
label.select.criteria
ggboxplot(expr,x="dataset",y=c("GATA3","PTEN","XBP1"),merge=T,color="dataset",palette="jco",ylab="Expression",label="bcr_patient_barcode",label.select =label.select.criteria,font.label=list(size=9,face="italic"),repel=T)
?ggboxplot

ggboxplot(expr,x="dataset",y="GATA3",title = "GATA3",ylab="Expression",color="dataset",palette = "jco",rotate=T)
library(RTCGA)
library(RTCGA.rnaseq)
expr=expressionsTCGA(BRCA.rnaseq,OV.rnaseq,LUSC.rnaseq,extract.cols = c("GATA3|2625", "PTEN|5728", "XBP1|7494","ESR1|2099", "MUC1|4582"))
expr
head(expr$bcr_patient_barcode)
nb_samples=table(expr$dataset)

expr %>% (substr(bcr_patient_barcode,14,15)=="01")

head(expr[substr(expr$bcr_patient_barcode,14,15)=="01",])
nb_samples
ggboxplot(expr,x="dataset",y="`PTEN|5728`",title="ESR1|2099",ylab="Expression",color="dataset",palette = "jco")
expr
sum(is.na(expr))
head(rnaseq.5genes.3cancers)

library(RTCGA.rnaseq)
library(dplyr)

expressionsTCGA(BRCA.rnaseq,OV.rnaseq,HNSC.rnaseq)%>%
  dplyr::rename(cohort=dataset)%>%
  filter(substr(bcr_patient_barcode,14,15)=="01")->BRCA.OV.HNSC.rnaseq.cancer
head(BRCA.OV.HNSC.rnaseq.cancer)
pcaTCGA(BRCA.OV.HNSC.rnaseq.cancer,"cohort")->pca_plot
expr
expr%>%rename(cohort=dataset)%>%head
plot(pca_plot)
expr[1:100,]%>%filter(substr(bcr_patient_barcode,14,15)=="01")->rnaseq.genes.cancers
rnaseq.5genes.3cancers
DT::datatable(rnaseq.genes.cancers)
pcaTCGA(rnaseq.genes.cancers[,c(1,2,4)],"dataset")->pca_plot1
dim(rnaseq.5genes.3cancers)
dim(expr)
rnaseq.genes.cance
table(expr$dataset)
table(is.na(rnaseq.genes.cancers))
expr=expressionsTCGA(BRCA.rnaseq,OV.rnaseq,LUSC.rnaseq,extract.cols = c("GATA3|2625", "PTEN|5728", "XBP1|7494","ESR1|2099", "MUC1|4582"))
?pcaTCGA
library(RTCGA.mutations)
library(survminer)

mu=mutationsTCGA(BRCA.mutations,OV.mutations)
bcrcode=mu$bcr_patient_barcode
table(substr(bcrcode,14,15))
dim(mu)
head(mu)
?subset
mu_p53=subset(mu,substr(bcr_patient_barcode,14,15)=="01"& Hugo_Symbol=="TP53")
head(mu_p53)
library(RTCGA.clinical)
library(dplyr)
sur_p53=survivalTCGA(BRCA.clinical,OV.clinical,extract.cols = "admin.disease_code")
head(sur_p53)
head(mu)
mu_p53$bcr_patient_barcode=substr(mu_p53$bcr_patient_barcode,1,12)
head(mu)
mu_p53=rename(mu_p53,disease=admin.disease_code)
?rename
head(sur_p53)
mer_p53=merge(sur_p53,mu_p53,by="bcr_patient_barcode",all.x=T)

head(mu_p53)
head(mer_p53)
table(is.na(mer_p53$Variant_Classification))


head(sur_p53)
dim(mu_p53)
?merge
dim(mer_p53)
head(mer_p53)
mer_p53=rename(mer_p53,TP53=Variant_Classification)
head(mer_p53)
mer_p53$TP53=ifelse(!is.na(mer_p53$TP53),"Mut","WILDorNOINFO")
head(mer_p53)

BRCA_OV.plot=subset(mer_p53,select=c("times","patient.vital_status","disease","TP53"))
head(BRCA_OV.plot)
km_plot=kmTCGA(BRCA_OV.plot,explanatory.names = c("TP53","disease"),break.time.by=400,xlim=c(0,2000),pval = T)
print(km_plot)
png("survival.jpg")
dev.off()
library(RTCGA.rnaseq)
exp=expressionsTCGA(ACC.rnaseq,BLCA.rnaseq,BRCA.rnaseq,OV.rnaseq,extract.cols = c("MET|4233","ZNF500|26048","ZNF501|115560"))
head(exp)
heatmapTCGA(exp,"dataset","`MET|4233`","`ZNF500|26048`")
h=rename(exp,cohort=dataset,MET=`MET|4233`)

head(h)
round(quantile(h$MET,prob=seq(0,1,0.25)),-2)
seq(0,1,0.25)
?cut
h$MET=cut(h$MET,round(quantile(h$MET,prob=seq(0,1,0.25)),-2),include.lowest = T,diag.lab=5)
head(h)
table(h$MET)
hh=h[,-1]
head(hh)
group=group_by(cohort,MET)
h=rename(h,ZNF500=`ZNF500|26048`,ZNF501=`ZNF501|115560`)
head(h)
aggdata=aggregate(hh,by=list(cohort,MET),FUN=median)
library(car)
agg=aggregate(mtcars,by=list(cyl,gear),FUN=mean)
library(stats)
?aggregate
remove.packages("stats")
install.packages("stats")
library(BiocInstaller)
biocLite("stats")
library(stats)
library(car)
install.packages("stats")
