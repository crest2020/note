library(swirl)
cars
which(cars[,2]==85)
mean(cars[,2])
names(cars)[2]
nrows(cars)
nrow(cars)
install.packages("downloader")
library(downloader) 
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleMiceWeights.csv"
filename <- "femaleMiceWeights.csv" 
download(url, destfile=filename)

filename
head(filename)
dat<-read.csv(filename)
head(dat)
dat[12,2]
dat$Bodyweight[11]
length(dat$Diet)
length(dat$Bodyweight)
dat
mean(dat$Bodyweight[13:24])
set.seed(1)
sample(13:24,1)
?sample
dat[16,2]
library(downloader)
url="https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/msleep_ggplot2.csv"
filename <- basename(url)
download(url,filename)
dat1<-read.csv(filename)
head(dat1)
install.packages("dplyr")
class(dat1)
dat1
names(dat1)
library(dplyr)
dat2<-filter(dat1,order=="Primates")
filter(dat, order=="Primates") %>% summarise( mean( sleep_total) )
filter(dat1, order=="Primates") %>% summarise( mean( sleep_total) )
?summarise
library(downloader) 
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename <- basename(url)
download(url, destfile=filename)
x <- unlist( read.csv(filename) )

set.seed(1)
for(i in 1:10000){mean1[i]<-mean(sample(x,5))}
sum(abs(mean1-mean(x))>1)

install.packages("gapminder")

hist(mean1)
library(gapminder)
data(gapminder)
head(gapminder)

life<-filter(gapminder,year==1952)
head(life)
a<-select(life,pop)
head(a)
plot(life$country,life$lifeExp)

sum(life$lifeExp<=40)/nrow(life)


dat1952<-gapminder[gapminder$year==1952,]
head(dat1952)
class(gapminder)
x<-dat1952$lifeExp
mean(x<=40)
mean(x<=60 )-mean(x<=40)
mean(x<=60 & x>40)
plot(ecdf(x))
prop<-function(q){mean(x<=q)}
prop(40)
qs = seq(from=min(x), to=max(x), length=20)
props = sapply(qs, prop)
plot(qs,props)
b<-read.csv("femaleControlsPopulation.csv")
head(b)
x<-unlist(b)


set.seed(1)
mean5<-vector()
mean50<-vector()
for(i in 1:1000){mean5[i]<-mean(sample(x,5))}
for(i in 1:1000){mean50[i]<-mean(sample(x,50))}

par(1,2)
hist(mean5)
hist(mean50)
mean(mean50<=25)-mean(mean50<=23)
mean( mean50 < 25 & mean50 > 23)
pnorm( (25-23.9) / 0.43)  - pnorm( (23-23.9) / 0.43) 

library(downloader) 
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/mice_pheno.csv"
filename <- basename(url)
download(url, destfile=filename)
dat <- read.csv(filename) 

dat <- na.omit( dat )
x<-filter(dat,Sex=="M",Diet=="chow")
head(x)
mean(x$Bodyweight)
library(rafalib)
popsd(x)

y<-filter(dat,Diet=="hf" & Sex=="M") %>%select(Bodyweight) %>% unlist()
mean(y)
install.packages("rafalib")
library(rafalib)
?popsd()
?sd()
?set.seed()
ph<-read.csv("mice_pheno.csv")
head(ph)
yy<-filter(dat,Diet=="chow" & Sex=="M")%>%select(Bodyweight)%>%unlist()
sd<-popsd(yy)
sd
mean(abs(yy-mean(yy))<sd)
mypar(1,1)
class(y)
qqnorm(yy);abline(0,1)

class(yy)
mypar(2,2)
y <- filter(dat, Sex=="M" & Diet=="chow") %>% select(Bodyweight) %>% unlist
z <- ( y - mean(y) ) / popsd(y)
qqnorm(z);abline(0,1)
y <- filter(dat, Sex=="F" & Diet=="chow") %>% select(Bodyweight) %>% unlist
z <- ( y - mean(y) ) / popsd(y)
qqnorm(z);abline(0,1)
y <- filter(dat, Sex=="M" & Diet=="hf") %>% select(Bodyweight) %>% unlist
z <- ( y - mean(y) ) / popsd(y)
qqnorm(z);abline(0,1)
y <- filter(dat, Sex=="F" & Diet=="hf") %>% select(Bodyweight) %>% unlist
z <- ( y - mean(y) ) / popsd(y)
qqnorm(z);abline(0,1)


y <- filter(dat, Sex=="M" & Diet=="chow") %>% select(Bodyweight) %>% unlist
avgs <- replicate(10000, mean( sample(y, 25)))
mypar(1,2)
hist(avgs)
qqnorm(avgs)
qqline(avgs)

mean(avgs)
popsd(avgs)
sd(avgs)


dat<-read.csv("femaleMiceWeights.csv")

set.seed(1)
replicate(1000,{x<-sample(1:6,100,replace=TRUE)
pp<-mean(x==6)})

?replicate
help(replicate)


set.seed(1)
pp<-replicate(10000,{x<-sample(1:6,100,replace=TRUE);mean(x==6)})
p=1/6
z<-(pp-p)/sqrt(p*(1-p)/100)
qqnorm(z)
abline(z)
qqline(z)
abline(0,1)
?rafalib

2*(1-pnorm(sqrt(12)*2/3.022541))

1 - pt(3,df=3)
1 - pt(3,df=15)
1 - pt(3,df=30)
2*(1 - pnorm(2.055174))

url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/babies.txt"
filename <- basename(url)
download(url, destfile=filename)
babies <- read.table("babies.txt", header=TRUE)

bwt.nonsmoke <- filter(babies, smoke==0) %>% select(bwt) %>% unlist 
bwt.smoke <- filter(babies, smoke==1) %>% select(bwt) %>% unlist

library(rafalib)
mean(bwt.nonsmoke)-mean(bwt.smoke)
popsd(bwt.nonsmoke)
popsd(bwt.smoke)

set.seed(1)
dat.ns<-sample(bwt.nonsmoke,25)
dat.s<-sample(bwt.smoke,25)
tval<-(mean(dat.ns)-mean(dat.s))/sqrt(sd(dat.ns)^2/25+sd(dat.ns)^2/25)

2*(1-pnorm(2.120904))
2*(1-pnorm(3))
2.576*sqrt((sd(dat.ns)^2/25+sd(dat.ns)^2/25))
tval



X.ns <- mean(dat.ns)
sd.ns <- sd(dat.ns)

X.s <- mean(dat.s)
sd.s <- sd(dat.s)

sd.diff <- sqrt(sd.ns^2/N+sd.s^2/N)
tval <- (X.ns - X.s)/sd.diff
tval

X.ns - X.s
sd.diff
2.326*sd.diff
N <- 25
qnorm(0.995)*sqrt( sd( dat.ns)^2/N + sd( dat.s)^2/N )

mn<-matrix(,10000,4)
set.seed(1)
N<-c(30,60,90,120)
for(i in 1:4){
mn[,i]<-replicate(10000,{

dat.ns <- sample(bwt.nonsmoke , N[i])
dat.s <- sample(bwt.smoke , N[i])
t.test(dat.s,dat.ns)$p.value})}

t.test(dat.s,dat.ns)

set.seed(1)
b<-rnorm(5)
tv<-mean(b)/sd(b)*sqrt(5)
tv
set.seed(1)
tv<-replicate(1000,{b<-rnorm(5)
mean(b)/sd(b)*sqrt(5)})

mean(tv>2)
1-pt(2,df=4)


library(rafalib)
mypar(3,2)

Ns<-seq(5,30,5)
B <- 1000
mypar(3,2)
LIM <- c(-4.5,4.5)
for(N in Ns){
  ts <- replicate(B, {
    X <- rnorm(N)
    sqrt(N)*mean(X)/sd(X)
  })
  ps <- seq(1/(B+1),1-1/(B+1),len=B)
  qqplot(qt(ps,df=N-1),ts,main=N,
         xlab="Theoretical",ylab="Observed",
         xlim=LIM, ylim=LIM)
  abline(0,1)
} 

set.seed(1)
N <- 5
B<- 1000

tstats <- replicate(B,{
  X <- rnorm(N)
  sqrt(N)*mean(X)/sd(X)
})
mean(tstats>2)


library(rafalib)
mypar(3,2)
lim<-c(-3.5,3.5)
B=1000
Ns<-seq(5,30,5)
for (N in Ns){
  tss<-replicate(B,{
    x1<-rnorm(N)
    x2<-rnorm(N)
    t.test(x1,x2,var.equal = TRUE)$statistic
  })
  ps<-seq(1/(B+1),1-1/(B+1),len=B)
  qqplot(qt(ps,df=2*N-2),tss,main=N,xlab="Therotical",ylab="Observed",xlim=lim,ylim=lim)
  abline(0,1)
}
?t.test



library(rafalib)
set.seed(1)
lim<-c(-3.5,3.5)
B=1000
N<-100

  tss<-replicate(B,{
    x<-rnorm(N)
    
    sqrt(N)*median(x)/sd(x)
  })
  ps<-seq(1/(B+1),1-1/(B+1),len=B)
  qqplot(qnorm(ps),tss,main=N,xlab="Therotical",ylab="Observed",xlim=lim,ylim=lim)
  abline(0,1)

  
  ?range
  
  
  set.seed(1)
  Ns <- seq(5,45,5)
  library(rafalib)
  mypar(3,3)
  for(N in Ns){
    medians <- replicate(10000, median ( rnorm(N) ) )
    title <- paste("N=",N,", avg=",round( mean(medians), 2) , ", sd/sqrt(N)=", round( sd(medians)/sqrt(N),2) )
    qqnorm(medians, main = title )
    qqline(medians)
  }
  ?round
  nss<-1/sqrt(Ns)
  nss
  
  library(dplyr)
  babies <- read.table("babies.txt", header=TRUE)
  bwt.nonsmoke <- filter(babies, smoke==0) %>% select(bwt) %>% unlist 
  bwt.smoke <- filter(babies, smoke==1) %>% select(bwt) %>% unlist
  
  
  N=10
  set.seed(1)
  nonsmokers <- sample(bwt.nonsmoke , N)
  smokers <- sample(bwt.smoke , N)
  obs <- mean(smokers) - mean(nonsmokers)
  set.seed(1)
 null<-replicate(1000,{
  dat <- c(smokers,nonsmokers)
  shuffle <- sample( dat )
  smokerstar<-shuffle[1:N]
  nonsmokerstar<-shuffle[(N+1):(2*N)]
  median(smokerstar)-median(nonsmokerstar)})
  
obs1<-median(smokers)-median(nonsmokers)
(sum(abs(null)>=abs(obs1))+1)/(length(null)+1)

d<-read.csv("assoctest.csv")
head(d)


head(InsectSprays)
tail(InsectSprays)
boxplot(InsectSprays$count ~ InsectSprays$spray)
install.packages("UsingR")
library(dplyr)
data(nym.2002, package="UsingR")

head(nym.2002)
mal<-filter(nym.2002,gender=="Male") %>%select(time)%>%unlist
fem<-filter(nym.2002,gender=="Female")%>%select(time)%>%unlist

boxplot(nym.2002$gender,nym.2002$time)
hist(mal)
hist(fem)
boxplot(mal,fem)

library(rafalib)
mypar(1,3)
hist(mal)
hist(fem)
boxplot(mal,fem)

mal1<-filter(nym.2002,gender=="Male") 
fem1<-filter(nym.2002,gender=="Female")
cor(mal1[,c("age","time")],method="pearson")
cor(fem1[,c("age","time")],method="pearson")
cor(fem1$age,fem1$time)
?scatterplot
scatterplots(nym.2002$age~nym.2002$time)
boxplot(nym.2002$time~nym.2002$age)

library(rafalib)
mypar(2,2)
plot(fem1$age, fem1$time)
plot(mal1$age, mal1$time)
group <- floor(fem1$age/5) * 5
boxplot(fem1$time~group)
group <- floor(mal1$age/5) * 5
boxplot(mal1$time~group)

group<-floor(fem1$age/5)*5

time <-sort(nym.2002$time)
min(time)/median(time)

max(time)/median(time)


par(mfrow=c(1,2))
plot(time/median(time), ylim=c(1/4,4))
abline(h=c(1/2,1,2))

plot(log2(time/median(time)),ylim=c(-2,2))
abline(h=-1:1)

nrow(nym.2002)
plot(time)
?reshape

data(ChickWeight)
head(ChickWeight)
plot( ChickWeight$Time, ChickWeight$weight, col=ChickWeight$Diet)

install.packages("reshape")
nrow(ChickWeight)
chick = reshape(ChickWeight, idvar=c("Chick","Diet"), timevar="Time",
                direction="wide")

head(chick)
head(ChickWeight,40)
chick = na.omit(chic
                
nrow(chick)
(sum(chick$weight.4)+3000)/46/mean(chick$weight.4)
median(c(chick$weight.4,3000))/median(chick$weight.4)
mad(c(chick$weight.4,3000))/mad(chick$weight.4)
plot(chick$weight.4,chick$weight.21)
cor(c(chick$weight.4,3000),c(chick$weight.21,3000))/cor(chick$weight.4,chick$weight.21)
x<-chick$weight.4[chick$Diet==1]
x
y<-chick$weight.4[chick$Diet==4]
t.test(x,y)
wilcox.test(x,y)
wilcox.test(c(x,200),y,exact=FALSE)


library(rafalib)
mypar(1,3)
boxplot(x,y)
boxplot(x,y+10)
boxplot(x,y+100)
t.test(x,y+10)$statistic-t.test(x,y+100)$statistic
x1<-c(1,2,3)
y1<-c(400,500,600)
wilcox.test(x1,y1)


datt<-read.csv("assoctest.csv")
sum(datt$allele)
sum(datt$case)
?table
datt1<-table(datt)
datt1
