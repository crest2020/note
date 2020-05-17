install.packages("UsingR")
library(UsingR)
data("father.son",package="UsingR")

head("father.son")
class("father.son")
dat<-read.table("father.son")
view("father.son")
data(father.son,package="UsingR")
class(father.son)
head(father.son)
nrow(father.son)
mean(father.son$sheight)

?round
round(2.51)
x<-father.son$sheight[round(father.son$fheight)==71]
mean(x)

X<- matrix(1:1000,100,10)
X[25,3]
matrix(1:60,20,3,byrow=TRUE)

x1<-matrix(1:10,10,1)
for(i in 2:5){
  x2<-cbind(x1,i*x1)
  x1<-x2
}
x1<-matrix(1:10,10,1)
x3<-cbind(x1,2*x1,3*x1,4*x1,5*x1)
sum(x3[7,])
?solve
a<-matrix(c(3,4,-5,1,2,2,2,-1,1,-1,5,-5,5,0,0,1),4,4,byrow=TRUE)
b<-matrix(c(10,5,7,4),4,1)
solve(a,b)

a <- matrix(1:12, nrow=4)
b <- matrix(1:15, nrow=3)
d<-a%*%b
d[3,2]
a[3,]*b[,2]
?var
?rnorm

beta<-c(5,2)
X <- matrix(c(1,1,1,1,0,0,1,1),nrow=4)
rownames(X) <- c("a","a","b","b")

fitted <- X %*% beta

fitted[ 1:2, ]

?lse
??lse
rnorm(25)

h0= 56.67
g=9.8
v0=0
n = 25
tt = seq(0,3.4,len=n)
set.seed(1)
gg<-replicate(100000,{
y <- h0 + v0*tt -0.5*g*tt^2+rnorm(n,sd=1)
X <- cbind(1,tt,tt^2)

A <- solve(crossprod(X))%*%t(X)
-2*(A%*%y)[3]})

sd(gg)
sqrt(mean((gg-mean(gg))^2))


library(UsingR)


X = father.son$fheight


Y = father.son$sheight


n = length(Y)
set.seed(1)
betahat<-replicate(10000,{
N =  50


index = sample(n,N)


sampledat = father.son[index,]


x = sampledat$fheight


y = sampledat$sheight


lm(y~x)$coef[2]})

mean( (Y - mean(Y))*(X-mean(X) ) )
sd(betahat)

lm(Y~X)$fitted.values




library(UsingR)

X = father.son$fheight

Y = father.son$sheight

n = length(Y)

N = 50

set.seed(1)

index = sample(n,N)

sampledat = father.son[index,]

x = sampledat$fheight

y = sampledat$sheight

betahat = lm(y~x)$coef
fit<-lm(y~x)
sum((fit$fitted.value-y)^2)
class(betahat)
betahat1<-matrix(c(betahat[1],betahat[2]),1,2)
r<-y-x%*%betahat1
RSS<-crossprod(r)
RSS

X1 <- cbind(rep(1,N), x)
XX<-solve(crossprod(X1))
XX[1,1]
sigma<-sum((fit$fitted.value-y)^2)/48
sqrt(sigma%*%diag(XX))
summary(fit)
XX[1,2]
nrow(XX)
ncol(XX)

?formula
?model.matrix
?factor
?split

n <- 10; nn <- 100
g <- factor(round(n * runif(n * nn)))
x <- rnorm(n * nn) + sqrt(as.numeric(g))
xg <- split(x, g)
boxplot(xg, col = "lavender", notch = TRUE, varwidth = TRUE)

?runif
?rnorm
qqplot(x2)
qqnorm(x2)
abline(0,1)

species <- factor(c("A","A","B","B"))
condition <- factor(c("control","treated","control","treated"))
fit1<-model.matrix(~ species + condition)
fit1
install.packages("contrast")
library(contrast)
atob<-contrast(fit1,list(species="B",condition="control"),list(species="A",condition="treated"))
atob
?contrast

url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/spider_wolff_gorb_2013.csv"
filename <- "spider_wolff_gorb_2013.csv"
library(downloader)
if (!file.exists(filename)) download(url, filename)

spider <- read.csv(filename, skip=1)

head(spider)

fit3<-lm(friction~type+leg,data=spider)
contrast(fit3,list(leg="L4"),list(leg="L2"))
spidersub<-subset(spider,leg=="L2"|leg=="L4")
spidersub
fit4<-lm(friction~type+leg,data=spider)

library(contrast)
contrast(fit4,list(leg="L4",type="pull"),list(leg="L2",type="pull"))
contrast(fit4,list(leg="L4",type="push"),list(leg="L2",type="push"))

contrast(fit4,list(leg="L4",type="pull"),list(leg="L2",type="push"))
contrast(fit4,list(leg="L4",type="push"),list(leg="L2",type="pull"))

spider$group<-factor(paste0(spider$leg,spider$type))
head(spider)

?paste0
?cov
summary(fit4)
X <- model.matrix(~ type + leg, data=spider)
(Sigma <- sum(fit4$residuals^2)/(nrow(X) - ncol(X)) * solve(t(X) %*% X))
summary(fit4)
C <- matrix(c(0,0,-1,0,1),1,5)
(C%*%Sigma%*%t(C))

spider$log2friction <- log2(spider$friction)
boxplot(log2friction ~ type*leg, data=spider)
fit5<-lm(log2friction~type*leg,data=spider)
summary(fit5)
anova(fit5)
head(spider)
nrow(spider)
l4<-spider[spider$leg=="L4",]
l4
fit6<-lm(log2friction~type,data=l4)
summary(fit6)
library(multcomp)
attach(cholesterol)
head(cholesterol)

table(trt)
table(response)

aggregate(response,by=list(trt),FUN=mean)
aggregate(response,by=list(trt),FUN=sd)

fit<-aov(response~trt)
summary(fit)
install.packages("gplots")
library(gplots)
plotmeans(response~trt,xlab="Treatment",ylab="Response",main="Mean Plot\nwith 95% CI")
detach(cholesterol)

library(contrast)
l2vsl1pull<-contrast(fit5,list(leg="L2",type="pull"),list(leg="L1",type="pull"))
l2vsl1pull

l2vsl1push<-contrast(fit5,list(leg="L2",type="push"),list(leg="L1",type="push"))
l2vsl1push

set.seed(1)
f.value <-replicate(1000,{
N <- 40
p <- 4
group <- factor(rep(1:p,each=N/p))
X <- model.matrix(~ group)

Y <- rnorm(N,mean=42,7)
mu0 <- mean(Y)
initial.ss <- sum((Y - mu0)^2)

s <- split(Y, group)
after.group.ss <- sum(sapply(s, function(x) sum((x - mean(x))^2)))

group.ss <- initial.ss - after.group.ss

group.ms <- group.ss / (p - 1)

X
after.group.ms <- after.group.ss / (N - p)
 group.ms / after.group.ms
})
mean(f.value)
hist(f.value, col="grey", border="white", breaks=50, freq=FALSE)
xs <- seq(from=0,to=6,length=100)
N<-40;p<-4
lines(xs, df(xs, df1 = p - 1, df2 = N - p), col="red")


group <- factor(rep(1:p,each=N/p))
X <- model.matrix(~ group)
X


sex <- factor(rep(c("female","male"),each=4))
trt <- factor(c("A","A","B","B","C","C","D","D"))
X <- model.matrix( ~ sex + trt)
qr(X)$rank
Y <- 1:8

#sum( ( Y - X_male* beta_male - X_D*beta_D - X_R*beta_R )^2 )

makeYstar <- function(a,b) {Y - X[,2] * a - X[,5] * b}

fitTheRest <- function(a,b) {
  Ystar <- makeYstar(a,b)
  Xrest <- X[,-c(2,5)]
  betarest <- solve(t(Xrest) %*% Xrest) %*% t(Xrest) %*% Ystar
  residuals <- Ystar - Xrest %*% betarest
  sum(residuals^2)
}

fitTheRest(1,2)
fitTheRest(8,-2)
fitTheRest(6,0)
fitTheRest(1,5)
X
X[,-c(2,5)]

expand.grid(1:3,1:3)
betas <- expand.grid(-2:8,-2:8)

rss <- apply(betas,1,function(x) fitTheRest(x[1],x[2]))
rss
library(rafalib)
themin=min(rss)

plot(betas[which(rss==themin),])



fit <- lm(friction ~ type + leg, data=spider)
betahat <- coef(fit)
Y <- matrix(spider$friction, ncol=1)

X <- model.matrix(~ type + leg, data=spider)
QR<-qr(X)
Q<-qr.Q(QR)
R<-qr.R(QR)
Q[1,1]
R[1,1]

crossprod(Q,Y)
solve.qr(QR,y)

chol2inv(QR$qr)
?chol2inv
crossprod(X)
