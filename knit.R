install.packages("Seurat")
getwd()
.libPaths()
help("[[")
example(merge)
merge
help(GRanges)
.packages()
head(installed.packages(),n=2)
library(GenomicRanges)
data(package="car")
wd=getwd()
wd
basename(wd)
strsplit(wd,'/')
example("strsplit")
unlist(strsplit("A text I want to display with spaces"," "))
?strsplit
c="abc"
d="c"
c-d
wd="D:/str/sjrioa/sjdifa/3.bam"
wd
substr(basename,"",wd)
basename=basename(wd)
example(substr)
?substring
rep("abcdef", 4)
x=c("12345678","abcdeughah")
substr(x,c(2,4),c(4,5,7))
substring(x,c(2,4),c(4,5,8))
bases=c("A","T","G","C")
bases
DNA=paste(sample(bases,12,replace = T),collapse = "")
DNA
DNA
substring(DNA,seq(1,10,3),seq(3,12,3))
strtrim(c("abcdef","abcdef","abcdef"),c(1,5,10))
substr("123456","5","")
sub(unlist(strsplit("123456","")),"5","")
unlist(strsplit("123456",""))
c=c("1", "2", "3" ,"4", "5", "6")
c
sub("123","2","")
?sub
?substr
?sub
sub(pattern = "2",replacement = "","123")
bases
wd
pat=basename(wd)
sub(pat,"",wd)

nchar("hello world")
substr("abcdef",2,4)
x="1234567890"
x
substr(x,3,3)
substr(x,5,7)
substr(x,4,4)="A"
x
substr(x,2,4)="TTD"
x
substring(x,5)
substring(x,5)='...'
x
substr(x,9,12)="ABCD"
x
wd
substr(basename(wd),3,5)
x='123456789'
substr(x,c(2,4),c(4,5,8))
substring(x,c(2,4),c(4,5,8))
x=c("123456789","fhgqaqhu")
substr(x,c(2,4),c(4,5,8))
substring(x,c(2,4),c(4,5,8))
bases=c("A","T","G","C")
bases
DNA=paste(sample(bases,12,replace = T),collapse = "")
DNA
substring(DNA,seq(1,10,by=3),seq(3,12,by=3))
substr("abcdef",1:6,1:6)
substr("12346",c(2,3),c(4,5,6))
substring("123455566",c(2,3),c(4,5,6))
x="hello world"
grep(pattern = "hello",x)
grep(pattern = "l",x)
grepl(pattern = "l",x)
x=c("hello","sys","hi")
x
grep(pattern = "l",x=x,value = T)
grep(pattern = )
txt=c("arm","foot","lefroo","bafoobar")
grep("foo",txt)
library(BiocManager)
BiocManager::install("DESeq2")

grep("foo",txt)
txt
grep("oo",txt)

grep("foo",txt,ignore.case = T,invert = F,value = T)
txt
toupper("abc")
tolower("ABC")
x=c("My","First","Second")
x
tolower(x)
toupper(x)
casefold(x)
x=rnorm(5)
x
a=c(1,2,3,4)
a
b=c("one","two","three","four")
b
c=c(T,T,F,F)
c
data.frame(a,b,c)
?data.frame
tmp=data.frame('letter'=letters[1:10],'number'=1:10,'value'=c('+',"-"))
tmp
tmp[[1]]
do.call("paste",c(tmp,sep=""))
number_add=list(101:110,1:10)
number_add
add=function(x,y){x+y}
add
do.call(add,number_add)
add(number_add[[1]],number_add[[2]])
x1=1:10
x2=11:20
x3=21:30
data.frame(x1,x2,x3)
do.call("data.frame",list(x1,x2,x3))
list <- list(matrix(1:25, ncol = 5), matrix(4:28, ncol = 5), matrix(21:45, ncol=5))
list
list.sum=do.call("sum",list)
list.sum
sum(list)
do.call(cbind,list)
do.call(rbind,list)
eval(2+3)
?eval
eval(2 ^ 2 ^ 3)
iris[,1:4]
do.call(kmeans,list(x=iris[,1:4],centers=3))
c=1:10   
c
m1=matrix(c<-(1:10),nrow=5,ncol = 6)
m1
apply(m1,1,sum)
movies <- c("SPYDERMAN","BATMAN","VERTIGO","CHINATOWN")
movies
class(movies)
movies_lower=unlist(lapply(movies, tolower))
str(movies_lower)
dt=cars
dt
plot(cars$speed,cars$dist)
lapply(dt,min)
apply(dt,2,min)
sapply(dt,max,simplify = F)
sapply(dt,mean)
hist(dt$speed)
data(iris)
head(iris)
tapply(iris$Sepal.Width, iris$Species, mean)
x=1:9
y=1

cbind(x,y)
cbind(y,x)
x=cbind(x1=3,x2=c(4:1,2:5))
x
myfun=function(x,c1,c2){c(sum(x[c1],1),mean(x[c2]))}
myfun
sum(x[x1])    
apply(x,1,myfun,c1='x1',c2=c('x1','x2'))
sum(iris$Sepal.Length,by=iris$Species)
fivenum(6)
m=matrix(c(1:10),nrow=2)
m
apply(m, 1, sum)
apply(m,2,sum)
v=1:5
v
ind=c('a','a','a','b','b')
ind
tapply(v, ind, sum)
v
tapply(v,ind,fivenum)
fivenum(c(4,5))
df=data.frame(a=c(1:5),b=c(6:10),ind=ind)        
df
tapply(df$a,df$ind,mean)
df
df
lst=list(a=c(1:5),b=c(6:10))
lst
lapply(lst,mean)
lapply(lst,fivenum)
sapply(lst,fivenum)
res=vapply(lst, function(x) c(min(x),max(x)), FUN.VALUE = c(min.=0,max.=0))
res
mapply(sum, list(a=1,b=2,c=4),list(a=10,b=20,d=30))
mapply(function(x,y)x^y,c(1:5),c(1:5))
mapply(function(x,y)c(x+y,x^y),c(1:5),c(1:5))
lst=list(a=list(aa=c(1:5),ab=c(6:10)),b=list(ba=c(1:10)))
lst
rapply(lst,sum,how='replace')
df <- data.frame(year=kronecker(2001:2003, rep(1,4)), loc=c('beijing','beijing','shanghai','shanghai'), type=rep(c('A','B'),6), sale=rep(1:12))
df
tapply(df$sale, df[,c('year','loc')], sum)
tapply(df$sale,df[,c('year','type')],sum)
data(iris)
head(iris)
tapply(iris$Sepal.Width,iris$Species,median)
library(dplyr)
iris%>%group_by(Species)%>%summarise(mean(Sepal.Width))
f=function(x){x^6}
do.call(f,list(c(1:4)))
do.call
?do.call
showMethods("do.call")
tmp <- expand.grid(letters[1:2], 1:3, c("+", "-"))
tmp
do.call(paste,c(tmp,sep=""),quote = T)
do.call(paste, list(as.name("A"), as.name("B")), quote = T)
x=array(1:24,2:4)
x
install.packages("rmarkdown")
.libPaths()
list.files()
list.dirs()
getwd()
library(rmarkdown)
render("1-example.Rmd")


x=c(33,55,11)
seq(5)
rep(c(1,3),2)
sqrt(-1)
sqrt(-1+0i)
c(1,3)%in%c(2,3,4,5,6)
match(c(1,3),c(2,3,4,3))
which((11:15)>12)
duplicated(c(1,2,3,3,4,8))
substr("JAN07",1,3)
substr(c("jgoqhig","ejfoqjf"),4,5)
substring(c("JAN07","MArgjei"),4)
x='10,8,7'
strsplit(x,',',fixed = T)
x <- '1, 3; 5'
gsub(';',',',x,fixed = T)
strsplit(gsub(';',',',x,fixed = T),',')
?gsub

x=c(1,4,6.25)
x[-2]
x[-c(1,3)]
x[x>3]
which(x>3)
sex <- c("男", "男", "女", "女", "男", "女", "女", "女", "女", "男")
sex
sex.color=c('男'='blue','女'='red')
sex.color
sex.color[sex]
unname(sex.color[sex])
unique(c(1,5,2,5))
match(5,c(1,5,2))
intersect(c(5,7),c(1,3,5,7))
typeof(1:3)
typeof(c(1,2,3))
typeof(c(1,2.1,4))
is.integer(c(1L,3L))
class(factor(c('F',"M")))
attributes(x)
x
s=101:200
s
attr(s,'author')='li xiaoming'
s
attr(s,'date')='2016-09-12'
s
str(s)
today()
x=as.POSIXct(c('1992-12-24','2020-01-12'))
difftime(x[2],x[1],units = "years")
sex=factor(sex)
sex
table(sex)
length(sex)
h=1:10
tapply(h, sex, mean)
library(forcats)
set.seed(1)
fac=sample(c('red','green','blue'),30,replace = T)
fac
fac=factor(fac,levels = c('red','green','blue'))
fac
x=round(100*(10+rt(30,2)))
x
res1=tapply(x, fac, sd)
res1
barplot(res1)
f=function(){
  x=seq(0,2*pi,length=50)
  y1=sin(x)
  y2=cos(x)
  plot(x,y1,type='l',lwd=2,col="red",xlab='x',ylab='')
  lines(x,y2,lwd=2,col='blue')
  abline(h=0,col="gray")
}

f()
fsub=function(x,y=0){
  cat("x=",x,"y=",y,"\n")
  x-y
}
fsub(9,6)
body(fsub)
body(cat)
body(GRanges)
formals(fsub)
environment(fsub)
fsub(10)
do.call(fsub,list(3,y=1))        
sapply(1:5,fsub,y=2:5)
example(sapply)
sapply(1:5,"-",2)
fib1=function(n){
  if(n==0)return(0)
  else if(n==1) return(1)
  else if(n>=2){
    Recall(n-1)+Recall(n-2)
  }
}
fib1(10)
for(i in 0:10) cat('i=',i,'x[i]=',fib1(i),'\n')
sapply(1:10,fib1)
gv=function(x){ifelse(abs(x)<6,x^2,1)}
gc
gv(10)
sapply(1:10,gv)
d=data.frame(x=c(1,7,2),y=c(3,5,8))
d
Map(function(x)sum(x^2),d)
lapply(d, function(x)max(x))
Reduce(sum,1:4)
set.seed(2)
x <- replicate(4, sample(1:5, 5, replace=TRUE), simplify=FALSE); x
Reduce(intersect,x)
integrate(sin,0,pi)
uniroot(function(x)x*(x-1)*(x+1),c(0,0.1))
y=2:10
x=3:11
lm(y~x)
summary(lm(y~x))
x=rnorm(1000)
hist(x)
library(dplyr)
library(knitr)
kable(res)
?tibble
d.stu <- tibble(
  sid=c(1,2,3,4,5,6),
  cid=c(1,2,1,2,1,2),
  sname=c("John", "Mary", "James", "Kitty", "Jasmine", "Kim"),
  sex=c("M", "F", "M", "F", "F", "M"))

d.stu
kable(d.stu)
x=rnorm(30,mean=100,sd=1)
x
hist(x)
hist(x,col=rainbow(15),main = "正态随机数",xlab='',ylab='频数',freq = F)
tmp.dens=density(x)
lines(tmp.dens,lwd=2,col="red")
qqnorm(x)
qqline(x,lwd=2,col="red")
