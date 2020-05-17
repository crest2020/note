
bases=c("A","T","G","C")
DNA=paste(sample(bases,12,replace = T),collapse = "")
DNA
DNA
substring(DNA,seq(1,10,3),seq(3,12,3))


install.packages("ggpubr")
library(ggpubr)
set.seed(1234)
df=data.frame(sex=factor(rep(c('f','m'),each=200)),weight=c(rnorm(200,55),rnorm(200,58)))
head(df)
example(rnorm)

?rnorm
head(df)
ggdensity(df,x="weight",add="mean",rug = T,color = "sex",fill='sex',palette = c('#00AFBB','#E7B800'))
gghistogram(df,x="weight",add="mean",rug=T,color = "sex",fill='sex',palette = c('#00AFBB',"#E7B800"))
data('ToothGrowth')
df1=ToothGrowth
df1
head(df1)
table(df1$dose)
p=ggboxplot(df1,x='dose',y='len',color = 'dose',palette = c('#00AFBB','#E7B800',"#FC4E07"),add='jitter',shape='dose')
p

my_comparisons=list(c("0.5",'1'),c('1','2'),c('0.5','2'))
ggviolin(df1,x='dose',y='len',fill='dose',palette = c('#00AFBB',"#E7B800","#FC4E07"),add="boxplot",add.params = list(fill="white"))+stat_compare_means(comparisons = my_comparisons,label = "p.signif")+stat_compare_means(label.y = 50)

p+stat_compare_means(comparisons = my_comparisons)+stat_compare_means(label.y = 50)
data("mtcars")
df2=mtcars
head(df2)
df2$cyl=factor(df2$cyl)
head(df2)
table(df2$cyl)
df2$name=rownames(df2)
colnames(df2)
head(df2[,c('name','wt','mpg','cyl')])
ggbarplot(df2,x='name',y='mpg',fill='cyl',color='white',palette = 'npg',sort.val = "desc",sort.by.groups = F,x.text.angle=60)


ggbarplot(df2,x='name',y='mpg',fill='cyl',color = 'white',palette = 'aaas',sort.val = 'asc',sort.by.groups = T,x.text.angle=60)

df2$mpg_z=(df2$mpg-mean(df2$mpg))/sd(df2$mpg)
df2$mpg_grp=factor(ifelse(df2$mpg_z<0,'low','high'),levels = c('low','high'))
table(df2$mpg_grp)
?split
?cut
ggbarplot(df2,x='name',y='mpg_z',fill="mpg_grp",color="white",palette = 'jco',sort.val = 'asc',sort.by.groups = F,x.text.angle=40,ylab="MPG z-score",xlab=F,lengend.title="MPG Group")
ggbarplot(df2,x='name',y="mpg_z",fill="mpg_grp",color="white",palette = 'jco',sort.val = "desc",sort.by.groups = F,x.text.angle=90,ylab="MPG z-score",xlab=F,legend.title='MPG Group',rotate=T,ggtheme = theme_minimal())
ggdotchart(df2,x='name',y='mpg',color='cyl',palette = c('#00AFBB','#E7B800',"#FC4E07"),sorting="ascending",add='segments',ggtheme = theme_pubr())

ggdotchart(df2,x='name',y='mpg',color = 'cyl',palette = c('#00AFBB','#E7B800',"#FC4E07"),sorting='descending',add='segments',rotate = T,group = 'cyl',dot.size = 6,label=round(df2$mpg),font.label = list(color='white',size=9,vjust=0.5),ggtheme = theme_pubr())
ggdotchart(df2,x='name',y='mpg_z',color='cyl',palette = c("#00AFBB","#ECB800","#FC4E07"),sorting = "descending",add='segments',add.params=list(color="lightgray",size=2),group = "cyl",dot.size = 6,label = round(df2$mpg_z,1),font.label = list(color="white",size=9,vjust=0.5),ggtheme = theme_pubr())+geom_hline(yintercept = 0,linetype=2,color="lightgray")
head(dfm)
ggdotchart(df2,x='name',y='mpg',color='cyl',palette = c('#00AFBB',"#E7B800","#FC4E07"),sorting = "descending",rotate = T,dot.size = 2,y.text.col=T,ggtheme = theme_pubr())+theme_cleveland()
ggboxplot(df2,y="mpg",x="cyl",color = "cyl",shape="cyl",add='jitter')
example("ggdensity")


library(ggplot2)
ggplot(diamonds,aes(x=carat,y=price,color=cut))+geom_point(alpha=0.1,size=1,shape=21,fill='white',stroke=2)+geom_smooth(method='glm')
ggplot(diamonds,aes(x=carat,y=price,color=cut))+geom_point()+scale_color_brewer(type='qual',palette = 8)
gg=ggplot(diamonds,aes(x=carat,y=price,color=cut))+geom_point()+scale_color_brewer(type="qual",palette = 8)+labs(title='scatterplot',x='carat',y='price')

gg1=gg+theme(plot.title = element_text(size=30,face = "bold"),axis.text.x = element_text(size=15),axis.text.y = element_text(size=15),axis.title.x = element_text(size=25),axis.title.y = element_text(size=25))
gg1

theme_set(theme_classic())
#plot
g=ggplot(mpg,aes(cty))
g+geom_density(aes(fill=factor(cyl)),alpha=0.8)+labs(title="Denstity plot",subtitle="City Mileage Grouped by number of cylinders",caption="source.mpg",x='city mileage',fill="# cylinders")
g=ggplot(mpg,aes(class,cty))
library(ggplot2)
library(ggpubr)
g+geom_violin()+labs(title="vionlin plot",subtitle = "city mileage vs class of vehivle",caption="source mpg",x='class of vehicle',y='city mileage')

img.file=system.file(file.path("images","background-image.png"),package = "ggpubr")
img=png::readPNG(img.file)
img
.libPaths()
ggplot(iris,aes(Species,Sepal.Length))+background_image(img)+geom_boxplot(aes(fill=Species),color='white')+fill_palette("jco")
list.files(pattern="\\.count$")
?list.files

ex=function(x){
  result=read.table(list.files(pattern = "\\.count$"))
  return(result)
}
name=  list.files(pattern = "\\.count$")

a=lapply(name, read.table)
b=sapply(name,read.table)
