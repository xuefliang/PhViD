n<-500
#通过控制饮食而恢复正常的
diet<-0.1
effect<-c(0,0.95)
names(effect)<-c('FA','FB')

#对病患总体指派因素F
set.seed(1)
f.chance<-runif(n)
f<-ifelse(f.chance<0.9,'FA','FB')
table(f)

#指派控制组和治疗组
set.seed(1)
group<-runif(n)
group<-ifelse(group<0.5,'control','drug')
table(group)

#在关于药物D的临床实验中，将背景相似，病情相似的500名患者分为治疗组（针对病情控制饮食，服用药物D）和控制组（针对病情控制饮食，服用安慰剂）。
#设饮食的控制可以使10%的病患恢复正常。
#而由于某种未知的隐藏因素F（比如基因），会影响药物D正常的发挥作用。
#设病患的90%为FA，药物D是不起作用的；
#而病患的10%的病患为FB，其中95%可以通过药物D的治疗恢复正常。

#指派治疗结果
set.seed(1)
diet.chance<-runif(n)
drug.chance<-runif(n)
#| 表示或
#是饮食控制10%或者药物起作用,将药物组效果*0或0.95
outcome<-((diet.chance<diet)|(drug.chance<effect[f]*(group=='drug')))
result <- effect[f]*(group=='drug')
trail<-data.frame(group=group,F=f,treatment=outcome)
summary(trail)

#分因素治疗效果
table(trail$F,trail$treatment)

#控制组和治疗组
#控制组
with(trail[group=='control',],table(F,treatment))
with(trail[group=='drug',],table(F,treatment))

#分组的治疗效果
treat.group<-with(trail,table(group,treatment))

chisq.test(treat.group)

library(ggplot2)
p<-ggplot(data.frame(x=c(0,1)),aes(x=x))
p+stat_function(fun=dbeta,args=list(0.1,0.9),colour='red')

betad.mean<-function(alpha,beta)
{alpha/(alpha+beta)}
betad.mode<-function(alpha,beta)
{(alpha+1)/(alpha+beta-2)}

#先验分布的情况
alpha<-0.1
beta<-0.9
#分控制组和治疗组分别讨论
false.control<-treat.group[1,1]
true.control<-treat.group[1,2]
false.drug<-treat.group[2,1]
true.drug<-treat.group[2,2]
#对控制组，后验分布的参数
alpha.control<-alpha+true.control
beta.control<-beta+false.control
#对治疗组，后验分布的参数
alpha.drug<-alpha+true.drug
beta.drug<-beta+false.drug

p<-ggplot(data.frame(x=c(0,.3)),aes(x=x))
p+stat_function(fun=dbeta,args=list(alpha.drug,beta.drug),colour='red')+
  stat_function(fun=dbeta,args=list(alpha.control,beta.control),colour='blue')+
  annotate("text",x=.03,y=20,label="control")+
  annotate("text",x=.23,y=15,label="drug")

#计算控制组的均值和众数（p的后验估计）
betad.mean(alpha.control,beta.control)
betad.mode(alpha.control,beta.control)

#计算治疗组的均值和众数(p的后验估计)
betad.mean(alpha.drug,beta.drug)
betad.mode(alpha.drug,beta.drug)