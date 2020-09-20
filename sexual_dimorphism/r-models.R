rm(list=ls())

#reading data
a<-read.table('raw_data.txt',header=T)
names(a)
a$height[a$height < 10]<-100*a$height[a$height < 10]

#add 0 for missing hormonal data
a$horm2<-rep(0,length(a$horm))
a$horm2[a$horm == '1']<- -0.5
a$horm2[a$horm == '2']<- 0.5

a$looking2<-rep(0,length(a$looking))
a$looking2[a$looking == '1']<- 1
a$looking2[a$looking == '2']<- 2
a$looking2[a$looking == '3']<- 3

#analyses for heterosexual women
a<-a[a$hetero == '1' & a$age<36 & a$weight>0,]

#converting the scan ratings to range [-4.5, 4.5] (from feminine to masculine) 
a$scan1b<-a$scan1-4.5
a$scan2b<-(9-a$scan2)-4.5
a$scan3b<-(9-a$scan3)-4.5
a$scan4b<-(9-a$scan4)-4.5
a$scan5b<-a$scan5-4.5
a$scan6b<-a$scan6-4.5
a$scan7b<-a$scan7-4.5
a$scan8b<-a$scan8-4.5
a$scan9b<-(9-a$scan9)-4.5
a$scan10b<-(9-a$scan10)-4.5

#grouping scans based on their initial masculinity level
s<-length(a$scan1)
group<-c(rep(2,s),rep(-1,s),rep(-2,s),rep(0,s),rep(2,s),
rep(-1,s),rep(1,s),rep(-2,s),rep(0,s),rep(1,s))
group2<-(group>0)-0.5
group2[group == '0']<- 0
group<-as.factor(group)
group2<-as.factor(group2)

id<-as.factor(rep(1:s,10))
attr<-rep(a$attr,10)
horm<-rep(a$horm2,10)

#in a steady relationship (1 = yes, 2=No)
a$steady2<-rep(0,length(a$horm))
a$steady2[a$steady == '1']<- -0.5
a$steady2[a$steady == '2']<- 0.5
steady<-rep(a$steady2,10)

# looking for (1 = steady relationship, 2 = short-term relationship/one night stand, 3 = no relationship
a$looking2<-rep(0,length(a$looking))
a$looking2[a$looking == '1']<- 1
a$looking2[a$looking == '2']<- 2
a$looking2[a$looking == '3']<- 3
looking<-rep(a$looking2,10)
looking<-as.factor(looking)

#age of menarche
men<-rep(a$mearche,10)

#masculinity score
scans<-c(a$scan1b,a$scan2b,a$scan3b,a$scan4b,a$scan5b,a$scan6b,
a$scan7b,a$scan8b,a$scan9b,a$scan10b)

#scan number
nr_scan<-as.factor(c(rep(1,s),rep(2,s),rep(3,s),rep(4,s),rep(5,s),
rep(6,s),rep(7,s),rep(8,s),rep(9,s),rep(10,s)))

#SOI score
a$soi1<-(a$number_sex1+a$number_sex2+a$number_sex3) #behaviour
a$soi2<-(a$SOI1+a$SOI2-a$SOI3) #attitude
a$soi3<-(a$SOI4+a$SOI5+a$SOI6) #desire
a$soi<-(a$soi1+a$soi2+a$soi3)/3 # average

soi1<-rep(a$soi1,10)
soi2<-rep(a$soi2,10)
soi3<-rep(a$soi3,10)
soi<-rep(a$soi,10)


library(lmerTest)
library(MuMIn)
library(nlme)

# null model
lmNull<-lmer(scans ~ 1+ (1|id)+(1|nr_scan))
summary(lmNull)
confint(lmNull)
anova(lmNull)

#second model with groups
lm2<-lmer(scans~ -1 + group+(1|id)+(1|nr_scan:group))
summary(lm2)
confint(lm2)
anova(lm2)

#Full Model
lm1c<-lmer(scans~-1+group+attr+horm+looking+soi+steady+men+(1|id)+(1|nr_scan:group))

summary(lm1c)
confint(lm1c)
anova(lm1c)
