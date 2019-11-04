# Import csv file, call it data
setwd("~/Desktop/PhD Research/Manuscripts/In prep/LTM")
data<-read.csv("LTM_all.csv",header=T)
str(data)
levels(data$Treatment)
levels(data$Time)
levels(data$Yearf)

dSOC<-read.csv("LTM_forSOC3.csv", header=T)
#remove november sampling
data<- subset(data, Month!= '11')

#remove 2014 since it has controls only
data<-subset(data, Year!='2014')

M <- table(data$Yearf, data$Treatment)
M

#separate 'before', 'during', and 'recovery'
during <- subset(data, Time== 'during')
recovery <- subset(data, Time== 'recovery')
before <- subset(data, Time== 'before')
library(lattice)
library(nlme)
library(lsmeans)
library(multcomp)
library(MuMIn)

#Begin by looking at pH and EC and chla by treatment during the experiment
plot(pH~Treatment, during)
plot(conductivity~Treatment, during)
plot(SOC~Treatment, during)
plot(Chl_a~Treatment, during)
plot(SJL~Treatment, during)
plot(SJL~Treatment, recovery)
plot(STD~Yearf,data)
plot(STL~Yearf,data)
plot(ETL~Yearf,data)
plot(SOC~Yearf, data)

### check nematode abundance in recovery
plot(total_live~Treatment,recovery)
## log transform to compare to walter's results
logL<-log(recovery$total_live+1)
recovery<-cbind(logL, recovery)
plot(logL~Treatment,recovery)

#now by pre-experiment
plot(pH~Treatment, before)
plot(SOC~Treatment, before)

#now during recovery
plot(pH~Treatment, recovery)
plot(conductivity~Treatment, recovery)
plot(SOC~Treatment, recovery)
plot(Chl_a~Treatment, recovery)

#analysis of variance for pH by year, treatment, and time 
A1<-aov(pH~ Treatment*Time*Yearf, data=data, na.action=na.exclude)
summary(A1)
qqnorm(residuals(A1))
qqline(residuals(A1))

A2<-aov(pH~ Treatment*Yearf, data=data, na.action=na.exclude)
summary(A2)
qqnorm(residuals(A2))
qqline(residuals(A2))

#data are pretty normally distributed

#because time was significant plot by pre, during, and post-experiment
plot(pH~Time, data)


#Yearf was significant, plot this too
plot(pH~Yearf, data)
#looks like 2004-05 is driving the differences in pH
#possibly a problem with the meter or methods that year, talked to ross about this and a whole order of magnitude drop and recovery would be very significant
#remove year 2004-05 pH data for analysis
pHdata<- subset(during, Yearf!= '2004-05')
plot(pH~Yearf, pHdata)
pHdataall<- subset(data, Yearf!='2004-05')
plot(pH~Yearf,pHdataall)

#mixed effects model for pH
#mixed effects model during
#did not include year as a fixed effect because only two years during experiment measured ph
mph<-lme(pH ~Treatment, random=~1|Block, pHdata, na.action=na.exclude)
summary(mph)
anova(lme(pH ~Treatment, random=~1|Block, pHdata, na.action=na.exclude))
qqnorm(residuals(mph))
qqline(residuals(mph))
m2<-residuals(mph)
shapiro.test(m2)
LS2<-lsmeans(mph, ~Treatment)
contrast(LS2, "pairwise")
#no differences

#test for differences by trt compared to just control
lsmeans(mph, list(trt.vs.ctrl ~ Treatment), adjust="tukey")
#no differences in pH compared to control

#mixed effects model for recovery
mphrec<-lme(pH ~Treatment, random=~1|Block, recovery, na.action=na.exclude)
summary(mphrec)
anova(mphrec)
qqnorm(residuals(mphrec))
qqline(residuals(mphrec))
mrec<-residuals(mphrec)
shapiro.test(mrec) #not normal, log is worse, leave alone

LS9<-lsmeans(mphrec, ~Treatment)
contrast(LS9, "pairwise")
#no differences

#test for differences by trt compared to just control
lsmeans(mphrec, list(trt.vs.ctrl ~ Treatment), adjust="tukey")
#no differences in pH compared to control

#compare pH across all years
allph<-lme(pH ~Treatment*Year, random=~1|Block, pHdataall, na.action=na.exclude)
summary(allph)
plot(lme(pH ~Treatment*Year, random=~1|Block, pHdataall, na.action=na.exclude))
anova(allph)
LS3<-lsmeans(allph, ~Year)
contrast(LS3, "pairwise", adjust="tukey")
plot(pH~Yearf,pHdataall)
aggregate(pH ~ Treatment, pHdataall, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))


#analysis of variance for salinity by year, treatment, and time 
C1<-aov(conductivity~ Treatment*Time*Yearf, data=data, na.action=na.exclude)
summary(C1)
qqnorm(residuals(C1))
qqline(residuals(C1))

C2<-aov(conductivity~ Treatment*Yearf, data=data, na.action=na.exclude)
summary(C2)
qqnorm(residuals(C2))
qqline(residuals(C2))
m3<-residuals(C2)
shapiro.test(m3)
#data are not normally distributed

#because Yearf was significant plot by year
plot(conductivity~Yearf, data)



#looks like conductivity fluctuates by year, highest in recent years
#new meter or methods change??


#mixed effects model for conductivity
#mixed effects model during
mc<-lme(conductivity ~Treatment*Year, random=~1|Block, during, na.action=na.exclude)
summary(mc)
anova(mc)
qqnorm(residuals(mc))
qqline(residuals(mc))
m4<-residuals(mc)
shapiro.test(m4)
#not normal
mc2<-lme(log(conductivity+1)  ~Treatment*Year, random=~1|Block, during, na.action=na.exclude)
summary(mc2)
anova(mc2)
qqnorm(residuals(mc2))
qqline(residuals(mc2))
m5<-residuals(mc2)
shapiro.test(m5)#normal


#test for differences by trt compared to just control
lsmeans(mc2, list(trt.vs.ctrl ~ Treatment), adjust="tukey")
#no differences in pH compared to control

#compare conductivity during and after experiment
allcon<-lme(log(conductivity+1)  ~Time, random=~1|Year/Block, data, na.action=na.exclude)
summary(allcon)
plot(allcon)
anova(allcon)
qqnorm(residuals(allcon))
qqline(residuals(allcon))
m6<-residuals(allcon)
shapiro.test(m6)# normal
LS7<-lsmeans(allcon, ~Time)
contrast(LS7, "pairwise")
#no difference in conductivity between during and recovery
plot(conductivity~Time,data)
aggregate(conductivity ~ Time, data, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))

#is it the chambers that cause differences in conductivity?
Tcon<-lme(log(conductivity+1)  ~temp_manip, random=~1|Year/Block, during, na.action=na.exclude)
summary(Tcon)
plot(Tcon)
anova(Tcon)
LS8<-lsmeans(Tcon, ~temp_manip)
contrast(LS8, "pairwise")
#no difference

#does moisture manipulation cause differences in conductivity?
Wcon<-lme(log(conductivity+1)  ~moist_manip, random=~1|Yearf/Block, during, na.action=na.exclude)
summary(Wcon)
plot(Wcon)
anova(Wcon)
r.squaredGLMM(Wcon) #2% explained by fixed effects, 62% by whole model
LS8<-lsmeans(Wcon, ~moist_manip)
contrast(LS8, "pairwise")
#significant difference!
plot(conductivity~moist_manip,during)
aggregate(conductivity ~ moist_manip, during, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))

#is this effect persistent during the legacy phase?
Wcon2<-lme(log(conductivity+1)  ~moist_manip, random=~1|Yearf/Block, recovery, na.action=na.exclude)
summary(Wcon2)
plot(Wcon2)
anova(Wcon2)
r.squaredGLMM(Wcon2)#2% by fixed, 76% by whole model
LS8<-lsmeans(Wcon2, ~moist_manip)
contrast(LS8, "pairwise")
#significant difference!
plot(conductivity~moist_manip,recovery)
aggregate(conductivity ~ moist_manip, recovery, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))

#mixed effects model for chla
#mixed effects model during
plot(Chl_a~Yearf, data)
mchl<-lme(Chl_a ~Treatment*Year, random=~1|Block, during, na.action=na.exclude)
summary(mchl)
anova(mchl) #significant effect of year
qqnorm(residuals(mchl))
qqline(residuals(mchl))
m1<-residuals(mchl)
shapiro.test(m1)
#not normal
mchl2<-lme(log(Chl_a+1)  ~Treatment*Year, random=~1|Block, during, na.action=na.exclude)
summary(mchl2)
anova(mchl2) #year significant (declines over time)
r.squaredGLMM(mchl2) #5% of variation explained by fixed, 5% by whole model
qqnorm(residuals(mchl2))
qqline(residuals(mchl2))
m8<-residuals(mchl2)
shapiro.test(m8)
#test for differences by trt compared to just control
lsmeans(mchl2, list(trt.vs.ctrl ~ Treatment), adjust="tukey")
#no differences in chla compared to control
plot(Chl_a~Yearf, during)

#mixed effects model recovery
mchlr<-lme(Chl_a ~Treatment, random=~1|Year/Block, recovery, na.action=na.exclude)
summary(mchlr)
anova(mchlr) 
r.squaredGLMM(mchlr)#3% of variation explained by fixed, 31% by whole model
qqnorm(residuals(mchlr))
qqline(residuals(mchlr))
shapiro.test(residuals(mchlr))
#not normal
mchlr2<-lme(sqrt(Chl_a)  ~Treatment, random=~1|Year/Block, recovery, na.action=na.exclude)
summary(mchlr2)
anova(mchlr2) 
r.squaredGLMM(mchlr2) #3% of variation explained by fixed, 31% by whole model
qqnorm(residuals(mchlr2))
qqline(residuals(mchlr2))
shapiro.test(residuals(mchl2))#worse, do not transform
#test for differences by trt compared to just control
lsmeans(mchl2, list(trt.vs.ctrl ~ Treatment), adjust="tukey")
#no differences in chla compared to control
plot(Chl_a~Yearf, during)

#compare chla during and after experiment
allchl<-lme(sqrt(Chl_a)  ~Time, random=~1|Year/Block, data, na.action=na.exclude)
summary(allchl)
plot(allchl)
anova(allchl)
qqnorm(residuals(allchl))
qqline(residuals(allchl))
m9<-residuals(allchl)
shapiro.test(m9)
# not normal
LS10<-lsmeans(allchl, ~Time)
contrast(LS10, "pairwise")
#no difference in chla between during and recovery
plot(Chl_a~Time,data)
aggregate(Chl_a ~ Time, data, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))

#effect of yearf on chla during experiment
mchly<-lme(Chl_a ~Yearf, random=~1|Block, during, na.action=na.exclude)
summary(mchly)
anova(mchly)
qqnorm(residuals(mchly))
qqline(residuals(mchly))
m10<-residuals(mchly)
shapiro.test(m10)
#not normal
mchly2<-lme(log(Chl_a+1)  ~Yearf, random=~1|Block, during, na.action=na.exclude)
summary(mchly2)
anova(mchly2)
qqnorm(residuals(mchly2))
qqline(residuals(mchly2))
m11<-residuals(mchly2)
shapiro.test(m11)
LS15<-lsmeans(mchly2, ~Yearf)
contrast(LS15, "pairwise")

#effect of yearf on chla in recovery 
mchlyr2<-lme(Chl_a ~Yearf, random=~1|Block, recovery, na.action=na.exclude)
summary(mchlyr2)
anova(mchlyr2)
qqnorm(residuals(mchlyr2))
qqline(residuals(mchlyr2))
m12<-residuals(mchlyr2)
shapiro.test(m12)
#not normal
mchlyr2<-lme(log(Chl_a+1)  ~Yearf, random=~1|Block, recovery, na.action=na.exclude)
summary(mchlyr2)
anova(mchlyr2)
qqnorm(residuals(mchlyr2))
qqline(residuals(mchlyr2))
m13<-residuals(mchlyr2)
shapiro.test(m13)
#still not normal
LS15<-lsmeans(mchlyr2, ~Yearf)
contrast(LS15, "pairwise")
#summarize chla across years
aggregate(Chl_a ~ Yearf, data, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))


#do the chambers that cause differences in chla?
Tchl<-lme(log(Chl_a+1)  ~temp_manip, random=~1|Yearf/Block, during, na.action=na.exclude)
summary(Tchl)
plot(Tchl)
anova(Tchl)
LS11<-lsmeans(Tchl, ~temp_manip)
contrast(LS11, "pairwise")
#no difference

#any lag effect of chambers?
Tchl2<-lme(log(Chl_a+1)  ~temp_manip, random=~1|Yearf/Block, recovery, na.action=na.exclude)
summary(Tchl2)
plot(Tchl2)
anova(Tchl2)
LS14<-lsmeans(Tchl2, ~temp_manip)
contrast(LS14, "pairwise")
#no difference

#does moisture manipulation cause differences in chla?
Wchl<-lme(log(Chl_a+1)  ~moist_manip, random=~1|Yearf/Block, during, na.action=na.exclude)
summary(Wchl)
plot(Wchl)
anova(Wchl)
LS12<-lsmeans(Wchl, ~moist_manip)
contrast(LS12, "pairwise")
#no difference
plot(Chl_a~moist_manip,during)
aggregate(Chl_a ~ moist_manip, during, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))

#is there a lag effect during the legacy phase?
Wchl2<-lme(log(Chl_a+1)  ~moist_manip, random=~1|Yearf/Block, recovery, na.action=na.exclude)
summary(Wchl2)
plot(Wchl2)
anova(Wchl2)
LS13<-lsmeans(Wchl2, ~moist_manip)
contrast(LS13, "pairwise")
#no difference
plot(Chl_a~moist_manip,recovery)
aggregate(Chl_a ~ moist_manip, recovery, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))


#does air or soil temperature affect chla?
temp<-read.csv("LTM temperatures.csv",header=T)
summer<-subset(temp, Month!= '11')
summer<-subset(summer, Month!= '10')
summer<-subset(summer, Month!= '9')
summer<-subset(summer, Month!= '8')
summer<-subset(summer, Month!= '7')
summer<-subset(summer, Month!= '6')
summer<-subset(summer, Month!= '5')
summer<-subset(summer, Month!= '4')
summer<-subset(summer, Month!= '3')
summer<-subset(summer, Month!= '2')
write.csv(summer, file= "LTM_summer_temp.csv")
aggregate(Open_surface ~ Month*Year, summer, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
summert<-read.csv("LTM_summer_temp.csv")
airT<-aggregate(Open_surface ~ Yearf, summert, function(x) c(M = mean(x)))
soilT<-aggregate(Open_5cm ~ Yearf, summert, function(x) c(M = mean(x)))
library(plyr)
data <- join(data, airT, by="Yearf")
data <- join(data,soilT, by="Yearf")
during <- join(during, airT, by="Yearf")
during <- join(during,soilT, by="Yearf")
dSOC <- join(dSOC, airT, by= "Yearf")
dSOC <- join(dSOC, soilT, by= "Yearf")

Chla<-aggregate(Chl_a ~ Yearf, data, function(x) c(M = mean(x)))
Chla
chla_t<-join(Chla,airT,by="Yearf")

ct<-aov(Chl_a ~ Open_surface*moisture, data=during, na.action=na.exclude)
summary(ct)
cs<-aov(Chl_a~ Open_5cm*moisture, during, na.action=na.exclude)
summary(cs)
plot(data$Open_surface, data$Chl_a)
plot(data$Open_5cm, sqrt(data$Chl_a))

data$tempcat = NA
data$tempcat [data$Open_surface > 1] = '0 to 3.5'
data$tempcat [data$Open_surface > 3.5] = '3.5 to 5'
data$tempcat [data$Open_surface > 5] = '5 to 6'
data$tempcat [data$Open_surface > 6] = '6+'

during$tempcat = NA
during$tempcat [during$Open_surface > 0] = '0 to 4'
during$tempcat [during$Open_surface > 4] = '4 to 6'
during$tempcat [during$Open_surface > 6] = '6+'

library(ggplot2)
ggplot(during, aes(x=moisture, y=Chl_a, group=tempcat, color=tempcat))+
  geom_point()+
  geom_smooth(method=lm)+
  labs(x = 'Soil Moisture', y = 'Chlorophyll a')

ggplot(during, aes(x=moisture, y=Chl_a, group=tempcat, color=tempcat))+
  geom_smooth(method=lm, se=FALSE)+
  labs(x = 'Soil Moisture', y = 'Chlorophyll a')

ggplot(dSOC, aes(x=Open_surface, y=SOC.mg.g, group= Treatment, color=Treatment))+
geom_point()+
  geom_smooth(method=lm, se=FALSE)+
  labs(x = 'Average summer air temperature', y = 'Soil Organic Carbon')


#AOV for SOC, use during data because there is currently no data during recover
B1<-aov(SOC~ Treatment*Yearf, data=during, na.action=na.exclude)
summary(B1)
qqnorm(residuals(B1))
qqline(residuals(B1))
shapiro.test(residuals(B1))
#not normally distributed,
#square root transformation 
sqrtSOC <- sqrt(during$SOC)
during<- cbind(sqrtSOC, during)
shapiro.test(sqrtSOC)
#p<0.05 proceed
B2<-aov(sqrtSOC~moist_manip, data=during, na.action=na.exclude)
summary(B2)
qqnorm(residuals(B2))
qqline(residuals(B2))
z<-residuals(B2)
shapiro.test(z)
lsmeans(B2, list(trt.vs.ctrl ~ moist_manip), adjust="tukey")
library(ggplot2)

##Add analysis of C in recovery
R1<-aov(sqrt(SOC)~ moisture, data=recovery, na.action=na.exclude)
summary(R1)
qqnorm(residuals(R1))
qqline(residuals(R1))
shapiro.test(residuals(R1))
lsmeans(R1, list(trt.vs.ctrl ~ moisture), adjust="tukey")



SOC_mean<-aggregate(SOC.mg.g ~ Yearf*Treatment, dSOC, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))  
write.csv(SOC_mean, 'SOC_mean.csv')
SOC_mean<-read.csv('SOC_mean.csv')
ggplot(d=SOC_mean, aes(x=Yearf, y=SOC.mg.g.M, color=Treatment, group=Treatment)) +
geom_errorbar(aes(ymin=SOC.mg.g.M-SOC.mg.g.SE, ymax=SOC.mg.g.M+SOC.mg.g.SE), width=.1)+
  geom_point() +
  geom_line()
  
SE<-aggregate(percent.of.added.C.remaining ~ Yearf*Treatment, dSOC, function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))  
write.csv(SE, 'SOC_SE.csv')
SE<-read.csv('SOC_SE.csv')
ggplot(d=SE, aes(x=Yearf, y=percent.of.added.C.remaining.M, color=Treatment, group=Treatment)) +
  geom_errorbar(aes(ymin=percent.of.added.C.remaining.M-percent.of.added.C.remaining.SE, ymax=percent.of.added.C.remaining.M+percent.of.added.C.remaining.SE), width=.1)+
  geom_point() +
geom_line()

Z1<-aov(SOC.mg.g~Treatment*Yearf, data=dSOC, na.action=na.exclude)
summary(Z1)
qqnorm(residuals(Z1))
zz<- residuals(Z1)
shapiro.test(zz)
lsmeans(Z1, list(trt.vs.ctrl ~ Treatment), adjust="tukey")

plot(sqrtSOC~Treatment, during)
plot(sqrtSOC~Yearf, during)
plot(SOC~Treatment, during)
plot(SOC~Yearf, during)
plot(SOC~Open_surface, during)

#AOV for conductivity
C1<-aov(conductivity~ Treatment*Yearf, data=data, na.action=na.exclude)
summary(C1)
qqnorm(residuals(C1))
qqline(residuals(C1))
#not normally distributed
# data are not normally distributed, so log transform
LogEC<- log(data$conductivity)
data<- cbind(LogEC, data)
shapiro.test(LogEC)

#p is >0.05, no need to transform

sqrtEC <- sqrt(data$conductivity)
data <- cbind (sqrtEC, data)
shapiro.test(sqrtEC)
#p is  <0.05, so not normal, check qqplot just in case
qqnorm(sqrtEC)
qqline(sqrtEC)
#better, not perfect. Try AOV again with sqrtEC values

#aov for during expt
C2<-aov(conductivity~Treatment*Yearf, data=during, na.action=na.exclude)
summary(C2)
qqnorm(residuals(C2))
qqline(residuals(C2))
lsmeans(C2, list(trt.vs.ctrl ~ Treatment), adjust="tukey")
plot(conductivity~Treatment, during)

#aov for overall but with sqrt transformation
C3<-aov(sqrtEC~ Treatment*Yearf, data=data, na.action=na.exclude)
summary(C3)
qqnorm(residuals(C3))
qqline(residuals(C3))
lsmeans(C3, list(trt.vs.ctrl ~ Treatment), adjust="tukey")
#lsmeans won't run, lets try again without 1997 
data4 <-subset(data, Year!='1997')

C4<-aov(sqrtEC~ Treatment*Yearf, data=data, na.action=na.exclude)
summary(C4)
qqnorm(residuals(C4))
qqline(residuals(C4))
lsmeans(C4, list(trt.vs.ctrl ~ Treatment), adjust="tukey")
#still won't give me anything... WHY?!

C5<- aov(conductivity~Treatment*Yearf, data=recovery, na.action=na.exclude)
summary(C5)
qqnorm(residuals(C5))
qqline(residuals(C5))

plot(sqrtEC~Yearf, data)
plot(sqrtEC~Treatment, data)
plot(conductivity~Yearf, during)
plot(conductivity~Treatment, during)

lsmeans(C4, list(trt.vs.ctrl ~ Treatment), adjust="tukey")

ggplot(data, aes(x=moisture, y=conductivity, color=Yearf)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)

#EC increasing in later years. But no early expt data to compare to

#AOV for chlorophyll
D1<-aov(Chl_a~ Treatment*Yearf, data=data, na.action=na.exclude)
summary(D1)
qqnorm(residuals(D1))
qqline(residuals(D1))
#data not normally distributed, try sqrt transformation
sqrtchla <- sqrt(data$Chl_a)
data <- cbind (sqrtchla, data)
shapiro.test(sqrtchla)
#p<0.05 so normality now ok, try AOV again
D2<- aov(sqrtchla~Treatment*Yearf, data=data, na.action=na.exclude)
summary(D2)
qqnorm(residuals(D2))
qqline(residuals(D2))

plot(Chl_a~Yearf,data)
plot(sqrtchla~Time,data)
plot(sqrtchla~Yearf,data)

D2<-aov(Chl_a~Treatment*Yearf, data=recovery, na.action=na.exclude)
summary(D2)

D3<-aov(Chl_a~Treatment*Yearf, data=during, na.action=na.exclude)
summary(D3)
qqnorm(residuals(D3))
qqline(residuals(D3))
lsmeans(D3, list(trt.vs.ctrl ~ Treatment), adjust="tukey")

finalyr <- subset(recovery, Year!=2006)
finalyr <- subset (finalyr, Year!=2011)
finalyr1 <- subset (finalyr, Year!=2014)
finalyr <- subset (finalyr1, Year!=2008)
yr2008 <-subset(finalyr1, Year!= 2015)

B3 <- aov(Chl_a~Treatment, data=finalyr, na.action=na.exclude)
summary(B3)
qqnorm(residuals(B3))
qqline(residuals(B3))
sqrtCHLA <- sqrt(finalyr$Chl_a)
finalyr <- cbind (sqrtCHLA, finalyr)
shapiro.test(sqrtCHLA)
B4 <- aov(sqrtCHLA~Treatment*Block, data=finalyr, na.action=na.exclude)
summary(B4)
lsmeans(B4, list(trt.vs.ctrl ~ Treatment), adjust="tukey")
plot(sqrtCHLA~Treatment, finalyr)

#check 2008 for any post treatment effects on chla
C1 <- aov(Chl_a~Treatment, data=yr2008, na.action=na.exclude)
summary(C1)
lsmeans(C1, list(trt.vs.ctrl ~ Treatment), adjust="tukey")

#AOV for soil moisture
A7 <- aov(moisture~Treatment*Yearf, data=data, na.action=na.exclude)
summary(A7)
qqnorm(residuals(A7))
qqline(residuals(A7))
#data not normal, try sqrt transformation
sqrtmoist <- sqrt(data$moisture)
data <- cbind (sqrtmoist, data)
shapiro.test(sqrtmoist)
#p<0.05 ok, try this

A8 <- aov(sqrtmoist~Treatment*Yearf, data=data, na.action=na.exclude)
summary(A8)
qqnorm(residuals(A8))
qqline(residuals(A8))

plot(sqrtmoist~Treatment,data)
plot(sqrtmoist~Yearf,data)

### test chla variability by comparing to P3
d<-read.csv("P3_chla.csv",header=T)
#make sure data are numeric
d<-transform(d,ug.chla.g.soil=as.character(ug.chla.g.soil))
d<-transform(d,ug.chla.g.soil=as.numeric(ug.chla.g.soil))
#drop data <0 (no such thing as a negative chla conc)
d<- d[d[,13] >0,]
d$Experiment= factor(d$Experiment, levels=c("LTM - 2016", "P3-2013", "P3- 2016", "P3- 2017", "ET SSLH - 2013","ET SSLH 2016", "ET SSLH - 2017", "BEE F6 2016", "BEE F6- 2017","BEE Bonney 2016", "BEE Bonney -2017", "BEE SSLH 2016", "BEE SSLH - 2017"))
plot(d$ug.chla.g.soil~d$Experiment,  ylim=c(0, 0.5), col= c("orange", "white", "khaki","white", "white","khaki", "white","khaki", "white","khaki", "white","khaki", "white"), xlab="", ylab="ug Chla/g soil", las=2, cex.axis=0.4, tck=-.01, )
plot(ug.chla.g.soil~Experiment, d, main="Chlorophyll a comparison",  ylim=c(0, 0.1))
