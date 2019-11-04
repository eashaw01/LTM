#the purpose of this script is to run mixed models assessing treatment and year effects on EC, pH, and chla
#block as a fixed effect
library(lattice)
library(nlme)
library(lsmeans)
library(multcomp)
library(MuMIn)
library(tidyverse)
library(ggplot2)

# Import csv file, call it data
#setwd("~/Desktop/PhD Research/Manuscripts/In prep/LTM")
#data<-read.csv("LTM_all2.csv",header=T)
setwd("~/Desktop")
data<-read.csv("LTM_all.csv",header=T)
str(data)
levels(data$Treatment)
levels(data$Time)
levels(data$Yearf)
data$Block<-as.factor(data$Block)

#remove november sampling and 2014 since it has controls only
november<-data %>%subset(Month=='11')
data14<-data %>%subset(Month!='11') #keep 2014 for controls
data<- data %>% subset(Month!= '11') %>% subset(Year!='2014') #remove 2014 for treatment analyses

M <- table(data$Yearf, data$Treatment)
M

#separate 'before', 'during', and 'recovery'
during <- subset(data, Time== 'during')
recovery <- subset(data, Time== 'recovery')
before <- subset(data, Time== 'before')
finalyr <- subset(data, Yearf=="2015-16")

#summarize soil chem variables for table
ChemTBL1 <- before %>% summarize(ph.m=mean(pH), ph.se=sd(pH)/sqrt(length(pH)),
                                EC.m=mean(conductivity), EC.se=sd(conductivity)/sqrt(length(conductivity)),
                                SOC.m=mean(SOC), SOC.se=sd(SOC)/sqrt(length(SOC)),
                                TN.m=mean(TN), TN.se=sd(TN)/sqrt(length(TN)),
                                TC.m=mean(TC), TC.se=sd(TC)/sqrt(length(TC)))
ChemTBL1_chla <- november %>% filter(!is.na(Chl_a)) %>% summarize(Chla.m=mean(Chl_a), Chla.se=sd(Chl_a)/sqrt(length(Chl_a)))
ChemTBL2_ph <- during %>% filter(!is.na(pH)) %>% filter(Yearf!="2004-05")%>% filter(Treatment=="C")%>%summarize(ph.m=mean(pH), ph.se=sd(pH)/sqrt(length(pH)))
ChemTBL2_EC <- during %>% filter(!is.na(conductivity)) %>% filter(Treatment=="C")%>%summarize(EC.m=mean(conductivity), EC.se=sd(conductivity)/sqrt(length(conductivity)))
ChemTBL2_SOC <- during %>% filter(!is.na(SOC2)) %>% filter(Treatment=="C")%>%summarize(SOC.m=mean(SOC2), SOC.se=sd(SOC2)/sqrt(length(SOC2)))
ChemTBL2_TN <- during %>% filter(!is.na(TN)) %>% filter(Treatment=="C")%>%summarize(TN.m=mean(TN), TN.se=sd(TN)/sqrt(length(TN)))
ChemTBL2_TC <- during %>% filter(!is.na(TC)) %>% filter(Treatment=="C")%>%summarize(TC.m=mean(TC), TC.se=sd(TC)/sqrt(length(TC)))
ChemTBL2_chla <- during %>% filter(!is.na(Chl_a)) %>% filter(Treatment=="C")%>%summarize(Chla.m=mean(Chl_a), Chla.se=sd(Chl_a)/sqrt(length(Chl_a)))

ChemTBL3_ph <- finalyr %>% filter(!is.na(pH)) %>%  filter(Treatment=="C")%>%summarize(ph.m=mean(pH), ph.se=sd(pH)/sqrt(length(pH)))
ChemTBL3_EC <- finalyr %>% filter(!is.na(conductivity)) %>% filter(Treatment=="C")%>%summarize(EC.m=mean(conductivity), EC.se=sd(conductivity)/sqrt(length(conductivity)))
ChemTBL3_SOC <- finalyr %>% filter(!is.na(SOC2)) %>% filter(Treatment=="C")%>%summarize(SOC.m=mean(SOC2), SOC.se=sd(SOC2)/sqrt(length(SOC2)))
ChemTBL3_TN <- finalyr %>% filter(!is.na(TN)) %>% filter(Treatment=="C")%>%summarize(TN.m=mean(TN), TN.se=sd(TN)/sqrt(length(TN)))
ChemTBL3_TC <- finalyr %>% filter(!is.na(TC)) %>% filter(Treatment=="C")%>%summarize(TC.m=mean(TC), TC.se=sd(TC)/sqrt(length(TC)))
ChemTBL3_chla <- finalyr %>% filter(!is.na(Chl_a)) %>% filter(Treatment=="C")%>%summarize(Chla.m=mean(Chl_a), Chla.se=sd(Chl_a)/sqrt(length(Chl_a)))

#mixed effects model for pH
#mixed effects model during
plot(pH~Yearf, data) #2004-05 has much lower pH, discussions with Ross suggest this likely isn't possible
pHdata<- during %>% subset(Yearf!='2004-05') #remove 2004-05 from analysis
mph<-lme(pH ~Treatment, random=~1|Year/Block, pHdata, na.action=na.exclude)
summary(mph)
anova(mph)
r.squaredGLMM(mph)#2% by fixed, 74% by whole model
qqnorm(residuals(mph))
qqline(residuals(mph))
m2<-residuals(mph)
shapiro.test(m2) #normal
LS2<-lsmeans(mph, ~Treatment)
contrast(LS2, "pairwise")
#no differences

#mixed effects model for pH
#mixed effects model recovery
plot(pH~Yearf, data)
mphr<-lme(log(pH+1) ~Treatment*Year, random=~1|Block, recovery, na.action=na.exclude)
summary(mphr)
anova(mphr)
r.squaredGLMM(mphr)#4% by fixed, 32% by whole model
qqnorm(residuals(mphr))
qqline(residuals(mphr))
shapiro.test(residuals(mphr)) #notnormal
LSr<-lsmeans(mphr, ~Treatment)
contrast(LSr, "pairwise")
#no differences

#does chamber or water addition have effects on ph?
plot(pH~moist_manip, data) #doesn't look like it
plot(pH~temp_manip, data)#doesn't look like it
mphtw<-lme(pH ~moist_manip+temp_manip, random=~1|Year/Block, during, na.action=na.exclude)
summary(mphtw)
anova(mphtw)
r.squaredGLMM(mphtw)#4% by fixed, 32% by whole model
qqnorm(residuals(mphtw))
qqline(residuals(mphtw))
shapiro.test(residuals(mphtw)) #normal
#no differences

#does chamber or water addition have effects on ph?
plot(pH~moist_manip, data) #doesn't look like it
plot(pH~temp_manip, data)#doesn't look like it
mphtwr<-lme(pH ~moist_manip+temp_manip, random=~1|Year/Block, recovery, na.action=na.exclude)
summary(mphtwr)
anova(mphtwr)
r.squaredGLMM(mphtwr)#4% by fixed, 32% by whole model
qqnorm(residuals(mphtwr))
qqline(residuals(mphtwr))
shapiro.test(residuals(mphtwr)) #not normal
#no differences

#mixed effects model for EC
#mixed effects model during
plot(conductivity~Yearf, data)
mEC<-lme(log(conductivity+1) ~Treatment*Year, random=~1|Block, during, na.action=na.exclude)
summary(mEC)
anova(mEC) #significant effect of year
r.squaredGLMM(mEC)#17% by fixed, 27% by whole model
qqnorm(residuals(mEC))
qqline(residuals(mEC))
shapiro.test(residuals(mEC)) #notnormal
LSec<-lsmeans(mEC, ~Treatment)
contrast(LSec, "pairwise")
#no differences in treatment


#mixed effects EC model recovery
plot(conductivity~Treatment, recovery)
mECr<-lme(log(conductivity+1) ~Treatment*Year, random=~1|Block, recovery, na.action=na.exclude)
summary(mECr)
anova(mECr)
r.squaredGLMM(mECr)#48% by fixed, 67% by whole model
qqnorm(residuals(mECr))
qqline(residuals(mECr))
shapiro.test(residuals(mECr)) #normal
LSec<-lsmeans(mECr, ~Treatment)
contrast(LSec, "pairwise", adjust="tukey")
lsmeans(mECr, list(trt.vs.ctrl ~ Treatment), adjust="tukey")
#no differences


#does chamber or water addition have effects on ec?
plot(conductivity~moist_manip, during) #maybe
plot(conductivity~temp_manip, during)#doesn't look like it
mectw<-lme(log(conductivity+1) ~moist_manip+temp_manip, random=~1|Year/Block, during, na.action=na.exclude)
summary(mectw)
anova(mectw)
r.squaredGLMM(mectw)#4% by fixed, 60% by whole model
qqnorm(residuals(mectw))
qqline(residuals(mectw))
shapiro.test(residuals(mectw)) #normal
plot(mectw)
LSecw1<-lsmeans(mectw, ~moist_manip)
contrast(LSecw1, "pairwise")
#water manipulation effects 
#plot it
EC_graph <- data %>% subset(Time!='before') %>% filter(conductivity>0)%>%
  group_by(Time, moist_manip) %>%
  summarize(meanec=mean(conductivity), seec=sd(conductivity)/sqrt(length(conductivity)))

ggplot(EC_graph, aes(x=moist_manip, y=meanec, color=moist_manip, shape=moist_manip))+
  geom_point()+
  facet_wrap(~Time)+
  theme_bw()+
  scale_shape_manual(name="Water Addition", values=c(16, 25))+
  scale_color_manual(name="Water Addition", values=c("black","blue"))+
  geom_errorbar(aes(ymin=meanec-seec, ymax=meanec+seec), width=.2,position=position_dodge(0.05))+
  labs(x="Moisture Treatment", y="Conductivity uS/cm")
  annotate("text", x = c(1.5), y=(50), label = c("*","*"), fontface = 3,size=12)

#does chamber or water addition have effects on EC?
plot(conductivity~moist_manip, recovery) #maybe
plot(conductivity~temp_manip, recovery) #doesn't look like it
mectwr<-lme(log(conductivity+1) ~moist_manip+temp_manip, random=~1|Year/Block, recovery, na.action=na.exclude)
summary(mectwr)
anova(mectwr)
r.squaredGLMM(mectwr)#2% by fixed, 67% by whole model
qqnorm(residuals(mectwr))
qqline(residuals(mectwr))
shapiro.test(residuals(mectwr)) # normal
LSecw<-lsmeans(mectwr, ~moist_manip)
contrast(LSecw, "pairwise", adjust="tukey")
#water manipulation effects 

#mixed effects model for chlorophyll a
#mixed effects model during
plot(Chl_a~Yearf, data)
mChla<-lme(sqrt(Chl_a) ~Treatment*Year, random=~1|Block, during, na.action=na.exclude)
summary(mChla)
anova(mChla)#year is significant
r.squaredGLMM(mChla)#5% by fixed, 5% by whole model
qqnorm(residuals(mChla))
qqline(residuals(mChla))
shapiro.test(residuals(mChla)) #close to normal
lsmeans(mChla, list(trt.vs.ctrl ~ Treatment), adjust="tukey")
#no treatment differences

#plot chla over time
Chla_graph <- data %>% filter(Chl_a != 'NA')%>%
  group_by(Yearf, Time) %>% 
  summarize(meanChla=mean(Chl_a), seChla=sd(Chl_a)/sqrt(length(Chl_a)))

ggplot(Chla_graph, aes(x=Yearf, y=meanChla, group=Time))+
  geom_point()+
  geom_rect(data=NULL,aes(xmin=0,xmax=5,ymin=-Inf,ymax=Inf),
            fill="gray")+
  geom_vline(xintercept =  4, linetype="dotted")+
  geom_point()+
  #geom_smooth(method = "lm", se=F)+
  theme_bw()+
  geom_errorbar(aes(ymin=meanChla-seChla, ymax=meanChla+seChla), width=.2,position=position_dodge(0.05))+
  labs(x="Year", y="Chlorophyll a (ug/g)")+
  annotate("text", x = c(3.92), y=(0.02), label = c("Flood year"), fontface = 3,size=4)+
  annotate("text", x = c(2), y=(0.025), label = c("Treatment"), fontface = 2,size=4)+
  annotate("text", x = c(6.5), y=(0.025), label = c("Post-Treatment"), fontface = 2,size=4)

#mixed effects chla model recovery
plot(Chl_a~Yearf, recovery)
mChlar<-lme(sqrt(Chl_a) ~Treatment*Year, random=~1|Block, recovery, na.action=na.exclude)
summary(mChlar)
anova(mChlar)#year is significant
r.squaredGLMM(mChlar)#11% by fixed, 13% by whole model
qqnorm(residuals(mChlar))
qqline(residuals(mChlar))
shapiro.test(residuals(mChlar)) #closer to normal
lsmeans(mChlar, list(trt.vs.ctrl ~ Treatment), adjust="tukey")
#no differences

#does chamber or water addition have effects on chla?
plot(Chl_a~moist_manip, during) #maybe
plot(Chl_a~temp_manip, during)#maybe
mchltw<-lme(sqrt(Chl_a) ~moist_manip+temp_manip, random=~1|Year/Block, during, na.action=na.exclude)
summary(mchltw)
anova(mchltw)
r.squaredGLMM(mchltw)#1% by fixed, 46% by whole model
qqnorm(residuals(mchltw))
qqline(residuals(mchltw))
shapiro.test(residuals(mchltw)) #not normal, but best transformation
plot(mchltw)
#no effects 

#does chamber or water addition have effects on ph?
plot(Chl_a~moist_manip, recovery) #maybe
plot(Chl_a~temp_manip, recovery) #doesn't look like it
mchltwr<-lme(sqrt(Chl_a) ~moist_manip+temp_manip, random=~1|Year/Block, recovery, na.action=na.exclude)
summary(mchltwr)
anova(mchltwr)
r.squaredGLMM(mchltwr)#2% by fixed, 67% by whole model
qqnorm(residuals(mchltwr))
qqline(residuals(mchltwr))
shapiro.test(residuals(mchltwr)) # better
plot(mchltwr)
#no effects 


#what about the change in soil C from initial (1993)?
#first create column of initial values
initial<- data %>% filter(Time=="before") %>% dplyr::select(Treatment, Block, SOC) %>% rename(SOCi=SOC)
data1<- data %>% filter(Time!="before")
idata<- data1 %>% full_join(initial, by=c("Treatment", "Block")) %>% mutate(dSOC=SOC-SOCi) %>% mutate(dSOCpc=(SOC-SOCi)/SOCi*100)
SOC_mean3 <- idata %>% group_by(Year, Treatment) %>%filter(!is.na(dSOC))%>%
  summarize(dSOC.M=mean(dSOC), dSOC.SE=sd(dSOC)/sqrt(length(dSOC)))
idata_2001<-idata %>% filter(Year=="2001")
idata_2015<-idata %>% filter(Year=="2015")#use this for recovery
idata_during<-idata %>% filter(Time=="during")

##NOW USE NEW VALUES FROM MATT!!!!!##
#what about the change in soil C from initial (1993)?
#first create column of initial values
initial<- data %>% filter(Time=="before") %>% dplyr::select(Treatment, Block, SOC2) %>% rename(SOCi=SOC2)
data1<- data %>% filter(Time!="before")
idata_2<- data1 %>% full_join(initial, by=c("Treatment", "Block")) %>% mutate(dSOC=SOC2-SOCi, na.rm=T) %>% mutate(dSOCpc=(SOC2-SOCi)/SOCi*100, na.rm=T)
SOC_mean4 <- idata_2 %>% group_by(Year, Treatment) %>% filter(!is.na(dSOC))%>%
  summarize(dSOC.M=mean(dSOC), dSOC.SE=sd(dSOC)/sqrt(length(dSOC)))
idata_2001_2<-idata_2 %>% filter(Year=="2001")
idata_2015_2<-idata_2 %>% filter(Year=="2015")#use this for recovery
idata_during_2<-idata_2 %>% filter(Time=="during")

#summarize percent change
SOC_mean_pc <- idata_2 %>% group_by(Year, Treatment) %>% filter(!is.na(dSOCpc))%>%
  summarize(dSOC.M=mean(dSOCpc), dSOC.SE=sd(dSOCpc)/sqrt(length(dSOCpc)))

#USING NEW VALUES FROM MATT check effects of treatments on change in SOC during experiment
idata_during_2<-filter(idata_during_2, dSOC!="NA")
disoc<-lme(dSOC ~Treatment, random=~1|Year/Block, idata_during_2, na.action=na.exclude)
summary(disoc)
plot(disoc)
anova(disoc)#treatment highly sign
qqnorm(residuals(disoc))
qqline(residuals(disoc))
shapiro.test(residuals(disoc))#close to normal
r.squaredGLMM(disoc)#25% of variation explained by fixed effects, 55% by whole model
diLS2<-lsmeans(disoc, ~Treatment)
contrast(diLS2, "pairwise", adjust="tukey")
lsmeans(disoc, list(trt.vs.ctrl ~ Treatment), adjust="tukey")
#WM is different from control

##NOW USE NEW VALUES FROM MATT!!!!!##
#check effects of treatments on change in SOC during experiment
disoc_2<-lme(dSOC  ~Treatment*Year, random=~1|Block, idata_during_2, na.action=na.exclude)
summary(disoc_2)
plot(disoc_2)
anova(disoc_2)#treatment highly sign
qqnorm(residuals(disoc_2))
qqline(residuals(disoc_2))
shapiro.test(residuals(disoc_2))#normal
r.squaredGLMM(disoc_2)#23% of variation explained by fixed effects, 42% by whole model
diLS2_2<-lsmeans(disoc_2, ~Treatment)
contrast(diLS2_2, "pairwise", adjust="tukey")
lsmeans(disoc_2, list(trt.vs.ctrl ~ Treatment), adjust="tukey")
#WM  diff from control
plot(dSOC~Treatment, idata_during_2)

#do water or temp treatments have effect when pooled?
msoctw<-lme(dSOC ~moist_manip+temp_manip, random=~1|Year/Block, idata_during_2, na.action=na.exclude)
summary(msoctw)
anova(msoctw)
r.squaredGLMM(msoctw)#1% by fixed, 46% by whole model
qqnorm(residuals(msoctw))
qqline(residuals(msoctw))
shapiro.test(residuals(msoctw)) #not normal, but best transformation
plot(msoctw)
#water addition has sign effects

ggplot()

disoc_2t<-lme(log(dSOC+1)  ~Treatment, random=~1|Year/Block, idata_during_2, na.action=na.exclude)
summary(disoc_2t)
plot(disoc_2t)
anova(disoc_2t)#treatment, year, & interaction highly sign
qqnorm(residuals(disoc_2t))
qqline(residuals(disoc_2t))
shapiro.test(residuals(disoc_2t))#normal
r.squaredGLMM(disoc_2t)#15% of variation explained by fixed effects, 66% by whole model
diLS2_2t<-lsmeans(disoc_2t, ~Treatment)
contrast(diLS2_2t, "pairwise", adjust="tukey")
lsmeans(disoc_2t, list(trt.vs.ctrl ~ Treatment), adjust="tukey")
#TWM and TWS are diff from control
plot(dSOC~Treatment, idata_during_2)


#check effects of treatments on change in SOC during flood year
isoc2<-lme(dSOC  ~Treatment, random=~1|Block, idata_2001, na.action=na.exclude)
summary(isoc2)
plot(isoc2)
anova(isoc2)#treatment highly sign
qqnorm(residuals(isoc2))
qqline(residuals(isoc2))
shapiro.test(residuals(isoc2))#normal
r.squaredGLMM(isoc2)#35% of variation explained by fixed effects, 55% by whole model
iLS2<-lsmeans(isoc2, ~Treatment)
contrast(iLS2, "pairwise", adjust="tukey")
lsmeans(isoc2, list(trt.vs.ctrl ~ Treatment), adjust="tukey")
#WM and WS are diff from control

SOC_mean5<-filter(SOC_mean3, Year=="2001")
ggplot(SOC_mean5, aes(x=Treatment, y=dSOC.M, color=Treatment, shape=Treatment)) + 
  geom_point(size=3) + 
  geom_hline(yintercept=0)+
  theme_bw() + 
  geom_errorbar(aes(ymax = dSOC.M+dSOC.SE, ymin = dSOC.M-dSOC.SE), width=.25, position=position_dodge(width=0.9)) + 
  labs(x="Treatment", y="Change in SOC (%) from 1993", title = "2001-2002", fill="Treatment") +
  theme(text = element_text(size=25))+
  scale_shape_manual(values=c(16, 17, 18, 19, 4, 25, 8, 15))+
  scale_color_manual(values=c("C"="black", "T"="red", "TWM"="firebrick", "TW"="firebrick", "TWS"="firebrick", "W"="blue", "WM"="navy", "WS"="purple"))+
  annotate("text", x = as.factor(unique(SOC_mean4$Treatment)),y=c(0,0,0,0,0,0,0.29, 0.23),
         label = c("","","", "","","","*","*"),fontface = 3,size=12)

#check effects of treatments on change in SOC at end of experiment
idata_2004<- idata_2 %>% filter(Year=="2004")
isoc4<-lme(dSOC  ~Treatment, random=~1|Block, idata_2004, na.action=na.exclude)
summary(isoc4)
plot(isoc4)
anova(isoc4)#treatment  sign
qqnorm(residuals(isoc4))
qqline(residuals(isoc4))
shapiro.test(residuals(isoc4))#normal
r.squaredGLMM(isoc4)#20% of variation explained by fixed effects, 38% by whole model
iLS4<-lsmeans(isoc4, ~Treatment)
contrast(iLS4, "pairwise", adjust="tukey")
lsmeans(isoc4, list(trt.vs.ctrl ~ Treatment), adjust="tukey")
#WM and WS are diff from control

SOC_new_graphic<-filter(SOC_mean4, Year=="2004")
ggplot(SOC_new_graphic, aes(x=Treatment, y=dSOC.M, color=Treatment, shape=Treatment)) + 
  geom_point(size=3) + 
  geom_hline(yintercept=0)+
  theme_bw() + 
  geom_errorbar(aes(ymax = dSOC.M+dSOC.SE, ymin = dSOC.M-dSOC.SE), width=.25, position=position_dodge(width=0.9)) + 
  labs(x="Treatment", y="Change in SOC (%) from 1993", title = "2004-2005", fill="Treatment") +
  theme(text = element_text(size=25))+
  scale_shape_manual(values=c(16, 17, 18, 19, 4, 25, 8, 15))+
  scale_color_manual(values=c("C"="black", "T"="red", "TWM"="firebrick", "TW"="firebrick", "TWS"="firebrick", "W"="blue", "WM"="navy", "WS"="purple"))+
  annotate("text", x = as.factor(unique(SOC_mean4$Treatment)),y=c(0,0,0,0,0,-0.1,0,0),
           label = c("","","", "","","*","",""),fontface = 3,size=12)

#look at 1995
idata_1995<- idata_2 %>% filter(Year=="1995")
isoc4<-lme(dSOC  ~Treatment, random=~1|Block, idata_1995, na.action=na.exclude)
summary(isoc4)
plot(isoc4)
anova(isoc4)#treatment  sign
qqnorm(residuals(isoc4))
qqline(residuals(isoc4))
shapiro.test(residuals(isoc4))#normal
r.squaredGLMM(isoc4)#20% of variation explained by fixed effects, 38% by whole model
iLS4<-lsmeans(isoc4, ~Treatment)
contrast(iLS4, "pairwise", adjust="tukey")
lsmeans(isoc4, list(trt.vs.ctrl ~ Treatment), adjust="tukey")
#none differ from control


SOC_new_graphic2<-filter(SOC_mean4, Year=="1995")
ggplot(SOC_new_graphic2, aes(x=Treatment, y=dSOC.M, color=Treatment, shape=Treatment)) + 
  geom_point(size=3) + 
  geom_hline(yintercept=0)+
  theme_bw() + 
  geom_errorbar(aes(ymax = dSOC.M+dSOC.SE, ymin = dSOC.M-dSOC.SE), width=.25, position=position_dodge(width=0.9)) + 
  labs(x="Treatment", y="Change in SOC (%) from 1993", title = "1995", fill="Treatment") +
  theme(text = element_text(size=25))+
  scale_shape_manual(values=c(16, 17, 18, 19, 4, 25, 8, 15))+
  scale_color_manual(values=c("C"="black", "T"="red", "TWM"="firebrick", "TW"="firebrick", "TWS"="firebrick", "W"="blue", "WM"="navy", "WS"="purple"))+
  annotate("text", x = as.factor(unique(SOC_mean4$Treatment)),y=c(0,0,0,0,0,-0.1,0,0),
           label = c("","","", "","","*","",""),fontface = 3,size=12)

#look at 1996
idata_1996<- idata_2 %>% filter(Year=="1996")
isoc4<-lme(dSOC  ~Treatment, random=~1|Block, idata_1996, na.action=na.exclude)
summary(isoc4)
plot(isoc4)
anova(isoc4)#treatment  sign
qqnorm(residuals(isoc4))
qqline(residuals(isoc4))
shapiro.test(residuals(isoc4))#normal
r.squaredGLMM(isoc4)#20% of variation explained by fixed effects, 38% by whole model
iLS4<-lsmeans(isoc4, ~Treatment)
contrast(iLS4, "pairwise", adjust="tukey")
lsmeans(isoc4, list(trt.vs.ctrl ~ Treatment), adjust="tukey")
#none diff from control

SOC_new_graphic3<-filter(SOC_mean4, Year=="1996")
ggplot(SOC_new_graphic3, aes(x=Treatment, y=dSOC.M, color=Treatment, shape=Treatment)) + 
  geom_point(size=3) + 
  geom_hline(yintercept=0)+
  theme_bw() + 
  geom_errorbar(aes(ymax = dSOC.M+dSOC.SE, ymin = dSOC.M-dSOC.SE), width=.25, position=position_dodge(width=0.9)) + 
  labs(x="Treatment", y="Change in SOC (%) from 1993", title = "1996", fill="Treatment") +
  theme(text = element_text(size=25))+
  scale_shape_manual(values=c(16, 17, 18, 19, 4, 25, 8, 15))+
  scale_color_manual(values=c("C"="black", "T"="red", "TWM"="firebrick", "TW"="firebrick", "TWS"="firebrick", "W"="blue", "WM"="navy", "WS"="purple"))+
  annotate("text", x = as.factor(unique(SOC_mean4$Treatment)),y=c(0,0,0,0,0,-0.1,0,0),
           label = c("","","", "","","*","",""),fontface = 3,size=12)

#look at 1997
idata_1997<- idata_2 %>% filter(Year=="1997")
isoc4<-lme(dSOC  ~Treatment, random=~1|Block, idata_1997, na.action=na.exclude)
summary(isoc4)
plot(isoc4)
anova(isoc4)#treatment  sign
qqnorm(residuals(isoc4))
qqline(residuals(isoc4))
shapiro.test(residuals(isoc4))#normal
r.squaredGLMM(isoc4)#20% of variation explained by fixed effects, 38% by whole model
iLS4<-lsmeans(isoc4, ~Treatment)
contrast(iLS4, "pairwise", adjust="tukey")
lsmeans(isoc4, list(trt.vs.ctrl ~ Treatment), adjust="tukey")
#none

SOC_new_graphic4<-filter(SOC_mean4, Year=="1997")
ggplot(SOC_new_graphic4, aes(x=Treatment, y=dSOC.M, color=Treatment, shape=Treatment)) + 
  geom_point(size=3) + 
  geom_hline(yintercept=0)+
  theme_bw() + 
  geom_errorbar(aes(ymax = dSOC.M+dSOC.SE, ymin = dSOC.M-dSOC.SE), width=.25, position=position_dodge(width=0.9)) + 
  labs(x="Treatment", y="Change in SOC (%) from 1993", title = "1997", fill="Treatment") +
  theme(text = element_text(size=25))+
  scale_shape_manual(values=c(16, 17, 18, 19, 4, 25, 8, 15))+
  scale_color_manual(values=c("C"="black", "T"="red", "TWM"="firebrick", "TW"="firebrick", "TWS"="firebrick", "W"="blue", "WM"="navy", "WS"="purple"))+
  annotate("text", x = as.factor(unique(SOC_mean4$Treatment)),y=c(0,0,0,0,0,-0.1,0,0),
           label = c("","","", "","","*","",""),fontface = 3,size=12)

#look at 1997
SOC_new_graphic4<-filter(SOC_mean4, Year=="1997")
ggplot(SOC_new_graphic4, aes(x=Treatment, y=dSOC.M, color=Treatment, shape=Treatment)) + 
  geom_point(size=3) + 
  geom_hline(yintercept=0)+
  theme_bw() + 
  geom_errorbar(aes(ymax = dSOC.M+dSOC.SE, ymin = dSOC.M-dSOC.SE), width=.25, position=position_dodge(width=0.9)) + 
  labs(x="Treatment", y="Change in SOC (%) from 1993", title = "1997", fill="Treatment") +
  theme(text = element_text(size=25))+
  scale_shape_manual(values=c(16, 17, 18, 19, 4, 25, 8, 15))+
  scale_color_manual(values=c("C"="black", "T"="red", "TWM"="firebrick", "TW"="firebrick", "TWS"="firebrick", "W"="blue", "WM"="navy", "WS"="purple"))+
  annotate("text", x = as.factor(unique(SOC_mean4$Treatment)),y=c(0,0,0,0,0,-0.1,0,0),
           label = c("","","", "","","*","",""),fontface = 3,size=12)

#look at 1994
idata_1994<- idata_2 %>% filter(Year=="1994")
isoc4<-lme(dSOC  ~Treatment, random=~1|Block, idata_1994, na.action=na.exclude)
summary(isoc4)
plot(isoc4)
anova(isoc4)#treatment  sign
qqnorm(residuals(isoc4))
qqline(residuals(isoc4))
shapiro.test(residuals(isoc4))#normal
r.squaredGLMM(isoc4)#20% of variation explained by fixed effects, 38% by whole model
iLS4<-lsmeans(isoc4, ~Treatment)
contrast(iLS4, "pairwise", adjust="tukey")
lsmeans(isoc4, list(trt.vs.ctrl ~ Treatment), adjust="tukey")
#none

#check effects of treatments on change in SOC for recovery
isoc<-lme(dSOC  ~Treatment, random=~1|Block, idata_2015, na.action=na.exclude)
summary(isoc)
plot(isoc)
anova(isoc)#treatment significant
qqnorm(residuals(isoc))
qqline(residuals(isoc))
shapiro.test(residuals(isoc))#normal
r.squaredGLMM(isoc) #11% of variation explained by fixed effects, 22% by whole model
iLS<-lsmeans(isoc, ~Treatment)
contrast(iLS, "pairwise", adjust="tukey")
lsmeans(isoc, list(trt.vs.ctrl ~ Treatment), adjust="tukey")
#W treatment diff from control

SOC_mean4_graph<-filter(SOC_mean3, Year=="2015")
ggplot(SOC_mean4_graph, aes(x=Treatment, y=dSOC.M, color=Treatment, shape=Treatment)) + 
  geom_point(size=3) + 
  theme_bw() + 
  geom_hline(yintercept=0)+
  geom_errorbar(aes(ymax = dSOC.M+dSOC.SE, ymin = dSOC.M-dSOC.SE), width=.25, position=position_dodge(width=0.9)) + 
  labs(x="Treatment", y="Change in SOC (%) from 1993", title="2015-2016", fill="Treatment") +
  theme(text = element_text(size=25))+
  ylim(-0.4, 0)+
  scale_shape_manual(values=c(16, 17, 18, 19, 4, 25, 8, 15))+
  scale_color_manual(values=c("C"="black", "T"="red", "TWM"="firebrick", "TW"="firebrick", "TWS"="firebrick", "W"="blue", "WM"="navy", "WS"="purple"))+
  annotate("text", x = as.factor(unique(SOC_mean4$Treatment)),y=c(0,0,0,0,0,-0.25,0, 0),
         label = c("","","", "","","","",""),fontface = 3,size=12)


#show change over time
ggplot(d=SOC_mean3, aes(x=Year, y=dSOC.M, color=Treatment, group=Treatment, shape=Treatment)) +
  geom_errorbar(aes(ymin=dSOC.M-dSOC.SE, ymax=dSOC.M+dSOC.SE), width=.1)+
  geom_point(size=3) +
  theme_bw()+
  geom_hline(yintercept=0)+
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  theme(text = element_text(size=25))+
  scale_x_continuous(breaks = round(seq(min(SOC_mean3$Year), max(SOC_mean3$Year), by = 1),1)) +
  labs(x="Year", y="Change in SOC (%) from initial 1993 measurement", title="Change in SOC over time by treatment")+
  scale_shape_manual(values=c(16, 17, 18, 19, 4, 25, 8, 15))+
  scale_color_manual(values=c("C"="black", "T"="red", "TWM"="orange", "TW"="orangered1", "TWS"="firebrick", "W"="blue", "WM"="navy", "WS"="purple"))

ggplot(d=SOC_mean4, aes(x=Year, y=dSOC.M, color=Treatment, group=Treatment, shape=Treatment)) +
  geom_rect(data=NULL,aes(xmin=-Inf,xmax=2005,ymin=-Inf,ymax=Inf),
            fill="lightgray", color="lightgray")+
  theme_bw()+
  scale_x_continuous(breaks = round(seq(min(SOC_mean3$Year), max(SOC_mean3$Year), by = 1),1)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  theme(text = element_text(size=18))+
  geom_hline(yintercept=0)+
  geom_vline(xintercept =  2001, linetype="dotted")+
  geom_errorbar(aes(ymin=dSOC.M-dSOC.SE, ymax=dSOC.M+dSOC.SE), width=.1, position = position_dodge(1))+
  geom_point(size=3, position = position_dodge(1)) +
  labs(x="Year", y="Change in SOC (mg/g) from pre-treatment", title="Change in SOC over time by treatment")+
  scale_shape_manual(values=c(1, 2, 0, 17, 15, 5, 19, 18), labels=c("U","T","TW","TWM", "TWS", "W", "WM", "WS"))+
  scale_color_manual(values=c("C"="black", "T"="red", "TWM"="orange", "TW"="orangered1", "TWS"="firebrick", "W"="blue", "WM"="navy", "WS"="purple"), labels=c("U","T","TW","TWM", "TWS", "W", "WM", "WS"))+
  annotate("text", x = c(2001), y=(0.32), label = c("Flood year"), fontface = 3,size=4)+
  annotate("text", x = c(1998), y=(0.35), label = c("Treatment"), fontface = 2,size=4)+
  annotate("text", x = c(2011), y=(0.35), label = c("Post-Treatment"), fontface = 2,size=4)

gcarbon<-ggplot(subset(SOC_mean4, Treatment %in% c("C", "TWM", "TWS", "WM", "WS")), aes(x=Year, y=dSOC.M, color=Treatment, group=Treatment, shape=Treatment)) +
  geom_rect(data=NULL,aes(xmin=-Inf,xmax=2005,ymin=-Inf,ymax=Inf),
            fill="lightgray", color="lightgray")+
  theme_bw()+
  scale_x_continuous(breaks = round(seq(min(SOC_mean3$Year), max(SOC_mean3$Year), by = 1),1)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  theme(text = element_text(size=12))+
  geom_hline(yintercept=0)+
  geom_vline(xintercept =  2001, linetype="dotted")+
  geom_errorbar(aes(ymin=dSOC.M-dSOC.SE, ymax=dSOC.M+dSOC.SE), width=.1, position = position_dodge(1))+
  geom_point(size=3, position = position_dodge(1)) +
  labs(x="Year", y="Change in SOC (mg/g) from pre-treatment", title="A. Change in SOC over time for carbon addition treatments")+
  scale_shape_manual(values=c(1, 17, 15, 19, 18), labels=c("U","TWM", "TWS", "WM", "WS"))+
  scale_color_manual(values=c("C"="black", "T"="red", "TWM"="orange", "TW"="orangered1", "TWS"="firebrick", "W"="blue", "WM"="navy", "WS"="purple"), labels=c("U","TWM", "TWS", "WM", "WS"))+
  annotate("text", x = c(2001), y=(0.32), label = c("Flood year"), fontface = 3,size=4)+
  annotate("text", x = c(1998), y=(0.35), label = c("Treatment"), fontface = 2,size=4)+
  annotate("text", x = c(2011), y=(0.35), label = c("Post-Treatment"), fontface = 2,size=4)

gnocarbon<-ggplot(subset(SOC_mean4, Treatment %in% c("C", "T", "TW", "W")), aes(x=Year, y=dSOC.M, color=Treatment, group=Treatment, shape=Treatment)) +
  geom_rect(data=NULL,aes(xmin=-Inf,xmax=2005,ymin=-Inf,ymax=Inf),
            fill="lightgray", color="lightgray")+
  theme_bw()+
  scale_x_continuous(breaks = round(seq(min(SOC_mean3$Year), max(SOC_mean3$Year), by = 1),1)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  theme(text = element_text(size=12))+
  geom_hline(yintercept=0)+
  geom_vline(xintercept =  2001, linetype="dotted")+
  geom_errorbar(aes(ymin=dSOC.M-dSOC.SE, ymax=dSOC.M+dSOC.SE), width=.1, position = position_dodge(1))+
  geom_point(size=3, position = position_dodge(1)) +
  labs(x="Year", y="Change in SOC (mg/g) from pre-treatment", title="B. Change in SOC over time for non-carbon treatments")+
  scale_shape_manual(values=c(1, 2, 0, 5), labels=c("U","T","TW", "W"))+
  scale_color_manual(values=c("C"="black", "T"="red", "TWM"="orange", "TW"="orangered1", "TWS"="firebrick", "W"="blue", "WM"="navy", "WS"="purple"), labels=c("U","T","TW", "W"))+
  annotate("text", x = c(2001), y=(0.1), label = c("Flood year"), fontface = 3,size=4)+
  annotate("text", x = c(1998), y=(0.12), label = c("Treatment"), fontface = 2,size=4)+
  annotate("text", x = c(2011), y=(0.12), label = c("Post-Treatment"), fontface = 2,size=4)
library(gridExtra)
grid.arrange(gcarbon,gnocarbon, ncol=1)

#check effects of temp and water manipulations
#does chamber or water addition have effects on soc?
plot(dSOC~moist_manip, idata_2001) #maybe
plot(dSOC~temp_manip, idata_2001) #maybe
msoctw<-lme(dSOC ~moist_manip+temp_manip, random=~1|Block, idata_2001, na.action=na.exclude)
summary(msoctw)
anova(msoctw)
r.squaredGLMM(msoctw)#8% by fixed, 28% by whole model
qqnorm(residuals(msoctw))
qqline(residuals(msoctw))
shapiro.test(residuals(msoctw)) #normal
plot(msoctw)
#no effects 

#does chamber or water addition have effects on soc in recovery?
plot(dSOC~moist_manip, idata_2015) #prob not
plot(dSOC~temp_manip, idata_2015) #maybe
msoctwr<-lme(dSOC ~moist_manip+temp_manip, random=~1|Block, idata_2015, na.action=na.exclude)
summary(msoctwr)
anova(msoctwr)
r.squaredGLMM(msoctwr)#11% by fixed, 31% by whole model
qqnorm(residuals(msoctwr))
qqline(residuals(msoctwr))
shapiro.test(residuals(msoctwr)) #close
plot(msoctwr)
riLS<-lsmeans(msoctwr, ~temp_manip)
contrast(riLS, "pairwise", adjust="tukey")
#temp manip effects

#plot it
SOC_graph <- idata %>% filter(dSOC != 'NA')%>%
  group_by(Time, temp_manip) %>%
  summarize(meansoc=mean(dSOC), sesoc=sd(dSOC)/sqrt(length(dSOC)))

ggplot(SOC_graph, aes(x=temp_manip, y=meansoc, color=temp_manip, shape=temp_manip))+
  geom_point()+
  facet_wrap(~Time)+
  theme_bw()+
  scale_shape_manual(name="Temperature Chamber", values=c(16, 17))+
  scale_color_manual(name="Temperature Chamber", values=c("black","red"))+
  geom_errorbar(aes(ymin=meansoc-sesoc, ymax=meansoc+sesoc), width=.2,position=position_dodge(0.05))+
  labs(x="Temperature Chamber Treatment", y="Change in SOC (%) from 1993")+
  annotate("text", x = c(1.5), y=(-0.2), label = c("","*"), fontface = 3,size=12)

#check effects of treatments on TN  during experiment
dTN<-lme(sqrt(TN)  ~Treatment, random=~1|Year/Block, during, na.action=na.exclude)
summary(dTN)
plot(dTN)
anova(dTN)#treatment not sign
qqnorm(residuals(dTN))
qqline(residuals(dTN))
shapiro.test(residuals(dTN))#not normal
r.squaredGLMM(dTN)#<1% of variation explained by fixed effects, 98% by whole model
lsmeans(dTN, list(trt.vs.ctrl ~ Treatment), adjust="tukey")
#no differences

#check effects of treatments on TN  post experiment
dTNr<-lme(TN  ~Treatment, random=~1|Block, finalyr, na.action=na.exclude)
summary(dTNr)
plot(dTNr)
anova(dTNr)#treatment not sign
qqnorm(residuals(dTNr))
qqline(residuals(dTNr))
shapiro.test(residuals(dTNr))# normal
r.squaredGLMM(dTNr)#22% of variation explained by fixed effects, 23% by whole model
rtnLS<-lsmeans(dTNr, ~Treatment)
contrast(rtnLS, "pairwise", adjust="tukey")
lsmeans(dTNr, list(trt.vs.ctrl ~ Treatment), adjust="tukey")
#W treatment diff from control
#no dif effects

#does chamber or water addition have effects on TN during expt?
plot(TN~moist_manip, during) #?
plot(TN~temp_manip, during) #prob
mTNtw<-lme(sqrt(TN) ~moist_manip+temp_manip, random=~1|Year/Block, during, na.action=na.exclude)
summary(mTNtw)
anova(mTNtw)
r.squaredGLMM(mTNtw)#<1% by fixed, 31% by whole model
qqnorm(residuals(mTNtw))
qqline(residuals(mTNtw))
shapiro.test(residuals(mTNtw)) #better
plot(mTNtw)
tnLS<-lsmeans(mTNtw, ~temp_manip*moist_manip)
contrast(tnLS, "pairwise", adjust="tukey")
#interactive effects
Npct1 <- during %>% filter(TN != 'NA')%>%
  group_by(temp_manip, moist_manip) %>%
  summarize(meanN=mean(TN), seN=sd(TN)/sqrt(length(TN)))


#does chamber or water addition have effects on TN after expt?
plot(TN~moist_manip, finalyr) #?
plot(TN~temp_manip, finalyr) #prob
mTNtwr<-lme(TN ~moist_manip+temp_manip, random=~1|Block, finalyr, na.action=na.exclude)
summary(mTNtwr)
anova(mTNtwr)
r.squaredGLMM(mTNtwr)#10% by fixed, 10% by whole model
qqnorm(residuals(mTNtwr))
qqline(residuals(mTNtwr))
shapiro.test(residuals(mTNtwr)) #normal
plot(mTNtwr)
tnLSr<-lsmeans(mTNtwr, ~temp_manip)
contrast(tnLSr, "pairwise", adjust="tukey")
#temp manipulation is significant

Npct <- data %>% filter(TN != 'NA')%>%
  group_by(Time, temp_manip) %>% filter(Time!="before")%>%
  summarize(meanN=mean(TN), seN=sd(TN)/sqrt(length(TN)))

ggplot(Npct, aes(x=temp_manip, y=meanN, color=temp_manip, shape=temp_manip))+
  geom_point()+
  facet_wrap(~Time)+
  theme_bw()+
  scale_shape_manual(name="Temperature Chamber", values=c(16, 17))+
  scale_color_manual(name="Temperature Chamber", values=c("black","red"))+
  geom_errorbar(aes(ymin=meanN-seN, ymax=meanN+seN), width=.2,position=position_dodge(0.05))+
  labs(x="Temperature Chamber Treatment", y="Total Nitrogen %")+
  annotate("text", x = c(1.5), y=(0.015), label = c("","*"), fontface = 3,size=12)

#note: cannot check effects on TC during expt bc no data
#check effects of treatments on TC  after experiment
dTC<-lme(TC  ~Treatment, random=~1|Year/Block, finalyr, na.action=na.exclude, weights=varIdent(form=~1|Year*Treatment))
summary(dTC)
plot(dTC)
anova(dTC)#treatment not sign
qqnorm(residuals(dTC))
qqline(residuals(dTC))
shapiro.test(residuals(dTC))#not normal, but close
r.squaredGLMM(dTC)#5% of variation explained by fixed effects, 14% by whole model
lsmeans(dTC, list(trt.vs.ctrl ~ Treatment), adjust="tukey")
#no differences



#does chamber or water addition have effects on TN after expt?
plot(TC~moist_manip, finalyr) #maybe
plot(TC~temp_manip, finalyr) #maybe
mTCtwr<-lme(TC ~moist_manip+temp_manip, random=~1|Block, finalyr, na.action=na.exclude)
summary(mTCtwr)
anova(mTCtwr)
r.squaredGLMM(mTCtwr)#10% by fixed, 10% by whole model
qqnorm(residuals(mTCtwr))
qqline(residuals(mTCtwr))
shapiro.test(residuals(mTCtwr)) #normal
plot(mTCtwr)
tnLSr<-lsmeans(mTCtwr, ~temp_manip)
contrast(tnLSr, "pairwise", adjust="tukey")
#temp manipulation is significant

Cpct <-data %>% filter (TC != 'NA')%>%
  group_by(Year, Treatment) %>%
  summarize(meanC=mean(TC), seN=sd(TC)/sqrt(length(TC)))


##Moisture (chamber and moisture manip effects)
#does chamber or water addition have effects on moisture during expt?
plot(moisture~moist_manip, during) #maybe
plot(moisture~temp_manip, during) #prob
mWtwr<-lme(log(moisture+1) ~moist_manip+temp_manip, random=~1|Year/Block, during, na.action=na.exclude, weights=varIdent(form=~1|Year))
summary(mWtwr)
anova(mWtwr)#both are significant
r.squaredGLMM(mWtwr)#9% by fixed, 73% by whole model
qqnorm(residuals(mWtwr))
qqline(residuals(mWtwr))
shapiro.test(residuals(mWtwr)) #better
plot(mWtwr)
WLS<-lsmeans(mWtwr, ~temp_manip+moist_manip)
contrast(WLS, "pairwise", adjust="tukey")
#both water and temp manipulation is significant
Wtemp <- data %>% filter(moisture != 'NA')%>%
  group_by(Time, temp_manip) %>% filter(Time!="before")%>%
  summarize(meanW=mean(moisture), seW=sd(moisture)/sqrt(length(moisture)))

ggplot(Wtemp, aes(x=temp_manip, y=meanW, color=temp_manip, shape=temp_manip))+
  geom_point()+
  facet_wrap(~Time)+
  theme_bw()+
  scale_shape_manual(name="Temperature Chamber", values=c(16, 17))+
  scale_color_manual(name="Temperature Chamber", values=c("black","red"))+
  geom_errorbar(aes(ymin=meanW-seW, ymax=meanW+seW), width=.2,position=position_dodge(0.05))+
  labs(x="Temperature Chamber Treatment", y="Gravimetric Moisture %")+
  annotate("text", x = c(1.5), y=(1.5), label = c("*",""), fontface = 3,size=12)

Wmoist <- data %>% filter(moisture != 'NA')%>%
  group_by(Time, moist_manip) %>% filter(Time!="before")%>%
  summarize(meanW=mean(moisture), seW=sd(moisture)/sqrt(length(moisture)))

ggplot(Wmoist, aes(x=moist_manip, y=meanW, color=moist_manip, shape=moist_manip))+
  geom_point()+
  facet_wrap(~Time)+
  theme_bw()+
  scale_shape_manual(name="Water Addition", values=c(16, 25))+
  scale_color_manual(name="Water Addition", values=c("black","blue"))+
  geom_errorbar(aes(ymin=meanW-seW, ymax=meanW+seW), width=.2,position=position_dodge(0.05))+
  labs(x="Moisture Treatment", y="Gravimetric Moisture %")+
  annotate("text", x = c(1.5), y=(1.3), label = c("*",""), fontface = 3,size=12)

#does chamber or water addition have effects on moisture during expt?
plot(moisture~moist_manip, recovery) #maybe
plot(moisture~temp_manip, recovery) #?
mWtw<-lme(log(moisture+1) ~moist_manip+temp_manip, random=~1|Year/Block, recovery, na.action=na.exclude, weights=varIdent(form=~1|Year))
summary(mWtw)
anova(mWtw)#both are significant
r.squaredGLMM(mWtw)#<1% by fixed, 19% by whole model
qqnorm(residuals(mWtw))
qqline(residuals(mWtw))
shapiro.test(residuals(mWtw)) #not normal at all - tried log, sqrt, weights=valIdent
plot(mWtw)
WLSr<-lsmeans(mWtw, ~temp_manip+moist_manip)
contrast(WLSr, "pairwise", adjust="tukey")
#neither water or temp manipulation is significant

#what about the change in Scottnema from initial (1993)?
#first create column of initial values
initial.s<- data %>% filter(Time=="before") %>% dplyr::select(Treatment, Block, STL) %>% rename(STLi=STL)
idata_2_s<- data1 %>% full_join(initial.s, by=c("Treatment", "Block")) %>% mutate(dSTL=STL-STLi, na.rm=T) %>% mutate(dSTLpc=(STL-STLi)/STLi*100, na.rm=T)
STL_mean <- idata_2_s %>% group_by(Year, Treatment) %>% filter(!is.na(dSTL))%>%
  summarize(dSTL.M=mean(dSTL), dSTL.SE=sd(dSTL)/sqrt(length(dSTL)))

ggplot(d=STL_mean, aes(x=Year, y=dSTL.M, color=Treatment, group=Treatment, shape=Treatment)) +
  geom_rect(data=NULL,aes(xmin=-Inf,xmax=2005,ymin=-Inf,ymax=Inf),
            fill="lightgray", color="lightgray")+
  theme_bw()+
  scale_x_continuous(breaks = round(seq(min(SOC_mean3$Year), max(SOC_mean3$Year), by = 1),1)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  theme(text = element_text(size=18))+
  geom_hline(yintercept=0)+
  geom_vline(xintercept =  2001, linetype="dotted")+
  geom_errorbar(aes(ymin=dSTL.M-dSTL.SE, ymax=dSTL.M+dSTL.SE), width=.1, position = position_dodge(1))+
  geom_point(size=3, position = position_dodge(1)) +
  labs(x="Year", y="Change in Scottnema (abundance/kg dry soil) from pre-treatment", title="Change in Scottnema populations over time by treatment")+
  scale_shape_manual(values=c(1, 2, 0, 17, 15, 5, 19, 18), labels=c("U","T","TW","TWM", "TWS", "W", "WM", "WS"))+
  scale_color_manual(values=c("C"="black", "T"="red", "TWM"="orange", "TW"="orangered1", "TWS"="firebrick", "W"="blue", "WM"="navy", "WS"="purple"), labels=c("U","T","TW","TWM", "TWS", "W", "WM", "WS"))+
  annotate("text", x = c(2001), y=(3300), label = c("Flood year"), fontface = 3,size=4)+
  annotate("text", x = c(1996), y=(2900), label = c("Treatment"), fontface = 2,size=4)+
  annotate("text", x = c(2011), y=(2900), label = c("Post-Treatment"), fontface = 2,size=4)

#what about the change in Eudorylaimus from initial (1993)?
#first create column of initial values
initial.e<- data %>% filter(Time=="before") %>% dplyr::select(Treatment, Block, ETL) %>% rename(ETLi=ETL)
idata_2_e<- data1 %>% full_join(initial.e, by=c("Treatment", "Block")) %>% mutate(dETL=ETL-ETLi, na.rm=T) %>% mutate(dETLpc=(ETL-ETLi)/ETLi*100, na.rm=T)
ETL_mean <- idata_2_e %>% group_by(Year, Treatment) %>% filter(!is.na(dETL))%>%
  summarize(dETL.M=mean(dETL), dETL.SE=sd(dETL)/sqrt(length(dETL)))

ggplot(d=ETL_mean, aes(x=Year, y=dETL.M, color=Treatment, group=Treatment, shape=Treatment)) +
  geom_rect(data=NULL,aes(xmin=-Inf,xmax=2005,ymin=-Inf,ymax=Inf),
            fill="lightgray", color="lightgray")+
  theme_classic()+
  scale_x_continuous(breaks = round(seq(min(SOC_mean3$Year), max(SOC_mean3$Year), by = 1),1)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  theme(text = element_text(size=18))+
  geom_hline(yintercept=0)+
  geom_vline(xintercept =  2001, linetype="dotted")+
  geom_errorbar(aes(ymin=dETL.M-dETL.SE, ymax=dETL.M+dETL.SE), width=.1, position = position_dodge(1))+
  geom_point(size=3, position = position_dodge(1)) +
  labs(x="Year", y="Change in Eudorylaimus (abundance/kg dry soil) from pre-treatment", title="Change in Eudorylaimus populations over time by treatment")+
  scale_shape_manual(values=c(1, 2, 0, 17, 15, 5, 19, 18), labels=c("U","T","TW","TWM", "TWS", "W", "WM", "WS"))+
  scale_color_manual(values=c("C"="black", "T"="red", "TWM"="orange", "TW"="orangered1", "TWS"="firebrick", "W"="blue", "WM"="navy", "WS"="purple"), labels=c("U","T","TW","TWM", "TWS", "W", "WM", "WS"))+
  annotate("text", x = c(2001), y=(210), label = c("Flood year"), fontface = 3,size=4)+
  annotate("text", x = c(1996), y=(180), label = c("Treatment"), fontface = 2,size=4)+
  annotate("text", x = c(2011), y=(180), label = c("Post-Treatment"), fontface = 2,size=4)

#what about the change in Plectus from initial (1993)?
#first create column of initial values
initial.p<- data %>% filter(Time=="before") %>% dplyr::select(Treatment, Block, PTL) %>% rename(PTLi=PTL)
idata_2_p<- data1 %>% full_join(initial.p, by=c("Treatment", "Block")) %>% mutate(dPTL=PTL-PTLi, na.rm=T) %>% mutate(dPTLpc=(PTL-PTLi)/PTLi*100, na.rm=T)
PTL_mean <- idata_2_p %>% group_by(Year, Treatment) %>% filter(!is.na(dPTL))%>%
  summarize(dPTL.M=mean(dPTL), dPTL.SE=sd(dPTL)/sqrt(length(dPTL)))

ggplot(d=PTL_mean, aes(x=Year, y=dPTL.M, color=Treatment, group=Treatment, shape=Treatment)) +
  geom_rect(data=NULL,aes(xmin=-Inf,xmax=2005,ymin=-Inf,ymax=Inf),
            fill="lightgray", color="lightgray")+
  theme_bw()+
  scale_x_continuous(breaks = round(seq(min(SOC_mean3$Year), max(SOC_mean3$Year), by = 1),1)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  theme(text = element_text(size=25))+
  geom_hline(yintercept=0)+
  geom_vline(xintercept =  2001, linetype="dotted")+
  geom_errorbar(aes(ymin=dPTL.M-dPTL.SE, ymax=dPTL.M+dPTL.SE), width=.1, position = position_dodge(1))+
  geom_point(size=3, position = position_dodge(1)) +
  labs(x="Year", y="Change in Plectus (abundance/kg dry soil) from pre-treatment", title="Change in Plectus populations over time by treatment")+
  scale_shape_manual(values=c(1, 2, 0, 17, 15, 5, 19, 18), labels=c("U","T","TW","TWM", "TWS", "W", "WM", "WS"))+
  scale_color_manual(values=c("C"="black", "T"="red", "TWM"="orange", "TW"="orangered1", "TWS"="firebrick", "W"="blue", "WM"="navy", "WS"="purple"), labels=c("U","T","TW","TWM", "TWS", "W", "WM", "WS"))+
  annotate("text", x = c(2001), y=(20), label = c("Flood year"), fontface = 3,size=4)+
  annotate("text", x = c(1996), y=(18), label = c("Treatment"), fontface = 2,size=4)+
  annotate("text", x = c(2011), y=(18), label = c("Post-Treatment"), fontface = 2,size=4)

###Calculate log response ratios and confidence intervals
library(SingleCaseES)
lnrr.s<-LRRd(A_data=idata_2_s$SLTi, B_data=idata_2_s$STL, condition=Treatment)

###Panel of control conditions
##a) temperature, b) moisture, c) soc, d) chla, e) STL, f) ETL
###Scottnema and Eudorylaimus
nema.c<- data14 %>% filter(Treatment=="C") %>% dplyr::select(Treatment, Block, Year, STL, ETL) %>% 
  group_by(Year) %>% summarize(STL.m=mean(STL), STL.SE=sd(STL)/sqrt(length(STL)), ETL.m=mean(ETL), ETL.SE=sd(ETL)/sqrt(length(ETL))) %>%
  gather(key=nematode, value=value, STL.m, ETL.m, STL.SE, ETL.SE) %>%
  separate(nematode, c("nematode", "x")) %>%
  spread(x, value) %>%
  mutate(time="Pre-flood", time=ifelse(Year>2001, "Post-flood", time))

stl.p1<-ggplot(d=subset(nema.c, nematode=="STL"), aes(x=Year, y=m, shape=nematode)) +
  #geom_rect(data=NULL,aes(xmin=-Inf,xmax=2005,ymin=-Inf,ymax=Inf),fill="lightgray", color="lightgray")+
  theme_classic()+
  geom_boxplot(aes(y=m, x=Year, group=time, color=time))+
  scale_x_continuous(breaks = round(seq(min(nema.c$Year), max(nema.c$Year), by = 1),1)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  theme(text = element_text(size=14), legend.position="none")+
  #geom_hline(yintercept=0)+
  geom_vline(xintercept =  2001, linetype="dotted")+
  geom_errorbar(aes(ymin=m-SE, ymax=m+SE), width=.1, position = position_dodge(1))+
  geom_point(size=3, position = position_dodge(1)) +
  labs(x="Year", y="Scottnema (#/kg dry soil)", title="e)")+
  #scale_shape_manual(values=c(1, 2, 0, 17, 15, 5, 19, 18), labels=c("U","T","TW","TWM", "TWS", "W", "WM", "WS"))+
  scale_color_manual(values=c("Pre-flood"="grey", "Post-flood"="grey"), labels=c("Pre-flood","Post-flood"))+
  annotate("text", x = c(2000.8), y=(1700), label = c("Flood year"), fontface = 3,size=4)
  #annotate("text", x = c(1996), y=(2900), label = c("Treatment"), fontface = 2,size=4)+
  #annotate("text", x = c(2011), y=(2900), label = c("Post-Treatment"), fontface = 2,size=4)

etl.p1<-ggplot(d=subset(nema.c, nematode=="ETL"), aes(x=Year, y=m, group=nematode, shape=nematode)) +
  #geom_rect(data=NULL,aes(xmin=-Inf,xmax=2005,ymin=-Inf,ymax=Inf),fill="lightgray", color="lightgray")+
  theme_classic()+
  geom_boxplot(aes(y=m, x=Year, group=time, color=time))+
  scale_x_continuous(breaks = round(seq(min(nema.c$Year), max(nema.c$Year), by = 1),1)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  theme(text = element_text(size=14), legend.position="none")+
  #geom_hline(yintercept=0)+
  geom_vline(xintercept =  2001, linetype="dotted")+
  geom_errorbar(aes(ymin=m-SE, ymax=m+SE), width=.1, position = position_dodge(1))+
  geom_point(size=3, position = position_dodge(1)) +
  labs(x="Year", y="Eudorylaimus (#/kg dry soil)", title="f)")+
  #scale_shape_manual(values=c(1, 2, 0, 17, 15, 5, 19, 18), labels=c("U","T","TW","TWM", "TWS", "W", "WM", "WS"))+
  scale_color_manual(values=c("Pre-flood"="grey", "Post-flood"="grey"), labels=c("Pre-flood","Post-flood"))+
  annotate("text", x = c(2000.8), y=(65), label = c("Flood year"), fontface = 3,size=4)

#soil moisture
sm.c<- data14 %>% filter(Treatment=="C") %>% dplyr::select(Treatment, Block, Year, moisture) %>% 
  filter(moisture != "NA") %>% group_by(Year) %>% summarize(m=mean(moisture), SE=sd(moisture)/sqrt(length(moisture)))%>%
  mutate(time="Pre-flood", time=ifelse(Year>2001, "Post-flood", time)) 

sm.p1<-ggplot(d=sm.c, aes(x=Year, y=m)) +
  #geom_rect(data=NULL,aes(xmin=-Inf,xmax=2005,ymin=-Inf,ymax=Inf),fill="lightgray", color="lightgray")+
  theme_classic()+
  geom_boxplot(aes(y=m, x=Year, group=time, color=time))+
  scale_x_continuous(breaks = round(seq(min(nema.c$Year), max(nema.c$Year), by = 1),1)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  theme(text = element_text(size=14), legend.position="none")+
  #geom_hline(yintercept=0)+
  geom_vline(xintercept =  2001, linetype="dotted")+
  geom_errorbar(aes(ymin=m-SE, ymax=m+SE), width=.1, position = position_dodge(1))+
  geom_point(size=3, position = position_dodge(1)) +
  labs(x="Year", y="Soil moisture (g/g)", title="b)")+
  #scale_shape_manual(values=c(1, 2, 0, 17, 15, 5, 19, 18), labels=c("U","T","TW","TWM", "TWS", "W", "WM", "WS"))+
  scale_color_manual(values=c("Pre-flood"="grey", "Post-flood"="grey"), labels=c("Pre-flood","Post-flood"))+
  annotate("text", x = c(2000.8), y=(3), label = c("Flood year"), fontface = 3,size=4)

#soil moisture
soc.c<- data14 %>% filter(Treatment=="C") %>% dplyr::select(Treatment, Block, Year, SOC2) %>% 
  filter(SOC2!="NA") %>% group_by(Year) %>% summarize(m=mean(SOC2), SE=sd(SOC2)/sqrt(length(SOC2)))%>%
  mutate(time="Pre-flood", time=ifelse(Year>2001, "Post-flood", time))

soc.p1<-ggplot(d=soc.c, aes(x=Year, y=m)) +
  #geom_rect(data=NULL,aes(xmin=-Inf,xmax=2005,ymin=-Inf,ymax=Inf),fill="lightgray", color="lightgray")+
  theme_classic()+
  geom_boxplot(aes(y=m, x=Year, group=time, color=time))+
  scale_x_continuous(breaks = round(seq(min(nema.c$Year), max(nema.c$Year), by = 1),1)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  theme(text = element_text(size=14), legend.position="none")+
  #geom_hline(yintercept=0)+
  geom_vline(xintercept =  2001, linetype="dotted")+
  geom_errorbar(aes(ymin=m-SE, ymax=m+SE), width=.1, position = position_dodge(1))+
  geom_point(size=3, position = position_dodge(1)) +
  labs(x="Year", y="SOC (mg/g)", title="d)")+
  #scale_shape_manual(values=c(1, 2, 0, 17, 15, 5, 19, 18), labels=c("U","T","TW","TWM", "TWS", "W", "WM", "WS"))+
  scale_color_manual(values=c("Pre-flood"="grey", "Post-flood"="grey"), labels=c("Pre-flood","Post-flood"))+
  annotate("text", x = c(2000.8), y=(0.3), label = c("Flood year"), fontface = 3,size=4)

#soil moisture
chla.c<- data14 %>% filter(Treatment=="C") %>% dplyr::select(Treatment, Block, Year, Chl_a) %>% 
  filter(Chl_a!="NA") %>% group_by(Year) %>% summarize(m=mean(Chl_a), SE=sd(Chl_a)/sqrt(length(Chl_a)))%>%
  mutate(time="Pre-flood", time=ifelse(Year>2001, "Post-flood", time))

chla.p1<-ggplot(d=chla.c, aes(x=Year, y=m)) +
  #geom_rect(data=NULL,aes(xmin=-Inf,xmax=2005,ymin=-Inf,ymax=Inf),fill="lightgray", color="lightgray")+
  theme_classic()+
  geom_boxplot(aes(y=m, x=Year, group=time, color=time))+
  scale_x_continuous(breaks = round(seq(min(nema.c$Year), max(nema.c$Year), by = 1),1)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  theme(text = element_text(size=14), legend.position="none")+
  #geom_hline(yintercept=0)+
  geom_vline(xintercept =  2001, linetype="dotted")+
  geom_errorbar(aes(ymin=m-SE, ymax=m+SE), width=.1, position = position_dodge(1))+
  geom_point(size=3, position = position_dodge(1)) +
  labs(x="Year", y="Chlorophyll-a (ug/g)", title="c)")+
  #scale_shape_manual(values=c(1, 2, 0, 17, 15, 5, 19, 18), labels=c("U","T","TW","TWM", "TWS", "W", "WM", "WS"))+
  scale_color_manual(values=c("Pre-flood"="grey", "Post-flood"="grey"), labels=c("Pre-flood","Post-flood"))+
  annotate("text", x = c(2000.8), y=(0.02), label = c("Flood year"), fontface = 3,size=4)
chla.p1

library(vegan)
sp<-data%>%dplyr::select(STL,ETL,PTL)
H <- diversity(sp)
data$Pie <- H/log(specnumber(sp))
piej.c<- data %>% filter(Treatment=="C") %>% dplyr::select(Treatment, Block, Year, Pie) %>% 
  filter(Pie!="NaN") %>% group_by(Year) %>% summarize(m=mean(Pie), SE=sd(Pie)/sqrt(length(Pie))) %>%
  mutate(time="Pre-flood", time=ifelse(Year>2001, "Post-flood", time))

j.p1<-ggplot(d=piej.c, aes(x=Year, y=m)) +
  #geom_rect(data=NULL,aes(xmin=-Inf,xmax=2005,ymin=-Inf,ymax=Inf),fill="lightgray", color="lightgray")+
  theme_classic()+
  geom_boxplot(aes(y=m, x=Year, group=time, color=time))+
  scale_x_continuous(breaks = round(seq(min(nema.c$Year), max(nema.c$Year), by = 1),1)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  theme(text = element_text(size=14), legend.position="none")+
  #geom_hline(yintercept=0)+
  geom_vline(xintercept =  2001, linetype="dotted")+
  geom_errorbar(aes(ymin=m-SE, ymax=m+SE), width=.1, position = position_dodge(1))+
  geom_point(size=3, position = position_dodge(1)) +
  labs(x="Year", y="Diversity (Pielou's J)", title="g)")+
  #scale_shape_manual(values=c(1, 2, 0, 17, 15, 5, 19, 18), labels=c("U","T","TW","TWM", "TWS", "W", "WM", "WS"))+
  scale_color_manual(values=c("Pre-flood"="grey", "Post-flood"="grey"), labels=c("Pre-flood","Post-flood"))+
  annotate("text", x = c(2000.8), y=(0.82), label = c("Flood year"), fontface = 3,size=4)

library(gridExtra)
t <- ggplot(d=piej.c, x=year)
grid.arrange(t, sm.p1 ,chla.p1,soc.p1, stl.p1, etl.p1, ncol=2)

###Panel of expt conditions
##a) soc, b) chla, c) STL, d) ETL
###Scottnema and Eudorylaimus

stl.p2<-ggplot(d=subset(STL_mean,Treatment!="C"), aes(x=Year, y=dSTL.M, color=Treatment, group=Treatment, shape=Treatment)) +
  geom_rect(data=NULL,aes(xmin=-Inf,xmax=2005,ymin=-Inf,ymax=Inf),
            fill="lightgray", color="lightgray")+
  theme_classic()+
  scale_x_continuous(breaks = round(seq(min(STL_mean$Year), max(STL_mean$Year), by = 1),1)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  theme(text = element_text(size=14), legend.position="none")+
  geom_hline(yintercept=0)+
  geom_vline(xintercept =  2001, linetype="dotted")+
  geom_errorbar(aes(ymin=dSTL.M-dSTL.SE, ymax=dSTL.M+dSTL.SE), width=.1, position = position_dodge(1))+
  geom_point(size=3, position = position_dodge(1)) +
  labs(x="Year", y="Change in Scottnema (#/kg dry soil)", title="c)")+
  scale_shape_manual(values=c( 2, 0, 17, 15, 5, 19, 18), labels=c("T","TW","TWM", "TWS", "W", "WM", "WS"))+
  scale_color_manual(values=c( "T"="red", "TWM"="orange", "TW"="orangered1", "TWS"="firebrick", "W"="blue", "WM"="navy", "WS"="purple"), labels=c("T","TW","TWM", "TWS", "W", "WM", "WS"))+
  annotate("text", x = c(2001), y=(3300), label = c("Flood year"), fontface = 3,size=4)+
  annotate("text", x = c(1996), y=(2900), label = c("Treatment"), fontface = 2,size=4)+
  annotate("text", x = c(2011), y=(2900), label = c("Post-Treatment"), fontface = 2,size=4)

etl.p2<-ggplot(d=subset(ETL_mean, Treatment!="C"), aes(x=Year, y=dETL.M, color=Treatment, group=Treatment, shape=Treatment)) +
  geom_rect(data=NULL,aes(xmin=-Inf,xmax=2005,ymin=-Inf,ymax=Inf),
            fill="lightgray", color="lightgray")+
  theme_classic()+
  scale_x_continuous(breaks = round(seq(min(STL_mean$Year), max(STL_mean$Year), by = 1),1)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  theme(text = element_text(size=14), legend.position="none")+
  geom_hline(yintercept=0)+
  geom_vline(xintercept =  2001, linetype="dotted")+
  geom_errorbar(aes(ymin=dETL.M-dETL.SE, ymax=dETL.M+dETL.SE), width=.1, position = position_dodge(1))+
  geom_point(size=3, position = position_dodge(1)) +
  labs(x="Year", y="Change in Eudorylaimus (#/kg dry soil)", title="d)")+
  scale_shape_manual(values=c( 2, 0, 17, 15, 5, 19, 18), labels=c("T","TW","TWM", "TWS", "W", "WM", "WS"))+
  scale_color_manual(values=c( "T"="red", "TWM"="orange", "TW"="orangered1", "TWS"="firebrick", "W"="blue", "WM"="navy", "WS"="purple"), labels=c("T","TW","TWM", "TWS", "W", "WM", "WS"))+
  annotate("text", x = c(2001), y=(210), label = c("Flood year"), fontface = 3,size=4)+
  annotate("text", x = c(1996), y=(180), label = c("Treatment"), fontface = 2,size=4)+
  annotate("text", x = c(2011), y=(180), label = c("Post-Treatment"), fontface = 2,size=4)

soc.p2<-ggplot(d=subset(SOC_mean4, Treatment!="C"), aes(x=Year, y=dSOC.M, color=Treatment, group=Treatment, shape=Treatment)) +
  geom_rect(data=NULL,aes(xmin=-Inf,xmax=2005,ymin=-Inf,ymax=Inf),
            fill="lightgray", color="lightgray")+
  theme_classic()+
  scale_x_continuous(breaks = round(seq(min(STL_mean$Year), max(STL_mean$Year), by = 1),1)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  theme(text = element_text(size=14), legend.position="none")+
  geom_hline(yintercept=0)+
  geom_vline(xintercept =  2001, linetype="dotted")+
  geom_errorbar(aes(ymin=dSOC.M-dSOC.SE, ymax=dSOC.M+dSOC.SE), width=.1, position = position_dodge(1))+
  geom_point(size=3, position = position_dodge(1)) +
  labs(x="Year", y="Change in SOC (mg/g)", title="b)")+
  scale_shape_manual(values=c(2, 0, 17, 15, 5, 19, 18), labels=c("T","TW","TWM", "TWS", "W", "WM", "WS"))+
  scale_color_manual(values=c("T"="red", "TWM"="orange", "TW"="orangered1", "TWS"="firebrick", "W"="blue", "WM"="navy", "WS"="purple"), labels=c("T","TW","TWM", "TWS", "W", "WM", "WS"))+
  annotate("text", x = c(2001), y=(0.32), label = c("Flood year"), fontface = 3,size=4)+
  annotate("text", x = c(1998), y=(0.35), label = c("Treatment"), fontface = 2,size=4)+
  annotate("text", x = c(2011), y=(0.35), label = c("Post-Treatment"), fontface = 2,size=4)

initial.chl<- data %>% filter(Year=="1995") %>% dplyr::select(Treatment, Block, Chl_a) %>% rename(chlai=Chl_a)
idata_2_chl<- data1 %>% full_join(initial.chl, by=c("Treatment", "Block")) %>% mutate(dchla=Chl_a-chlai, na.rm=T) %>% mutate(dchlapc=(Chl_a-chlai)/chlai*100, na.rm=T)
CHL_mean <- idata_2_chl %>% group_by(Year, Treatment) %>% filter(!is.na(dchla))%>%
  summarize(dchla.M=mean(dchla), dchla.SE=sd(dchla)/sqrt(length(dchla)))

chl.p2<-ggplot(d=subset(CHL_mean, Treatment!="C"), aes(x=Year, y=dchla.M, color=Treatment, group=Treatment, shape=Treatment)) +
  geom_rect(data=NULL,aes(xmin=-Inf,xmax=2005,ymin=-Inf,ymax=Inf),
            fill="lightgray", color="lightgray")+
  theme_classic()+
  scale_x_continuous(breaks = round(seq(min(STL_mean$Year), max(STL_mean$Year), by = 1),1)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  theme(text = element_text(size=14), legend.position="none")+
  geom_hline(yintercept=0)+
  geom_vline(xintercept =  2001, linetype="dotted")+
  geom_errorbar(aes(ymin=dchla.M-dchla.SE, ymax=dchla.M+dchla.SE), width=.1, position = position_dodge(1))+
  geom_point(size=3, position = position_dodge(1)) +
  labs(x="Year", y="Change in Chla (ug/g)", title="a)")+
  scale_shape_manual(values=c(2, 0, 17, 15, 5, 19, 18), labels=c("T","TW","TWM", "TWS", "W", "WM", "WS"))+
  scale_color_manual(values=c("T"="red", "TWM"="orange", "TW"="orangered1", "TWS"="firebrick", "W"="blue", "WM"="navy", "WS"="purple"), labels=c("T","TW","TWM", "TWS", "W", "WM", "WS"))+
  annotate("text", x = c(2000.8), y=(0.025), label = c("Flood year"), fontface = 3,size=4)+
  annotate("text", x = c(1998), y=(0.03), label = c("Treatment"), fontface = 2,size=4)+
  annotate("text", x = c(2011), y=(0.03), label = c("Post-Treatment"), fontface = 2,size=4)

#Pielou's J
initial.piej<- data %>% filter(Time=="before") %>% dplyr::select(Treatment, Block, Pie) %>% rename(Ji=Pie)
idata_2_j<- data %>% full_join(initial.piej, by=c("Treatment", "Block")) %>% mutate(dJ=Pie-Ji, na.rm=T) %>% mutate(dJpc=(Pie-Ji)/Ji*100, na.rm=T)
dJ_mean <- idata_2_j %>% group_by(Year, Treatment) %>% filter(!is.na(dJ))%>%
  summarize(J.M=mean(dJ), J.SE=sd(dJ)/sqrt(length(dJ)))

j.p2<-ggplot(d=subset(dJ_mean, Treatment!="C"), aes(x=Year, y=J.M, color=Treatment, group=Treatment, shape=Treatment)) +
  geom_rect(data=NULL,aes(xmin=-Inf,xmax=2005,ymin=-Inf,ymax=Inf),
            fill="lightgray", color="lightgray")+
  theme_classic()+
  scale_x_continuous(breaks = round(seq(min(STL_mean$Year), max(STL_mean$Year), by = 1),1)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  theme(text = element_text(size=14), legend.position="none")+
  geom_hline(yintercept=0)+
  geom_vline(xintercept =  2001, linetype="dotted")+
  geom_errorbar(aes(ymin=J.M-J.SE, ymax=J.M+J.SE), width=.1, position = position_dodge(1))+
  geom_point(size=3, position = position_dodge(1)) +
  labs(x="Year", y="Diversity (Pielou's J)", title="e)")+
  scale_shape_manual(values=c(2, 0, 17, 15, 5, 19, 18), labels=c("T","TW","TWM", "TWS", "W", "WM", "WS"))+
  scale_color_manual(values=c("T"="red", "TWM"="orange", "TW"="orangered1", "TWS"="firebrick", "W"="blue", "WM"="navy", "WS"="purple"), labels=c("T","TW","TWM", "TWS", "W", "WM", "WS"))+
  annotate("text", x = c(2000.8), y=(0.65), label = c("Flood year"), fontface = 3,size=4)+
  annotate("text", x = c(1998), y=(0.7), label = c("Treatment"), fontface = 2,size=4)+
  annotate("text", x = c(2011), y=(0.7), label = c("Post-Treatment"), fontface = 2,size=4)

grid.arrange(chl.p2,soc.p2, stl.p2, etl.p2,j.p2, ncol=2)

