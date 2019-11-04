library(tidyverse)
library(readr)
library(vegan)
library(MASS)
library(dplyr)
library(RVAideMemoire) #for posthoc tests on permanova

# Import csv file, call it data
setwd("~/Desktop/")
data<-read.csv("LTM_all.csv",header=T)
str(data)
levels(data$Treatment)
levels(data$Time)
levels(data$Yearf)

#remove november sampling and 2014 since it has controls only
november<-data %>%subset(Month=='11')
data<- data %>% subset(Month!= '11') %>% subset(Year!='2014')

M <- table(data$Yearf, data$Treatment)
M

#separate 'before', 'during', and 'recovery'
during <- subset(data, Time== 'during')
recovery <- subset(data, Time== 'recovery')
before <- subset(data, Time== 'before')
finalyr <- subset(data, Yearf=="2015-16")
data1<- data %>% filter(Time!="before")

#dat2<-dat %>% mutate(cover=cover+0.00001) #a work around for data with lots of zeros
dat<- data %>% dplyr::select(Block, Treatment, Year, SML,SMD, SFL,SFD, SJL,SJD, EML,EMD, EFL, EFD, EJL,EJD, PTL, PTD) 
dat[is.na(dat)] <- 0 #replace NAs with 0 (species not counted in plots have NAs when wide dataset created)

#Import environmental data - NADP pinnacles deposition data and NOAA San Jose temperature data
env <- data %>% dplyr::select(Block, Treatment, Year, moisture, pH, TN, TC, SOC2, conductivity)
env.all <- env %>% dplyr::select(pH, TN, TC, SOC2, conductivity)
env.all<- as.data.frame(env.all[-c(642, 683, 695), ])
dat.drop <- as.data.frame(dat[-c(642, 683, 695), ])


data$ID <- seq.int(nrow(data))
plotnames<-data[,58]
cover.Bio<- dat %>% dplyr::select(-c(1:3)) #wide data with ID columns removed, only species/cover for NMDS
rownames(cover.Bio)<-plotnames

#check for empty rows
cover.Biodrop1<-cover.Bio[rowSums(cover.Bio[, (1:14)]) ==0, ] #no empty rows, next step not needed
cover.Biodrop<-cover.Bio[rowSums(cover.Bio[, (1:14)])  >0, ]#remove empty rows


#if needed, relativize by row or column or calculate presence/absence
cover.rowsums <- rowSums(cover.Biodrop [1:14])
cover.relrow <- data.frame(cover.Biodrop /cover.rowsums)
#cover.colmax<-sapply(cover.Bio ,max)
#cover.relcolmax <- data.frame(sweep(cover.Bio ,2,cover.colmax,'/'))
#cover.pa <- cover.Bio %>% mutate_each(funs(ifelse(.>0,1,0)), 1:57)


######################
#2. NMS
#see also section 2.1 of vegan tutorial: 
#http://cc.oulu.fi/~jarioksa/opetus/metodi/vegantutor.pdf
######################

#make bray-curtis dissimilarity matrix
spp.bcd <- vegdist(cover.relrow)

#quick run to check out NMS ordination
spp.mds0 <-isoMDS(spp.bcd) #runs nms only once
spp.mds0  #by default 2 dimensions returned, stress is 6.4, converged
ordiplot(spp.mds0) #ugly

#prefer to run multiple NMS ordinations
#with different starting configurations, and choose the best
#this function does that, and in addition does several other steps too including: 
#standardizing the data (though fuction call below turns this off with autotransform=F)
#calculating distance matrix (default bray-curtis)
#running NMDS with random starts
#rotation of axes to maximize variance of site scores on axis 1
#calculate species scores based on weighted averaging
#help(metaMDS)
spp.mds<-metaMDS(cover.relrow, trace = TRUE, autotransform=T, trymax=100, k=6) #runs several with different starting configurations
#trace= TRUE will give output for step by step what its doing
#default is 2 dimensions, can put k=4 for 4 dimensions
spp.mds #solution did not converge after 100 tries
summary(spp.mds)

#quick plot of results
stressplot(spp.mds, spp.bcd) #stressplot to show fit
ordiplot(spp.mds)

#overlay environmental variables on full data
envvec.nms<-envfit(spp.mds,env.all, na.rm=TRUE)
envvec.nms
plot(spp.mds)
plot(envvec.nms) #add vectors to previous ordination

#store scores in new dataframe
spscores1<-scores(spp.mds,display="sites",choices=1)
spscores2<-scores(spp.mds,display="sites",choices=2)
tplots<-dat.drop[,2]
tplot_levels<-levels(tplots)
spscoresall<-data.frame(tplots,spscores1,spscores2)

#make nicer plot colored based on treatment, shapes on pre/post flood
#help(ordiplot)
#first, set colors and shapes
cols1<- dat.drop %>% dplyr::select(Treatment) %>% mutate(color = "black", 
                                                   color = ifelse(Treatment == "T","red", 
                                                                  ifelse(Treatment=="TWM","orange", 
                                                                         ifelse(Treatment=="TW","orangered1",
                                                                                ifelse(Treatment=="TWS","firebrick",
                                                                                       ifelse(Treatment== "W","blue",
                                                                                              ifelse(Treatment== "WM","navy",
                                                                                                     ifelse(Treatment=="WS","purple",color)))))))) #colors based on thermal group
Lcols <- rep(c("Black", "Red", "Orange", "Orangered1", "Firebrick","Blue", "Navy", "Purple")) #colors for the legend
shapes <- dat.drop %>% dplyr::select(Year) %>%
  mutate(shape = 16, shape = ifelse(Year == "2001", 8, 
                                   ifelse(Year<"2001" & Year<"2006", 16,
                                          ifelse(Year>="2006", 15,shape)))) #shapes based on year 
shapes<-shapes %>% mutate(time="Treatment", 
                          time= ifelse(Year==2001, "Flood",
                                       ifelse(Year<"2001" & Year<"2006", "Treatment",
                                              ifelse(Year>="2006", "Post-Treatment",time)))) 
Lshapes <-rep(c(8,15, 16))#shapes for legend
#make the plot
bio.plot <- ordiplot(spp.mds, choices=c(1,2), type = "none")   #Set up the plot
points(spscoresall$NMDS1,spscoresall$NMDS2,col=cols1$color,pch=shapes$shape) 
plot(envvec.nms, col="green")
text(spp.mds, display = "species", cex=0.5, col="grey30") #label species
legend("bottomright",legend=levels(as.factor(cols1$Treatment)), col=Lcols, pch=15, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
legend("bottomleft",legend=levels(as.factor(shapes$time)), col="black", pch=Lshapes, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)

#make nicer plot colored based on treatment, shapes on pre/post flood
#help(ordiplot)
#first, set colors and shapes
cols2<- dat.drop %>% dplyr::select(Year) %>% mutate(color = "black", 
                                                         color = ifelse(Year == "1994","red", 
                                                                        ifelse(Year=="1995","orange", 
                                                                               ifelse(Year=="1996","orangered1",
                                                                                      ifelse(Year=="1997","firebrick",
                                                                                             ifelse(Year== "2001","blue",
                                                                                                    ifelse(Year== "2004","navy",
                                                                                                           ifelse(Year=="2006","purple",
                                                                                                                  ifelse(Year=="2008", "pink",
                                                                                                                         ifelse(Year=="2011", "magenta", 
                                                                                                                                ifelse(Year=="2015", "grey",color))))))))))) #colors based on year 
#Lcols <- rep(c("Black", "Red", "Orange", "Orangered1", "Firebrick","Blue", "Navy", "Purple")) #colors for the legend
#shapes <- dat.drop %>% dplyr::select(Year) %>%
  mutate(shape = 16, shape = ifelse(Year == "2001", 8, 
                                    ifelse(Year<"2001" & Year<"2006", 16,
                                           ifelse(Year>="2006", 15,shape)))) #shapes based on year 
#shapes<-shapes %>% mutate(time="Treatment", 
                          time= ifelse(Year==2001, "Flood",
                                       ifelse(Year<"2001" & Year<"2006", "Treatment",
                                              ifelse(Year>="2006", "Post-Treatment",time)))) 
#Lshapes <-rep(c(8,15, 16))#shapes for legend
#make the plot
bio.plot <- ordiplot(spp.mds, choices=c(1,2), type = "none")   #Set up the plot
points(spscoresall$NMDS1,spscoresall$NMDS2,col=cols2$color,pch=shapes$shape) 
plot(envvec.nms, col="green")
text(spp.mds, display = "species", cex=0.5, col="grey30") #label species
legend("bottomright",legend=levels(as.factor(cols1$Treatment)), col=Lcols, pch=15, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
legend("bottomleft",legend=levels(as.factor(shapes$time)), col="black", pch=Lshapes, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
