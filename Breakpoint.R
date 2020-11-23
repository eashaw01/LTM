library(nlme)
library(lsmeans)
library(tidyverse)
library(ggplot2)
library(segmented)
library(gridExtra)

#breakpoint analysis for 1993-2015 data for Scottnema populations, Eudorylaimus populations, Soil Organic Carbon, and Chlorophyll a
## Set ggplot2 theme
theme_set(theme_bw())
theme_update( panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              strip.background = element_blank(),
              text = element_text(size = 16),
              strip.text= element_text(size = 14), 
              axis.text = element_text(size = 12))


# Import csv file, call it data
#setwd("~/Desktop/PhD Research/Manuscripts/In prep/LTM")
#data<-read.csv("LTM_all2.csv",header=T)
setwd("~/Desktop/LTM")
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


##breakpoint analysis for Scottnema
C.STL<-data %>% subset(Treatment=="C")%>%group_by(Year) %>% filter(STL!="NA")%>%summarise(mean=mean(STL))

#linear model of STL over time
stl.lm <- lm(mean ~ Year, data = C.STL)
summary(stl.lm)

# Extract te coefficients from the overall model
stl.coef <- coef(stl.lm)
stl.coef

p <- ggplot(C.STL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
         y = "Scottnema total live")
p

p <- p + geom_abline(intercept = stl.coef[1], 
                     slope = stl.coef[2], 
                     aes(colour = "overall"))
p

stl.lm3 <- lm(mean ~ poly(Year, 3), data = C.STL)

p + geom_smooth(method = "lm",
                formula = y ~ poly(x, degree = 3), 
                se = FALSE, colour = "orange")
my.seg <- segmented(stl.lm, 
                    seg.Z = ~ Year, 
                    psi = list(Year = c(2004)))
summary(my.seg)
# get the breakpoints
my.seg$psi
slope(my.seg)
# get the fitted data
my.fitted <- fitted(my.seg)
my.model <- data.frame(Year = C.STL$Year, STL = my.fitted)
ggplot(my.model, aes(x = Year, y = STL)) + geom_line()

pSTL<- ggplot(C.STL, aes(x = Year, y = mean)) +
  geom_segment(aes(x=1993, y=2090, xend=2004.5, yend=2090), color = "grey0") + 
  geom_point(aes(x=1993, y=2090), color="grey0")+
  geom_point(aes(x=2004.5, y=2090), color="grey0")+
  geom_segment(aes(x=2005, y=2090, xend=2015, yend=2090), color = "grey60") + 
  geom_point(aes(x=2005, y=2090), color="grey60")+
  geom_point(aes(x=2015, y=2090), color="grey60")+
  geom_segment(aes(x=2001, y=165, xend=2001, yend=0), arrow=arrow(length = unit(0.03, "npc")), color = "grey0") + 
  annotate("text", x = c(2001), y=(220), label = c("Flood year"), fontface = 3,size=4)+
  annotate("text", x = c(1999), y=(2200), label = c("Treatment"), fontface = 2,size=4)+
  annotate("text", x = c(2010), y=(2200), label = c("Post-Treatment"), fontface = 2,size=4, color="grey60")+
  labs(x = "Year",y = "Scottnema (#/kg dry soil)")+ 
  geom_point(color="black", alpha=0.3)+
  geom_line(data = my.model, aes(x = Year, y = STL), colour = "black", cex=1.5)+
  ylim(0,2200)
pSTL

#Add STL for WM mannitol treatment
M.STL<-data %>% subset(Treatment=="WM")%>%group_by(Year) %>% filter(STL!="NA")%>%summarise(mean=mean(STL))

#linear model of STL over time
stl.lm.m <- lm(mean ~ Year, data = M.STL)
summary(stl.lm.m)

# Extract te coefficients from the overall model
stl.coef.m <- coef(stl.lm.m)
stl.coef.m

p.2 <- ggplot(M.STL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Scottnema total live")
p.2 

p.2  <- p.2  + geom_abline(intercept = stl.coef.m[1], 
                     slope = stl.coef.m[2], 
                     aes(colour = "overall"))
p.2 

stl.lm3 <- lm(mean ~ poly(Year, 3), data = M.STL)

p.2  + geom_smooth(method = "lm",
                formula = y ~ poly(x, degree = 3), 
                se = FALSE, colour = "orange")
my.seg.m <- segmented(stl.lm.m, 
                    seg.Z = ~ Year, 
                    psi = list(Year = c(2004)))
summary(my.seg.m)
# get the breakpoints
my.seg.m$psi
slope(my.seg.m)
# get the fitted data
my.fitted.m <- fitted(my.seg.m)
my.model.m <- data.frame(Year = M.STL$Year, STL = my.fitted.m)
ggplot(my.model.m, aes(x = Year, y = STL)) + geom_line()

p.2 + geom_line(data = my.model.m, aes(x = Year, y = STL), colour = "tomato")

##add to original plot with controls

pSTL <- pSTL + geom_line(data = my.model.m, aes(x = Year, y = STL), colour = "steelblue4", cex=1.5)
pSTL <- pSTL +  geom_point(data=M.STL, aes(x = Year, y = mean), color="steelblue4", alpha=0.3)
pSTL

#Add STL for WS sucrose treatment
S.STL<-data %>% subset(Treatment=="WS")%>%group_by(Year) %>% filter(STL!="NA")%>%summarise(mean=mean(STL))

#linear model of STL over time
stl.lm.s <- lm(mean ~ Year, data = S.STL)
summary(stl.lm.s)

# Extract te coefficients from the overall model
stl.coef.s <- coef(stl.lm.s)
stl.coef.s

p.3 <- ggplot(S.STL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Scottnema total live")
p.3 

p.3  <- p.3  + geom_abline(intercept = stl.coef.s[1], 
                           slope = stl.coef.s[2], 
                           aes(colour = "overall"))
p.3 

stl.lm3.s <- lm(mean ~ poly(Year, 3), data = S.STL)

p.3  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.s <- segmented(stl.lm.s, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(2001)))
summary(my.seg.s)
# get the breakpoints
my.seg.s$psi
slope(my.seg.s)
# get the fitted data
my.fitted.s <- fitted(my.seg.s)
my.model.s <- data.frame(Year = S.STL$Year, STL = my.fitted.s)
ggplot(my.model.s, aes(x = Year, y = STL)) + geom_line()

p.3 + geom_line(data = my.model.s, aes(x = Year, y = STL), colour = "purple")

##add to original plot with controls

pSTL <- pSTL + geom_line(data = my.model.s, aes(x = Year, y = STL), colour = "orchid3", cex=1.5)
pSTL <- pSTL +  geom_point(data=S.STL, aes(x = Year, y = mean), color="orchid3", alpha=0.3)
pSTL

#Add STL for W water treatment
W.STL<-data %>% subset(Treatment=="W")%>%group_by(Year) %>% filter(STL!="NA")%>%summarise(mean=mean(STL))

#linear model of STL over time
stl.lm.w <- lm(mean ~ Year, data = W.STL)
summary(stl.lm.w)

# Extract te coefficients from the overall model
stl.coef.w <- coef(stl.lm.w)
stl.coef.w

p.4 <- ggplot(W.STL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Scottnema total live")
p.4 

p.4  <- p.4  + geom_abline(intercept = stl.coef.w[1], 
                           slope = stl.coef.w[2], 
                           aes(colour = "overall"))
p.4 

stl.lm3.w <- lm(mean ~ poly(Year, 3), data = W.STL)

p.4  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.w <- segmented(stl.lm.w, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(2001)))
summary(my.seg.w)
# get the breakpoints
my.seg.w$psi
slope(my.seg.w)
# get the fitted data
my.fitted.w <- fitted(my.seg.w)
my.model.w <- data.frame(Year = W.STL$Year, STL = my.fitted.w)
ggplot(my.model.w, aes(x = Year, y = STL)) + geom_line()

p.4 + geom_line(data = my.model.w, aes(x = Year, y = STL), colour = "purple")

##add to original plot with controls

pSTL <- pSTL + geom_line(data = my.model.w, aes(x = Year, y = STL), colour = "deepskyblue", cex=1.5)
pSTL <- pSTL +  geom_point(data=W.STL, aes(x = Year, y = mean), color="deepskyblue", alpha=0.3)
pSTL
pSTL1<-pSTL

#Add STL for T temp treatment
T.STL<-data %>% subset(Treatment=="T")%>%group_by(Year) %>% filter(STL!="NA")%>%summarise(mean=mean(STL))

#linear model of STL over time
stl.lm.t <- lm(mean ~ Year, data = T.STL)
summary(stl.lm.t)

# Extract the coefficients from the overall model
stl.coef.t <- coef(stl.lm.t)
stl.coef.t

p.5 <- ggplot(T.STL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Scottnema total live")
p.5 

p.5  <- p.5  + geom_abline(intercept = stl.coef.t[1], 
                           slope = stl.coef.t[2], 
                           aes(colour = "overall"))
p.5 

stl.lm3.t <- lm(mean ~ poly(Year, 3), data = T.STL)

p.5  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.t <- segmented(stl.lm.t, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(1995)))
summary(my.seg.t)
# get the breakpoints
my.seg.t$psi
slope(my.seg.t)
# get the fitted data
my.fitted.t <- fitted(my.seg.t)
my.model.t <- data.frame(Year = T.STL$Year, STL = my.fitted.t)
ggplot(my.model.t, aes(x = Year, y = STL)) + geom_line()

p.5 + geom_line(data = my.model.t, aes(x = Year, y = STL), colour = "red")

##add to original plot with controls

pSTL <- pSTL + geom_line(data = my.model.t, aes(x = Year, y = STL), colour = "red")
pSTL

#Add STL for TW temp/water treatment
TW.STL<-data %>% subset(Treatment=="TW")%>%group_by(Year) %>% filter(STL!="NA")%>%summarise(mean=mean(STL))

#linear model of STL over time
stl.lm.tw <- lm(mean ~ Year, data = TW.STL)
summary(stl.lm.tw)

# Extract the coefficients from the overall model
stl.coef.tw <- coef(stl.lm.tw)
stl.coef.tw

p.6 <- ggplot(TW.STL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Scottnema total live")
p.6 

p.6  <- p.6  + geom_abline(intercept = stl.coef.tw[1], 
                           slope = stl.coef.tw[2], 
                           aes(colour = "overall"))
p.6 

stl.lm3.tw <- lm(mean ~ poly(Year, 3), data = TW.STL)

p.6  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.tw <- segmented(stl.lm.tw, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(2001)))
summary(my.seg.tw)
# get the breakpoints
my.seg.tw$psi
slope(my.seg.tw)
# get the fitted data
my.fitted.tw <- fitted(my.seg.tw)
my.model.tw <- data.frame(Year = TW.STL$Year, STL = my.fitted.tw)
ggplot(my.model.tw, aes(x = Year, y = STL)) + geom_line()

p.6 + geom_line(data = my.model.tw, aes(x = Year, y = STL), colour = "orangered1")

##add to original plot with controls

pSTL <- pSTL + geom_line(data = my.model.tw, aes(x = Year, y = STL), colour = "orangered1")
pSTL

#Add STL for TWM temp/water/mannitol treatment
TWM.STL<-data %>% subset(Treatment=="TWM")%>%group_by(Year) %>% filter(STL!="NA")%>%summarise(mean=mean(STL))

#linear model of STL over time
stl.lm.twm <- lm(mean ~ Year, data = TWM.STL)
summary(stl.lm.twm)

# Extract the coefficients from the overall model
stl.coef.twm <- coef(stl.lm.twm)
stl.coef.twm

p.7 <- ggplot(TWM.STL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Scottnema total live")
p.7

p.7  <- p.7  + geom_abline(intercept = stl.coef.twm[1], 
                           slope = stl.coef.twm[2], 
                           aes(colour = "overall"))
p.7 

stl.lm3.twm <- lm(mean ~ poly(Year, 3), data = TWM.STL)

p.7  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.twm <- segmented(stl.lm.twm, 
                       seg.Z = ~ Year, 
                       psi = list(Year = c(2001)))
summary(my.seg.twm)
# get the breakpoints
my.seg.twm$psi
slope(my.seg.twm)
# get the fitted data
my.fitted.twm <- fitted(my.seg.twm)
my.model.twm <- data.frame(Year = TWM.STL$Year, STL = my.fitted.twm)
ggplot(my.model.twm, aes(x = Year, y = STL)) + geom_line()

p.7 + geom_line(data = my.model.twm, aes(x = Year, y = STL), colour = "orange")

##add to original plot with controls

pSTL <- pSTL + geom_line(data = my.model.twm, aes(x = Year, y = STL), colour = "orange")
pSTL

#Add STL for TWM temp/water/mannitol treatment
TWS.STL<-data %>% subset(Treatment=="TWS")%>%group_by(Year) %>% filter(STL!="NA")%>%summarise(mean=mean(STL))

#linear model of STL over time
stl.lm.tws <- lm(mean ~ Year, data = TWS.STL)
summary(stl.lm.tws)

# Extract the coefficients from the overall model
stl.coef.tws <- coef(stl.lm.tws)
stl.coef.tws

p.8 <- ggplot(TWS.STL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Scottnema total live")
p.8 

p.8   <- p.8  + geom_abline(intercept = stl.coef.tws[1], 
                           slope = stl.coef.tws[2], 
                           aes(colour = "overall"))
p.8 

stl.lm3.tws <- lm(mean ~ poly(Year, 3), data = TWS.STL)

p.8  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.tws <- segmented(stl.lm.tws, 
                        seg.Z = ~ Year, 
                        psi = list(Year = c(2001)))
summary(my.seg.tws)
# get the breakpoints
my.seg.tws$psi
slope(my.seg.tws)
# get the fitted data
my.fitted.tws <- fitted(my.seg.tws)
my.model.tws <- data.frame(Year = TWS.STL$Year, STL = my.fitted.tws)
ggplot(my.model.tws, aes(x = Year, y = STL)) + geom_line()

p.8 + geom_line(data = my.model.tws, aes(x = Year, y = STL), colour = "firebrick")

##add to original plot with controls

pSTL <- pSTL + geom_line(data = my.model.tws, aes(x = Year, y = STL), colour = "firebrick")
pSTL

################
#scottnema juveniles

C.SJL<-data %>% subset(Treatment=="C")%>%group_by(Year) %>% filter(SJL!="NA")%>%summarise(mean=mean(SJL))

#linear model of SJL over time
SJL.lm <- lm(mean ~ Year, data = C.SJL)
summary(SJL.lm)

# Extract te coefficients from the overall model
SJL.coef <- coef(SJL.lm)
SJL.coef

p <- ggplot(C.SJL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Scottnema total live")
p

p <- p + geom_abline(intercept = SJL.coef[1], 
                     slope = SJL.coef[2], 
                     aes(colour = "overall"))
p

SJL.lm3 <- lm(mean ~ poly(Year, 3), data = C.SJL)

p + geom_smooth(method = "lm",
                formula = y ~ poly(x, degree = 3), 
                se = FALSE, colour = "orange")
my.seg <- segmented(SJL.lm, 
                    seg.Z = ~ Year, 
                    psi = list(Year = c(2004)))
summary(my.seg)
# get the breakpoints
my.seg$psi
slope(my.seg)
# get the fitted data
my.fitted <- fitted(my.seg)
my.model <- data.frame(Year = C.SJL$Year, SJL = my.fitted)
ggplot(my.model, aes(x = Year, y = SJL)) + geom_line()

pSJL<- ggplot(C.SJL, aes(x = Year, y = mean)) + 
  labs(x = "Year",y = "Scottnema juveniles live")+ 
  geom_line(data = my.model, aes(x = Year, y = SJL), colour = "black")+
  ylim(0,1800)
pSJL

#Add SJL for WM mannitol treatment
M.SJL<-data %>% subset(Treatment=="WM")%>%group_by(Year) %>% filter(SJL!="NA")%>%summarise(mean=mean(SJL))

#linear model of SJL over time
SJL.lm.m <- lm(mean ~ Year, data = M.SJL)
summary(SJL.lm.m)

# Extract te coefficients from the overall model
SJL.coef.m <- coef(SJL.lm.m)
SJL.coef.m

p.2 <- ggplot(M.SJL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Scottnema total live")
p.2 

p.2  <- p.2  + geom_abline(intercept = SJL.coef.m[1], 
                           slope = SJL.coef.m[2], 
                           aes(colour = "overall"))
p.2 

SJL.lm3 <- lm(mean ~ poly(Year, 3), data = M.SJL)

p.2  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.m <- segmented(SJL.lm.m, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(2004)))
summary(my.seg.m)
# get the breakpoints
my.seg.m$psi
slope(my.seg.m)
# get the fitted data
my.fitted.m <- fitted(my.seg.m)
my.model.m <- data.frame(Year = M.SJL$Year, SJL = my.fitted.m)
ggplot(my.model.m, aes(x = Year, y = SJL)) + geom_line()

p.2 + geom_line(data = my.model.m, aes(x = Year, y = SJL), colour = "tomato")

##add to original plot with controls

pSJL <- pSJL + geom_line(data = my.model.m, aes(x = Year, y = SJL), colour = "navy")
pSJL

#Add SJL for WS sucrose treatment
S.SJL<-data %>% subset(Treatment=="WS")%>%group_by(Year) %>% filter(SJL!="NA")%>%summarise(mean=mean(SJL))

#linear model of SJL over time
SJL.lm.s <- lm(mean ~ Year, data = S.SJL)
summary(SJL.lm.s)

# Extract te coefficients from the overall model
SJL.coef.s <- coef(SJL.lm.s)
SJL.coef.s

p.3 <- ggplot(S.SJL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Scottnema total live")
p.3 

p.3  <- p.3  + geom_abline(intercept = SJL.coef.s[1], 
                           slope = SJL.coef.s[2], 
                           aes(colour = "overall"))
p.3 

SJL.lm3.s <- lm(mean ~ poly(Year, 3), data = S.SJL)

p.3  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.s <- segmented(SJL.lm.s, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(2001)))
summary(my.seg.s)
# get the breakpoints
my.seg.s$psi
slope(my.seg.s)
# get the fitted data
my.fitted.s <- fitted(my.seg.s)
my.model.s <- data.frame(Year = S.SJL$Year, SJL = my.fitted.s)
ggplot(my.model.s, aes(x = Year, y = SJL)) + geom_line()

p.3 + geom_line(data = my.model.s, aes(x = Year, y = SJL), colour = "purple")

##add to original plot with controls

pSJL <- pSJL + geom_line(data = my.model.s, aes(x = Year, y = SJL), colour = "purple")
pSJL

#Add SJL for W water treatment
W.SJL<-data %>% subset(Treatment=="W")%>%group_by(Year) %>% filter(SJL!="NA")%>%summarise(mean=mean(SJL))

#linear model of SJL over time
SJL.lm.w <- lm(mean ~ Year, data = W.SJL)
summary(SJL.lm.w)

# Extract te coefficients from the overall model
SJL.coef.w <- coef(SJL.lm.w)
SJL.coef.w

p.4 <- ggplot(W.SJL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Scottnema total live")
p.4 

p.4  <- p.4  + geom_abline(intercept = SJL.coef.w[1], 
                           slope = SJL.coef.w[2], 
                           aes(colour = "overall"))
p.4 

SJL.lm3.w <- lm(mean ~ poly(Year, 3), data = W.SJL)

p.4  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.w <- segmented(SJL.lm.w, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(2001)))
summary(my.seg.w)
# get the breakpoints
my.seg.w$psi
slope(my.seg.w)
# get the fitted data
my.fitted.w <- fitted(my.seg.w)
my.model.w <- data.frame(Year = W.SJL$Year, SJL = my.fitted.w)
ggplot(my.model.w, aes(x = Year, y = SJL)) + geom_line()

p.4 + geom_line(data = my.model.w, aes(x = Year, y = SJL), colour = "purple")

##add to original plot with controls

pSJL <- pSJL + geom_line(data = my.model.w, aes(x = Year, y = SJL), colour = "blue")
pSJL

#Add SJL for T temp treatment
T.SJL<-data %>% subset(Treatment=="T")%>%group_by(Year) %>% filter(SJL!="NA")%>%summarise(mean=mean(SJL))

#linear model of SJL over time
SJL.lm.t <- lm(mean ~ Year, data = T.SJL)
summary(SJL.lm.t)

# Extract the coefficients from the overall model
SJL.coef.t <- coef(SJL.lm.t)
SJL.coef.t

p.5 <- ggplot(T.SJL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Scottnema total live")
p.5 

p.5  <- p.5  + geom_abline(intercept = SJL.coef.t[1], 
                           slope = SJL.coef.t[2], 
                           aes(colour = "overall"))
p.5 

SJL.lm3.t <- lm(mean ~ poly(Year, 3), data = T.SJL)

p.5  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.t <- segmented(SJL.lm.t, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(2001)))
summary(my.seg.t)
# get the breakpoints
my.seg.t$psi
slope(my.seg.t)
# get the fitted data
my.fitted.t <- fitted(my.seg.t)
my.model.t <- data.frame(Year = T.SJL$Year, SJL = my.fitted.t)
ggplot(my.model.t, aes(x = Year, y = SJL)) + geom_line()

p.5 + geom_line(data = my.model.t, aes(x = Year, y = SJL), colour = "red")

##add to original plot with controls

pSJL <- pSJL + geom_line(data = my.model.t, aes(x = Year, y = SJL), colour = "red")
pSJL

#Add SJL for TW temp/water treatment
TW.SJL<-data %>% subset(Treatment=="TW")%>%group_by(Year) %>% filter(SJL!="NA")%>%summarise(mean=mean(SJL))

#linear model of SJL over time
SJL.lm.tw <- lm(mean ~ Year, data = TW.SJL)
summary(SJL.lm.tw)

# Extract the coefficients from the overall model
SJL.coef.tw <- coef(SJL.lm.tw)
SJL.coef.tw

p.6 <- ggplot(TW.SJL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Scottnema total live")
p.6 

p.6  <- p.6  + geom_abline(intercept = SJL.coef.tw[1], 
                           slope = SJL.coef.tw[2], 
                           aes(colour = "overall"))
p.6 

SJL.lm3.tw <- lm(mean ~ poly(Year, 3), data = TW.SJL)

p.6  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.tw <- segmented(SJL.lm.tw, 
                       seg.Z = ~ Year, 
                       psi = list(Year = c(2001)))
summary(my.seg.tw)
# get the breakpoints
my.seg.tw$psi
slope(my.seg.tw)
# get the fitted data
my.fitted.tw <- fitted(my.seg.tw)
my.model.tw <- data.frame(Year = TW.SJL$Year, SJL = my.fitted.tw)
ggplot(my.model.tw, aes(x = Year, y = SJL)) + geom_line()

p.6 + geom_line(data = my.model.tw, aes(x = Year, y = SJL), colour = "orangered1")

##add to original plot with controls

pSJL <- pSJL + geom_line(data = my.model.tw, aes(x = Year, y = SJL), colour = "orangered1")
pSJL

#Add SJL for TWM temp/water/mannitol treatment
TWM.SJL<-data %>% subset(Treatment=="TWM")%>%group_by(Year) %>% filter(SJL!="NA")%>%summarise(mean=mean(SJL))

#linear model of SJL over time
SJL.lm.twm <- lm(mean ~ Year, data = TWM.SJL)
summary(SJL.lm.twm)

# Extract the coefficients from the overall model
SJL.coef.twm <- coef(SJL.lm.twm)
SJL.coef.twm

p.7 <- ggplot(TWM.SJL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Scottnema total live")
p.7

p.7  <- p.7  + geom_abline(intercept = SJL.coef.twm[1], 
                           slope = SJL.coef.twm[2], 
                           aes(colour = "overall"))
p.7 

SJL.lm3.twm <- lm(mean ~ poly(Year, 3), data = TWM.SJL)

p.7  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.twm <- segmented(SJL.lm.twm, 
                        seg.Z = ~ Year, 
                        psi = list(Year = c(2001)))
summary(my.seg.twm)
# get the breakpoints
my.seg.twm$psi
slope(my.seg.twm)
# get the fitted data
my.fitted.twm <- fitted(my.seg.twm)
my.model.twm <- data.frame(Year = TWM.SJL$Year, SJL = my.fitted.twm)
ggplot(my.model.twm, aes(x = Year, y = SJL)) + geom_line()

p.7 + geom_line(data = my.model.twm, aes(x = Year, y = SJL), colour = "orange")

##add to original plot with controls

pSJL <- pSJL + geom_line(data = my.model.twm, aes(x = Year, y = SJL), colour = "orange")
pSJL

#Add SJL for TWM temp/water/mannitol treatment
TWS.SJL<-data %>% subset(Treatment=="TWS")%>%group_by(Year) %>% filter(SJL!="NA")%>%summarise(mean=mean(SJL))

#linear model of SJL over time
SJL.lm.tws <- lm(mean ~ Year, data = TWS.SJL)
summary(SJL.lm.tws)

# Extract the coefficients from the overall model
SJL.coef.tws <- coef(SJL.lm.tws)
SJL.coef.tws

p.8 <- ggplot(TWS.SJL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Scottnema total live")
p.8 

p.8   <- p.8  + geom_abline(intercept = SJL.coef.tws[1], 
                            slope = SJL.coef.tws[2], 
                            aes(colour = "overall"))
p.8 

SJL.lm3.tws <- lm(mean ~ poly(Year, 3), data = TWS.SJL)

p.8  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.tws <- segmented(SJL.lm.tws, 
                        seg.Z = ~ Year, 
                        psi = list(Year = c(2001)))
summary(my.seg.tws)
# get the breakpoints
my.seg.tws$psi
slope(my.seg.tws)
# get the fitted data
my.fitted.tws <- fitted(my.seg.tws)
my.model.tws <- data.frame(Year = TWS.SJL$Year, SJL = my.fitted.tws)
ggplot(my.model.tws, aes(x = Year, y = SJL)) + geom_line()

p.8 + geom_line(data = my.model.tws, aes(x = Year, y = SJL), colour = "firebrick")

##add to original plot with controls

pSJL <- pSJL + geom_line(data = my.model.tws, aes(x = Year, y = SJL), colour = "firebrick")
pSJL

################
#scottnema adults

C.SAL<-data %>% subset(Treatment=="C")%>%group_by(Year) %>% filter(STL!="NA")%>%summarise(mean=mean(SML+SFL))

#linear model of SAL over time
SAL.lm <- lm(mean ~ Year, data = C.SAL)
summary(SAL.lm)

# Extract te coefficients from the overall model
SAL.coef <- coef(SAL.lm)
SAL.coef

p <- ggplot(C.SAL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Scottnema total live")
p

p <- p + geom_abline(intercept = SAL.coef[1], 
                     slope = SAL.coef[2], 
                     aes(colour = "overall"))
p

SAL.lm3 <- lm(mean ~ poly(Year, 3), data = C.SAL)

p + geom_smooth(method = "lm",
                formula = y ~ poly(x, degree = 3), 
                se = FALSE, colour = "orange")
my.seg <- segmented(SAL.lm, 
                    seg.Z = ~ Year, 
                    psi = list(Year = c(2004)))
summary(my.seg)
# get the breakpoints
my.seg$psi
slope(my.seg)
# get the fitted data
my.fitted <- fitted(my.seg)
my.model <- data.frame(Year = C.SAL$Year, SAL = my.fitted)
ggplot(my.model, aes(x = Year, y = SAL)) + geom_line()

pSAL<- ggplot(C.SAL, aes(x = Year, y = mean)) + 
  labs(x = "Year",y = "Scottnema adults live")+ 
  geom_line(data = my.model, aes(x = Year, y = SAL), colour = "black")+
  ylim(0,1800)
pSAL

#Add SAL for WM mannitol treatment
M.SAL<-data %>% subset(Treatment=="WM")%>%group_by(Year) %>% filter(STL!="NA")%>%summarise(mean=mean(SML+SFL))

#linear model of SAL over time
SAL.lm.m <- lm(mean ~ Year, data = M.SAL)
summary(SAL.lm.m)

# Extract te coefficients from the overall model
SAL.coef.m <- coef(SAL.lm.m)
SAL.coef.m

p.2 <- ggplot(M.SAL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Scottnema total live")
p.2 

p.2  <- p.2  + geom_abline(intercept = SAL.coef.m[1], 
                           slope = SAL.coef.m[2], 
                           aes(colour = "overall"))
p.2 

SAL.lm3 <- lm(mean ~ poly(Year, 3), data = M.SAL)

p.2  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.m <- segmented(SAL.lm.m, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(2004)))
summary(my.seg.m)
# get the breakpoints
my.seg.m$psi
slope(my.seg.m)
# get the fitted data
my.fitted.m <- fitted(my.seg.m)
my.model.m <- data.frame(Year = M.SAL$Year, SAL = my.fitted.m)
ggplot(my.model.m, aes(x = Year, y = SAL)) + geom_line()

p.2 + geom_line(data = my.model.m, aes(x = Year, y = SAL), colour = "tomato")

##add to original plot with controls

pSAL <- pSAL + geom_line(data = my.model.m, aes(x = Year, y = SAL), colour = "navy")
pSAL

#Add SAL for WS sucrose treatment
S.SAL<-data %>% subset(Treatment=="WS")%>%group_by(Year) %>% filter(STL!="NA")%>%summarise(mean=mean(SFL+SML))

#linear model of SAL over time
SAL.lm.s <- lm(mean ~ Year, data = S.SAL)
summary(SAL.lm.s)

# Extract te coefficients from the overall model
SAL.coef.s <- coef(SAL.lm.s)
SAL.coef.s

p.3 <- ggplot(S.SAL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Scottnema total live")
p.3 

p.3  <- p.3  + geom_abline(intercept = SAL.coef.s[1], 
                           slope = SAL.coef.s[2], 
                           aes(colour = "overall"))
p.3 

SAL.lm3.s <- lm(mean ~ poly(Year, 3), data = S.SAL)

p.3  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.s <- segmented(SAL.lm.s, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(2001)))
summary(my.seg.s)
# get the breakpoints
my.seg.s$psi
slope(my.seg.s)
# get the fitted data
my.fitted.s <- fitted(my.seg.s)
my.model.s <- data.frame(Year = S.SAL$Year, SAL = my.fitted.s)
ggplot(my.model.s, aes(x = Year, y = SAL)) + geom_line()

p.3 + geom_line(data = my.model.s, aes(x = Year, y = SAL), colour = "purple")

##add to original plot with controls

pSAL <- pSAL + geom_line(data = my.model.s, aes(x = Year, y = SAL), colour = "purple")
pSAL

#Add SAL for W water treatment
W.SAL<-data %>% subset(Treatment=="W")%>%group_by(Year) %>% filter(STL!="NA")%>%summarise(mean=mean(SFL+SML))

#linear model of SAL over time
SAL.lm.w <- lm(mean ~ Year, data = W.SAL)
summary(SAL.lm.w)

# Extract te coefficients from the overall model
SAL.coef.w <- coef(SAL.lm.w)
SAL.coef.w

p.4 <- ggplot(W.SAL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Scottnema total live")
p.4 

p.4  <- p.4  + geom_abline(intercept = SAL.coef.w[1], 
                           slope = SAL.coef.w[2], 
                           aes(colour = "overall"))
p.4 

SAL.lm3.w <- lm(mean ~ poly(Year, 3), data = W.SAL)

p.4  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.w <- segmented(SAL.lm.w, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(2001)))
summary(my.seg.w)
# get the breakpoints
my.seg.w$psi
slope(my.seg.w)
# get the fitted data
my.fitted.w <- fitted(my.seg.w)
my.model.w <- data.frame(Year = W.SAL$Year, SAL = my.fitted.w)
ggplot(my.model.w, aes(x = Year, y = SAL)) + geom_line()

p.4 + geom_line(data = my.model.w, aes(x = Year, y = SAL), colour = "purple")

##add to original plot with controls

pSAL <- pSAL + geom_line(data = my.model.w, aes(x = Year, y = SAL), colour = "blue")
pSAL

#Add SAL for T temp treatment
T.SAL<-data %>% subset(Treatment=="T")%>%group_by(Year) %>% filter(STL!="NA")%>%summarise(mean=mean(SFL+SML))

#linear model of SAL over time
SAL.lm.t <- lm(mean ~ Year, data = T.SAL)
summary(SAL.lm.t)

# Extract the coefficients from the overall model
SAL.coef.t <- coef(SAL.lm.t)
SAL.coef.t

p.5 <- ggplot(T.SAL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Scottnema total live")
p.5 

p.5  <- p.5  + geom_abline(intercept = SAL.coef.t[1], 
                           slope = SAL.coef.t[2], 
                           aes(colour = "overall"))
p.5 

SAL.lm3.t <- lm(mean ~ poly(Year, 3), data = T.SAL)

p.5  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.t <- segmented(SAL.lm.t, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(2001)))
summary(my.seg.t)
# get the breakpoints
my.seg.t$psi
slope(my.seg.t)
# get the fitted data
my.fitted.t <- fitted(my.seg.t)
my.model.t <- data.frame(Year = T.SAL$Year, SAL = my.fitted.t)
ggplot(my.model.t, aes(x = Year, y = SAL)) + geom_line()

p.5 + geom_line(data = my.model.t, aes(x = Year, y = SAL), colour = "red")

##add to original plot with controls

pSAL <- pSAL + geom_line(data = my.model.t, aes(x = Year, y = SAL), colour = "red")
pSAL

#Add SAL for TW temp/water treatment
TW.SAL<-data %>% subset(Treatment=="TW")%>%group_by(Year) %>% filter(STL!="NA")%>%summarise(mean=mean(SFL+SML))

#linear model of SAL over time
SAL.lm.tw <- lm(mean ~ Year, data = TW.SAL)
summary(SAL.lm.tw)

# Extract the coefficients from the overall model
SAL.coef.tw <- coef(SAL.lm.tw)
SAL.coef.tw

p.6 <- ggplot(TW.SAL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Scottnema total live")
p.6 

p.6  <- p.6  + geom_abline(intercept = SAL.coef.tw[1], 
                           slope = SAL.coef.tw[2], 
                           aes(colour = "overall"))
p.6 

SAL.lm3.tw <- lm(mean ~ poly(Year, 3), data = TW.SAL)

p.6  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.tw <- segmented(SAL.lm.tw, 
                       seg.Z = ~ Year, 
                       psi = list(Year = c(2001)))
summary(my.seg.tw)
# get the breakpoints
my.seg.tw$psi
slope(my.seg.tw)
# get the fitted data
my.fitted.tw <- fitted(my.seg.tw)
my.model.tw <- data.frame(Year = TW.SAL$Year, SAL = my.fitted.tw)
ggplot(my.model.tw, aes(x = Year, y = SAL)) + geom_line()

p.6 + geom_line(data = my.model.tw, aes(x = Year, y = SAL), colour = "orangered1")

##add to original plot with controls

pSAL <- pSAL + geom_line(data = my.model.tw, aes(x = Year, y = SAL), colour = "orangered1")
pSAL

#Add SAL for TWM temp/water/mannitol treatment
TWM.SAL<-data %>% subset(Treatment=="TWM")%>%group_by(Year) %>% filter(STL!="NA")%>%summarise(mean=mean(SFL+SML))

#linear model of SAL over time
SAL.lm.twm <- lm(mean ~ Year, data = TWM.SAL)
summary(SAL.lm.twm)

# Extract the coefficients from the overall model
SAL.coef.twm <- coef(SAL.lm.twm)
SAL.coef.twm

p.7 <- ggplot(TWM.SAL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Scottnema total live")
p.7

p.7  <- p.7  + geom_abline(intercept = SAL.coef.twm[1], 
                           slope = SAL.coef.twm[2], 
                           aes(colour = "overall"))
p.7 

SAL.lm3.twm <- lm(mean ~ poly(Year, 3), data = TWM.SAL)

p.7  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.twm <- segmented(SAL.lm.twm, 
                        seg.Z = ~ Year, 
                        psi = list(Year = c(2001)))
summary(my.seg.twm)
# get the breakpoints
my.seg.twm$psi
slope(my.seg.twm)
# get the fitted data
my.fitted.twm <- fitted(my.seg.twm)
my.model.twm <- data.frame(Year = TWM.SAL$Year, SAL = my.fitted.twm)
ggplot(my.model.twm, aes(x = Year, y = SAL)) + geom_line()

p.7 + geom_line(data = my.model.twm, aes(x = Year, y = SAL), colour = "orange")

##add to original plot with controls

pSAL <- pSAL + geom_line(data = my.model.twm, aes(x = Year, y = SAL), colour = "orange")
pSAL

#Add SAL for TWM temp/water/mannitol treatment
TWS.SAL<-data %>% subset(Treatment=="TWS")%>%group_by(Year) %>% filter(STL!="NA")%>%summarise(mean=mean(SFL+SML))

#linear model of SAL over time
SAL.lm.tws <- lm(mean ~ Year, data = TWS.SAL)
summary(SAL.lm.tws)

# Extract the coefficients from the overall model
SAL.coef.tws <- coef(SAL.lm.tws)
SAL.coef.tws

p.8 <- ggplot(TWS.SAL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Scottnema total live")
p.8 

p.8   <- p.8  + geom_abline(intercept = SAL.coef.tws[1], 
                            slope = SAL.coef.tws[2], 
                            aes(colour = "overall"))
p.8 

SAL.lm3.tws <- lm(mean ~ poly(Year, 3), data = TWS.SAL)

p.8  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.tws <- segmented(SAL.lm.tws, 
                        seg.Z = ~ Year, 
                        psi = list(Year = c(2001)))
summary(my.seg.tws)
# get the breakpoints
my.seg.tws$psi
slope(my.seg.tws)
# get the fitted data
my.fitted.tws <- fitted(my.seg.tws)
my.model.tws <- data.frame(Year = TWS.SAL$Year, SAL = my.fitted.tws)
ggplot(my.model.tws, aes(x = Year, y = SAL)) + geom_line()

p.8 + geom_line(data = my.model.tws, aes(x = Year, y = SAL), colour = "firebrick")

##add to original plot with controls

pSAL <- pSAL + geom_line(data = my.model.tws, aes(x = Year, y = SAL), colour = "firebrick")
pSAL

##put scottnema plots together
legend<-ggplot(subset(data, STL!="NA"), aes(y=STL, x=Year, color=Treatment))+
  scale_color_manual(labels=c("U","T", "TW", "TWM", "TWS","W","WM","WS"),values=c("C"="black", "T"="red", "TWM"="orange", "TW"="orangered1", "TWS"="firebrick", "W"="blue", "WM"="navy", "WS"="purple"))+
  geom_point()
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(legend)
  
grid.arrange(pSTL,pSAL,pSJL,mylegend, ncol=4)


##breakpoint analysis for Eudorylaimus
C.ETL<-data %>% subset(Treatment=="C")%>%group_by(Year) %>% filter(ETL!="NA")%>%summarise(mean=mean(ETL))

#linear model of ETL over time
ETL.lm <- lm(mean ~ Year, data = C.ETL)
summary(ETL.lm)

# Extract te coefficients from the overall model
ETL.coef <- coef(ETL.lm)
ETL.coef

p <- ggplot(C.ETL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Eudorylaimus total live")
p

p <- p + geom_abline(intercept = ETL.coef[1], 
                     slope = ETL.coef[2], 
                     aes(colour = "overall"))
p

ETL.lm3 <- lm(mean ~ poly(Year, 3), data = C.ETL)

p + geom_smooth(method = "lm",
                formula = y ~ poly(x, degree = 3), 
                se = FALSE, colour = "orange")
my.seg <- segmented(ETL.lm, 
                    seg.Z = ~ Year, 
                    psi = list(Year = c(2004)))
summary(my.seg)
# get the breakpoints
my.seg$psi
slope(my.seg)
# get the fitted data
my.fitted <- fitted(my.seg)
my.model <- data.frame(Year = C.ETL$Year, ETL = my.fitted)
ggplot(my.model, aes(x = Year, y = ETL)) + geom_line()

pETL<- ggplot(C.ETL, aes(x = Year, y = mean)) + 
  geom_point(alpha=0.3)+
  geom_segment(aes(x=1993, y=114, xend=2004.5, yend=114), color = "grey0") + 
  geom_point(aes(x=1993, y=114), color="grey0")+
  geom_point(aes(x=2004.5, y=114), color="grey0")+
  geom_segment(aes(x=2005, y=114, xend=2015, yend=114), color = "grey60") + 
  geom_point(aes(x=2005, y=114), color="grey60")+
  geom_point(aes(x=2015, y=114), color="grey60")+
  geom_segment(aes(x=2001, y=9, xend=2001, yend=0), arrow=arrow(length = unit(0.03, "npc")), color = "grey0") + 
  annotate("text", x = c(2001), y=(12), label = c("Flood year"), fontface = 3,size=4)+
  annotate("text", x = c(1999), y=(120), label = c("Treatment"), fontface = 2,size=4)+
  annotate("text", x = c(2010), y=(120), label = c("Post-Treatment"), fontface = 2,size=4, color="grey60")+
  labs(x = "Year",y = "Eudorylaimus (#/kg dry soil)")+ 
  geom_line(data = my.model, aes(x = Year, y = ETL), colour = "black", cex=1.5)+
  ylim(0,120)
pETL

#Add ETL for WM mannitol treatment
M.ETL<-data %>% subset(Treatment=="WM")%>%group_by(Year) %>% filter(ETL!="NA")%>%summarise(mean=mean(ETL))

#linear model of ETL over time
ETL.lm.m <- lm(mean ~ Year, data = M.ETL)
summary(ETL.lm.m)

# Extract te coefficients from the overall model
ETL.coef.m <- coef(ETL.lm.m)
ETL.coef.m

p.2 <- ggplot(M.ETL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Eudorylaimus total live")
p.2 

p.2  <- p.2  + geom_abline(intercept = ETL.coef.m[1], 
                           slope = ETL.coef.m[2], 
                           aes(colour = "overall"))
p.2 

ETL.lm3 <- lm(mean ~ poly(Year, 3), data = M.ETL)

p.2  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.m <- segmented(ETL.lm.m, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(2009)))
summary(my.seg.m)
# get the breakpoints
my.seg.m$psi
slope(my.seg.m)
# get the fitted data
my.fitted.m <- fitted(my.seg.m)
my.model.m <- data.frame(Year = M.ETL$Year, ETL = my.fitted.m)
ggplot(my.model.m, aes(x = Year, y = ETL)) + geom_line()

p.2 + geom_line(data = my.model.m, aes(x = Year, y = ETL), colour = "tomato")

##add to original plot with controls

pETL <- pETL + geom_line(data = my.model.m, aes(x = Year, y = ETL),colour = "steelblue4", cex=1.5)
pETL <- pETL +  geom_point(data=M.ETL, aes(x = Year, y = mean), color="steelblue4", alpha=0.3)
pETL

#Add ETL for WS sucrose treatment
S.ETL<-data %>% subset(Treatment=="WS")%>%group_by(Year) %>% filter(ETL!="NA")%>%summarise(mean=mean(ETL))

#linear model of ETL over time
ETL.lm.s <- lm(mean ~ Year, data = S.ETL)
summary(ETL.lm.s)

# Extract te coefficients from the overall model
ETL.coef.s <- coef(ETL.lm.s)
ETL.coef.s

p.3 <- ggplot(S.ETL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Eudorylaimus total live")
p.3 

p.3  <- p.3  + geom_abline(intercept = ETL.coef.s[1], 
                           slope = ETL.coef.s[2], 
                           aes(colour = "overall"))
p.3 

ETL.lm3.s <- lm(mean ~ poly(Year, 3), data = S.ETL)

p.3  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.s <- segmented(ETL.lm.s, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(2001)))
summary(my.seg.s)
# get the breakpoints
my.seg.s$psi
slope(my.seg.s)
# get the fitted data
my.fitted.s <- fitted(my.seg.s)
my.model.s <- data.frame(Year = S.ETL$Year, ETL = my.fitted.s)
ggplot(my.model.s, aes(x = Year, y = ETL)) + geom_line()

p.3 + geom_line(data = my.model.s, aes(x = Year, y = ETL), colour = "purple")

##add to original plot with controls
pETL <- pETL + geom_line(data = my.model.s, aes(x = Year, y = ETL),colour = "orchid3", cex=1.5)
pETL <- pETL +  geom_point(data=S.ETL, aes(x = Year, y = mean), color="orchid3", alpha=0.3)
pETL

#Add ETL for W water treatment
W.ETL<-data %>% subset(Treatment=="W")%>%group_by(Year) %>% filter(ETL!="NA")%>%summarise(mean=mean(ETL))

#linear model of ETL over time
ETL.lm.w <- lm(mean ~ Year, data = W.ETL)
summary(ETL.lm.w)

# Extract te coefficients from the overall model
ETL.coef.w <- coef(ETL.lm.w)
ETL.coef.w

p.4 <- ggplot(W.ETL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Eudorylaimus total live")
p.4 

p.4  <- p.4  + geom_abline(intercept = ETL.coef.w[1], 
                           slope = ETL.coef.w[2], 
                           aes(colour = "overall"))
p.4 

ETL.lm3.w <- lm(mean ~ poly(Year, 3), data = W.ETL)

p.4  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.w <- segmented(ETL.lm.w, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(1995)))
summary(my.seg.w)
# get the breakpoints
my.seg.w$psi
slope(my.seg.w)
# get the fitted data
my.fitted.w <- fitted(my.seg.w)
my.model.w <- data.frame(Year = W.ETL$Year, ETL = my.fitted.w)
ggplot(my.model.w, aes(x = Year, y = ETL)) + geom_line()

p.4 + geom_line(data = my.model.w, aes(x = Year, y = ETL), colour = "purple")

##add to original plot with controls
pETL <- pETL + geom_line(data = my.model.w, aes(x = Year, y = ETL),colour = "deepskyblue", cex=1.5)
pETL <- pETL +  geom_point(data=W.ETL, aes(x = Year, y = mean), color="deepskyblue", alpha=0.3)
pETL
pETL1<-pETL

#Add ETL for T temp treatment
T.ETL<-data %>% subset(Treatment=="T")%>%group_by(Year) %>% filter(ETL!="NA")%>%summarise(mean=mean(ETL))

#linear model of ETL over time
ETL.lm.t <- lm(mean ~ Year, data = T.ETL)
summary(ETL.lm.t)

# Extract the coefficients from the overall model
ETL.coef.t <- coef(ETL.lm.t)
ETL.coef.t

p.5 <- ggplot(T.ETL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Eudorylaimus total live")
p.5 

p.5  <- p.5  + geom_abline(intercept = ETL.coef.t[1], 
                           slope = ETL.coef.t[2], 
                           aes(colour = "overall"))
p.5 

ETL.lm3.t <- lm(mean ~ poly(Year, 3), data = T.ETL)

p.5  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.t <- segmented(ETL.lm.t, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(2001)))
summary(my.seg.t)
# get the breakpoints
my.seg.t$psi
slope(my.seg.t)
# get the fitted data
my.fitted.t <- fitted(my.seg.t)
my.model.t <- data.frame(Year = T.ETL$Year, ETL = my.fitted.t)
ggplot(my.model.t, aes(x = Year, y = ETL)) + geom_line()

p.5 + geom_line(data = my.model.t, aes(x = Year, y = ETL), colour = "red")

##add to original plot with controls

pETL <- pETL + geom_line(data = my.model.t, aes(x = Year, y = ETL), colour = "red")
pETL

#Add ETL for TW temp/water treatment
TW.ETL<-data %>% subset(Treatment=="TW")%>%group_by(Year) %>% filter(ETL!="NA")%>%summarise(mean=mean(ETL))

#linear model of ETL over time
ETL.lm.tw <- lm(mean ~ Year, data = TW.ETL)
summary(ETL.lm.tw)

# Extract the coefficients from the overall model
ETL.coef.tw <- coef(ETL.lm.tw)
ETL.coef.tw

p.6 <- ggplot(TW.ETL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Eudorylaimus total live")
p.6 

p.6  <- p.6  + geom_abline(intercept = ETL.coef.tw[1], 
                           slope = ETL.coef.tw[2], 
                           aes(colour = "overall"))
p.6 

ETL.lm3.tw <- lm(mean ~ poly(Year, 3), data = TW.ETL)

p.6  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.tw <- segmented(ETL.lm.tw, 
                       seg.Z = ~ Year, 
                       psi = list(Year = c(2001)))
summary(my.seg.tw)
# get the breakpoints
my.seg.tw$psi
slope(my.seg.tw)
# get the fitted data
my.fitted.tw <- fitted(my.seg.tw)
my.model.tw <- data.frame(Year = TW.ETL$Year, ETL = my.fitted.tw)
ggplot(my.model.tw, aes(x = Year, y = ETL)) + geom_line()

p.6 + geom_line(data = my.model.tw, aes(x = Year, y = ETL), colour = "orangered1")

##add to original plot with controls

pETL <- pETL + geom_line(data = my.model.tw, aes(x = Year, y = ETL), colour = "orangered1")
pETL

#Add ETL for TWM temp/water/mannitol treatment
TWM.ETL<-data %>% subset(Treatment=="TWM")%>%group_by(Year) %>% filter(ETL!="NA")%>%summarise(mean=mean(ETL))

#linear model of ETL over time
ETL.lm.twm <- lm(mean ~ Year, data = TWM.ETL)
summary(ETL.lm.twm)

# Extract the coefficients from the overall model
ETL.coef.twm <- coef(ETL.lm.twm)
ETL.coef.twm

p.7 <- ggplot(TWM.ETL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Eudorylaimus total live")
p.7

p.7  <- p.7  + geom_abline(intercept = ETL.coef.twm[1], 
                           slope = ETL.coef.twm[2], 
                           aes(colour = "overall"))
p.7 

ETL.lm3.twm <- lm(mean ~ poly(Year, 3), data = TWM.ETL)

p.7  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.twm <- segmented(ETL.lm.twm, 
                        seg.Z = ~ Year, 
                        psi = list(Year = c(2001)))
summary(my.seg.twm)
# get the breakpoints
my.seg.twm$psi
slope(my.seg.twm)
# get the fitted data
my.fitted.twm <- fitted(my.seg.twm)
my.model.twm <- data.frame(Year = TWM.ETL$Year, ETL = my.fitted.twm)
ggplot(my.model.twm, aes(x = Year, y = ETL)) + geom_line()

p.7 + geom_line(data = my.model.twm, aes(x = Year, y = ETL), colour = "orange")

##add to original plot with controls

pETL <- pETL + geom_line(data = my.model.twm, aes(x = Year, y = ETL), colour = "orange")
pETL

#Add ETL for TWM temp/water/mannitol treatment
TWS.ETL<-data %>% subset(Treatment=="TWS")%>%group_by(Year) %>% filter(ETL!="NA")%>%summarise(mean=mean(ETL))

#linear model of ETL over time
ETL.lm.tws <- lm(mean ~ Year, data = TWS.ETL)
summary(ETL.lm.tws)

# Extract the coefficients from the overall model
ETL.coef.tws <- coef(ETL.lm.tws)
ETL.coef.tws

p.8 <- ggplot(TWS.ETL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Eudorylaimus total live")
p.8 

p.8   <- p.8  + geom_abline(intercept = ETL.coef.tws[1], 
                            slope = ETL.coef.tws[2], 
                            aes(colour = "overall"))
p.8 

ETL.lm3.tws <- lm(mean ~ poly(Year, 3), data = TWS.ETL)

p.8  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.tws <- segmented(ETL.lm.tws, 
                        seg.Z = ~ Year, 
                        psi = list(Year = c(2001)))
summary(my.seg.tws)
# get the breakpoints
my.seg.tws$psi
slope(my.seg.tws)
# get the fitted data
my.fitted.tws <- fitted(my.seg.tws)
my.model.tws <- data.frame(Year = TWS.ETL$Year, ETL = my.fitted.tws)
ggplot(my.model.tws, aes(x = Year, y = ETL)) + geom_line()

p.8 + geom_line(data = my.model.tws, aes(x = Year, y = ETL), colour = "firebrick")

##add to original plot with controls

pETL <- pETL + geom_line(data = my.model.tws, aes(x = Year, y = ETL), colour = "firebrick")
pETL

################
#Eudorylaimus juveniles

C.EJL<-data %>% subset(Treatment=="C")%>%group_by(Year) %>% filter(EJL!="NA")%>%summarise(mean=mean(EJL))

#linear model of EJL over time
EJL.lm <- lm(mean ~ Year, data = C.EJL)
summary(EJL.lm)

# Extract te coefficients from the overall model
EJL.coef <- coef(EJL.lm)
EJL.coef

p <- ggplot(C.EJL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Eudorylaimus total live")
p

p <- p + geom_abline(intercept = EJL.coef[1], 
                     slope = EJL.coef[2], 
                     aes(colour = "overall"))
p

EJL.lm3 <- lm(mean ~ poly(Year, 3), data = C.EJL)

p + geom_smooth(method = "lm",
                formula = y ~ poly(x, degree = 3), 
                se = FALSE, colour = "orange")
my.seg <- segmented(EJL.lm, 
                    seg.Z = ~ Year, 
                    psi = list(Year = c(2004)))
summary(my.seg)
# get the breakpoints
my.seg$psi
slope(my.seg)
# get the fitted data
my.fitted <- fitted(my.seg)
my.model <- data.frame(Year = C.EJL$Year, EJL = my.fitted)
ggplot(my.model, aes(x = Year, y = EJL)) + geom_line()

pEJL<- ggplot(C.EJL, aes(x = Year, y = mean)) + 
  labs(x = "Year",y = "Eudorylaimus juveniles live")+ 
  geom_line(data = my.model, aes(x = Year, y = EJL), colour = "black")+
  ylim(0,100)
pEJL

#Add EJL for WM mannitol treatment
M.EJL<-data %>% subset(Treatment=="WM")%>%group_by(Year) %>% filter(EJL!="NA")%>%summarise(mean=mean(EJL))

#linear model of EJL over time
EJL.lm.m <- lm(mean ~ Year, data = M.EJL)
summary(EJL.lm.m)

# Extract te coefficients from the overall model
EJL.coef.m <- coef(EJL.lm.m)
EJL.coef.m

p.2 <- ggplot(M.EJL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Eudorylaimus total live")
p.2 

p.2  <- p.2  + geom_abline(intercept = EJL.coef.m[1], 
                           slope = EJL.coef.m[2], 
                           aes(colour = "overall"))
p.2 

EJL.lm3 <- lm(mean ~ poly(Year, 3), data = M.EJL)

p.2  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.m <- segmented(EJL.lm.m, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(2004)))
summary(my.seg.m)
# get the breakpoints
my.seg.m$psi
slope(my.seg.m)
# get the fitted data
my.fitted.m <- fitted(my.seg.m)
my.model.m <- data.frame(Year = M.EJL$Year, EJL = my.fitted.m)
ggplot(my.model.m, aes(x = Year, y = EJL)) + geom_line()

p.2 + geom_line(data = my.model.m, aes(x = Year, y = EJL), colour = "tomato")

##add to original plot with controls

pEJL <- pEJL + geom_line(data = my.model.m, aes(x = Year, y = EJL), colour = "navy")
pEJL

#Add EJL for WS sucrose treatment
S.EJL<-data %>% subset(Treatment=="WS")%>%group_by(Year) %>% filter(EJL!="NA")%>%summarise(mean=mean(EJL))

#linear model of EJL over time
EJL.lm.s <- lm(mean ~ Year, data = S.EJL)
summary(EJL.lm.s)

# Extract te coefficients from the overall model
EJL.coef.s <- coef(EJL.lm.s)
EJL.coef.s

p.3 <- ggplot(S.EJL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Eudorylaimus total live")
p.3 

p.3  <- p.3  + geom_abline(intercept = EJL.coef.s[1], 
                           slope = EJL.coef.s[2], 
                           aes(colour = "overall"))
p.3 

EJL.lm3.s <- lm(mean ~ poly(Year, 3), data = S.EJL)

p.3  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.s <- segmented(EJL.lm.s, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(2001)))
summary(my.seg.s)
# get the breakpoints
my.seg.s$psi
slope(my.seg.s)
# get the fitted data
my.fitted.s <- fitted(my.seg.s)
my.model.s <- data.frame(Year = S.EJL$Year, EJL = my.fitted.s)
ggplot(my.model.s, aes(x = Year, y = EJL)) + geom_line()

p.3 + geom_line(data = my.model.s, aes(x = Year, y = EJL), colour = "purple")

##add to original plot with controls

pEJL <- pEJL + geom_line(data = my.model.s, aes(x = Year, y = EJL), colour = "purple")
pEJL

#Add EJL for W water treatment
W.EJL<-data %>% subset(Treatment=="W")%>%group_by(Year) %>% filter(EJL!="NA")%>%summarise(mean=mean(EJL))

#linear model of EJL over time
EJL.lm.w <- lm(mean ~ Year, data = W.EJL)
summary(EJL.lm.w)

# Extract te coefficients from the overall model
EJL.coef.w <- coef(EJL.lm.w)
EJL.coef.w

p.4 <- ggplot(W.EJL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Eudorylaimus total live")
p.4 

p.4  <- p.4  + geom_abline(intercept = EJL.coef.w[1], 
                           slope = EJL.coef.w[2], 
                           aes(colour = "overall"))
p.4 

EJL.lm3.w <- lm(mean ~ poly(Year, 3), data = W.EJL)

p.4  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.w <- segmented(EJL.lm.w, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(2001)))
summary(my.seg.w)
# get the breakpoints
my.seg.w$psi
slope(my.seg.w)
# get the fitted data
my.fitted.w <- fitted(my.seg.w)
my.model.w <- data.frame(Year = W.EJL$Year, EJL = my.fitted.w)
ggplot(my.model.w, aes(x = Year, y = EJL)) + geom_line()

p.4 + geom_line(data = my.model.w, aes(x = Year, y = EJL), colour = "purple")

##add to original plot with controls

pEJL <- pEJL + geom_line(data = my.model.w, aes(x = Year, y = EJL), colour = "blue")
pEJL

#Add EJL for T temp treatment
T.EJL<-data %>% subset(Treatment=="T")%>%group_by(Year) %>% filter(EJL!="NA")%>%summarise(mean=mean(EJL))

#linear model of EJL over time
EJL.lm.t <- lm(mean ~ Year, data = T.EJL)
summary(EJL.lm.t)

# Extract the coefficients from the overall model
EJL.coef.t <- coef(EJL.lm.t)
EJL.coef.t

p.5 <- ggplot(T.EJL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Eudorylaimus total live")
p.5 

p.5  <- p.5  + geom_abline(intercept = EJL.coef.t[1], 
                           slope = EJL.coef.t[2], 
                           aes(colour = "overall"))
p.5 

EJL.lm3.t <- lm(mean ~ poly(Year, 3), data = T.EJL)

p.5  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.t <- segmented(EJL.lm.t, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(2001)))
summary(my.seg.t)
# get the breakpoints
my.seg.t$psi
slope(my.seg.t)
# get the fitted data
my.fitted.t <- fitted(my.seg.t)
my.model.t <- data.frame(Year = T.EJL$Year, EJL = my.fitted.t)
ggplot(my.model.t, aes(x = Year, y = EJL)) + geom_line()

p.5 + geom_line(data = my.model.t, aes(x = Year, y = EJL), colour = "red")

##add to original plot with controls

pEJL <- pEJL + geom_line(data = my.model.t, aes(x = Year, y = EJL), colour = "red")
pEJL

#Add EJL for TW temp/water treatment
TW.EJL<-data %>% subset(Treatment=="TW")%>%group_by(Year) %>% filter(EJL!="NA")%>%summarise(mean=mean(EJL))

#linear model of EJL over time
EJL.lm.tw <- lm(mean ~ Year, data = TW.EJL)
summary(EJL.lm.tw)

# Extract the coefficients from the overall model
EJL.coef.tw <- coef(EJL.lm.tw)
EJL.coef.tw

p.6 <- ggplot(TW.EJL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Eudorylaimus total live")
p.6 

p.6  <- p.6  + geom_abline(intercept = EJL.coef.tw[1], 
                           slope = EJL.coef.tw[2], 
                           aes(colour = "overall"))
p.6 

EJL.lm3.tw <- lm(mean ~ poly(Year, 3), data = TW.EJL)

p.6  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.tw <- segmented(EJL.lm.tw, 
                       seg.Z = ~ Year, 
                       psi = list(Year = c(2001)))
summary(my.seg.tw)
# get the breakpoints
my.seg.tw$psi
slope(my.seg.tw)
# get the fitted data
my.fitted.tw <- fitted(my.seg.tw)
my.model.tw <- data.frame(Year = TW.EJL$Year, EJL = my.fitted.tw)
ggplot(my.model.tw, aes(x = Year, y = EJL)) + geom_line()

p.6 + geom_line(data = my.model.tw, aes(x = Year, y = EJL), colour = "orangered1")

##add to original plot with controls

pEJL <- pEJL + geom_line(data = my.model.tw, aes(x = Year, y = EJL), colour = "orangered1")
pEJL

#Add EJL for TWM temp/water/mannitol treatment
TWM.EJL<-data %>% subset(Treatment=="TWM")%>%group_by(Year) %>% filter(EJL!="NA")%>%summarise(mean=mean(EJL))

#linear model of EJL over time
EJL.lm.twm <- lm(mean ~ Year, data = TWM.EJL)
summary(EJL.lm.twm)

# Extract the coefficients from the overall model
EJL.coef.twm <- coef(EJL.lm.twm)
EJL.coef.twm

p.7 <- ggplot(TWM.EJL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Eudorylaimus total live")
p.7

p.7  <- p.7  + geom_abline(intercept = EJL.coef.twm[1], 
                           slope = EJL.coef.twm[2], 
                           aes(colour = "overall"))
p.7 

EJL.lm3.twm <- lm(mean ~ poly(Year, 3), data = TWM.EJL)

p.7  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.twm <- segmented(EJL.lm.twm, 
                        seg.Z = ~ Year, 
                        psi = list(Year = c(2001)))
summary(my.seg.twm)
# get the breakpoints
my.seg.twm$psi
slope(my.seg.twm)
# get the fitted data
my.fitted.twm <- fitted(my.seg.twm)
my.model.twm <- data.frame(Year = TWM.EJL$Year, EJL = my.fitted.twm)
ggplot(my.model.twm, aes(x = Year, y = EJL)) + geom_line()

p.7 + geom_line(data = my.model.twm, aes(x = Year, y = EJL), colour = "orange")

##add to original plot with controls

pEJL <- pEJL + geom_line(data = my.model.twm, aes(x = Year, y = EJL), colour = "orange")
pEJL

#Add EJL for TWM temp/water/mannitol treatment
TWS.EJL<-data %>% subset(Treatment=="TWS")%>%group_by(Year) %>% filter(EJL!="NA")%>%summarise(mean=mean(EJL))

#linear model of EJL over time
EJL.lm.tws <- lm(mean ~ Year, data = TWS.EJL)
summary(EJL.lm.tws)

# Extract the coefficients from the overall model
EJL.coef.tws <- coef(EJL.lm.tws)
EJL.coef.tws

p.8 <- ggplot(TWS.EJL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Eudorylaimus total live")
p.8 

p.8   <- p.8  + geom_abline(intercept = EJL.coef.tws[1], 
                            slope = EJL.coef.tws[2], 
                            aes(colour = "overall"))
p.8 

EJL.lm3.tws <- lm(mean ~ poly(Year, 3), data = TWS.EJL)

p.8  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.tws <- segmented(EJL.lm.tws, 
                        seg.Z = ~ Year, 
                        psi = list(Year = c(2001)))
summary(my.seg.tws)
# get the breakpoints
my.seg.tws$psi
slope(my.seg.tws)
# get the fitted data
my.fitted.tws <- fitted(my.seg.tws)
my.model.tws <- data.frame(Year = TWS.EJL$Year, EJL = my.fitted.tws)
ggplot(my.model.tws, aes(x = Year, y = EJL)) + geom_line()

p.8 + geom_line(data = my.model.tws, aes(x = Year, y = EJL), colour = "firebrick")

##add to original plot with controls

pEJL <- pEJL + geom_line(data = my.model.tws, aes(x = Year, y = EJL), colour = "firebrick")
pEJL

################
#Eudorylaimus juveniles

C.EAL<-data %>% subset(Treatment=="C")%>%group_by(Year) %>% filter(ETL!="NA")%>%summarise(mean=mean(EML+EFL))

#linear model of EAL over time
EAL.lm <- lm(mean ~ Year, data = C.EAL)
summary(EAL.lm)

# Extract te coefficients from the overall model
EAL.coef <- coef(EAL.lm)
EAL.coef

p <- ggplot(C.EAL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Eudorylaimus total live")
p

p <- p + geom_abline(intercept = EAL.coef[1], 
                     slope = EAL.coef[2], 
                     aes(colour = "overall"))
p

EAL.lm3 <- lm(mean ~ poly(Year, 3), data = C.EAL)

p + geom_smooth(method = "lm",
                formula = y ~ poly(x, degree = 3), 
                se = FALSE, colour = "orange")
my.seg <- segmented(EAL.lm, 
                    seg.Z = ~ Year, 
                    psi = list(Year = c(2004)))
summary(my.seg)
# get the breakpoints
my.seg$psi
slope(my.seg)
# get the fitted data
my.fitted <- fitted(my.seg)
my.model <- data.frame(Year = C.EAL$Year, EAL = my.fitted)
ggplot(my.model, aes(x = Year, y = EAL)) + geom_line()

pEAL<- ggplot(C.EAL, aes(x = Year, y = mean)) + 
  labs(x = "Year",y = "Eudorylaimus adults live")+ 
  geom_line(data = my.model, aes(x = Year, y = EAL), colour = "black")+
  ylim(0,100)
pEAL

#Add EAL for WM mannitol treatment
M.EAL<-data %>% subset(Treatment=="WM")%>%group_by(Year) %>% filter(ETL!="NA")%>%summarise(mean=mean(EML+EFL))

#linear model of EAL over time
EAL.lm.m <- lm(mean ~ Year, data = M.EAL)
summary(EAL.lm.m)

# Extract te coefficients from the overall model
EAL.coef.m <- coef(EAL.lm.m)
EAL.coef.m

p.2 <- ggplot(M.EAL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Eudorylaimus total live")
p.2 

p.2  <- p.2  + geom_abline(intercept = EAL.coef.m[1], 
                           slope = EAL.coef.m[2], 
                           aes(colour = "overall"))
p.2 

EAL.lm3 <- lm(mean ~ poly(Year, 3), data = M.EAL)

p.2  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.m <- segmented(EAL.lm.m, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(2004)))
summary(my.seg.m)
# get the breakpoints
my.seg.m$psi
slope(my.seg.m)
# get the fitted data
my.fitted.m <- fitted(my.seg.m)
my.model.m <- data.frame(Year = M.EAL$Year, EAL = my.fitted.m)
ggplot(my.model.m, aes(x = Year, y = EAL)) + geom_line()

p.2 + geom_line(data = my.model.m, aes(x = Year, y = EAL), colour = "tomato")

##add to original plot with controls

pEAL <- pEAL + geom_line(data = my.model.m, aes(x = Year, y = EAL), colour = "navy")
pEAL

#Add EAL for WS sucrose treatment
S.EAL<-data %>% subset(Treatment=="WS")%>%group_by(Year) %>% filter(ETL!="NA")%>%summarise(mean=mean(EFL+EML))

#linear model of EAL over time
EAL.lm.s <- lm(mean ~ Year, data = S.EAL)
summary(EAL.lm.s)

# Extract te coefficients from the overall model
EAL.coef.s <- coef(EAL.lm.s)
EAL.coef.s

p.3 <- ggplot(S.EAL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Eudorylaimus total live")
p.3 

p.3  <- p.3  + geom_abline(intercept = EAL.coef.s[1], 
                           slope = EAL.coef.s[2], 
                           aes(colour = "overall"))
p.3 

EAL.lm3.s <- lm(mean ~ poly(Year, 3), data = S.EAL)

p.3  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.s <- segmented(EAL.lm.s, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(2001)))
summary(my.seg.s)
# get the breakpoints
my.seg.s$psi
slope(my.seg.s)
# get the fitted data
my.fitted.s <- fitted(my.seg.s)
my.model.s <- data.frame(Year = S.EAL$Year, EAL = my.fitted.s)
ggplot(my.model.s, aes(x = Year, y = EAL)) + geom_line()

p.3 + geom_line(data = my.model.s, aes(x = Year, y = EAL), colour = "purple")

##add to original plot with controls

pEAL <- pEAL + geom_line(data = my.model.s, aes(x = Year, y = EAL), colour = "purple")
pEAL

#Add EAL for W water treatment
W.EAL<-data %>% subset(Treatment=="W")%>%group_by(Year) %>% filter(ETL!="NA")%>%summarise(mean=mean(EFL+EML))

#linear model of EAL over time
EAL.lm.w <- lm(mean ~ Year, data = W.EAL)
summary(EAL.lm.w)

# Extract te coefficients from the overall model
EAL.coef.w <- coef(EAL.lm.w)
EAL.coef.w

p.4 <- ggplot(W.EAL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Eudorylaimus total live")
p.4 

p.4  <- p.4  + geom_abline(intercept = EAL.coef.w[1], 
                           slope = EAL.coef.w[2], 
                           aes(colour = "overall"))
p.4 

EAL.lm3.w <- lm(mean ~ poly(Year, 3), data = W.EAL)

p.4  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.w <- segmented(EAL.lm.w, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(2001)))
summary(my.seg.w)
# get the breakpoints
my.seg.w$psi
slope(my.seg.w)
# get the fitted data
my.fitted.w <- fitted(my.seg.w)
my.model.w <- data.frame(Year = W.EAL$Year, EAL = my.fitted.w)
ggplot(my.model.w, aes(x = Year, y = EAL)) + geom_line()

p.4 + geom_line(data = my.model.w, aes(x = Year, y = EAL), colour = "purple")

##add to original plot with controls

pEAL <- pEAL + geom_line(data = my.model.w, aes(x = Year, y = EAL), colour = "blue")
pEAL

#Add EAL for T temp treatment
T.EAL<-data %>% subset(Treatment=="T")%>%group_by(Year) %>% filter(ETL!="NA")%>%summarise(mean=mean(EFL+EML))

#linear model of EAL over time
EAL.lm.t <- lm(mean ~ Year, data = T.EAL)
summary(EAL.lm.t)

# Extract the coefficients from the overall model
EAL.coef.t <- coef(EAL.lm.t)
EAL.coef.t

p.5 <- ggplot(T.EAL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Eudorylaimus total live")
p.5 

p.5  <- p.5  + geom_abline(intercept = EAL.coef.t[1], 
                           slope = EAL.coef.t[2], 
                           aes(colour = "overall"))
p.5 

EAL.lm3.t <- lm(mean ~ poly(Year, 3), data = T.EAL)

p.5  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.t <- segmented(EAL.lm.t, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(2001)))
summary(my.seg.t)
# get the breakpoints
my.seg.t$psi
slope(my.seg.t)
# get the fitted data
my.fitted.t <- fitted(my.seg.t)
my.model.t <- data.frame(Year = T.EAL$Year, EAL = my.fitted.t)
ggplot(my.model.t, aes(x = Year, y = EAL)) + geom_line()

p.5 + geom_line(data = my.model.t, aes(x = Year, y = EAL), colour = "red")

##add to original plot with controls

pEAL <- pEAL + geom_line(data = my.model.t, aes(x = Year, y = EAL), colour = "red")
pEAL

#Add EAL for TW temp/water treatment
TW.EAL<-data %>% subset(Treatment=="TW")%>%group_by(Year) %>% filter(ETL!="NA")%>%summarise(mean=mean(EFL+EML))

#linear model of EAL over time
EAL.lm.tw <- lm(mean ~ Year, data = TW.EAL)
summary(EAL.lm.tw)

# Extract the coefficients from the overall model
EAL.coef.tw <- coef(EAL.lm.tw)
EAL.coef.tw

p.6 <- ggplot(TW.EAL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Eudorylaimus total live")
p.6 

p.6  <- p.6  + geom_abline(intercept = EAL.coef.tw[1], 
                           slope = EAL.coef.tw[2], 
                           aes(colour = "overall"))
p.6 

EAL.lm3.tw <- lm(mean ~ poly(Year, 3), data = TW.EAL)

p.6  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.tw <- segmented(EAL.lm.tw, 
                       seg.Z = ~ Year, 
                       psi = list(Year = c(2001)))
summary(my.seg.tw)
# get the breakpoints
my.seg.tw$psi
slope(my.seg.tw)
# get the fitted data
my.fitted.tw <- fitted(my.seg.tw)
my.model.tw <- data.frame(Year = TW.EAL$Year, EAL = my.fitted.tw)
ggplot(my.model.tw, aes(x = Year, y = EAL)) + geom_line()

p.6 + geom_line(data = my.model.tw, aes(x = Year, y = EAL), colour = "orangered1")

##add to original plot with controls

pEAL <- pEAL + geom_line(data = my.model.tw, aes(x = Year, y = EAL), colour = "orangered1")
pEAL

#Add EAL for TWM temp/water/mannitol treatment
TWM.EAL<-data %>% subset(Treatment=="TWM")%>%group_by(Year) %>% filter(ETL!="NA")%>%summarise(mean=mean(EFL+EML))

#linear model of EAL over time
EAL.lm.twm <- lm(mean ~ Year, data = TWM.EAL)
summary(EAL.lm.twm)

# Extract the coefficients from the overall model
EAL.coef.twm <- coef(EAL.lm.twm)
EAL.coef.twm

p.7 <- ggplot(TWM.EAL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Eudorylaimus total live")
p.7

p.7  <- p.7  + geom_abline(intercept = EAL.coef.twm[1], 
                           slope = EAL.coef.twm[2], 
                           aes(colour = "overall"))
p.7 

EAL.lm3.twm <- lm(mean ~ poly(Year, 3), data = TWM.EAL)

p.7  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.twm <- segmented(EAL.lm.twm, 
                        seg.Z = ~ Year, 
                        psi = list(Year = c(2001)))
summary(my.seg.twm)
# get the breakpoints
my.seg.twm$psi
slope(my.seg.twm)
# get the fitted data
my.fitted.twm <- fitted(my.seg.twm)
my.model.twm <- data.frame(Year = TWM.EAL$Year, EAL = my.fitted.twm)
ggplot(my.model.twm, aes(x = Year, y = EAL)) + geom_line()

p.7 + geom_line(data = my.model.twm, aes(x = Year, y = EAL), colour = "orange")

##add to original plot with controls

pEAL <- pEAL + geom_line(data = my.model.twm, aes(x = Year, y = EAL), colour = "orange")
pEAL

#Add EAL for TWM temp/water/mannitol treatment
TWS.EAL<-data %>% subset(Treatment=="TWS")%>%group_by(Year) %>% filter(ETL!="NA")%>%summarise(mean=mean(EFL+EML))

#linear model of EAL over time
EAL.lm.tws <- lm(mean ~ Year, data = TWS.EAL)
summary(EAL.lm.tws)

# Extract the coefficients from the overall model
EAL.coef.tws <- coef(EAL.lm.tws)
EAL.coef.tws

p.8 <- ggplot(TWS.EAL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Eudorylaimus total live")
p.8 

p.8   <- p.8  + geom_abline(intercept = EAL.coef.tws[1], 
                            slope = EAL.coef.tws[2], 
                            aes(colour = "overall"))
p.8 

EAL.lm3.tws <- lm(mean ~ poly(Year, 3), data = TWS.EAL)

p.8  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.tws <- segmented(EAL.lm.tws, 
                        seg.Z = ~ Year, 
                        psi = list(Year = c(2001)))
summary(my.seg.tws)
# get the breakpoints
my.seg.tws$psi
slope(my.seg.tws)
# get the fitted data
my.fitted.tws <- fitted(my.seg.tws)
my.model.tws <- data.frame(Year = TWS.EAL$Year, EAL = my.fitted.tws)
ggplot(my.model.tws, aes(x = Year, y = EAL)) + geom_line()

p.8 + geom_line(data = my.model.tws, aes(x = Year, y = EAL), colour = "firebrick")

##add to original plot with controls

pEAL <- pEAL + geom_line(data = my.model.tws, aes(x = Year, y = EAL), colour = "firebrick")
pEAL

##put Eudorylaimus plots together
legend<-ggplot(subset(data, ETL!="NA"), aes(y=ETL, x=Year, color=Treatment))+
  scale_color_manual(labels=c("U","T", "TW", "TWM", "TWS","W","WM","WS"),values=c("C"="black", "T"="red", "TWM"="orange", "TW"="orangered1", "TWS"="firebrick", "W"="blue", "WM"="navy", "WS"="purple"))+
  geom_point()
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(legend)

grid.arrange(pETL,pEAL,pEJL,mylegend, ncol=4)



#############
###SOC######
C.SOC<-data %>% subset(Treatment=="C")%>%group_by(Year) %>% filter(SOC2!="NA")%>%summarise(mean=mean(SOC2))

#linear model of SOC over time
SOC.lm <- lm(mean ~ Year, data = C.SOC)
summary(SOC.lm)

# Extract te coefficients from the overall model
SOC.coef <- coef(SOC.lm)
SOC.coef

p <- ggplot(C.SOC, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "SOC (mg/g)")
p

p <- p + geom_abline(intercept = SOC.coef[1], 
                     slope = SOC.coef[2], 
                     aes(colour = "overall"))
p

SOC.lm3 <- lm(mean ~ poly(Year, 3), data = C.SOC)

p + geom_smooth(method = "lm",
                formula = y ~ poly(x, degree = 3), 
                se = FALSE, colour = "orange")
my.seg <- segmented(SOC.lm, 
                    seg.Z = ~ Year, 
                    psi = list(Year = c(1994)))
summary(my.seg)
# get the breakpoints
my.seg$psi
slope(my.seg)
# get the fitted data
my.fitted <- fitted(my.seg)
my.model <- data.frame(Year = C.SOC$Year, SOC = my.fitted)
ggplot(my.model, aes(x = Year, y = SOC)) + geom_line()

pSOC<- ggplot(C.SOC, aes(x = Year, y = mean)) + 
  geom_point(color="black", alpha=0.3)+
  labs(x = "Year",y = "Soil Organic Carbon (mg/g)")+ 
  geom_segment(aes(x=1993, y=0.7125, xend=2004.5, yend=0.7125), color = "grey0") + 
  geom_point(aes(x=1993, y=0.7125), color="grey0")+
  geom_point(aes(x=2004.5, y=0.7125), color="grey0")+
  geom_segment(aes(x=2005, y=0.7125, xend=2015, yend=0.7125), color = "grey60")+  
  geom_point(aes(x=2005, y=0.7125), color="grey60")+
  geom_point(aes(x=2015, y=0.7125), color="grey60")+
  geom_segment(aes(x=2001, y=0.05625, xend=2001, yend=0), arrow=arrow(length = unit(0.03, "npc")), color = "grey0") + 
  annotate("text", x = c(2001), y=(0.075), label = c("Flood year"), fontface = 3,size=4)+
  annotate("text", x = c(1999), y=(0.75), label = c("Treatment"), fontface = 2,size=4)+
  annotate("text", x = c(2010), y=(0.75), label = c("Post-Treatment"), fontface = 2,size=4, color="grey60")+
  geom_line(data = my.model, aes(x = Year, y = SOC), colour = "black", cex=1.5)+
  ylim(0,0.75)
pSOC

#Add SOC for WM mannitol treatment
M.SOC<-data %>% subset(Treatment=="WM")%>%group_by(Year) %>% filter(SOC!="NA")%>%summarise(mean=mean(SOC))

#linear model of SOC over time
SOC.lm.m <- lm(mean ~ Year, data = M.SOC)
summary(SOC.lm.m)

# Extract te coefficients from the overall model
SOC.coef.m <- coef(SOC.lm.m)
SOC.coef.m

p.2 <- ggplot(M.SOC, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "SOC (mg/g)")
p.2 

p.2  <- p.2  + geom_abline(intercept = SOC.coef.m[1], 
                           slope = SOC.coef.m[2], 
                           aes(colour = "overall"))
p.2 

SOC.lm3 <- lm(mean ~ poly(Year, 3), data = M.SOC)

p.2  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.m <- segmented(SOC.lm.m, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(1995)))
summary(my.seg.m)
# get the breakpoints
my.seg.m$psi
slope(my.seg.m)
# get the fitted data
my.fitted.m <- fitted(my.seg.m)
my.model.m <- data.frame(Year = M.SOC$Year, SOC = my.fitted.m)
ggplot(my.model.m, aes(x = Year, y = SOC)) + geom_line()

p.2 + geom_line(data = my.model.m, aes(x = Year, y = SOC), colour = "tomato")

##add to original plot with controls
pSOC <- pSOC + geom_line(data = my.model.m, aes(x = Year, y = SOC),colour = "steelblue4", cex=1.5)
pSOC <- pSOC +  geom_point(data=M.SOC, aes(x = Year, y = mean), color="steelblue4", alpha=0.3)
pSOC

#Add SOC for WS sucrose treatment
S.SOC<-data %>% subset(Treatment=="WS")%>%group_by(Year) %>% filter(SOC!="NA")%>%summarise(mean=mean(SOC))

#linear model of SOC over time
SOC.lm.s <- lm(mean ~ Year, data = S.SOC)
summary(SOC.lm.s)

# Extract te coefficients from the overall model
SOC.coef.s <- coef(SOC.lm.s)
SOC.coef.s

p.3 <- ggplot(S.SOC, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "SOC (mg/g)")
p.3 

p.3  <- p.3  + geom_abline(intercept = SOC.coef.s[1], 
                           slope = SOC.coef.s[2], 
                           aes(colour = "overall"))
p.3 

SOC.lm3.s <- lm(mean ~ poly(Year, 3), data = S.SOC)

p.3  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.s <- segmented(SOC.lm.s, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(1995)))
summary(my.seg.s)
# get the breakpoints
my.seg.s$psi
slope(my.seg.s)
# get the fitted data
my.fitted.s <- fitted(my.seg.s)
my.model.s <- data.frame(Year = S.SOC$Year, SOC = my.fitted.s)
ggplot(my.model.s, aes(x = Year, y = SOC)) + geom_line()

p.3 + geom_line(data = my.model.s, aes(x = Year, y = SOC), colour = "purple")

##add to original plot with controls

pSOC <- pSOC + geom_line(data = my.model.s, aes(x = Year, y = SOC),colour = "orchid3", cex=1.5)
pSOC <- pSOC +  geom_point(data=S.SOC, aes(x = Year, y = mean), color="orchid3", alpha=0.3)
pSOC

#Add SOC for W water treatment
W.SOC<-data %>% subset(Treatment=="W")%>%group_by(Year) %>% filter(SOC!="NA")%>%summarise(mean=mean(SOC))

#linear model of SOC over time
SOC.lm.w <- lm(mean ~ Year, data = W.SOC)
summary(SOC.lm.w)

# Extract te coefficients from the overall model
SOC.coef.w <- coef(SOC.lm.w)
SOC.coef.w

p.4 <- ggplot(W.SOC, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "SOC (mg/g)")
p.4 

p.4  <- p.4  + geom_abline(intercept = SOC.coef.w[1], 
                           slope = SOC.coef.w[2], 
                           aes(colour = "overall"))
p.4 

SOC.lm3.w <- lm(mean ~ poly(Year, 3), data = W.SOC)

p.4  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.w <- segmented(SOC.lm.w, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(1995)))
summary(my.seg.w)
# get the breakpoints
my.seg.w$psi
slope(my.seg.w)
# get the fitted data
my.fitted.w <- fitted(my.seg.w)
my.model.w <- data.frame(Year = W.SOC$Year, SOC = my.fitted.w)
ggplot(my.model.w, aes(x = Year, y = SOC)) + geom_line()

p.4 + geom_line(data = my.model.w, aes(x = Year, y = SOC), colour = "blue")

##add to original plot with controls

pSOC <- pSOC + geom_abline(intercept = SOC.coef.w[1], 
                           slope = SOC.coef.w[2], 
                           aes(colour = "blue"), color="blue")
pSOC

pSOC <- pSOC + geom_line(data = my.model.w, aes(x = Year, y = SOC),colour = "deepskyblue", cex=1.5)
pSOC <- pSOC +  geom_point(data=W.SOC, aes(x = Year, y = mean), color="deepskyblue", alpha=0.3)
pSOC
pSOC1<-pSOC

#panel for manuscript
#first a dummy plot to extract a legend from
legend.dat<-data %>% subset(Treatment=="C"|Treatment=="WM"|Treatment=="WS"|Treatment=="W")
l.plot<-ggplot(legend.dat, aes(y=SOC2, x=Year, color=Treatment))+
  geom_point(pch=16, cex=4)+
  scale_color_manual(values=c('black', 'steelblue', 'orchid3', 'deepskyblue'), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c('Untreated', 'Mannitol', 'Sucrose', 'Water'))

#code to extract legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(l.plot)

library(ggpubr)
bp_fig <- ggarrange(pSOC1, pSTL1, pETL1,mylegend,
                      ncol =4, nrow =1, 
                      align = "v",labels = c("a)", "b)", "c)", ""))
bp_fig

#Add SOC for T temp treatment
T.SOC<-data %>% subset(Treatment=="T")%>%group_by(Year) %>% filter(SOC!="NA")%>%summarise(mean=mean(SOC))

#linear model of SOC over time
SOC.lm.t <- lm(mean ~ Year, data = T.SOC)
summary(SOC.lm.t)

# Extract the coefficients from the overall model
SOC.coef.t <- coef(SOC.lm.t)
SOC.coef.t

p.5 <- ggplot(T.SOC, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "SOC (mg/g)")
p.5 

p.5  <- p.5  + geom_abline(intercept = SOC.coef.t[1], 
                           slope = SOC.coef.t[2], 
                           aes(colour = "overall"))
p.5 

SOC.lm3.t <- lm(mean ~ poly(Year, 3), data = T.SOC)

p.5  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.t <- segmented(SOC.lm.t, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(1995)))
summary(my.seg.t)
# get the breakpoints
my.seg.t$psi
slope(my.seg.t)
# get the fitted data
my.fitted.t <- fitted(my.seg.t)
my.model.t <- data.frame(Year = T.SOC$Year, SOC = my.fitted.t)
ggplot(my.model.t, aes(x = Year, y = SOC)) + geom_line()

p.5 + geom_line(data = my.model.t, aes(x = Year, y = SOC), colour = "red")

##add to original plot with controls

pSOC <- pSOC + geom_line(data = my.model.t, aes(x = Year, y = SOC), colour = "red")
pSOC

#Add SOC for TW temp/water treatment
TW.SOC<-data %>% subset(Treatment=="TW")%>%group_by(Year) %>% filter(SOC!="NA")%>%summarise(mean=mean(SOC))

#linear model of SOC over time
SOC.lm.tw <- lm(mean ~ Year, data = TW.SOC)
summary(SOC.lm.tw)

# Extract the coefficients from the overall model
SOC.coef.tw <- coef(SOC.lm.tw)
SOC.coef.tw

p.6 <- ggplot(TW.SOC, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "SOC (mg/g)")
p.6 

p.6  <- p.6  + geom_abline(intercept = SOC.coef.tw[1], 
                           slope = SOC.coef.tw[2], 
                           aes(colour = "overall"))
p.6 

SOC.lm3.tw <- lm(mean ~ poly(Year, 3), data = TW.SOC)

p.6  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.tw <- segmented(SOC.lm.tw, 
                       seg.Z = ~ Year, 
                       psi = list(Year = c(1995)))
summary(my.seg.tw)
# get the breakpoints
my.seg.tw$psi
slope(my.seg.tw)
# get the fitted data
my.fitted.tw <- fitted(my.seg.tw)
my.model.tw <- data.frame(Year = TW.SOC$Year, SOC = my.fitted.tw)
ggplot(my.model.tw, aes(x = Year, y = SOC)) + geom_line()

p.6 + geom_line(data = my.model.tw, aes(x = Year, y = SOC), colour = "orangered1")

##add to original plot with controls

pSOC <- pSOC + geom_line(data = my.model.tw, aes(x = Year, y = SOC), colour = "orangered1")
pSOC

#Add SOC for TWM temp/water/mannitol treatment
TWM.SOC<-data %>% subset(Treatment=="TWM")%>%group_by(Year) %>% filter(SOC!="NA")%>%summarise(mean=mean(SOC))

#linear model of SOC over time
SOC.lm.twm <- lm(mean ~ Year, data = TWM.SOC)
summary(SOC.lm.twm)

# Extract the coefficients from the overall model
SOC.coef.twm <- coef(SOC.lm.twm)
SOC.coef.twm

p.7 <- ggplot(TWM.SOC, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "SOC (mg/g)")
p.7

p.7  <- p.7  + geom_abline(intercept = SOC.coef.twm[1], 
                           slope = SOC.coef.twm[2], 
                           aes(colour = "overall"))
p.7 

SOC.lm3.twm <- lm(mean ~ poly(Year, 3), data = TWM.SOC)

p.7  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.twm <- segmented(SOC.lm.twm, 
                        seg.Z = ~ Year, 
                        psi = list(Year = c(1999)))
summary(my.seg.twm)
# get the breakpoints
my.seg.twm$psi
slope(my.seg.twm)
# get the fitted data
my.fitted.twm <- fitted(my.seg.twm)
my.model.twm <- data.frame(Year = TWM.SOC$Year, SOC = my.fitted.twm)
ggplot(my.model.twm, aes(x = Year, y = SOC)) + geom_line()

p.7 + geom_line(data = my.model.twm, aes(x = Year, y = SOC), colour = "orange")

##add to original plot with controls

pSOC <- pSOC + geom_line(data = my.model.twm, aes(x = Year, y = SOC), colour = "orange")
pSOC

#Add SOC for TWM temp/water/mannitol treatment
TWS.SOC<-data %>% subset(Treatment=="TWS")%>%group_by(Year) %>% filter(SOC!="NA")%>%summarise(mean=mean(SOC))

#linear model of SOC over time
SOC.lm.tws <- lm(mean ~ Year, data = TWS.SOC)
summary(SOC.lm.tws)

# Extract the coefficients from the overall model
SOC.coef.tws <- coef(SOC.lm.tws)
SOC.coef.tws

p.8 <- ggplot(TWS.SOC, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "SOC (mg/g)")
p.8 

p.8   <- p.8  + geom_abline(intercept = SOC.coef.tws[1], 
                            slope = SOC.coef.tws[2], 
                            aes(colour = "overall"))
p.8 

SOC.lm3.tws <- lm(mean ~ poly(Year, 3), data = TWS.SOC)

p.8  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.tws <- segmented(SOC.lm.tws, 
                        seg.Z = ~ Year, 
                        psi = list(Year = c(1995)))
summary(my.seg.tws)
# get the breakpoints
my.seg.tws$psi
slope(my.seg.tws)
# get the fitted data
my.fitted.tws <- fitted(my.seg.tws)
my.model.tws <- data.frame(Year = TWS.SOC$Year, SOC = my.fitted.tws)
ggplot(my.model.tws, aes(x = Year, y = SOC)) + geom_line()

p.8 + geom_line(data = my.model.tws, aes(x = Year, y = SOC), colour = "firebrick")

##add to original plot with controls

pSOC <- pSOC + geom_line(data = my.model.tws, aes(x = Year, y = SOC), colour = "firebrick")
pSOC

grid.arrange(pSOC,mylegend, ncol=2)


#############
###Scottnema x SOC######
STL.SOC<-data %>% filter(Treatment=="C"|Treatment=="WM")%>%group_by(Year,Treatment) %>% filter(SOC2!="NA", STL !="NA")%>%summarise(SOC=mean(SOC2), STL=mean(STL))

#linear model of SOC over time
SOC.lm.n <- lm(STL ~ SOC, data = STL.SOC)
summary(SOC.lm.n)

# Extract te coefficients from the overall model
SOC.coef.n <- coef(SOC.lm.n)
SOC.coef.n

p <- ggplot(STL.SOC, aes(x = SOC, y = STL)) + geom_point()+
  labs(x = "SOC (mg/g)",
       y = "Scottnema total live")
p

p <- p + geom_abline(intercept = SOC.coef.n[1], 
                     slope = SOC.coef.n[2], 
                     aes(colour = "overall"))
p

SOC.lm3.n <- lm(mean ~ poly(Year, 3), data = STL.SOC)

p + geom_smooth(method = "lm",
                formula = y ~ poly(x, degree = 3), 
                se = FALSE, colour = "orange")
my.seg <- segmented(SOC.lm.n, 
                    seg.Z = ~ SOC, 
                    psi = list(SOC = c(0.3)))
summary(my.seg)
# get the breakpoints
my.seg$psi
slope(my.seg)
# get the fitted data
my.fitted <- fitted(my.seg)
my.model <- data.frame(SOC = STL.SOC$SOC, STL = my.fitted)
ggplot(my.model, aes(x = SOC, y = STL)) + geom_line()

pSTL.SOC<- ggplot(STL.SOC, aes(x = SOC, y = STL)) + 
  labs(x = "Soil Organic Carbon (mg/g)",y ="Scottnema total live" )+ 
  geom_line(data = my.model, aes(x = SOC, y = STL), colour = "black")+
  ylim(0,3000)
pSTL.SOC

#############
###CHL######
C.CHL<-data %>% subset(Treatment=="C")%>%group_by(Year) %>% filter(Chl_a!="NA")%>%summarise(mean=mean(Chl_a))

#linear model of CHL over time
CHL.lm <- lm(mean ~ Year, data = C.CHL)
summary(CHL.lm)

# Extract te coefficients from the overall model
CHL.coef <- coef(CHL.lm)
CHL.coef

p <- ggplot(C.CHL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Chla ug/g")
p

p <- p + geom_abline(intercept = CHL.coef[1], 
                     slope = CHL.coef[2], 
                     aes(colour = "overall"))
p

CHL.lm3 <- lm(mean ~ poly(Year, 3), data = C.CHL)

p + geom_smooth(method = "lm",
                formula = y ~ poly(x, degree = 3), 
                se = FALSE, colour = "orange")
my.seg <- segmented(CHL.lm, 
                    seg.Z = ~ Year, 
                    psi = list(Year = c(2001)))
summary(my.seg)
# get the breakpoints
my.seg$psi
slope(my.seg)
# get the fitted data
my.fitted <- fitted(my.seg)
my.model <- data.frame(Year = C.CHL$Year, CHL = my.fitted)
ggplot(my.model, aes(x = Year, y = CHL)) + geom_line()

pCHL<- ggplot(C.CHL, aes(x = Year, y = mean)) + 
  geom_point(color="black", alpha=0.3)+
  labs(x = "Year",y = "Chlorophyll-a (ug/g)")+ 
  geom_segment(aes(x=1993, y=0.0285, xend=2004.5, yend=0.0285), color = "grey0") + 
  geom_point(aes(x=1993, y=0.0285), color="grey0")+
  geom_point(aes(x=2004.5, y=0.0285), color="grey0")+
  geom_segment(aes(x=2005, y=0.0285, xend=2015, yend=0.0285), color = "grey60")+  
  geom_point(aes(x=2005, y=0.0285), color="grey60")+
  geom_point(aes(x=2015, y=0.0285), color="grey60")+
  geom_segment(aes(x=2001, y=-0.003, xend=2001, yend=-0.005), arrow=arrow(length = unit(0.03, "npc")), color = "grey0") + 
  annotate("text", x = c(2001), y=(-0.002), label = c("Flood year"), fontface = 3,size=4)+
  annotate("text", x = c(1999), y=(0.03), label = c("Treatment"), fontface = 2,size=4)+
  annotate("text", x = c(2010), y=(0.03), label = c("Post-Treatment"), fontface = 2,size=4, color="grey60")+
  geom_line(data = my.model, aes(x = Year, y = CHL), colour = "black", cex=1.5)+
  ylim(-0.005,0.03)
pCHL

#Add CHL for WM mannitol treatment
M.CHL<-data %>% subset(Treatment=="WM")%>%group_by(Year) %>% filter(Chl_a!="NA")%>%summarise(mean=mean(Chl_a))

#linear model of CHL over time
CHL.lm.m <- lm(mean ~ Year, data = M.CHL)
summary(CHL.lm.m)

# Extract te coefficients from the overall model
CHL.coef.m <- coef(CHL.lm.m)
CHL.coef.m

p.2 <- ggplot(M.CHL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Chla (ug/g)")
p.2 

p.2  <- p.2  + geom_abline(intercept = CHL.coef.m[1], 
                           slope = CHL.coef.m[2], 
                           aes(colour = "overall"))
p.2 

CHL.lm3 <- lm(mean ~ poly(Year, 3), data = M.CHL)

p.2  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.m <- segmented(CHL.lm.m, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(1999)))
summary(my.seg.m)
# get the breakpoints
my.seg.m$psi
slope(my.seg.m)
# get the fitted data
my.fitted.m <- fitted(my.seg.m)
my.model.m <- data.frame(Year = M.CHL$Year, CHL = my.fitted.m)
ggplot(my.model.m, aes(x = Year, y = CHL)) + geom_line()

p.2 + geom_line(data = my.model.m, aes(x = Year, y = CHL), colour = "tomato")

##add to original plot with controls
pCHL <- pCHL + geom_line(data = my.model.m, aes(x = Year, y = CHL),colour = "steelblue4", cex=1.5)
pCHL <- pCHL +  geom_point(data=M.CHL, aes(x = Year, y = mean), color="steelblue4", alpha=0.3)
pCHL

#Add CHL for WS sucrose treatment
S.CHL<-data %>% subset(Treatment=="WS")%>%group_by(Year) %>% filter(Chl_a!="NA")%>%summarise(mean=mean(Chl_a))

#linear model of CHL over time
CHL.lm.s <- lm(mean ~ Year, data = S.CHL)
summary(CHL.lm.s)

# Extract te coefficients from the overall model
CHL.coef.s <- coef(CHL.lm.s)
CHL.coef.s

p.3 <- ggplot(S.CHL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "Chla (ug/g)")
p.3 

p.3  <- p.3  + geom_abline(intercept = CHL.coef.s[1], 
                           slope = CHL.coef.s[2], 
                           aes(colour = "overall"))
p.3 

CHL.lm3.s <- lm(mean ~ poly(Year, 3), data = S.CHL)

p.3  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.s <- segmented(CHL.lm.s, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(1999)))
summary(my.seg.s)
# get the breakpoints
my.seg.s$psi
slope(my.seg.s)
# get the fitted data
my.fitted.s <- fitted(my.seg.s)
my.model.s <- data.frame(Year = S.CHL$Year, CHL = my.fitted.s)
ggplot(my.model.s, aes(x = Year, y = CHL)) + geom_line()

p.3 + geom_line(data = my.model.s, aes(x = Year, y = CHL), colour = "purple")

##add to original plot with controls

pCHL <- pCHL + geom_line(data = my.model.s, aes(x = Year, y = CHL),colour = "orchid3", cex=1.5)
pCHL <- pCHL +  geom_point(data=S.CHL, aes(x = Year, y = mean), color="orchid3", alpha=0.3)
pCHL

#Add CHL for W water treatment
W.CHL<-data %>% subset(Treatment=="W")%>%group_by(Year) %>% filter(Chl_a!="NA")%>%summarise(mean=mean(Chl_a))

#linear model of CHL over time
CHL.lm.w <- lm(mean ~ Year, data = W.CHL)
summary(CHL.lm.w)

# Extract te coefficients from the overall model
CHL.coef.w <- coef(CHL.lm.w)
CHL.coef.w

p.4 <- ggplot(W.CHL, aes(x = Year, y = mean)) + geom_point()+
  labs(x = "Year",
       y = "CHL (ug/g)")
p.4 

p.4  <- p.4  + geom_abline(intercept = CHL.coef.w[1], 
                           slope = CHL.coef.w[2], 
                           aes(colour = "overall"))
p.4 

CHL.lm3.w <- lm(mean ~ poly(Year, 3), data = W.CHL)

p.4  + geom_smooth(method = "lm",
                   formula = y ~ poly(x, degree = 3), 
                   se = FALSE, colour = "orange")
my.seg.w <- segmented(CHL.lm.w, 
                      seg.Z = ~ Year, 
                      psi = list(Year = c(1999)))
summary(my.seg.w)
# get the breakpoints
my.seg.w$psi
slope(my.seg.w)
# get the fitted data
my.fitted.w <- fitted(my.seg.w)
my.model.w <- data.frame(Year = W.CHL$Year, CHL = my.fitted.w)
ggplot(my.model.w, aes(x = Year, y = CHL)) + geom_line()

p.4 + geom_line(data = my.model.w, aes(x = Year, y = CHL), colour = "blue")

##add to original plot with controls


pCHL <- pCHL + geom_line(data = my.model.w, aes(x = Year, y = CHL),colour = "deepskyblue", cex=1.5)
pCHL <- pCHL +  geom_point(data=W.CHL, aes(x = Year, y = mean), color="deepskyblue", alpha=0.3)
pCHL
pCHL1<-pCHL

#run time series to get psummer temperature plot
bp_fig_v2 <- ggarrange(psummer, pSOC1, pCHL1, pSTL1, pETL1, mylegend,
                    ncol =3, nrow =2, 
                    labels = c("a)", "b)", "c)", "d)", "e)", ""))
bp_fig_v2

