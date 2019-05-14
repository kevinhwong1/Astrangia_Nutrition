#Title: A. Poculata Nutrition Graphing
#Author: KH Wong
#Date Last Modified: 20190501
#See Readme file for details 

library(dbplyr)
library(tidyverse)
library(readr)
library(stringr)
library(Rmisc)

# Import Data
cs_raw <- read.csv("../Data/CS_20190428.csv")
protein <- read.csv("../Data/Protein.Calculation.csv")
well.meta <- read.csv("../Data/CS_Well_metadata.csv")
vial.meta <- read.csv("../Data/Vial.Homogenate.Data.csv")
coral.meta <- read.csv("../Data/Master_Fragment_Apoc2019.csv")

PAM <- read.csv("../Data/PAM_Apoc2019.csv")
PAM.T2 <- PAM %>%
  filter(Timepoint == "T2")

chl.Host <- read.csv("../Data/mean.chl.host.csv")
CS.Host <- read.csv("../Data/CS.host.calculated.csv")
color <- read.csv("../Data/Color_Score.csv")

### Comparing CS with Fo

PAM <- read.csv("../Data/PAM_Apoc2019.csv")
PAM.T2 <- PAM %>%
  filter(Timepoint == "T2")

CS.Host.PAM <- merge(mean.cs.host.plug, PAM.T2, by = "Plug")

CS.PAM.H <- ggplot(CS.Host.PAM, aes(x=Fo, y=cs.activity.mUmg)) +
  geom_point() +   
  #  geom_smooth(method=lm) +
  xlab("Fo") + ylab("Host CS Activity mUmg")+ #Axis titles
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.title.x = element_blank())+
  theme(legend.position = "none") +
  facet_wrap(~Exp.Cond)

ggsave(file = "../Output/CS.Fo.pdf", CS.PAM.H)

NL <- CS.Host.PAM %>%
  filter(Exp.Cond == "No.Light")
NL.lm <- lm (cs.activity.mUmg ~ Fo, data = NL)
NL.lmsummary <- summary(NL.lm)

L <- CS.Host.PAM %>%
  filter(Exp.Cond == "Light")
L.lm <- lm (cs.activity.mUmg ~ Fo, data = L)
L.lmsummary <- summary(L.lm)

LF <- CS.Host.PAM %>%
  filter(Exp.Cond == "Light.Feed")
LF.lm <- lm (cs.activity.mUmg ~ Fo, data = LF)
LF.lmsummary <- summary(LF.lm)


CS.PAM.H.all <- ggplot(CS.Host.PAM, aes(x=Fo, y=cs.activity.mUmg, color = Exp.Cond)) +
  geom_point() +   
  #  geom_smooth(method=lm) +
  xlab("Fo") + ylab("Host CS Activity mUmg")+ #Axis titles
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.title.x = element_blank())+
  theme(legend.position = "none") 


# comparing to chl 
chl.host <- read.csv("../Data/mean.chl.host.csv")

CS.Chl.host.plug <- merge(mean.cs.host.plug, chl.host, by = "Plug")
CS.Chl.host.exp <- merge(CS.Chl.host.plug, coral.meta, by = "Plug")

CS.Chl.H <- ggplot(CS.Chl.host.exp, aes(x=ChlA.ug.pro.mg, y=cs.activity.mUmg)) +
  geom_point() +   
  geom_smooth(method=lm) +
  xlab("Chlorophyll A ug/mg protein") + ylab("Host CS Activity mUmg")+ #Axis titles
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.title.x = element_blank())+
  theme(legend.position = "none") +
  facet_wrap(~Exp.Cond)

NL.cs.chl <- CS.Chl.host.exp %>%
  filter(Exp.Cond == "No.Light")
NL.cs.chl.lm <- lm (cs.activity.mUmg ~ ChlA.ug.pro.mg, data = NL.cs.chl)
NL.cs.chl.lmsummary <- summary(NL.cs.chl.lm)

L.cs.chl <- CS.Chl.host.exp %>%
  filter(Exp.Cond == "Light")
L.cs.chl.lm <- lm (cs.activity.mUmg ~ ChlA.ug.pro.mg, data = L.cs.chl)
L.cs.chl.lmsummary <- summary(L.cs.chl.lm)

LF.cs.chl <- CS.Chl.host.exp %>%
  filter(Exp.Cond == "Light.Feed")
LF.cs.chl.lm <- lm (cs.activity.mUmg ~ ChlA.ug.pro.mg, data = LF.cs.chl)
LF.cs.chl.lmsummary <- summary(LF.cs.chl.lm)


#chl by treatment
mean.Chl<- summarySE(CS.Chl.host.exp, measurevar="ChlA.ug.pro.mg", groupvars="Exp.Cond")
Treat.Chl.H <- ggplot(mean.Chl, aes(x=Exp.Cond, y=ChlA.ug.pro.mg)) +
  geom_point() +
  geom_errorbar(aes(ymin=ChlA.ug.pro.mg-se, ymax=ChlA.ug.pro.mg+se), width=.1, color="black") + #Error bars
#  geom_smooth(method=lm) +
  xlab("Treatment") + ylab("Chlorophyll A ug/mg proteing")+ #Axis titles
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.title.x = element_blank())+
  theme(legend.position = "none") 

ggsave(file = "../Output/Chl_Treatment.pdf", Treat.Chl.H)

## Correlating PAM and Chl

PAM.Chl <- merge(PAM.T2, chl.host, by = "Plug")
PAM.Chl.Plot <- ggplot(PAM.Chl, aes(x=ChlA.ug.pro.mg, y=Fo)) +
  geom_point() +   
  geom_smooth(method=lm) +
  xlab("Chlorophyll A ug/mg protein") + ylab("Fo")+ #Axis titles
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.title.x = element_blank())+
  theme(legend.position = "none") +
  facet_wrap(~Exp.Cond)

ggsave(file = "../Output/Fo_Chl_Treatment.pdf", PAM.Chl.Plot)


## Correlating CS and Color score

color.score <- read.csv("../Data/Color_Score.csv")
colnames(color.score)[colnames(color.score) == "PLUG.ID"] <- "Plug"

CS.color.T2 <- merge(CS.Host, color.score, by = "Plug") %>%
  filter(Timepoint == "t2")

CS.color.H.T2 <- ggplot(CS.color.T2, aes(x=Bleaching.Score, y=cs.activity.mUmg)) +
  geom_point() +   
  geom_smooth(method=lm) +
  xlab("Bleaching Score") + ylab("Host CS Activity mUmg")+ #Axis titles
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.title.x = element_blank())+
  theme(legend.position = "none") +
  facet_wrap(~Treatment)

ggsave(file = "../Output/CS_Color_Treatment.pdf", CS.color.H.T2)

## Correlating SAM and CS

CS.SAM.host <- merge(SAM.mean.plug, CS.Host, by = "Plug")

SAM.CS.H <- ggplot(CS.SAM.host, aes(x=cs.activity.mUmg, y=SAM.Pro)) +
  geom_point() +   
  geom_smooth(method=lm) +
  xlab("Host CS Activity mUmg") + ylab("ug SAM/mg Protein")+ #Axis titles
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.title.x = element_blank())+
  theme(legend.position = "none") +
  facet_wrap(~Exp.Cond)

ggsave(file = "../Output/SAM_CS_Treatment.pdf", SAM.CS.H)
