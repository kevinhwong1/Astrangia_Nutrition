#Title: PAM analysis
#Author: Kevin Wong
#Date last modified: 20190426

library(dbplyr)
library(tidyverse)
library(readr)
library(stringr)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(Rmisc)

#Import data 
PAM.raw <- read.csv("../Data/PAM_Apoc2019.csv")

#Filter NAs
PAM.filter <- PAM.raw %>% 
  filter(Y != "NA")

#Timepoints by treatment
PAM.means <- summarySE(PAM.filter, measurevar="Y", groupvars=c("Exp.Cond", "Timepoint"))

#plotting rxn norms 
pd <- position_dodge(0.1) # moves object .05 to the left and right
PAM.Treat<- ggplot(PAM.means, aes(x=Timepoint, y=Y)) + 
  geom_line(position=pd, color="black")+
  geom_errorbar(aes(ymin=Y-se, ymax=Y+se), width=.1, position=pd, color="black") + #Error bars
#  ylim(5,8)+
  xlab("Timepoint") + ylab("Fv/Fm")+ #Axis titles
  geom_point(color ="black", pch=21, size=8)+
#  scale_fill_manual(legend.title, values=c("dodgerblue3", "tomato1"))+ #colour modification
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(colour = "black", size=1), legend.position = "none") +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.title.x = element_blank()) +
  facet_wrap(~ Exp.Cond)

ggsave(file = "../Output/PAM.pdf", PAM.Treat)


#Making NAs = zero

PAM.zero <- PAM.raw %>%
  replace_na(list(Y = 0))

#Timepoints by treatment
PAM.Z.means <- summarySE(PAM.zero, measurevar="Y", groupvars=c("Exp.Cond", "Timepoint"))

pd <- position_dodge(0.1) # moves object .05 to the left and right
PAM.Z.Treat<- ggplot(PAM.Z.means, aes(x=Timepoint, y=Y)) + 
  geom_line(position=pd, color="black")+
  geom_errorbar(aes(ymin=Y-se, ymax=Y+se), width=.1, position=pd, color="black") + #Error bars
  #  ylim(5,8)+
  xlab("Timepoint") + ylab("Fv/Fm")+ #Axis titles
  geom_point(color ="black", pch=21, size=8)+
  #  scale_fill_manual(legend.title, values=c("dodgerblue3", "tomato1"))+ #colour modification
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(colour = "black", size=1), legend.position = "none") +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.title.x = element_blank()) +
  facet_wrap(~ Exp.Cond)

ggsave(file = "../Output/PAM.Z.pdf", PAM.Z.Treat)

#Taking difference between timepoints of each colony

PAM.Z.T1 <- PAM.zero %>%
  filter(Timepoint == "T1")

PAM.Z.T2 <- PAM.zero %>%
  filter(Timepoint == "T2") %>%
  select(Plug,Y)

PAM.Z.Diff <- merge(PAM.Z.T1, PAM.Z.T2, by = "Plug")
PAM.Z.Diff$delta.Y <- PAM.Z.Diff$Y.y - PAM.Z.Diff$Y.x
PAM.Z.Diff.means <- summarySE(PAM.Z.Diff, measurevar="delta.Y", groupvars=("Exp.Cond"))

pd <- position_dodge(0.1) # moves object .05 to the left and right
PAM.Z.Delta<- ggplot(PAM.Z.Diff.means, aes(x=Exp.Cond, y=delta.Y)) + 
#  geom_line(position=pd, color="black")+
  geom_errorbar(aes(ymin=delta.Y-se, ymax=delta.Y+se), width=.1, position=pd, color="black") + #Error bars
  #  ylim(5,8)+
  xlab("Timepoint") + ylab("delta Fv/Fm")+ #Axis titles
  geom_point(color ="black", pch=21, size=8)+
  #  scale_fill_manual(legend.title, values=c("dodgerblue3", "tomato1"))+ #colour modification
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(colour = "black", size=1), legend.position = "none") +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.title.x = element_blank()) 

ggsave(file = "../Output/PAM.Z.delta.pdf", PAM.Z.Delta)
