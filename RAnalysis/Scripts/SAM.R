#Title: A. Poculata Nutrition SAM analysis
#Author: KH Wong
#Date Last Modified: 20190429
#See Readme file for details 

library(dbplyr)
library(tidyverse)
library(readr)
library(stringr)
library(Rmisc)

#Import data

SAM_raw <- read.csv("../Data/SAM_Apoc_20190513.csv")
SAM_meta <- read.csv("../Data/SAM_Well_Meta.csv")
protein <- read.csv("../Data/Protein.Calculation.csv")
vial.meta <- read.csv("../Data/Vial.Homogenate.Data.csv")
coral.meta <- read.csv("../Data/Master_Fragment_Apoc2019.csv")

#Merging data frames
SAM_data <- merge(SAM_meta, SAM_raw, by = "Well")

#Log concentration values
SAM_data$log_conc <- log(SAM_data$Concentration)

#subsetting standard
SAM.standard <- SAM_data %>% 
  filter(Sample.Type == "Standard") 

#Plotting
SAM.plot.S1<- ggplot(data = SAM.standard, aes(x=log_conc, y=Absorbance))+
  ylab("Absorbance (nm)")+ xlab("Log Concentration") + 
  geom_point()+
  geom_smooth(method = "lm") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Creating model 
SAM.lmstandard1 <- lm (log_conc ~ Absorbance, data = SAM.standard)
SAM.lmsummary1 <- summary(SAM.lmstandard1)

#Interpolating samples values from model
SAM.Sample1 <- SAM_data %>% 
  filter(Sample.Type == "Sample") 
SAM.Sample1$log_conc_calc <- predict(SAM.lmstandard1, newdata = SAM.Sample1) #using model to get concentration

#Inverse log transformation
SAM.Sample1$conc_calc <- exp(SAM.Sample1$log_conc_calc)

#Formatting protein dataset
protein.host <- protein %>%
  separate(Vial, c("Vial", "Species"), "-") %>%
  filter(Species == "H")

#merging protein/SAM with meta
pro.meta <- merge(protein.host, vial.meta, by = "Vial") %>%
  select(Concentration.ug.uL, Plug, Homo.Vol.ml)
SAM.metad <- merge(SAM.Sample1, vial.meta, by = "Vial") %>%
  select(conc_calc, Plug)
SAM.Pro <- merge(SAM.metad, pro.meta, by = "Plug")
SAM.Pro.coral <- merge(SAM.Pro, coral.meta, by = "Plug")

#Calculating ngSAM/mg protein (ug/ul = mg/ml)
SAM.Pro.coral$SAM.Pro <- SAM.Pro.coral$conc_calc/SAM.Pro.coral$Concentration.ug.uL

#Meaning values per plug
SAM.mean.plug <- summarySE(SAM.Pro.coral, measurevar="SAM.Pro", groupvars=c("Plug", "Exp.Cond"))

#Comparing to Chl
chl.host <- read.csv("../Data/mean.chl.host.csv")
SAM.Chl <- merge(SAM.mean.plug, chl.host, by = "Plug")

SAM.Chl.H <- ggplot(SAM.Chl, aes(x=ChlA.ug.pro.mg, y=SAM.Pro)) +
  geom_point() +   
  geom_smooth(method=lm) +
  xlab("Chlorophyll A ug/mg protein") + ylab("ug SAM/mg Protein")+ #Axis titles
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.title.x = element_blank())+
  theme(legend.position = "none") +
  facet_wrap(~Exp.Cond)

ggsave(file = "../Output/SAM_Chl_Treatment.pdf", SAM.Chl.H)


# Mean SAM per treatment
mean.SAM.exp <- summarySE(SAM.mean.plug, measurevar="SAM.Pro", groupvars="Exp.Cond")

mean.cs.H <- ggplot(mean.SAM.exp, aes(x=Exp.Cond, y=SAM.Pro)) + 
  geom_point() +
  geom_errorbar(aes(ymin=SAM.Pro-se, ymax=SAM.Pro+se), width=.1, color="black") + #Error bars
  xlab("Treatment") + ylab("ugSAM/mgProtein")+ #Axis titles
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.title.x = element_blank())+
  theme(legend.position = "none")
ggsave(file = "../Output/SAM_Treatment.pdf", mean.cs.H)
