#Title: A. Poculata Nutrition citrate synthase analysis
#Author: KH Wong
#Date Last Modified: 20190429
#See Readme file for details 

#Load Libaries
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

# Calculating change in blanks
cs_raw$Blank <- cs_raw$'B.Run4' - cs_raw$'B.Run2'

# Calculating change in OA addition 
cs_raw$OA <- cs_raw$`O.Run4` - cs_raw$`O.Run2`

### Calculation of CS Specific activity (U mg-1) from Hawkins et al. 2016 ###
#((OA - Blank) X Vrxn X Dilution Factor) / 13.6 X Path Length X Vsample X Protein Content 

# Merging protein data and selecting columns of interest
Well.meta.pro <- merge(well.meta, protein, by = "Vial")

# making a unique column for run and well
cs_raw$run.well <- paste(cs_raw$Well, cs_raw$Run, sep = "-")
Well.meta.pro$run.well <- paste(Well.meta.pro$Well, Well.meta.pro$Run, sep = "-")

# merging dataframes 
cs.raw.pro <- merge(cs_raw, Well.meta.pro, by = "run.well") 

#Separting -H and -S to attach metadata
cs.calc <- cs.raw.pro %>%
  separate(Vial, c("Vial", "Species"), "-") %>%
  select(Blank, OA, Vial, Species, DLF, Concentration.ug.uL)

# Caclulating CS activity U/mg
cs.calc$cs.activity.Umg <- ((cs.calc$OA - cs.calc$Blank) * 0.2 * cs.calc$DLF)/(13.6 * 0.552 * 0.02 * cs.calc$Concentration.ug.uL)

# Converting to mU/mg to compare to Hawkins et al. 2016
cs.calc$cs.activity.mUmg <- cs.calc$cs.activity.Umg * 1000

# Merging vial metadata 
cs.calc.2 <- merge(cs.calc, vial.meta, by = "Vial")  %>%
  select(Vial, Species, cs.activity.mUmg, Plug)

#Merging coral metadata
cs.calc.3 <- merge(cs.calc.2, coral.meta, by = "Plug")

#Filtering only host CS
cs.calc.host <- cs.calc.3 %>%
  filter(Species == "H")

# Summarizing
mean.cs.host <- summarySE(cs.calc.host, measurevar="cs.activity.mUmg", groupvars="Exp.Cond")

mean.cs.host.plug <- summarySE(cs.calc.host, measurevar="cs.activity.mUmg", groupvars="Plug")

mean.cs.H <- ggplot(mean.cs.host, aes(x=Exp.Cond, y=cs.activity.mUmg)) + 
  geom_point() +
  geom_errorbar(aes(ymin=cs.activity.mUmg-se, ymax=cs.activity.mUmg+se), width=.1, color="black") + #Error bars
  xlab("Treatment") + ylab("CS Activity mUmg")+ #Axis titles
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.title.x = element_blank())+
  theme(legend.position = "none")

ggsave(file = "../Output/Mean.CS.H.pdf", mean.cs.H)
write.csv(mean.cs.host.plug, file = "../Data/CS.host.calculated.csv")

