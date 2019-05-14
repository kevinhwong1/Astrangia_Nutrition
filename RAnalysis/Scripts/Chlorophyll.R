#Title: Astrangia Nutrition 2019 Cholorphyll A analysis
#Author: KH Wong
#Date Last Modified: 20190430
#See Readme file for details 

library(dbplyr)
library(tidyverse)
library(readr)
library(stringr)

#Load data 

chla.run1 <- read.csv("../Data/Chl_Astrangia_20190429_Run1.csv")
chla.run2 <- read.csv("../Data/Chl_Astrangia_20190429_Run2.csv")
chla.run3 <- read.csv("../Data/Chl_Astrangia_20190429_Run3.csv")
chla.run4 <- read.csv("../Data/Chl_Astrangia_20190429_Run4.csv")
chla.run5 <- read.csv("../Data/Chl_Astrangia_20190429_Run5.csv")
chla.run6 <- read.csv("../Data/Chl_Astrangia_20190429_Run6.csv")
chla.run7 <- read.csv("../Data/Chl_Astrangia_20190429_Run7.csv")
chla.run8 <- read.csv("../Data/Chl_Astrangia_20190429_Run8.csv")
chla.run9 <- read.csv("../Data/Chl_Astrangia_20190429_Run9.csv")
chla.run10 <- read.csv("../Data/Chl_Astrangia_20190429_Run10.csv")
chla.run11 <- read.csv("../Data/Chl_Astrangia_20190429_Run11.csv")
chla.run12 <- read.csv("../Data/Chl_Astrangia_20190429_Run12.csv")

well.meta <- read.csv("../Data/Chl_Astrangia_20190429.csv")
vial.meta <- read.csv("../Data/Vial.Homogenate.Data.csv")
coral.meta <- read.csv("../Data/Master_Fragment_Apoc2019.csv")
protein <- read.csv("../Data/Protein.Calculation.csv")

# Adding Run number to datasets

chla.run1$Run <- 1
chla.run2$Run <- 2
chla.run3$Run <- 3
chla.run4$Run <- 4
chla.run5$Run <- 5
chla.run6$Run <- 6
chla.run7$Run <- 7
chla.run8$Run <- 8
chla.run9$Run <- 9
chla.run10$Run <- 10
chla.run11$Run <- 11
chla.run12$Run <- 12

# Combining Datasets
chla.all <- rbind(chla.run1, chla.run2, chla.run3, chla.run4, chla.run5, chla.run6, chla.run7, chla.run8, chla.run9, chla.run10,
                  chla.run11, chla.run12)

# Make a unique column 
chla.all$run.well <- paste(chla.all$Run, chla.all$Well, sep = "-")
well.meta$run.well <- paste(well.meta$Run, well.meta$Well, sep = "-")

# Attaching vial metadata
chla.meta <- merge(chla.all, well.meta, by = "run.well" )

#Removing Blank values
chla.data <- chla.meta %>%
  filter(Vial != "NA")

# Attaching colony metadata
chla.data.vial <- merge(chla.data, vial.meta, by = "Vial")
chla.data.vial <- chla.data.vial %>%
  select(`Chl.630`, `Chl.663`, `Chl.750`, Plug, Homo.Vol.ml, Vial)

#Subtracting 750 (blank) from 630 and 633 values
chla.data.vial$abs.630.corr <- chla.data.vial$`Chl.630` - chla.data.vial$`Chl.750`
chla.data.vial$abs.663.corr <- chla.data.vial$`Chl.663` - chla.data.vial$`Chl.750`

#Chlorphyll A concentration equation 
chla.data.vial$chlA.ug.sample <- 11.43*chla.data.vial$abs.663.corr - 0.64*chla.data.vial$abs.630.corr

#Protein cleanup
protein.host <- protein %>%
  separate(Vial, c("Vial", "Species"), "-") %>%
  filter(Species == "H")
protein.host.meta <- merge(protein.host, vial.meta, by = "Vial") %>%
  select(Plug, Concentration.ug.uL)

#adding protein data
chla.data.all <- merge(chla.data.vial, protein.host.meta, by = "Plug") 

#Standardization to ug Chla/ul
chla.data.all$ChlA.ug.uL <- (chla.data.all$chlA.ug.sample * (1000/500) * (chla.data.all$Homo.Vol.ml/1000))

#Standardization to ug Chla/ug protein
chla.data.all$ChlA.ug.pro.ug <- chla.data.all$ChlA.ug.uL / chla.data.all$Concentration.ug.uL

#Standardization to ug Chla/mg protein
chla.data.all$Concentration.mg.uL <- chla.data.all$Concentration.ug.uL * 0.001
chla.data.all$ChlA.ug.pro.mg <- chla.data.all$ChlA.ug.uL / chla.data.all$Concentration.mg.uL

# Summarizing 
mean.ChlA.ug.pro.mg <- summarySE(chla.data.all, measurevar="ChlA.ug.pro.mg", groupvars="Plug")

write.csv(mean.ChlA.ug.pro.mg, file = "../Data/mean.chl.host.csv")
