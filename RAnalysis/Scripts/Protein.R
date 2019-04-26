#Title:Total protien for A. Poculata Nutrition citrate synthase 
#Author: KH Wong
#Date Last Modified: 20190426
#See Readme file for details 

#Load Libaries
library(dbplyr)
library(tidyverse)
library(readr)
library(stringr)
library(Rmisc)

# Import Data
Protein1 <- read.csv("../Data/Protein_Astrangia_20190424_Run1.csv")
Protein2 <- read.csv("../Data/Protein_Astrangia_20190424_Run2.csv")
Protein3 <- read.csv("../Data/Protein_Astrangia_20190425_Run3.csv")
Protein4 <- read.csv("../Data/Protein_Astrangia_20190425_Run4.csv")

well.meta <- read.csv("../Data/Astrangia_Nutrition_Protein_Wells.csv")

# Adding Run numbers into datasets
Protein1$Run <- 1
Protein2$Run <- 2
Protein3$Run <- 3
Protein4$Run <- 4

# Make a unique column for merging
Protein1$run.well <- paste(Protein1$Run, Protein1$Well, sep = "-")
Protein2$run.well <- paste(Protein2$Run, Protein2$Well, sep = "-")
Protein3$run.well <- paste(Protein3$Run, Protein3$Well, sep = "-")
Protein4$run.well <- paste(Protein4$Run, Protein4$Well, sep = "-")

well.meta$run.well <- paste(well.meta$Run, well.meta$Well, sep = "-")

# Merge with metadata
Protein.Run1 <- merge(Protein1, well.meta, by = "run.well")
Protein.Run2 <- merge(Protein2, well.meta, by = "run.well")
Protein.Run3 <- merge(Protein3, well.meta, by = "run.well")
Protein.Run4 <- merge(Protein4, well.meta, by = "run.well")

#Making absorbance reads numeric 
Protein.Run1$X562 <- as.numeric(as.character(Protein.Run1$X562))
Protein.Run2$X562 <- as.numeric(as.character(Protein.Run2$X562))
Protein.Run3$X562 <- as.numeric(as.character(Protein.Run3$X562))
Protein.Run4$X562 <- as.numeric(as.character(Protein.Run4$X562))

# Subtract blanks means for each run (Runs 1 and 3 have standards and blanks)
Protein.standardblank1 <- Protein.Run1 %>% 
  filter(Sample.Type == "Blank") %>%
  summarise(blk.avg = mean(X562))

Protein.standardblank3 <- Protein.Run3 %>% 
  filter(Sample.Type == "Blank") %>%
  summarise(blk.avg = mean(X562))

# Correct for blanks
Protein.Run1$abs.corr <- Protein.Run1$X562 - Protein.standardblank1$blk.avg
Protein.Run2$abs.corr <- Protein.Run2$X562 - Protein.standardblank1$blk.avg

Protein.Run3$abs.corr <- Protein.Run3$X562 - Protein.standardblank3$blk.avg
Protein.Run4$abs.corr <- Protein.Run4$X562 - Protein.standardblank3$blk.avg

# Running standards
target <- c("Standard", "Blank")
Protein.standard1 <- Protein.Run1 %>% 
  filter(Sample.Type %in% target) 

Protein.plot.S1<- ggplot(data = Protein.standard1, aes(x=Concentration, y=abs.corr))+
  ylab("Absorbance (nm)")+ xlab("Concentration") + 
  geom_point()+
  geom_smooth(method = "lm") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

Protein.lmstandard1 <- lm (Concentration ~ abs.corr, data = Protein.standard1)
Protein.lmsummary1 <- summary(Protein.lmstandard1)


Protein.standard3 <- Protein.Run3 %>% 
  filter(Sample.Type %in% target) 

Protein.plot.S3<- ggplot(data = Protein.standard3, aes(x=Concentration, y=abs.corr))+
  ylab("Absorbance (nm)")+ xlab("Concentration") + 
  geom_point()+
  geom_smooth(method = "lm") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

Protein.lmstandard3 <- lm (Concentration ~ abs.corr, data = Protein.standard3)
Protein.lmsummary3 <- summary(Protein.lmstandard3)

# Extrapolate concentration values

Protein.Sample1 <- Protein.Run1 %>% #subsetting Samples
  filter(Sample.Type == "Sample") 
Protein.Sample1$ConcentrationS <- predict(Protein.lmstandard1, newdata = Protein.Sample1) #using model to get concentration

Protein.Sample2 <- Protein.Run2 %>% #subsetting Samples
  filter(Sample.Type == "Sample") 
Protein.Sample2$ConcentrationS <- predict(Protein.lmstandard1, newdata = Protein.Sample2) #using model to get concentration

Protein.Sample3 <- Protein.Run3 %>% #subsetting Samples
  filter(Sample.Type == "Sample") 
Protein.Sample3$ConcentrationS <- predict(Protein.lmstandard3, newdata = Protein.Sample3) #using model to get concentration

Protein.Sample4 <- Protein.Run4 %>% #subsetting Samples
  filter(Sample.Type == "Sample") 
Protein.Sample4$ConcentrationS <- predict(Protein.lmstandard3, newdata = Protein.Sample4) #using model to get concentration

Protein.all <- rbind(Protein.Sample1,Protein.Sample2,Protein.Sample3,Protein.Sample4)
Protein.all <- Protein.all %>%
  select("Vial", "ConcentrationS")

# Convert to mg/mL
Protein.all$Concentration.ug.uL <- Protein.all$ConcentrationS/1000

# Summarizing
mean.protein <- summarySE(Protein.all, measurevar="Concentration.ug.uL", groupvars="Vial")

# Determining ug/ul protein in 20 uL
mean.protein$ug.ul.20 <- mean.protein$Concentration.ug.uL*20

mean.protein <- mean.protein %>%
  select(Vial, Concentration.ug.uL, ug.ul.20)

write.csv(mean.protein, file = "../Data/Protein.Calculation.csv")
