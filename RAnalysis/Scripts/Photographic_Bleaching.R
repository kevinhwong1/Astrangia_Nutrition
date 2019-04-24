#Title: Bleaching estimates from pictures
#Project: Astrangia Nutrition
#Author: HM Putnam
#Edited by: HM Putnam
#Date Last Modified: 20190424
#See Readme file for details

rm(list=ls()) #clears workspace 

if ("vegan" %in% rownames(installed.packages()) == 'FALSE') install.packages('vegan') 
if ("ggpubr" %in% rownames(installed.packages()) == 'FALSE') install.packages('ggpubr') 
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra') 
if ("plyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('plyr') 
if ("dplyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('dplyr') 
if ("lsmeans" %in% rownames(installed.packages()) == 'FALSE') install.packages('lsmeans') 
if ("multcompView" %in% rownames(installed.packages()) == 'FALSE') install.packages('multcompView') 


#Read in required libraries
##### Include Versions of libraries
library("vegan")
library("ggpubr")
library("gridExtra")
library("plyr") 
library("lsmeans")
library("multcompView")
library("dplyr")

#Required Data files

# Set Working Directory:
setwd("~/MyProjects/Astrangia_Nutrition/RAnalysis/") #set working

Data <- read.csv("Data/Color_Data.csv", header=T, sep=",", na.string="NA") #read in data file
#Sample.Info <- read.csv("Data/Master_Fragment_Sheet.csv", header=T, sep=",", na.string="NA") #read in data file 
#Sample.Info <- Sample.Info[,c(7,5,9)]
#Tank.Info <- read.csv("Data/Tank_to_Treatment.csv", header=T, sep=",", na.string="NA") #read in data file 


#A <- Sample.Info$PLUG.ID
#B <- Data$PLUG.ID
#bad.id.tanks <- which(B %in% A ==FALSE)
#bad.id.tanks

#Data <- merge(Data, Sample.Info, by="PLUG.ID")
#Data <-na.omit(Data)
#Data <- merge(Data, Tank.Info, by="Tank")

Data$Red.Norm.Coral <- Data$Red.Coral/Data$Red.Standard #normalize to color standard
Data$Green.Norm.Coral <- Data$Green.Coral/Data$Green.Standard #normalize to color standard
Data$Blue.Norm.Coral <- Data$Blue.Coral/Data$Blue.Standard #normalize to color standard

par(mfrow=c(1,3))
boxplot(Red.Norm.Coral ~ Treatment*Time, data = Data, lwd = 1, ylab = 'Red Norm Coral', las=2, cex=0.8) #plot boxplot of PC1 color score by Genotype and timepoint
boxplot(Green.Norm.Coral ~ Treatment*Time, data = Data, lwd = 1, ylab = 'Green Norm Coral', las=2, cex=0.8) #plot boxplot of PC1 color score by Genotype and timepoint
boxplot(Blue.Norm.Coral ~ Treatment*Time, data = Data, lwd = 1, ylab = 'Blue Norm Coral', las=2, cex=0.8) #plot boxplot of PC1 color score by Genotype and timepoint

# xxx <- subset(Data, Timepoint=="Time5")
# xxx <- subset(xxx, Species.x=="Pocillopora")
# xxx <- subset(xxx, Treatment.x=="ATHC")
# par(mfrow=c(1,3))
# plot(xxx$Red.Norm.Coral ~ xxx$Tank)
# plot(xxx$Green.Norm.Coral ~ xxx$Tank)
# plot(xxx$Blue.Norm.Coral ~ xxx$Tank)

means <- Data %>%
  group_by(PLUG.ID, Time) %>% 
  summarise_each(funs(mean))


boxplot(Red.Norm.Coral ~ PLUG.ID*Time, data = means, lwd = 1, ylab = 'Red Norm Coral', las=2, cex=0.8) #plot boxplot of PC1 color score by Genotype and timepoint
boxplot(Green.Norm.Coral ~ PLUG.ID*Time, data = means, lwd = 1, ylab = 'Green Norm Coral', las=2, cex=0.8) #plot boxplot of PC1 color score by Genotype and timepoint
boxplot(Blue.Norm.Coral ~ PLUG.ID*Time, data = means, lwd = 1, ylab = 'Blue Norm Coral', las=2, cex=0.8) #plot boxplot of PC1 color score by Genotype and timepoint

blch.scor <- as.matrix(means[,c(13,14,15)]) #create matrix
rownames(blch.scor) <- means$PLUG.ID #name columns in dataframe

dist <- vegdist(blch.scor, method="euclidean") #calculate distance matrix of color scores

PCA.color <- princomp(dist) #run principal components Analysis
summary(PCA.color) # view variance explained by PCs

Blch <- as.data.frame(PCA.color$scores[,1]) #extract PC1
Blch$PLUG.ID <- rownames(blch.scor)
#Blch <- merge(Blch, Data, by="PLUG.ID")
Blch  <- cbind(Blch, means$Time, means$Treatment.1) #make a dataframe of PC1 and experiment factors
colnames(Blch) <- c("Bleaching.Score", "PLUG.ID", "Timepoint", "Treatment")
Blch$Group <- paste(Blch$Timepoint, Blch$Treatment)

#write.table(Blch, file = "~/MyProjects/Holobiont_Integration/RAnalysis/Output/Bleaching_Score.csv", append = FALSE, quote = TRUE, sep = ",",
            #eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            #col.names = TRUE, qmethod = c("escape", "double"),
            #fileEncoding = "")

x<- Blch #assign values to test data frame to look for outliers
par(mar=c(10,4,2,2)) #bottom, left, top and right margins respectively
boxplot(Bleaching.Score ~ Group, data = x, lwd = 1, ylab = 'PC1Color', las=2, cex=0.8) #plot boxplot of PC1 color score by Genotype and timepoint

##### Repeat investigation until min >-25
#min <- which(x$Bleaching.Score==min(x$Bleaching.Score))
#x[min,]
#x<- x[-min,]
#####
wordList <- c("t0")
Blch <- subset(Blch, !(Timepoint %in% wordList))

Blch$Trt <- ifelse(grepl(0, Blch$Treatment), "Dark", Blch$Treatment)
Blch$Trt2 <- ifelse(grepl(1, Blch$Treatment), "Light", Blch$Trt)
Blch$Trt3 <- ifelse(grepl(2, Blch$Treatment), "Light.Fed", Blch$Trt2)
Blch$Treatment <- Blch$Trt3

Blch.data <- Blch[,-c(5:8)]

write.table(Blch.data,"Data/Bleaching_Score.csv",sep=",", row.names=FALSE)

pdf("Output/Photographic_Bleaching.pdf")
par(mar=c(10,4,2,2)) #bottom, left, top and right margins respectively
boxplot(Bleaching.Score ~ Group, data = Blch, lwd = 1, ylab = 'PC1Color', las=2) #plot boxplot of PC1 color score by Genotype and timepoint
stripchart(Bleaching.Score ~ Group, vertical = TRUE, data = Blch, 
           method = "jitter", add = TRUE, pch = 20, col = 'blue', cex=0.2) #include all datapoints in blue overlaid on boxplots
text(x= 0.5, y= min(Blch$Bleaching.Score), labels= "pale") #add text to indicate dark and pale on graphic
text(x= 0.5, y= max(Blch$Bleaching.Score)+2, labels= "dark") #add text to indicate dark and pale on graphic
dev.off()

mod1 <- aov(sqrt(Bleaching.Score+50) ~ Treatment*Timepoint, data=Blch) #run an ANOVA by Genotype
hist(residuals(mod1)) #look at normality of data
boxplot(residuals(mod1)) #look at normality of data
summary(mod1)


marginal = lsmeans(mod1, ~ Treatment*Timepoint)

CLD = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="tukey")

CLD <- CLD[order(CLD$Timepoint, CLD$Treatment),]
CLD

All.Means <- ddply(Blch, c('Timepoint', 'Treatment'), summarize,
                    mean= mean(Bleaching.Score, na.rm=T), #mean pnet
                    N = sum(!is.na(Bleaching.Score)), # sample size
                    se = sd(Bleaching.Score, na.rm=T)/sqrt(N)) #SE
All.Means
#All.Means$se[is.na(All.Means$se)] <- 0
All.Means$Group <- paste(All.Means$Timepoint, All.Means$Treatment)
#All.Means$SpGroup <- paste(All.Means$Treatment, All.Means$Species)

cols <- c("black", "orange", "green")
All.Means$Timepoint <- factor(All.Means$Timepoint, levels = c("t1", "t2"))

#rec <- data.frame(x1=c(1,3,1,5,4))


Fig.All <- ggplot(All.Means, aes(x=Timepoint, y=mean, group=Treatment)) + 
  geom_line(aes(linetype= Treatment, colour=Treatment, group=Treatment), position = position_dodge(width = 0.1), alpha=0.5) + # colour, group both depend on cond2
  geom_errorbar(aes(ymin=All.Means$mean-All.Means$se, ymax=All.Means$mean+All.Means$se), colour="black", width=0, size=0.5, position = position_dodge(width = 0.1)) +
  geom_point(aes(colour=Treatment, shape=Treatment), size = 2, position = position_dodge(width = 0.1)) +
  scale_colour_manual(values=cols) +
  #annotate("text", x=43, y=1.85, label = "a", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=132, y=2.15, label = "b", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=c(43,43), y=c(2.45,2.2), label = "ab", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=c(141,141), y=c(3.15,3.4), label = "c", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=35, y=0.5, label = ".", size = 1) + #add text to the graphic for posthoc letters
  xlab("Timepoint") +
  ylab(expression(paste("Bleaching Score"))) +
  ylim(-3,4) +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.position= "right") + #remove legend background
  ggtitle("Photographic Bleaching") +
  theme(plot.title = element_text(face = 'bold', 
                                  size = 10, 
                                  hjust = 0))
Fig.All



ggsave(file="Output/Photo_Bch.pdf", Fig.All, width = 4, height = 3, units = c("in"))
