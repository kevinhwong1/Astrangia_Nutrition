#Title: Buoyant Weight Data
#Project: NSF BSF
#Author: HM Putnam 
#Edited by: HM Putnam
#Date Last Modified: 20181228
#See Readme file for details

rm(list=ls()) #clears workspace 

#Read in required libraries
##### Include Versions of libraries
library('plyr')
library('ggplot2')
library("sciplot")
library("plotrix")
library("reshape2")
library("nlme") #mixed model, repeated measures ANOVA
library("lsmeans") #mixed model posthoc  statistical comparisons
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra') 
library("gridExtra")

#Required Data files


#Timepoints
#Time0 - 20180920, 20180921
#Time1 - 20180927, 20180928
#Time2 - 201801004, 20181005
#Time3 - 201801018, 201801019
#Time4 - 201801101, 201801102
#Time5 - 201801115, 201801116
#Time6 - 201801129
#Time7 - 201801213
#Time8 - 201801227
#Time9 - 20190110


# Set Working Directory:
setwd("~/MyProjects/Holobiont_Integration/RAnalysis/") #set working

##### Standards ##### 
CalData<-read.csv('Data/Buoyant_Weight/Bouyant_weight_calibration_curve.csv', header=T, sep=",")
#Data with Standard weights in grams in deionized water, saltwater, and the temp that they were measured at

#plot relationship between temp and Standard in fresh and salt water
plot(CalData$Temp, CalData$StandardSalt, col = 'red', ylab = 'Salt Weight (g)', xlab = 'Temp C')
par(new=TRUE)
plot(CalData$Temp, CalData$StandardFresh, col = 'blue',xaxt="n",yaxt="n",xlab="",ylab="")
axis(4)
mtext("'Fresh Weight (g)'",side=4,line=3)
legend("topleft",col=c("red","blue"),lty=1,legend=c("Salt","Fresh"), bty = "n")

#linear regression between temp and Standard
#Salt water standard measures
StandardSaltModel <- lm(CalData$StandardSalt~CalData$Temp)
summary(StandardSaltModel)
summary(StandardSaltModel)$r.squared
StandardSaltModel$coef

#Fresh water standard measures
StandardFreshModel <- lm(CalData$StandardFresh~CalData$Temp)
summary(StandardFreshModel)
summary(StandardFreshModel)$r.squared
StandardFreshModel$coef



##### load coral weight data #####
BW.data <-read.csv('Data/Buoyant_Weight/Buoyant_Weight.csv', header=T, na.strings = "NA") 
BW <- as.numeric(as.character(BW.data$Mass.g)) #assign mass data to BW for function
Temp <- BW.data$Temp.C #assign temp to temp for function

##### BW Function #####
BWCalc <- function(StandardAir= 31.348, Temp, BW, CoralDensity = 2.93){
  #Parameters---------------------------------------------------------------
  # StandardAir is the weight of the Standard in air (Default set at the grams weighed on top of the balance)
  # Temp is the temperature of the water during the coral measurement
  # BW is the buoywant weight of the coral
  # CoralDensity <- 2.03 #g cm-3, set aragonite density for Montipora from literature 
  # Montipora Arag = 2.03 g cm-3 average from table 2 in Anthony 2003 Functional Ecology 17:246-259
  # You can change the density to literatrue values for the specific coral species of interest in the function. 
  
  #Calculation------------------------------------------------------------
  #Step 1: correct the Standard weights for temperature
  # StandardFresh is the weight of the Standard in  fresh water at temp x
  # StandardSalt is the weight of the Standard in salt water at temp x
  
  # This is based on a calibration curve for Standards weighed in fresh and salt water at many temps
StandardFresh <- StandardFreshModel$coef[-1]*Temp + StandardFreshModel$coef[1] 
StandardSalt <- StandardSaltModel$coef[-1]*Temp + StandardSaltModel$coef[1] 
  
  # Step 2: Use weight in air and freshwater of the glass Standard to calculate
  # the density of the Standard (Spencer Davies eq. 4)
FreshDensity <- 1 #Fresh water has a density of 1 g/cm3
StandardDensity <- FreshDensity/(1-(StandardFresh/StandardAir)) 
  
  # Step 3: Calculate the density of seawater using the density of the Standard
  # (Spencer Davies eq. 3)
SWDensity <- StandardDensity*(1-(StandardSalt/StandardAir))
  
  # Step 4: Calculate the dry weight of the coral (Spencer Davies eq. 1)
CoralWeight <- BW/(1-(SWDensity/CoralDensity))
  
  return(CoralWeight) #returns coral dry weights in g
}
  
BW.data$Dry.Weigh.g <- BWCalc(BW=BW, Temp=Temp) #use function to calculate dry weight
hist(BW.data$Dry.Weigh.g, las=2) #examine data
boxplot(BW.data$Dry.Weigh.g ~BW.data$Tank, las=2) #examine data

##### Calculating Rates #####

#Time0 - 20180920, 20180921
#Time1 - 20180927, 20180928
#Time2 - 201801004, 20181005
#Time3 - 201801018, 201801019
#Time4 - 201801101, 201801102
#Time5 - 201801115, 201801116
#Time6 - 201801129
#Time7 - 201801213
#Time8 - 201801227
#Time9 - 20190110

#Dry Weights
data.0.curve <- subset(BW.data, Date==20180920 | Date==20180921)
data.1.curve <- subset(BW.data, Date==20180927 | Date==20180928)
data.2.curve <- subset(BW.data, Date==20181004 | Date==20181005)
data.3.curve <- subset(BW.data, Date==20181018 | Date==20181019)
data.4.curve <- subset(BW.data, Date==20181101 | Date==20181102)
data.5.curve <- subset(BW.data, Date==20181115 | Date==20181116)
data.6.curve <- subset(BW.data, Date==20181129)
data.7.curve <- subset(BW.data, Date==20181213)
data.8.curve <- subset(BW.data, Date==20181227)
data.9.curve <- subset(BW.data, Date==20190110)

boxplot(data.0.curve$Dry.Weigh.g ~data.0.curve$Tank, las=2) #examine data
boxplot(data.1.curve$Dry.Weigh.g ~data.1.curve$Tank, las=2) #examine data
boxplot(data.2.curve$Dry.Weigh.g ~data.2.curve$Tank, las=2) #examine data
boxplot(data.3.curve$Dry.Weigh.g ~data.3.curve$Tank, las=2) #examine data
boxplot(data.4.curve$Dry.Weigh.g ~data.4.curve$Tank, las=2) #examine data
boxplot(data.5.curve$Dry.Weigh.g ~data.5.curve$Tank, las=2) #examine data
boxplot(data.6.curve$Dry.Weigh.g ~data.6.curve$Tank, las=2) #examine data
boxplot(data.7.curve$Dry.Weigh.g ~data.7.curve$Tank, las=2) #examine data
boxplot(data.8.curve$Dry.Weigh.g ~data.8.curve$Tank, las=2) #examine data
boxplot(data.9.curve$Dry.Weigh.g ~data.9.curve$Tank, las=2) #examine data


par(mfrow=c(2,4))
boxplot(data.0.curve$Dry.Weigh.g ~data.0.curve$Tank, main="Time0", ylim=c(0,15), las=2) #examine data
boxplot(data.1.curve$Dry.Weigh.g ~data.1.curve$Tank, main="Time1", ylim=c(0,15), las=2) #examine data
boxplot(data.2.curve$Dry.Weigh.g ~data.2.curve$Tank, main="Time2", ylim=c(0,15), las=2) #examine data
boxplot(data.3.curve$Dry.Weigh.g ~data.3.curve$Tank, main="Time3", ylim=c(0,15), las=2) #examine data
boxplot(data.4.curve$Dry.Weigh.g ~data.4.curve$Tank, main="Time4", ylim=c(0,15), las=2) #examine data
boxplot(data.5.curve$Dry.Weigh.g ~data.5.curve$Tank, main="Time5", ylim=c(0,15), las=2) #examine data
boxplot(data.6.curve$Dry.Weigh.g ~data.6.curve$Tank, main="Time6", ylim=c(0,15), las=2) #examine data
boxplot(data.7.curve$Dry.Weigh.g ~data.7.curve$Tank, main="Time7", ylim=c(0,15), las=2) #examine data
boxplot(data.8.curve$Dry.Weigh.g ~data.8.curve$Tank, main="Time8", ylim=c(0,15), las=2) #examine data
boxplot(data.9.curve$Dry.Weigh.g ~data.9.curve$Tank, main="Time9", ylim=c(0,15), las=2) #examine data

######Compare with single point standard #####
# #load coral weight data
# BWdata <-read.csv('Plob_Buoyant_Weight.csv', header=T, na.strings = "NA") 
# BW.std.data <- read.csv("BW_standards.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
# Arag <- 2.93 # density of aragonite g cm-3
# BW.std.data$FW.Dens <- ((-0.000005*BW.std.data$Temp*BW.std.data$Temp)+(0.000007*BW.std.data$Temp)+1.0001) #g cm-3
# BW.std.data$Ref.Dens <- ((BW.std.data$Ref.Dry.g*BW.std.data$FW.Dens)/(BW.std.data$Ref.Dry.g-BW.std.data$Ref.FW.g)) #g cm-3
# BW.std.data$SW.Dens <- ((BW.std.data$Ref.Dens*(BW.std.data$Ref.Dry.g-BW.std.data$Ref.SW.g))/(BW.std.data$Ref.Dry.g))  #g cm-3
# 
# data.0 <- subset(BWdata, Date==20180127)
# data.1 <- subset(BWdata, Date==20180317)
# data.2 <- subset(BWdata, Date==20180507)
# data.3 <- subset(BWdata, Date==20180827)
# data.0$dry.weight <- ((as.numeric(as.character(data.0$Mass.g)))/(1-(BW.std.data[1,9]/Arag)))
# data.1$dry.weight <- ((as.numeric(as.character(data.1$Mass.g)))/(1-(BW.std.data[2,9]/Arag)))
# data.2$dry.weight <- ((as.numeric(as.character(data.2$Mass.g)))/(1-(BW.std.data[3,9]/Arag)))
# data.3$dry.weight <- ((as.numeric(as.character(data.3$Mass.g)))/(1-(BW.std.data[4,9]/Arag)))
# 
# par(mfrow=c(2,2))
# plot(data.0.curve$Dry.Weigh.g,data.0$dry.weight)
# plot(data.1.curve$Dry.Weigh.g,data.1$dry.weight)
# plot(data.2.curve$Dry.Weigh.g,data.2$dry.weight)
# plot(data.3.curve$Dry.Weigh.g,data.3$dry.weight)
# 
# #can use standard curve data for all time points

#####

#####NEED TO DEAL WITH THOSE REGLUED #####
#If else statements

data <- reshape(BW.data, timevar = "TimePoint", drop = c("QC", "Mass.g","Temp.C", "Notes"), idvar=c("PLUG.ID","Species","Tank"), direction="wide")
T2T<- read.csv("Data/Tank_to_Treatment.csv", header=T, sep=",")
colnames(T2T) <- c("Tank", "Tank.Num", "Treatment")
data <- merge(data, T2T, by="Tank")

#calculating time between measurements
data$Days.Time1 <- as.Date(as.character(data$Date.Format.Time1), format="%m/%d/%y")-
  as.Date(as.character(data$Date.Format.Time0), format="%m/%d/%y")

data$Days.Time2 <- as.Date(as.character(data$Date.Format.Time2), format="%m/%d/%y")-
  as.Date(as.character(data$Date.Format.Time1), format="%m/%d/%y")

data$Days.Time3 <- as.Date(as.character(data$Date.Format.Time3), format="%m/%d/%y")-
  as.Date(as.character(data$Date.Format.Time2), format="%m/%d/%y")

data$Days.Time4 <- as.Date(as.character(data$Date.Format.Time4), format="%m/%d/%y")-
  as.Date(as.character(data$Date.Format.Time3), format="%m/%d/%y")

data$Days.Time5 <- as.Date(as.character(data$Date.Format.Time5), format="%m/%d/%y")-
  as.Date(as.character(data$Date.Format.Time4), format="%m/%d/%y")

data$Days.Time6 <- as.Date(as.character(data$Date.Format.Time6), format="%m/%d/%y")-
  as.Date(as.character(data$Date.Format.Time5), format="%m/%d/%y")

data$Days.Time7 <- as.Date(as.character(data$Date.Format.Time7), format="%m/%d/%y")-
  as.Date(as.character(data$Date.Format.Time6), format="%m/%d/%y")

data$Days.Time8 <- as.Date(as.character(data$Date.Format.Time8), format="%m/%d/%y")-
  as.Date(as.character(data$Date.Format.Time7), format="%m/%d/%y")

data$Days.Time9 <- as.Date(as.character(data$Date.Format.Time9), format="%m/%d/%y")-
  as.Date(as.character(data$Date.Format.Time8), format="%m/%d/%y")

#Growth Normalized to dry weight the time prior

par(mfrow=c(2,4))
data$time1.growth <- (((data$Dry.Weigh.g.Time1 - data$Dry.Weigh.g.Time0)/(data$Dry.Weigh.g.Time0))/as.numeric(data$Days.Time1))*100 #calculate growth rate 
boxplot(data$time1.growth ~ data$Species, ylim=c(0,3), las=2)

data$time2.growth <- (((data$Dry.Weigh.g.Time2 - data$Dry.Weigh.g.Time1)/(data$Dry.Weigh.g.Time1))/as.numeric(data$Days.Time2))*100 #calculate growth rate 
boxplot(data$time2.growth ~ data$Species, ylim=c(0,3), las=2)

data$time3.growth <- (((data$Dry.Weigh.g.Time3 - data$Dry.Weigh.g.Time2)/(data$Dry.Weigh.g.Time2))/as.numeric(data$Days.Time3))*100 #calculate growth rate 
boxplot(data$time3.growth ~ data$Species, ylim=c(0,3), las=2)

data$time4.growth <- (((data$Dry.Weigh.g.Time4 - data$Dry.Weigh.g.Time3)/(data$Dry.Weigh.g.Time3))/as.numeric(data$Days.Time4))*100 #calculate growth rate 
boxplot(data$time4.growth ~ data$Species, ylim=c(0,3), las=2)

data$time5.growth <- (((data$Dry.Weigh.g.Time5 - data$Dry.Weigh.g.Time4)/(data$Dry.Weigh.g.Time4))/as.numeric(data$Days.Time5))*100 #calculate growth rate 
boxplot(data$time5.growth ~ data$Species, ylim=c(0,3), las=2)

data$time6.growth <- (((data$Dry.Weigh.g.Time6 - data$Dry.Weigh.g.Time5)/(data$Dry.Weigh.g.Time5))/as.numeric(data$Days.Time6))*100 #calculate growth rate 
boxplot(data$time6.growth ~ data$Species, ylim=c(0,3), las=2)

data$time7.growth <- (((data$Dry.Weigh.g.Time7 - data$Dry.Weigh.g.Time6)/(data$Dry.Weigh.g.Time6))/as.numeric(data$Days.Time7))*100 #calculate growth rate 
boxplot(data$time7.growth ~ data$Species, ylim=c(0,3), las=2)

data$time8.growth <- (((data$Dry.Weigh.g.Time8 - data$Dry.Weigh.g.Time7)/(data$Dry.Weigh.g.Time7))/as.numeric(data$Days.Time8))*100 #calculate growth rate 
boxplot(data$time8.growth ~ data$Species, ylim=c(0,3), las=2)

data$time9.growth <- (((data$Dry.Weigh.g.Time9 - data$Dry.Weigh.g.Time8)/(data$Dry.Weigh.g.Time8))/as.numeric(data$Days.Time9))*100 #calculate growth rate 
boxplot(data$time9.growth ~ data$Species, ylim=c(0,3), las=2)


range(na.omit(data$time1.growth))
#outs <- which(data$time1.growth==min((na.omit(data$time1.growth))))
#data <- data[-outs,]

range(na.omit(data$time2.growth))
#outs <- which(data$time2.growth==min((na.omit(data$time2.growth))))
#data <- data[-outs,]
#range(na.omit(data$time2.growth))

range(na.omit(data$time3.growth))
#outs <- which(data$time3.growth==min((na.omit(data$time3.growth))))
#data <- data[-outs,]
#range(na.omit(data$time3.growth))

range(na.omit(data$time4.growth))
#outs <- which(data$time4.growth==min((na.omit(data$time4.growth))))
#data <- data[-outs,]
#range(na.omit(data$time4.growth))

range(na.omit(data$time5.growth))
#outs <- which(data$time5.growth==min((na.omit(data$time5.growth))))
#data <- data[-outs,]
#range(na.omit(data$time5.growth))

range(na.omit(data$time6.growth))
#outs <- which(data$time6.growth==min((na.omit(data$time6.growth))))
#data <- data[-outs,]
#range(na.omit(data$time6.growth))

range(na.omit(data$time7.growth))
#outs <- which(data$time7.growth==min((na.omit(data$time7.growth))))
#data <- data[-outs,]
#range(na.omit(data$time7.growth))

range(na.omit(data$time8.growth))
#outs <- which(data$time8.growth==min((na.omit(data$time8.growth))))
#data <- data[-outs,]
#range(na.omit(data$time8.growth))

range(na.omit(data$time9.growth))
#outs <- which(data$time9.growth==min((na.omit(data$time9.growth))))
#data <- data[-outs,]
#range(na.omit(data$time9.growth))


trun.data <- data[,c(1,2,3,75:84)]

trun.data <- melt(trun.data)
trun.data <- subset(trun.data, variable!="PLUG.ID")

All.Means <- ddply(trun.data, c('variable','Species', 'Treatment'), summarize,
                   mean= mean(value, na.rm=T), #mean pnet
                   N = sum(!is.na(value)), # sample size
                   se = sd(value, na.rm=T)/sqrt(N)) #SE
All.Means

#All.Means <- na.omit(All.Means)
colnames(All.Means) <- c("Timepoint", "Species", "Treatment", "mean", "N", "se")

All.Means$Group <- paste(All.Means$Timepoint, All.Means$Treatment, All.Means$Species)
All.Means$SpGroup <- paste(All.Means$Treatment, All.Means$Species)

cols <- c("lightblue", "blue", "pink", "red")

Fig.All <- ggplot(All.Means, aes(x=Timepoint, y=mean, group=SpGroup)) + 
  geom_line(aes(linetype= Species, colour=Treatment, group=SpGroup), position = position_dodge(width = 0.1), alpha=0.5) + # colour, group both depend on cond2
  geom_errorbar(aes(ymin=All.Means$mean-All.Means$se, ymax=All.Means$mean+All.Means$se), colour="darkgray", width=0, size=0.5, position = position_dodge(width = 0.1)) +
  geom_point(aes(colour=Treatment, shape=Species), size = 1, position = position_dodge(width = 0.1)) +
  scale_colour_manual(values=cols) +
  #annotate("text", x=43, y=1.85, label = "a", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=132, y=2.15, label = "b", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=c(43,43), y=c(2.45,2.2), label = "ab", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=c(141,141), y=c(3.15,3.4), label = "c", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=35, y=0.5, label = ".", size = 1) + #add text to the graphic for posthoc letters
  xlab("Timepoint") +
  ylab(expression(paste("Growth % per day"))) +
  ylim(-0.5,1.5) +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank())+  #Set the plot background
  ggtitle("Growth") +
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0))
Fig.All

Gro.Figs <- arrangeGrob(Fig.All, ncol=1)
ggsave(file="Output/Growth_percent_BW.pdf", Gro.Figs, width = 4, height = 3, units = c("in"))


