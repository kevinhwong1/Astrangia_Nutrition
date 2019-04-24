#Title: Bleaching estimates from pictures
#Project: Astrangia Nutrition
#Author: HM Putnam
#Edited by: HM Putnam
#Date Last Modified: 20190424
#See Readme file for details

if ("dplyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('dplyr') 
if ("plotrix" %in% rownames(installed.packages()) == 'FALSE') install.packages('plotrix') 

library("dplyr")
library(plotrix)

rm(list=ls()) #clears workspace 


setwd("~/MyProjects/Astrangia_Nutrition/RAnalysis/") #set working directory

Calib.Data <-read.table("Data/20190403.csv", header=TRUE, sep=",", na.string="NA", as.is=TRUE) #reads in the data files
model <-lm(mVTris ~ TTris, data=Calib.Data) #runs a linear regression of mV as a function of temperature
coe <- coef(model) #extracts the coeffecients
R2<-summary(model)$r.squared
plot(mVTris ~ TTris, data=Calib.Data)
abline(lm(mVTris ~ TTris, data=Calib.Data))
legend('topleft', legend = bquote(R^2 == .(format(R2, digits = 3))), bty='n')
coe <- coef(model)
Slope <- coe[2]
Intercept <- coe[1]
#constants for use in pH calculation 
R <- 8.31447215 #gas constant in J mol-1 K-1 
F <-96485.339924 #Faraday constant in coulombs mol-1

#read in probe measurements of pH, temperature, and salinity from tanks
daily <- read.csv("Data/Daily_Measurements.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
min(daily$Temp.C)
max(daily$Temp.C)


mvTris <- daily$Temp.C*Slope+Intercept #calculate the mV of the tris standard using the temperature mv relationships in the measured standard curves 
STris<-35 #salinity of the Tris
phTris<- (11911.08-18.2499*STris-0.039336*STris^2)*(1/(daily$Temp.C+273.15))-366.27059+ 0.53993607*STris+0.00016329*STris^2+(64.52243-0.084041*STris)*log(daily$Temp.C+273.15)-0.11149858*(daily$Temp.C+273.15) #calculate the pH of the tris (Dickson A. G., Sabine C. L. and Christian J. R., SOP 6a)
daily$pH.Total<-phTris+(mvTris/1000-daily$pH.mV/1000)/(R*(daily$Temp.C+273.15)*log(10)/F) #calculate the pH on the total scale (Dickson A. G., Sabine C. L. and Christian J. R., SOP 6a)

write.table(daily,"Data/Total_pH.csv",sep=",", row.names=FALSE)


daily <- daily[70:nrow(daily),]

means <- daily %>%
  group_by(Treatment) %>% 
  summarise_each(funs(mean,std.error))

cols <- c("black", "orange", "green")

Fig.Temp <- ggplot(means, aes(x=Treatment, y=Temp.C_mean, group=Treatment)) + 
  geom_errorbar(aes(ymin=means$Temp.C_mean-means$Temp.C_std.error, ymax=means$Temp.C_mean+means$Temp.C_std.error), colour="black", width=0, size=0.5, position = position_dodge(width = 0.1)) +
  geom_point(aes(colour=Treatment, shape=Treatment), size = 2, position = position_dodge(width = 0.1)) +
  scale_colour_manual(values=cols) +
  #annotate("text", x=43, y=1.85, label = "a", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=132, y=2.15, label = "b", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=c(43,43), y=c(2.45,2.2), label = "ab", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=c(141,141), y=c(3.15,3.4), label = "c", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=35, y=0.5, label = ".", size = 1) + #add text to the graphic for posthoc letters
  xlab("Treatment") +
  ylab(expression(paste("Temperature °C"))) +
  ylim(18,20) +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.position= "right") + #remove legend background
  theme(plot.title = element_text(face = 'bold', 
                                  size = 10, 
                                  hjust = 0))
Fig.Temp

Fig.pH <- ggplot(means, aes(x=Treatment, y=pH.Total_mean, group=Treatment)) + 
  geom_errorbar(aes(ymin=means$pH.Total_mean-means$pH.Total_std.error, ymax=means$pH.Total_mean+means$pH.Total_std.error), colour="black", width=0, size=0.5, position = position_dodge(width = 0.1)) +
  geom_point(aes(colour=Treatment, shape=Treatment), size = 2, position = position_dodge(width = 0.1)) +
  scale_colour_manual(values=cols) +
  #annotate("text", x=43, y=1.85, label = "a", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=132, y=2.15, label = "b", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=c(43,43), y=c(2.45,2.2), label = "ab", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=c(141,141), y=c(3.15,3.4), label = "c", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=35, y=0.5, label = ".", size = 1) + #add text to the graphic for posthoc letters
  xlab("Treatment") +
  ylab(expression(paste("pH Total Scale"))) +
  ylim(7.5,8.1) +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.position= "right") + #remove legend background
  theme(plot.title = element_text(face = 'bold', 
                                  size = 10, 
                                  hjust = 0))
Fig.pH

Fig.Sal <- ggplot(means, aes(x=Treatment, y=Sal.psu_mean, group=Treatment)) + 
  geom_errorbar(aes(ymin=means$Sal.psu_mean-means$Sal.psu_std.error, ymax=means$Sal.psu_mean+means$Sal.psu_std.error), colour="black", width=0, size=0.5, position = position_dodge(width = 0.1)) +
  geom_point(aes(colour=Treatment, shape=Treatment), size = 2, position = position_dodge(width = 0.1)) +
  scale_colour_manual(values=cols) +
  #annotate("text", x=43, y=1.85, label = "a", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=132, y=2.15, label = "b", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=c(43,43), y=c(2.45,2.2), label = "ab", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=c(141,141), y=c(3.15,3.4), label = "c", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=35, y=0.5, label = ".", size = 1) + #add text to the graphic for posthoc letters
  xlab("Treatment") +
  ylab(expression(paste("Salinity psu"))) +
  ylim(28,31) +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.position= "right") + #remove legend background
  theme(plot.title = element_text(face = 'bold', 
                                  size = 10, 
                                  hjust = 0))
Fig.Sal

Fig.Light <- ggplot(means, aes(x=Treatment, y=Light.PAR_mean, group=Treatment)) + 
  geom_errorbar(aes(ymin=means$Light.PAR_mean-means$Light.PAR_std.error, ymax=means$Light.PAR_mean+means$Light.PAR_std.error), colour="black", width=0, size=0.5, position = position_dodge(width = 0.1)) +
  geom_point(aes(colour=Treatment, shape=Treatment), size = 2, position = position_dodge(width = 0.1)) +
  scale_colour_manual(values=cols) +
  #annotate("text", x=43, y=1.85, label = "a", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=132, y=2.15, label = "b", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=c(43,43), y=c(2.45,2.2), label = "ab", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=c(141,141), y=c(3.15,3.4), label = "c", size = 3) + #add text to the graphic for posthoc letters
  #annotate("text", x=35, y=0.5, label = ".", size = 1) + #add text to the graphic for posthoc letters
  xlab("Treatment") +
  ylab(expression(paste("Light µmol photons m-2 s-1"))) +
  #ylim(28,31) +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.position= "right") + #remove legend background
  theme(plot.title = element_text(face = 'bold', 
                                  size = 10, 
                                  hjust = 0))
Fig.Light

ggsave(file="Output/Photo_Bch.pdf", Fig.All, width = 4, height = 3, units = c("in"))

Physical.Data <- arrangeGrob(Fig.Temp, Fig.pH, Fig.Sal, Fig.Light, ncol=2)
ggsave(file="Output/Tank.Conditions.pdf", Physical.Data, width =12, height = 6, units = c("in"))

