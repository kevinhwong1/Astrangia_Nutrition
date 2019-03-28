## Dark Adapted Yield Test
## 20190328 

library(dbplyr)
library(tidyverse)
library(readr)
library(stringr)

#Import Data
D.Yield.test <- read.csv("../Data/Dark_Adapt_Yield_Test.csv")

ggplot(data=D.Yield.test, aes(x=Time.in.Dark, y=Fo, group=Coral.ID)) +
  geom_line()+
  geom_point()

ggplot(data=D.Yield.test, aes(x=Time.in.Dark, y=Fm, group=Coral.ID)) +
  geom_line()+
  geom_point()

ggplot(data=D.Yield.test, aes(x=Time.in.Dark, y=Y, group=Coral.ID)) +
  geom_line()+
  geom_point()
