

library(readr)
library(readxl)
library(ggplot2)
library(tidyr)
library(plyr)
library(dplyr)
library(purrr)




data <- read_csv("MRR_all_collection_4_1.csv")

##exploring the data

#get total trapped adults in each round, we should expect a decay pattern by distance for marked but no pattern for unmarked
data1=data%>%group_by(round,marking_status)%>%summarise(adults1=sum(n_adults))

data2=data%>%group_by(round,marking_status,distance)%>%summarise(adults2=sum(n_adults))



data3=merge(data1,data2,by.y=c("round","marking_status"),all=T)

data3$probability=data3$adults2/data3$adults1

ggplot(data=data3,aes(x=distance,y=probability))+
  geom_point()+
  geom_smooth(se=F)+
  facet_wrap(~marking_status,scales = "free")




#####calculate MDT
data$annulus=data$distance
data$annulus=factor(data$annulus)
levels(data$annulus)
data$annulus=revalue(data$annulus, c("50"="1", "100"="2","150"="3","200"="4","250"="5","300"="6","400"="7"))

data1=data%>%group_by(round,day_collection,annulus,marking_status)%>%summarise(adults=sum(n_adults))



MDT=read_excel("MDT.xlsx")

data2=merge(data1[data1$marking_status=="marked",],MDT,by.y="annulus",all=T)

data2$ER=data2$adults*data2$CF/data2$number_traps


data3=data2%>%group_by(round,annulus,annulus_distance)%>%summarise(sum_ER=sum(ER))
data3$M=data3$annulus_distance*data3$sum_ER



data4=data3%>%group_by(round)%>%summarise(K=sum(sum_ER),N=sum(M))
data4$MDT=data4$N/data4$K








