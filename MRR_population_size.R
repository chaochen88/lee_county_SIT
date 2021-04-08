library(tidyverse)
library(knitr) # kable()
library(cowplot) # plot_grid() and theme_cowplot()
library(ggpubr) #stat_compare_means
library(maptools) #getKMLcoodenates()
library(scales)

#### Equations ####
# select data bank #
data.bank.mrr<-function(data, guide){
  a<-data%>%filter(round==guide, gender!="female", marking_status=="unmarked")%>%group_by(date)%>%summarise(unmarked=sum(n_adults))
  b<-data%>%filter(round==guide,gender!="female", marking_status=="marked")%>%group_by(date)%>%summarise(marked=sum(n_adults))
  mrr_data<-full_join(a,b,by="date")
  mrr_data[is.na(mrr_data)]<-0
  mrr_data$captured<-mrr_data$unmarked+mrr_data$marked
  mrr_data
}

# Lincoln Index #
lincoln<-function(data, released, plus=FALSE){
  if(plus == FALSE){
    marked<-data$marked
    unmarked<-data$unmarked
    captured<-data$captured
    n<-(released*unmarked)/marked
    result<-data.frame(days=seq(1,length(n)),lincoln=n)
    result
  }else{
    marked<-data$marked+1
    unmarked<-data$unmarked+1
    captured<-marked+unmarked
    n<-(released*unmarked)/marked
    result<-data.frame(days=seq(1,length(n)),lincoln=n)
    result
  }
}

# PDS&ALE #
pds.ale<-function(data){
  model<-lm(log10(marked+1)~date, data=data)
  PDS<-as.numeric(10^coef(model)[2])
  ALE<-1/-log(PDS)
  result<-data.frame(PDS=PDS, ALE=ALE)
  result
}

# Fisher-Ford #
fisher.ford<-function(data,released,plus=FALSE){
  model<-lm(log10(marked+1)~date, data=data)
  PDS<-as.numeric(10^coef(model)[2])
  if(plus == FALSE){
    marked<-data$marked
    captured<-data$captured
    n<-((PDS*captured*released)/marked)-released
    result<-data.frame(days=seq(1,length(n)),fisher=n)
    result
  }else{
    marked<-data$marked+1
    captured<-data$captured+1
    n<-((PDS*captured*released)/marked)-released
    result<-data.frame(days=seq(1,length(n)),fisher=n)
    result
  }
}

# Bailey # for low recaptures (?)
bailey<-function(data,released){
  marked<-data$marked+1
  captured<-data$captured+1
  model<-lm(log10(marked+1)~date, data=data)
  PDS<-as.numeric(10^coef(model)[2])
  n<-((PDS*captured*(released+1))/marked)-1-released
  result<-data.frame(days=seq(1,length(n)),bailey=n)
  result
}

# Cianci # No spatial component included in the model
cianci<-function(data, total_released, lambda, plus=FALSE) {
  if(plus == FALSE){
    marked<-data$marked
    unmarked<-data$unmarked
    captured<-marked+unmarked
    days<-as.numeric(factor(levels(factor(data$date))))
    pi<-marked/captured
    pi_1<-unmarked/captured
    beta<-log(lambda)
    a<-log10(pi/pi_1)
    b<-exp(beta*days)
    c<-log10(total_released)
    n<-10^(-a-b+c)
    result<-data.frame(days=seq(1,length(n)), cianci=n)
    result
  }else{
    marked<-data$marked+1
    unmarked<-data$unmarked+1
    captured<-marked+unmarked
    days<-as.numeric(factor(levels(factor(data$date))))
    pi<-marked/captured
    pi_1<-unmarked/captured
    beta<-log(lambda)
    a<-log10(pi/pi_1)
    b<-exp(beta*days)
    c<-log10(total_released)
    n<-10^(-a-b+c)
    result<-data.frame(days=seq(1,length(n)), cianci=n)
    result
  }
}

# Number released #
total.released<-function(data){
  result<-data%>%group_by(round)%>%
    summarise(total_dead=sum(no.pupae,no.adult),
              estimated=sum(estimates)/length(unique(evaluation)),
              survival=1-(total_dead/estimated),
              released=estimated*survival)
  result
}

##### Captiva Island Polygon ####
captiva_coordinates<-getKMLcoordinates(kmlfile = "./Field Activities/MRR/MRR outline area.kml",
                                       ignoreAltitude = T)
captiva_coord<-Polygon(captiva_coordinates)
captiva_poly<-Polygons(list(captiva_coord), ID = "Captiva Island")
captiva_poly1<-SpatialPolygons(list(captiva_poly),proj4string=CRS("+init=epsg:4326"))

captiva_coord_1<-as.data.frame(captiva_coordinates)
colnames(captiva_coord_1)<-c("long","lat")

#### Unique MRR features ####
release_point<-data.frame(long=-82.19055556, lat=26.52370278)
release_multipoints<-read.csv("./Field Activities/MRR/multipoints_lat_long_table.csv")
trap_position<-read.csv("./Field Activities/MRR/mrr_trap_location.csv")

#### Loading MRR Data ####
mrr_data<-read.csv("./Field Activities/MRR/MRR_all_collection.csv")
mrr_data$date<-as.Date(mrr_data$date, format='%Y-%m-%d')

mrr_release<-read.csv("./Field Activities/MRR/MRR_all_release.csv")

##### MRR Releases #####
m_released<-total.released(mrr_release)

#### MRR-A ####
mrr_data_A<-data.bank.mrr(mrr_data, "A")

#Probability Daily Survival - PDS#
#Average life expectancy#
pds.ale(mrr_data_A)

# Lincoln #
lincoln_mrr_A<-lincoln(mrr_data_A, m_released%>%filter(round=="A")%>%.$released)
mean(lincoln(mrr_data_A, m_released%>%filter(round=="A")%>%.$released)%>%.$lincoln)
lincoln_mrr_A

# Fisher-Ford #
fisher_ford_mrr_A<-fisher.ford(mrr_data_A, m_released%>%filter(round=="A")%>%.$released)
mean(fisher.ford(mrr_data_A, m_released%>%filter(round=="A")%>%.$released)%>%.$fisher)
fisher_ford_mrr_A

# Bailey #
bailey_A<-bailey(mrr_data_A, m_released%>%filter(round=="A")%>%.$released)
mean(bailey(mrr_data_A, m_released%>%filter(round=="A")%>%.$released)%>%.$bailey)
bailey_A

# Cianci #
cianci_A<-cianci(mrr_data_A, m_released%>%filter(round=="A")%>%.$released, 0.33)
cianci_A

pop_size_A<-merge(lincoln_mrr_A, fisher_ford_mrr_A)%>%merge(bailey_A)%>%merge(cianci_A)
pop_size_A$round<-"A"
summary(pop_size_A)

#### MRR-B ####
mrr_data_B<-data.bank.mrr(mrr_data, "B")

#Probability Daily Survival - PDS#
#Average life expectancy#
pds.ale(mrr_data_B)

# Lincoln #
lincoln_mrr_B<-lincoln(mrr_data_B, m_released%>%filter(round=="B")%>%.$released)
mean(lincoln(mrr_data_B, m_released%>%filter(round=="B")%>%.$released)%>%.$lincoln)
lincoln_mrr_B

# Fisher-Ford #
fisher_ford_mrr_B<-fisher.ford(mrr_data_B, m_released%>%filter(round=="B")%>%.$released)
mean(fisher.ford(mrr_data_B, m_released%>%filter(round=="B")%>%.$released)%>%.$fisher)
fisher_ford_mrr_B

# Bailey #
bailey_B<-bailey(mrr_data_B, m_released%>%filter(round=="B")%>%.$released)
mean(bailey(mrr_data_B, m_released%>%filter(round=="B")%>%.$released)%>%.$bailey)
bailey_B

# Cianci #
cianci_B<-cianci(mrr_data_B, m_released%>%filter(round=="B")%>%.$released, 0.33)
cianci_B

pop_size_B<-merge(lincoln_mrr_B, fisher_ford_mrr_B)%>%merge(bailey_B)%>%merge(cianci_B)
pop_size_B$round<-"B"
summary(pop_size_B)

#### MRR-C1 ####
mrr_data_C1<-data.bank.mrr(mrr_data, "C1")

#Probability Daily Survival - PDS#
#Average life expectancy#
pds.ale(mrr_data_C1)

# Lincoln #
lincoln_mrr_C1<-lincoln(mrr_data_C1, m_released%>%filter(round=="C1")%>%.$released)
mean(lincoln(mrr_data_C1, m_released%>%filter(round=="C1")%>%.$released)%>%.$lincoln)
lincoln_mrr_C1

# Fisher-Ford #
fisher_ford_mrr_C1<-fisher.ford(mrr_data_C1, m_released%>%filter(round=="C1")%>%.$released)
mean(fisher.ford(mrr_data_C1, m_released%>%filter(round=="C1")%>%.$released)%>%.$fisher)
fisher_ford_mrr_C1

# Bailey #
bailey_C1<-bailey(mrr_data_C1, m_released%>%filter(round=="C1")%>%.$released)
mean(bailey(mrr_data_C1, m_released%>%filter(round=="C1")%>%.$released)%>%.$bailey)
bailey_C1

# Cianci #
cianci_C1<-cianci(mrr_data_C1, m_released%>%filter(round=="C1")%>%.$released, 0.33)
mean(cianci(mrr_data_C1, m_released%>%filter(round=="C1")%>%.$released, 0.33)%>%.$cianci)
cianci_C1

pop_size_C1<-merge(lincoln_mrr_C1, fisher_ford_mrr_C1)%>%merge(bailey_C1)%>%merge(cianci_C1)
pop_size_C1$round<-"C1"
summary(pop_size_C1)

#### MRR-C2 ####
mrr_data_C2<-data.bank.mrr(mrr_data, "C2")

#Probability Daily Survival - PDS#
#Average life expectancy#
pds.ale(mrr_data_C2)

# Lincoln #
lincoln_mrr_C2<-lincoln(mrr_data_C2, m_released%>%filter(round=="C2")%>%.$released)
mean(lincoln(mrr_data_C2, m_released%>%filter(round=="C2")%>%.$released)%>%.$lincoln)
lincoln_mrr_C2

# Fisher-Ford #
fisher_ford_mrr_C2<-fisher.ford(mrr_data_C2, m_released%>%filter(round=="C2")%>%.$released)
mean(fisher.ford(mrr_data_C2, m_released%>%filter(round=="C2")%>%.$released)%>%.$fisher)
fisher_ford_mrr_C2

# Bailey #
bailey_C2<-bailey(mrr_data_C2, m_released%>%filter(round=="C2")%>%.$released)
mean(bailey(mrr_data_C2, m_released%>%filter(round=="C2")%>%.$released)%>%.$bailey)
bailey_C2

# Cianci #
cianci_C2<-cianci(mrr_data_C2, m_released%>%filter(round=="C2")%>%.$released, 0.33)
mean(cianci(mrr_data_C2, m_released%>%filter(round=="C2")%>%.$released, 0.33)%>%.$cianci)
cianci_C2

pop_size_C2<-merge(lincoln_mrr_C2, fisher_ford_mrr_C2)%>%merge(bailey_C2)%>%merge(cianci_C2)
pop_size_C2$round<-"C2"
summary(pop_size_C2)

#### MRR-D ####
mrr_data_D<-data.bank.mrr(mrr_data, "D")

#Probability Daily Survival - PDS#
#Average life expectancy#
pds.ale(mrr_data_D)

# Lincoln #
lincoln_mrr_D<-lincoln(mrr_data_D, m_released%>%filter(round=="D")%>%.$released, plu=TRUE)
mean(lincoln(mrr_data_D, m_released%>%filter(round=="D")%>%.$released, plus=TRUE)%>%.$lincoln)
lincoln_mrr_D

# Fisher-Ford #
fisher_ford_mrr_D<-fisher.ford(mrr_data_D, m_released%>%filter(round=="D")%>%.$released,plus = TRUE)
mean(fisher.ford(mrr_data_D, m_released%>%filter(round=="D")%>%.$released, plus = TRUE)%>%.$fisher)
fisher_ford_mrr_D

# Bailey #
bailey_D<-bailey(mrr_data_D, m_released%>%filter(round=="D")%>%.$released)
mean(bailey(mrr_data_D, m_released%>%filter(round=="D")%>%.$released)%>%.$bailey)
bailey_D

# Cianci #
cianci_D<-cianci(mrr_data_D, m_released%>%filter(round=="D")%>%.$released, 0.33, plus=TRUE)
mean(cianci(mrr_data_D, m_released%>%filter(round=="D")%>%.$released, 0.33, plus=TRUE)%>%.$cianci)
cianci_D

pop_size_D<-merge(lincoln_mrr_D, fisher_ford_mrr_D)%>%merge(bailey_D)%>%merge(cianci_D)
pop_size_D$round<-"D"
summary(pop_size_D)

#### MRR-E ####
mrr_data_E<-data.bank.mrr(mrr_data, "E")

#Probability Daily Survival - PDS#
#Average life expectancy#
pds.ale(mrr_data_E)

# Lincoln #
lincoln_mrr_E<-lincoln(mrr_data_E, m_released%>%filter(round=="E")%>%.$released)
mean(lincoln(mrr_data_E, m_released%>%filter(round=="E")%>%.$released)%>%.$lincoln)
lincoln_mrr_E

# Fisher-Ford #
fisher_ford_mrr_E<-fisher.ford(mrr_data_E, m_released%>%filter(round=="E")%>%.$released)
mean(fisher.ford(mrr_data_E, m_released%>%filter(round=="E")%>%.$released)%>%.$fisher)
fisher_ford_mrr_E

# Bailey #
bailey_E<-bailey(mrr_data_E, m_released%>%filter(round=="E")%>%.$released)
mean(bailey(mrr_data_E, m_released%>%filter(round=="E")%>%.$released)%>%.$bailey)
bailey_E

# Cianci #
cianci_E<-cianci(mrr_data_E, m_released%>%filter(round=="E")%>%.$released, 0.33)
mean(cianci(mrr_data_E, m_released%>%filter(round=="E")%>%.$released, 0.33)%>%.$cianci)
cianci_E

pop_size_E<-merge(lincoln_mrr_E, fisher_ford_mrr_E)%>%merge(bailey_E)%>%merge(cianci_E)
pop_size_E$round<-"E"
summary(pop_size_E)

#### MRR-F ####
mrr_data_F<-data.bank.mrr(mrr_data, "F")

#Probability Daily Survival - PDS#
#Average life expectancy#
pds.ale(mrr_data_F)

# Lincoln #
lincoln_mrr_F<-lincoln(mrr_data_F, m_released%>%filter(round=="F")%>%.$released, plus=TRUE)
mean(lincoln(mrr_data_F, m_released%>%filter(round=="F")%>%.$released, plus=TRUE)%>%.$lincoln)
lincoln_mrr_F

# Fisher-Ford #
fisher_ford_mrr_F<-fisher.ford(mrr_data_F, m_released%>%filter(round=="F")%>%.$released, plus=TRUE)
mean(fisher.ford(mrr_data_F, m_released%>%filter(round=="F")%>%.$released, plus=TRUE)%>%.$fisher)
fisher_ford_mrr_F

# Bailey #
bailey_F<-bailey(mrr_data_F, m_released%>%filter(round=="F")%>%.$released)
mean(bailey(mrr_data_F, m_released%>%filter(round=="F")%>%.$released)%>%.$bailey)
bailey_F

# Cianci #
cianci_F<-cianci(mrr_data_F, m_released%>%filter(round=="F")%>%.$released, 0.33, plus=TRUE)
mean(cianci(mrr_data_F, m_released%>%filter(round=="F")%>%.$released, 0.33, plus=TRUE)%>%.$cianci)
cianci_F

pop_size_F<-merge(lincoln_mrr_F, fisher_ford_mrr_F)%>%merge(bailey_F)%>%merge(cianci_F)
pop_size_F$round<-"F"
summary(pop_size_F)

#### MRR Summary ####
summary_table<-rbind(pop_size_A,pop_size_B,pop_size_C1,pop_size_C2,pop_size_D,pop_size_E,pop_size_F)
summary_table%>%ggplot()+
  geom_point(aes(x=days,y=lincoln), color="red")+
  geom_point(aes(x=days,y=fisher), color="blue")+
  geom_point(aes(x=days,y=bailey), color="green")+
  geom_point(aes(x=days,y=cianci), color="purple")+
  stat_summary(aes(x=days,y=lincoln), fun=mean, geom="line", size=1, color="red")+
  stat_summary(aes(x=days,y=fisher), fun=mean, geom="line", size=1, color="blue")+
  stat_summary(aes(x=days,y=bailey), fun=mean, geom="line", size=1, color="green")+
  stat_summary(aes(x=days,y=cianci), fun=mean, geom="line", size=1, color="purple")+
  scale_y_continuous(name="Population size")+
  theme_bw()+theme(strip.background = element_blank())+
  facet_wrap(.~round)