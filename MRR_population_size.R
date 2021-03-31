library(tidyverse)
library(knitr) # kable()
library(cowplot) # plot_grid() and theme_cowplot()
library(ggpubr) #stat_compare_means
library(geosphere) #calculate de distance between the release point and trap distHaversine() function
library(maptools) #getKMLcoodenates()
library(raster) #pointDistance()
library(gifski) #gifski()
library(scales)

#Cianci - no spatial component included into the model
cianci<-function(marked, unmarked, total_released, lambda, days) {
  pi<-marked/(marked+unmarked)
  pi_1<-unmarked/(marked+unmarked)
  beta<-log(lambda)
  a<-log10(pi/pi_1)
  b<-exp(beta*days)
  c<-log10(total_released)
  d<-10^(-a-b+c)
  mean(d)
}

##### Captiva Island Polygon ####
captiva_coordinates<-getKMLcoordinates(kmlfile = "./MRR outline area.kml",
                                       ignoreAltitude = T)
captiva_coord<-Polygon(captiva_coordinates)
captiva_poly<-Polygons(list(captiva_coord), ID = "Captiva Island")
captiva_poly1<-SpatialPolygons(list(captiva_poly),proj4string=CRS("+init=epsg:4326"))

captiva_coord_1<-as.data.frame(captiva_coordinates)
colnames(captiva_coord_1)<-c("long","lat")

#### unique MRR features ####
release_point<-data.frame(long=-82.19055556, lat=26.52370278)
release_multipoints<-read.csv("./multipoints_lat_long_table.csv")
trap_position<-read.csv("./mrr_trap_location.csv")

mrr_data<-read.csv("./MRR_all_collection.csv")
mrr_data$date<-as.Date(mrr_data$date, format='%Y-%m-%d')

mrr_release<-read.csv("./MRR_all_release.csv")
MRR_A<-mrr_release%>%filter(round=="A")%>%summarise(total=sum(no.pupae,no.adult),
                                                    estimated=sum(estimates)/2,
                                                    survival=1-(total/estimated),
                                                    released=estimated*survival)
MRR_B<-mrr_release%>%filter(round=="B")%>%summarise(total=sum(no.pupae,no.adult),
                                                    estimated=sum(estimates),
                                                    survival=1-(total/estimated),
                                                    released=estimated*survival)
MRR_C1<-mrr_release%>%filter(round=="C1")%>%summarise(total=sum(no.pupae,no.adult),
                                                      estimated=sum(estimates)/2,
                                                      survival=1-(total/estimated),
                                                      released=estimated*survival)
MRR_C2<-mrr_release%>%filter(round=="C2")%>%summarise(total=sum(no.pupae,no.adult),
                                                      estimated=sum(estimates)/2,
                                                      survival=1-(total/estimated),
                                                      released=estimated*survival)
MRR_D<-mrr_release%>%filter(round=="D")%>%summarise(total=sum(no.pupae,no.adult),
                                                    estimated=sum(estimates)/2,
                                                    survival=1-(total/estimated),
                                                    released=estimated*survival)
MRR_E<-mrr_release%>%filter(round=="E")%>%summarise(total=sum(no.pupae,no.adult),
                                                    estimated=sum(estimates)/2,
                                                    survival=1-(total/estimated),
                                                    released=estimated*survival)
MRR_F<-mrr_release%>%filter(round=="F")%>%summarise(total=sum(no.pupae,no.adult),
                                                    estimated=sum(estimates)/2,
                                                    survival=1-(total/estimated),
                                                    released=estimated*survival)

mrr_release_summary<-data.frame(round=c("A","B","C1","C2","D","E","F"),
                                rbind(MRR_A, MRR_B, MRR_C1, MRR_C2, MRR_D, MRR_E, MRR_F))

##### MRR #####
mrr_data_A<-mrr_data%>%filter(round=="A")
mrr_data_B<-mrr_data%>%filter(round=="B")
mrr_data_C1<-mrr_data%>%filter(round=="C1")
mrr_data_C2<-mrr_data%>%filter(round=="C2")
mrr_data_D<-mrr_data%>%filter(round=="D")
mrr_data_E<-mrr_data%>%filter(round=="E")
mrr_data_F<-mrr_data%>%filter(round=="F")

#### Population Size ####
## MRR-A ##
mrr_data_A%>%group_by(marking_status,gender)%>%summarise(total=sum(n_adults))
mrr_data_A%>%filter(marking_status=="marked", gender=="male")%>%
  group_by(date)%>%summarise(total=sum(n_adults),
                             recapture=mean(total/MRR_A$released)*100)
pop_size_A<-mrr_data_A
pop_size_A$t<-as.integer(as.factor(pop_size_A$date))
#Probability Daily Survival - PDS#
mrr_pop_size_model_A<-lm(log10(n_adults+1)~date,
                         data=pop_size_A%>%filter(gender=="male", marking_status=="marked"))

PDS_mrr_A<-as.numeric(10^coef(mrr_pop_size_model_A)[2])

#Average life expectancy#
ALE_mrr_A<-1/-log(PDS_mrr_A)

# Total mosquitoes captured
unmarked_mrr_A<-pop_size_A%>%filter(gender=="male",marking_status=="unmarked")%>%
  group_by(t)%>%
  summarise(n_adults=sum(n_adults))

marked_mrr_A<-pop_size_A%>%filter(gender=="male",marking_status=="marked")%>%
  group_by(t)%>%
  summarise(n_adults=sum(n_adults))

# Lincoln #
lincoln_mrr_A<-round((MRR_A$released*sum(unmarked_mrr_A$n_adults))/sum(marked_mrr_A$n_adults),0)

# Fisher-Ford #
fisher_ford_mrr_A<-round(mean(((PDS_mrr_A*MRR_A$released*unmarked_mrr_A$n_adults)/
                           marked_mrr_A$n_adults)-(PDS_mrr_A*MRR_A$released)),0)
fisher_ford_mrr_A

# Bailey #
bailey_A<-round(mean(((PDS_mrr_A*(MRR_A$released+1)*(unmarked_mrr_A$n_adults+1))/
                        (marked_mrr_A$n_adults+1))-(PDS_mrr_A*MRR_A$released)),0)
bailey_A

# Cianci #
cianci(marked_mrr_A$n_adults,unmarked_mrr_A$n_adults,MRR_A$released,
       lambda = 0.33, days = marked_mrr_A$t)

## MRR-B ##
mrr_data_B%>%group_by(marking_status,gender)%>%summarise(total=sum(n_adults))
mrr_data_B%>%filter(gender=="male")%>%
  group_by(date, marking_status)%>%summarise(total=sum(n_adults),
                             recapture=mean(total/MRR_B$released)*100)
pop_size_B<-mrr_data_B
pop_size_B$t<-as.integer(as.factor(pop_size_B$date))
#Probability Daily Survival - PDS#
mrr_pop_size_model_B<-lm(log10(n_adults+1)~date,
                         data=pop_size_B%>%filter(gender=="male", marking_status=="marked"))

PDS_mrr_B<-as.numeric(10^coef(mrr_pop_size_model_B)[2])

#Average life expectancy#
ALE_mrr_B<-1/-log(PDS_mrr_B)

# Total mosquitoes captured
unmarked_mrr_B<-pop_size_B%>%filter(gender=="male",marking_status=="unmarked")%>%
  group_by(t)%>%
  summarise(n_adults=sum(n_adults))

marked_mrr_B<-pop_size_B%>%filter(gender=="male",marking_status=="marked")%>%
  group_by(t)%>%
  summarise(n_adults=sum(n_adults))

# Lincoln #
lincoln_mrr_B<-round((MRR_B$released*sum(unmarked_mrr_B$n_adults))/sum(marked_mrr_B$n_adults),0)

# Fisher-Ford #
fisher_ford_mrr_B<-round(mean(((PDS_mrr_B*MRR_B$released*unmarked_mrr_B$n_adults)/
                                 marked_mrr_B$n_adults)-(PDS_mrr_B*MRR_B$released)),0)
fisher_ford_mrr_B

# Bailey #
bailey_B<-round(mean(((PDS_mrr_B*(MRR_B$released+1)*(unmarked_mrr_B$n_adults+1))/
                        (marked_mrr_B$n_adults+1))-(PDS_mrr_B*MRR_B$released)),0)
bailey_B

# Cianci #
cianci(marked_mrr_B$n_adults,unmarked_mrr_B$n_adults,MRR_B$released,
       lambda = 0.33, days = marked_mrr_B$t)

## MRR-C1 ##
mrr_data_C1%>%group_by(marking_status,gender)%>%summarise(total=sum(n_adults))
mrr_data_C1%>%filter(marking_status=="marked", gender=="male")%>%
  group_by(date)%>%summarise(total=sum(n_adults),
                             recapture=mean(total/MRR_C1$released)*100)
pop_size_C1<-mrr_data_C1
pop_size_C1$t<-as.integer(as.factor(pop_size_C1$date))
#Probability Daily Survival - PDS#
mrr_pop_size_model_C1<-lm(log10(n_adults+1)~date,
                         data=pop_size_C1%>%filter(gender=="male", marking_status=="marked"))

PDS_mrr_C1<-as.numeric(10^coef(mrr_pop_size_model_C1)[2])

#Average life expectancy#
ALE_mrr_C1<-1/-log(PDS_mrr_C1)

# Total mosquitoes captured
unmarked_mrr_C1<-pop_size_C1%>%filter(gender=="male",marking_status=="unmarked")%>%
  group_by(t)%>%
  summarise(n_adults=sum(n_adults))

marked_mrr_C1<-pop_size_C1%>%filter(gender=="male",marking_status=="marked")%>%
  group_by(t)%>%
  summarise(n_adults=sum(n_adults))

# Lincoln #
lincoln_mrr_C1<-round((MRR_C1$released*sum(unmarked_mrr_C1$n_adults))/sum(marked_mrr_C1$n_adults),0)

# Fisher-Ford #
fisher_ford_mrr_C1<-round(mean(((PDS_mrr_C1*MRR_C1$released*unmarked_mrr_C1$n_adults)/
                                 marked_mrr_C1$n_adults)-(PDS_mrr_C1*MRR_C1$released)),0)
fisher_ford_mrr_C1

round(mean(((PDS_mrr_C1^11*MRR_C1$released*sum(unmarked_mrr_C1$n_adults))/
              sum(marked_mrr_C1$n_adults))-(PDS_mrr_C1^11*MRR_C1$released)),0)

# Bailey #
bailey_C1<-round(mean(((PDS_mrr_C1*(MRR_C1$released+1)*(unmarked_mrr_C1$n_adults+1))/
                        (marked_mrr_C1$n_adults+1))-(PDS_mrr_C1*MRR_C1$released)),0)
bailey_C1

cianci(marked_mrr_C1$n_adults,unmarked_mrr_C1$n_adults,MRR_C1$released,
       lambda = 0.33, days = marked_mrr_C1$t)

## MRR-C2 ##
mrr_data_C2%>%group_by(marking_status,gender)%>%summarise(total=sum(n_adults))
mrr_data_C2%>%filter(marking_status=="unmarked", gender=="male")%>%
  group_by(date)%>%summarise(total=sum(n_adults),
                             recapture=mean(total/MRR_C2$released)*100)
pop_size_C2<-mrr_data_C2
pop_size_C2$t<-as.integer(as.factor(pop_size_C2$date))
#Probability Daily Survival - PDS#
mrr_pop_size_model_C2<-lm(log10(n_adults+1)~date,
                         data=pop_size_C2%>%filter(gender=="male", marking_status=="marked"))

PDS_mrr_C2<-as.numeric(10^coef(mrr_pop_size_model_C2)[2])

#Average life expectancy#
ALE_mrr_C2<-1/-log(PDS_mrr_C2)

# Total mosquitoes captured
unmarked_mrr_C2<-pop_size_C2%>%filter(gender=="male",marking_status=="unmarked")%>%
  group_by(t)%>%
  summarise(n_adults=sum(n_adults))

marked_mrr_C2<-pop_size_C2%>%filter(gender=="male",marking_status=="marked")%>%
  group_by(t)%>%
  summarise(n_adults=sum(n_adults))

# Lincoln #
lincoln_mrr_C2<-round((MRR_C2$released*sum(unmarked_mrr_C2$n_adults))/sum(marked_mrr_C2$n_adults),0)

# Fisher-Ford #
fisher_ford_mrr_C2<-round(mean(((PDS_mrr_C2*MRR_C2$released*unmarked_mrr_C2$n_adults)/
                                 marked_mrr_C2$n_adults)-(PDS_mrr_C2*MRR_C2$released)),0)
fisher_ford_mrr_C2

# Bailey #
bailey_C2<-round(mean(((PDS_mrr_C2*(MRR_C2$released+1)*(unmarked_mrr_C2$n_adults+1))/
                        (marked_mrr_C2$n_adults+1))-(PDS_mrr_C2*MRR_C2$released)),0)
bailey_C2

cianci(marked_mrr_C2$n_adults,unmarked_mrr_C2$n_adults,MRR_C2$released,
       lambda = 0.33, days = marked_mrr_C2$t)

## MRR-D ##
mrr_data_D%>%group_by(marking_status,gender)%>%summarise(total=sum(n_adults))
mrr_data_D%>%filter(gender=="male")%>%
  group_by(date, marking_status)%>%summarise(total=sum(n_adults),
                             recapture=mean(total/MRR_D$released)*100)
pop_size_D<-mrr_data_D
pop_size_D$t<-as.integer(as.factor(pop_size_D$date))
#Probability Daily Survival - PDS#
mrr_pop_size_model_D<-lm(log10(n_adults+1)~date,
                         data=pop_size_D%>%filter(gender=="male", marking_status=="marked"))

PDS_mrr_D<-as.numeric(10^coef(mrr_pop_size_model_D)[2])

#Average life expectancy#
ALE_mrr_D<-1/-log(PDS_mrr_D)

# Total mosquitoes captured
unmarked_mrr_D<-pop_size_D%>%filter(gender=="male",marking_status=="unmarked", !date %in% c(as.Date("2020-08-19"), 
                                                                                            as.Date("2020-08-20")))%>%
  group_by(t)%>%
  summarise(n_adults=sum(n_adults))

marked_mrr_D<-pop_size_D%>%filter(gender=="male",marking_status=="marked")%>%
  group_by(t)%>%
  summarise(n_adults=sum(n_adults))

# Lincoln #
lincoln_mrr_D<-round((MRR_D$released*sum(unmarked_mrr_D$n_adults))/sum(marked_mrr_D$n_adults),0)

# Fisher-Ford #
fisher_ford_mrr_D<-round(mean(((PDS_mrr_D*MRR_D$released*unmarked_mrr_D$n_adults)/
                                 marked_mrr_D$n_adults)-(PDS_mrr_D*MRR_D$released)),0)
fisher_ford_mrr_D

# Bailey #
bailey_D<-round(mean(((PDS_mrr_D*(MRR_D$released+1)*(unmarked_mrr_D$n_adults+1))/
                        (marked_mrr_D$n_adults+1))-(PDS_mrr_D*MRR_D$released)),0)
bailey_D

cianci(marked_mrr_D$n_adults,unmarked_mrr_D$n_adults,MRR_D$released,
       lambda = 0.33, days = marked_mrr_D$t)

## MRR-E ##
mrr_data_E%>%group_by(marking_status,gender)%>%summarise(total=sum(n_adults))
mrr_data_E%>%filter(marking_status=="marked", gender=="male")%>%
  group_by(date)%>%summarise(total=sum(n_adults),
                             recapture=mean(total/MRR_E$released)*100)
pop_size_E<-mrr_data_E
pop_size_E$t<-as.integer(as.factor(pop_size_E$date))
#Probability Daily Survival - PDS#
mrr_pop_size_model_E<-lm(log10(n_adults+1)~date,
                          data=pop_size_E%>%filter(gender=="male", marking_status=="marked"))

PDS_mrr_E<-as.numeric(10^coef(mrr_pop_size_model_E)[2])

#Average life expectancy#
ALE_mrr_E<-1/-log(PDS_mrr_E)

# Total mosquitoes captured
unmarked_mrr_E<-pop_size_E%>%filter(gender=="male",marking_status=="unmarked")%>%
  group_by(t)%>%
  summarise(n_adults=sum(n_adults))

marked_mrr_E<-pop_size_E%>%filter(gender=="male",marking_status=="marked")%>%
  group_by(t)%>%
  summarise(n_adults=sum(n_adults))

# Lincoln #
lincoln_mrr_E<-round((MRR_E$released*sum(unmarked_mrr_E$n_adults))/sum(marked_mrr_E$n_adults),0)
lincoln_mrr_E

# Fisher-Ford #
fisher_ford_mrr_E<-round(mean(((PDS_mrr_E*MRR_E$released*unmarked_mrr_E$n_adults)/
                                  marked_mrr_E$n_adults)-(PDS_mrr_E*MRR_E$released)),0)
fisher_ford_mrr_E

# Bailey #
bailey_E<-round(mean(((PDS_mrr_E*(MRR_E$released+1)*(unmarked_mrr_E$n_adults+1))/
                         (marked_mrr_E$n_adults+1))-(PDS_mrr_E*MRR_E$released)),0)
bailey_E

cianci(marked_mrr_E$n_adults,unmarked_mrr_E$n_adults,MRR_E$released,
       lambda = 0.33, days = marked_mrr_E$t)

## MRR-F ##
mrr_data_F%>%group_by(marking_status,gender)%>%summarise(total=sum(n_adults))
mrr_data_F%>%filter(marking_status=="marked", gender=="male")%>%
  group_by(date)%>%summarise(total=sum(n_adults),
                             recapture=mean(total/MRR_F$released)*100)

mrr_data_F%>%filter(gender=="male")%>%
  group_by(marking_status, date)%>%summarise(total=sum(n_adults),
                             recapture=mean(total/MRR_F$released)*100)

pop_size_F<-mrr_data_F
pop_size_F$t<-as.integer(as.factor(pop_size_F$date))
#Probability Daily Survival - PDS#
mrr_pop_size_model_F<-lm(log10(n_adults+1)~date,
                         data=pop_size_F%>%filter(gender=="male", marking_status=="marked"))

PDS_mrr_F<-as.numeric(10^coef(mrr_pop_size_model_F)[2])

#Average life expectancy#
ALE_mrr_F<-1/-log(PDS_mrr_F)

# Total mosquitoes captured
unmarked_mrr_F<-pop_size_F%>%filter(gender=="male",marking_status=="unmarked")%>%
  group_by(t)%>%
  summarise(n_adults=sum(n_adults))

marked_mrr_F<-pop_size_F%>%filter(gender=="male",marking_status=="marked")%>%
  group_by(t)%>%
  summarise(n_adults=sum(n_adults))

status_marking<-merge(marked_mrr_F,unmarked_mrr_F, by="t", all.x=TRUE, all.y=TRUE)
status_marking[is.na(status_marking)]<-0
colnames(status_marking)<-c("t", "marked", "unmarked")

# Lincoln #
lincoln_mrr_F<-round((MRR_F$released*sum(unmarked_mrr_F$n_adults))/sum(marked_mrr_F$n_adults),0)
lincoln_mrr_F

# Fisher-Ford #
fisher_ford_mrr_F<-round(mean(((PDS_mrr_F*MRR_F$released*(status_marking$unmarked))/
                                 (status_marking$marked))-(PDS_mrr_F*MRR_F$released)),0)
fisher_ford_mrr_F

# Bailey #
bailey_F<-round(mean(((PDS_mrr_F*(MRR_F$released)*(status_marking$unmarked+1))/
                        (status_marking$marked+1))-(PDS_mrr_F*MRR_F$released)),0)
bailey_F

cianci(marked_mrr_F$n_adults,unmarked_mrr_F$n_adults,MRR_F$released,
       lambda = 0.33, days = marked_mrr_F$t)

#### MRR Summary ####
summary_table<-data.frame(MRR=rep(c("A", "B","C1","C2", "D", "E", "F"),
                                     each=3),
                          method=rep(c("Lincoln Index", "Fisher-Ford", "Bailey"),times=7),
                          size=c(lincoln_mrr_A,fisher_ford_mrr_A,bailey_A,
                                 lincoln_mrr_B,fisher_ford_mrr_B,bailey_B,
                                 lincoln_mrr_C1,fisher_ford_mrr_C1,bailey_C1,
                                 lincoln_mrr_C2,fisher_ford_mrr_C2,bailey_C2,
                                 lincoln_mrr_D,fisher_ford_mrr_D,bailey_D,
                                 lincoln_mrr_E,fisher_ford_mrr_E,bailey_E,
                                 lincoln_mrr_F,fisher_ford_mrr_F,bailey_F))
summary_table1<-data.frame(MRR_ID=c("MRR A", "MRR B","MRR C.1","MRR C.2","MRR D","MRR E","MRR F"),
                           PDS=round(c(PDS_mrr_A,PDS_mrr_B, PDS_mrr_C1,PDS_mrr_C2, PDS_mrr_D,PDS_mrr_E,PDS_mrr_F),2),
                           ALE=round(c(ALE_mrr_A,ALE_mrr_B, ALE_mrr_C1,ALE_mrr_C2, ALE_mrr_D,ALE_mrr_E,ALE_mrr_F),2),
                           Recapture=round(c(sum(mrr_1_data%>%
                                                   filter(marking_status=="marked",
                                                          gender=="male")%>%
                                                   .$n_adults)/MRR_I_released*100,
                                             sum(mrr_2_data%>%
                                                   filter(marking_status=="marked", gender=="male")%>%
                                                   .$n_adults)/MRR_II_released*100,
                                             sum(mrr_3.1_data%>%
                                                   filter(marking_status=="marked", gender=="male")%>%
                                                   .$n_adults)/MRR_III.1_released*100,
                                             sum(mrr_3.2_data%>%
                                                   filter(marking_status=="marked", gender=="male")%>%
                                                   .$n_adults)/MRR_III.2_released*100),2))

