library(tidyverse)
library(knitr) # kable()
library(cowplot) # plot_grid() and theme_cowplot()
library(ggpubr) #stat_compare_means
library(geosphere) #calculate de distance between the release point and trap distHaversine() function
library(maptools) #getKMLcoodenates()
library(raster) #pointDistance()
library(gifski) #gifski()
library(scales)
library(gganimate)

### FORMULAS ###
plot_by_day<-function(a,b,c,day="", titulo="", cores=c("red","orange")){
  data_1<-a
  data_2<-b%>%filter(date==day)
  data_3<-c
  ggplot()+
    geom_path(data = data_1, aes(x=long, y = lat), col="black")+
    geom_point(data=data_2, aes(x=long, y=lat, col=marking_status, size=n_adults), 
               position = position_jitter(width = .00001, height = .00001, seed = 1))+
    scale_size_continuous(name = "Recaptured", limits = c(0,max(data_2$n_adults)))+
    scale_color_manual(name="", labels=c("Marked", "Unmarked"),
                       values=cores)+
    geom_point(aes(x=data_3[[1]],y=data_3[[2]]), col="black", fill="black", size=3, shape=25)+
    theme_minimal(base_size = 20)+ggtitle(day)+
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          plot.title=element_text(margin=margin(b=5), size = 14, hjust = 0.15),
          plot.subtitle=element_text(margin=margin(t=5, b = 5)))
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
mrr_data$release<-as.Date(mrr_data$release, format='%Y-%m-%d')
mrr_data$date_release<-mrr_data$date-mrr_data$release

mrr_release<-read.csv("./MRR_all_release.csv")
mrr_release$date<-as.Date(mrr_release$date, format='%Y-%m-%d')

MRR_A<-mrr_release%>%group_by(date)%>%filter(round=="A")%>%summarise(total=sum(no.pupae,no.adult),
                                                    estimated=sum(estimates)/2,
                                                    survival=1-(total/estimated),
                                                    released=estimated*survival)
MRR_B<-mrr_release%>%group_by(date)%>%filter(round=="B")%>%summarise(total=sum(no.pupae,no.adult),
                                                    estimated=sum(estimates),
                                                    survival=1-(total/estimated),
                                                    released=estimated*survival)
MRR_C1<-mrr_release%>%group_by(date)%>%filter(round=="C1")%>%summarise(total=sum(no.pupae,no.adult),
                                                      estimated=sum(estimates)/2,
                                                      survival=1-(total/estimated),
                                                      released=estimated*survival)
MRR_C2<-mrr_release%>%group_by(date)%>%filter(round=="C2")%>%summarise(total=sum(no.pupae,no.adult),
                                                      estimated=sum(estimates)/2,
                                                      survival=1-(total/estimated),
                                                      released=estimated*survival)
MRR_D<-mrr_release%>%group_by(date)%>%filter(round=="D")%>%summarise(total=sum(no.pupae,no.adult),
                                                    estimated=sum(estimates)/2,
                                                    survival=1-(total/estimated),
                                                    released=estimated*survival)
MRR_E<-mrr_release%>%group_by(date)%>%filter(round=="E")%>%summarise(total=sum(no.pupae,no.adult),
                                                    estimated=sum(estimates)/2,
                                                    survival=1-(total/estimated),
                                                    released=estimated*survival)
MRR_F<-mrr_release%>%group_by(date)%>%filter(round=="F")%>%summarise(total=sum(no.pupae,no.adult),
                                                    estimated=sum(estimates)/2,
                                                    survival=1-(total/estimated),
                                                    released=estimated*survival)

mrr_release_summary<-data.frame(round=c("A","B","C1","C2","D","E","F"),
                                rbind(MRR_A, MRR_B, MRR_C1, MRR_C2, MRR_D, MRR_E, MRR_F))
mrr_release_summary$date<-as.Date(mrr_release_summary$date, format="%Y-%m-%d")
mrr_release_summary$facet_year<-as.numeric(format(mrr_release_summary$date,format="%Y"))
mrr_release_summary$facet_week<-as.numeric(format(mrr_release_summary$date,format="%U"))

##### MRR-A #####
mrr_data_A<-mrr_data%>%filter(round=="A")
mrr_data_A$ring<-as.factor(mrr_data_A$distance)

mrr_table_A<-mrr_data_A%>%filter(gender=="male", marking_status=="marked")%>%
  #group_by(date)%>%
  summarise(mean_recap=mean(n_adults),
            sem=sd(n_adults)/sqrt(length(n_adults)),
            sum_total=sum(n_adults),
            recap_rate=(sum(n_adults)/MRR_A$released))

kable(mrr_table_A)

mrr_table_A1<-mrr_data_A%>%filter(gender=="male")%>%
  group_by(marking_status)%>%
  summarise(mean_recap=mean(n_adults),
            sem=sd(n_adults)/sqrt(length(n_adults)),
            sum_total=sum(n_adults),
            recap_rate=(sum(n_adults)/MRR_A$released))

kable(mrr_table_A1)

distance_11<-pointDistance(release_point,cbind(mrr_data_A$long, mrr_data_A$lat),lonlat=TRUE)
distance_12<-distHaversine(release_point, cbind(mrr_data_A$long, mrr_data_A$lat))
distance_13<-distVincentyEllipsoid(release_point, cbind(mrr_data_A$long, mrr_data_A$lat))

mrr_data_A$recap_rate<-mrr_data_A$n_adults/MRR_A$released
mrr_data_A$distance_pnt_dist<-distance_11
mrr_data_A$distance_havesine<-distance_12
mrr_data_A$distance_vinellip<-distance_13
mrr_data_A$mean_geral_distan<-(distance_11+distance_12+distance_13)/3

mrr_data_A%>%filter(marking_status=="marked")%>%
  group_by(gender,ring)%>%
  summarise(sum_adults=sum(n_adults),
            recapture=sum_adults/MRR_A$released,
            mean_dist=mean(mean_geral_distan),
            min_dist=min(mean_geral_distan),
            max_dist=max(mean_geral_distan))

distancia_A<-mrr_data_A%>%filter(marking_status=="marked")%>%
  group_by(gender)%>%
  summarise(sum_adults=sum(n_adults),
            recapture=sum_adults/MRR_A$released*100,
            mean_dist=mean(mean_geral_distan),
            median=median(mean_geral_distan),
            min_dist=min(mean_geral_distan),
            max_dist=max(mean_geral_distan))

mrr_data_A_male<-mrr_data_A%>%filter(gender=="male")

day_11<-plot_by_day(captiva_coord_1, mrr_data_A_male, release_point, "2019-11-06",
                    titulo = "MRR-I", cores = c("#ffa500", "#A9A9A9"))
day_12<-plot_by_day(captiva_coord_1, mrr_data_A_male, release_point, "2019-11-07",
                    titulo = "MRR-I", cores = c("#ffa500", "#A9A9A9"))
day_13<-plot_by_day(captiva_coord_1, mrr_data_A_male, release_point, "2019-11-08",
                    titulo = "MRR-I", cores = c("#ffa500", "#A9A9A9"))
day_14<-plot_by_day(captiva_coord_1, mrr_data_A_male, release_point, "2019-11-09",
                    titulo = "MRR-I", cores = c("#ffa500", "#A9A9A9"))
day_15<-plot_by_day(captiva_coord_1, mrr_data_A_male, release_point, "2019-11-10",
                    titulo = "MRR-I", cores = c("#ffa500", "#A9A9A9"))
day_16<-plot_by_day(captiva_coord_1, mrr_data_A_male, release_point, "2019-11-11",
                    titulo = "MRR-I", cores = c("#ffa500", "#A9A9A9"))

mrr_data_A_plot<-ggplot()+ #this is not properly working using macOS Mojave
  geom_path(data = captiva_coord_1, aes(x=long, y = lat), col="black")+
  geom_point(data= release_point, aes(x=long,y=lat), fill="black",
             col="black", size=3, shape=25)+
  geom_point(data=mrr_data_A_male,
             aes(x=long, y=lat, col=marking_status, size=n_adults), 
             position = position_jitter(width = .00001, height = .00001, seed = 1))+
  labs(title = 'Colletion: {closest_state}') +
  transition_states(mrr_data_A_male$date,
    transition_length = 0.5)+
  ease_aes("sine-in-out")+
  scale_size_continuous(name = "Recaptured", 
                        limits = c(0,max(mrr_data_A$n_adults)))+
  scale_color_manual(name="", labels=c("Marked", "Unmarked"),values=c("#ffa500", "#A9A9A9"))+
  theme_minimal(base_size = 20)+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.position = "bottom",
        plot.title=element_text(margin=margin(b=5), size = 14, hjust = 0.15),
        plot.subtitle=element_text(margin=margin(t=5, b = 5)))

#animate(mrr_data_A_plot, fps = 10, width = 750, height = 750)
anim_save("./mrr_A.gif")

###
MRR_A_plot<-ggarrange(day_11,day_12,day_13,day_14,day_15,day_16, nrow=2, ncol=3, align="hv",
                      legend="bottom", labels="AUTO", common.legend = TRUE)
MRR_A_plot
#ggsave("./Field Activities/MRR/MRR/MRR - I/MRR_A_plot.pdf", MRR_A_plot,
#       dpi=1000, width=36, height=32, units="cm")

###
MRR_A_plot1<-ggplot(mrr_data_A_male, aes(x=date, y=n_adults, col=marking_status))+
  geom_point(position = position_jitter(width=0.1,height=0.1,seed=1), alpha=0.5)+
  ggtitle("MRR-A")+
  scale_color_manual(name="Marking status", values = c("#ffa500", "#A9A9A9"),
                     label=c("Marked","Unmarked"))+
  scale_y_continuous(trans = "log10", name = "Recaptured males (log10)")+
  scale_x_date(name = "Recapture date", date_breaks = "1 day",date_labels = "%m/%d")+
  geom_smooth(method = "glm", method.args=list(family="poisson"), se = TRUE, alpha=0.25)+
  theme_bw(base_size = 20)+theme(legend.position="bottom")
MRR_A_plot1
#ggsave("./Field Activities/MRR/MRR/MRR_A_plot1.png",MRR_A_plot1,
#       dpi = 1000, width = 28, height = 18, units = "cm")

###
MRR_A_plot2<-ggplot(mrr_data_A%>%filter(marking_status=="marked", gender=="male"), 
                    aes(x=date, y=mean_geral_distan, col=n_adults))+
  geom_point(position = position_jitter(width=0.2,height=0.2,seed=1))+
  labs(x="Recapture date", y="Distance from release point (m)")+ggtitle("MRR-A")+
  scale_color_gradient(low = "#ffc04d", high = "#ff2500", na.value = NA, 
                       name="Recaptured")+
  scale_x_date(name = "Recapture date", date_breaks = "1 day",date_labels = "%m/%d")+
  geom_smooth(method = "glm", method.args=list(family="poisson"),
              se = TRUE, col="#ffa500")+
  theme_bw(base_size = 20)+
  theme(legend.position = c(0.7, 0.90),legend.background = element_rect(fill=NA))+
  guides(color = guide_colourbar(barwidth = 8, barheight = 1, direction = "horizontal"))
MRR_A_plot2
#ggsave("./Field Activities/MRR/MRR/MRR - I/MRR_A_plot2.pdf", MRR_A_plot2,
#       dpi=1000, width=28, height=18, units="cm")

###
MRR_data_A_plot3<-mrr_data_A%>%filter(marking_status=="marked")%>%
  ggplot()+
  geom_path(data = captiva_coord_1, aes(x=long, y = lat), col="black")+
  geom_point(data=release_point, aes(x=long, y=lat),
             col="black",fill="black", size=0.85, shape=4)+
  geom_point(data=mrr_data_A%>%filter(gender=="male"),
             aes(x=long, y=lat, size=n_adults, color=marking_status),
             position = position_jitter(width = .00001, height = .00001, seed = 1), alpha=0.35)+
  scale_color_manual(name="", 
                     values = c("#a9a9a9","#ffa500"),
                     label=c("Unmarked", "Marked"))+
  scale_size_area(limits = c(0,max(mrr_data_A%>%#filter(marking_status=="marked")%>%
                                     .$n_adults)), max_size = 10,
                  breaks = seq(0,max(mrr_data_A%>%#filter(marking_status=="marked")%>%
                                       .$n_adults),by=10),
                  name="")+
  theme_minimal(base_size = 20)+ggtitle("MRR-A Singlepoint release")+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        plot.title=element_text(margin=margin(b=5), size = 14, hjust = 0.15),
        plot.subtitle=element_text(margin=margin(t=5, b = 5)))+
  facet_wrap(.~date)

MRR_data_A_plot3
#ggsave("./Field Activities/MRR/MRR/MRR-A_plot3.png",MRR_data_A_plot3, dpi = 1000, width = 24, 
#       height = 18, units = "cm")

###
MRR_A_plot4<-ggplot(mrr_data_A%>%filter(marking_status=="marked", gender=="male"),
                    aes(x=date, y=recap_rate, col=mean_geral_distan))+
  geom_point(position = position_jitter(width=0.25,seed=1))+
  ggtitle("MRR-A")+
  scale_color_gradient(low = "#ffc04d", high = "#ff2500", na.value = NA, name="Dispersion (m)")+
  scale_y_continuous(name = "Recaptured males", limits=c(0,0.005),
                     labels = scales::percent_format(accuracy = .02))+
  scale_x_date(name = "Recapture date", date_breaks = "1 day",date_labels = "%m/%d")+
  geom_smooth(method = "glm", method.args=list(family="binomial"), se = FALSE, col="orange")+
  theme_bw(base_size = 20)+theme(legend.position = c(0.7, 0.90),
                                 legend.background = element_rect(fill=NA))+
  guides(color = guide_colourbar(barwidth = 10, barheight = 1, direction = "horizontal"))
MRR_A_plot4
#ggsave("./Field Activities/MRR/MRR/MRR - I/MRR_A_plot4.pdf", MRR_A_plot4, 
#       dpi=1000, width=28, height=18, units="cm")

##### MRR-B #####
mrr_data_B<-mrr_data%>%filter(round=="B")
mrr_data_B$ring<-as.factor(mrr_data_B$distance)

mrr_table_B<-mrr_data_B%>%filter(gender=="male",marking_status=="marked")%>%
  #group_by(date)%>%
  summarise(mean_recap=mean(n_adults),
            sem=sd(n_adults)/sqrt(length(n_adults)),
            sum_total=sum(n_adults),
            recap=(sum(n_adults)/MRR_B$released))

kable(mrr_table_B)

mrr_table_B1<-mrr_data_B%>%filter(gender=="male")%>%
  group_by(marking_status)%>%
  summarise(mean_recap=mean(n_adults),
            sem=sd(n_adults)/sqrt(length(n_adults)),
            sum_total=sum(n_adults),
            recap_rate=(sum(n_adults)/MRR_B$released))

kable(mrr_table_B1)

distance_21<-pointDistance(release_point,cbind(mrr_data_B$long, mrr_data_B$lat),lonlat=TRUE)
distance_22<-distHaversine(release_point, cbind(mrr_data_B$long, mrr_data_B$lat))
distance_23<-distVincentyEllipsoid(release_point, cbind(mrr_data_B$long, mrr_data_B$lat))

mrr_data_B$recap_rate<-mrr_data_B$n_adults/MRR_B$released
mrr_data_B$distance_pnt_dist<-distance_21
mrr_data_B$distance_havesine<-distance_22
mrr_data_B$distance_vinellip<-distance_23
mrr_data_B$mean_geral_distan<-(distance_21+distance_22+distance_23)/3

distancia_B<-mrr_data_B%>%filter(marking_status=="marked")%>%
  group_by(gender)%>%
  summarise(sum_adults=sum(n_adults),
            recapture=sum_adults/MRR_B$released,
            mean_dist=mean(mean_geral_distan),
            median=median(mean_geral_distan),
            min_dist=min(mean_geral_distan),
            max_dist=max(mean_geral_distan))

mrr_data_B%>%filter(marking_status=="marked")%>%
  group_by(date)%>%
  summarise(sum_adults=sum(n_adults),
            recapture=sum_adults/MRR_B$released*100,
            mean_dist=mean(mean_geral_distan),
            min_dist=min(mean_geral_distan),
            max_dist=max(mean_geral_distan))

day_21<-plot_by_day(captiva_coord_1, mrr_data_B, release_point, "2020-03-09",
                    titulo = "MRR-II", cores = c("#00b300", "#A9A9A9"))
day_22<-plot_by_day(captiva_coord_1, mrr_data_B, release_point, "2020-03-10",
                    titulo = "MRR-II", cores = c("#00b300", "#A9A9A9"))
day_23<-plot_by_day(captiva_coord_1, mrr_data_B, release_point, "2020-03-11",
                    titulo = "MRR-II", cores = c("#00b300", "#A9A9A9"))
day_24<-plot_by_day(captiva_coord_1, mrr_data_B, release_point, "2020-03-12",
                    titulo = "MRR-II", cores = c("#00b300", "#A9A9A9"))
day_25<-plot_by_day(captiva_coord_1, mrr_data_B, release_point, "2020-03-13",
                    titulo = "MRR-II", cores = c("#00b300", "#A9A9A9"))
day_26<-plot_by_day(captiva_coord_1, mrr_data_B, release_point, "2020-03-14",
                    titulo = "MRR-II", cores = c("#00b300", "#A9A9A9"))
day_27<-plot_by_day(captiva_coord_1, mrr_data_B, release_point, "2020-03-15",
                    titulo = "MRR-II", cores = c("#00b300", "#A9A9A9"))
day_28<-plot_by_day(captiva_coord_1, mrr_data_B, release_point, "2020-03-16",
                    titulo = "MRR-II", cores = c("#00b300", "#A9A9A9"))
day_29<-plot_by_day(captiva_coord_1, mrr_data_B, release_point, "2020-03-17",
                    titulo = "MRR-II", cores = c("#00b300", "#A9A9A9"))

mrr_data_B_male<-mrr_data_B%>%filter(gender=="male")

mrr_data_B_plot<-ggplot()+ # this is not working using macOS Mojave
  geom_path(data = captiva_coord_1, aes(x=long, y = lat), col="black")+
  geom_point(data= release_point, aes(x=long,y=lat), fill="black",
             col="black", size=3, shape=25)+
  geom_point(data=mrr_data_B_male,
             aes(x=long, y=lat, col=marking_status, size=n_adults), 
             position = position_jitter(width = .00001, height = .00001, seed = 1))+
  labs(title = 'Colletion: {closest_state}') +
  transition_states(mrr_data_B_male$date,
                    transition_length = 0.5)+
  ease_aes("sine-in-out")+
  scale_size_continuous(name = "Recaptured", 
                        limits = c(0,max(mrr_data_B_male$n_adults)))+
  scale_color_manual(name="", labels=c("Marked", "Unmarked"),values=c("#00b300", "#A9A9A9"))+
  theme_minimal(base_size = 20)+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.position = "bottom",
        plot.title=element_text(margin=margin(b=5), size = 14, hjust = 0.15),
        plot.subtitle=element_text(margin=margin(t=5, b = 5)))

#animate(mrr_data_B_plot, fps = 10, width = 750, height = 750)
#anim_save("./Field Activities/MRR/MRR/MRR - II/mrr_2.gif")

###
MRR_data_B_plot<-ggarrange(day_21,day_22,day_23,day_24,day_25,day_26,day_27,day_28,day_29,
                       nrow=2, ncol=5, align="hv",
                      legend="bottom", labels="AUTO", common.legend = TRUE)
MRR_data_B_plot
#ggsave("./Field Activities/MRR/MRR/MRR - II/MRR-II_plot.pdf", MRR_II_plot, 
#       dpi=1000, width=56, height=32, units="cm")

###
MRR_data_B_plot1<-ggplot(mrr_data_B, aes(x=date, y=n_adults, col=marking_status))+
  geom_point(position = position_jitter(width=0.1,height=0.1,seed=1), alpha=0.5)+
  ggtitle("MRR-B")+
  scale_color_manual(name="Marking status", values = c("#00b300", "#A9A9A9"))+
  scale_y_continuous(name="Recaptured males (log10)",trans = "log10")+
  scale_x_date(name = "Recapture date", date_breaks = "2 day",date_labels = "%m/%d")+
  geom_smooth(method = "glm", method.args=list(family="poisson"), se = TRUE, alpha=0.25)+
  theme_bw(base_size = 20)+theme(legend.position="bottom")
MRR_data_B_plot1
#ggsave("./Field Activities/MRR/MRR/MRR_B_plot1.png",MRR_data_B_plot1, dpi = 1000,
#       width = 28, height = 18, units = "cm")

###
MRR_data_B_plot2<-ggplot(mrr_data_B%>%filter(marking_status=="marked"),
                     aes(x=date, y=mean_geral_distan, col=n_adults))+
  geom_point(position = position_jitter(width=0.2,height=0.2,seed=1))+
  labs(y="Distance from release point (m)")+ggtitle("MRR-B")+
  scale_color_gradient(low = "#00ff00", high = "#006700", na.value = NA, name="Recaptured")+
  scale_x_date(name = "Recapture date", date_breaks = "2 day",date_labels = "%m/%d")+
  geom_smooth(method = "glm", method.args=list(family="poisson"),
              se = TRUE, col="#00b300")+
  theme_bw(base_size = 20)+
  theme(legend.position = c(0.6, 0.90),legend.background = element_rect(fill=NA))+
  guides(color = guide_colourbar(barwidth = 8, barheight = 1, direction = "horizontal"))
MRR_data_B_plot2
#ggsave("./Field Activities/MRR/MRR/MRR - II/MRR_B_plot2.pdf", MRR_data_B_plot2, dpi=1000, width=28, height=18, units="cm")

###
MRR_data_B_plot3<-mrr_data_B%>%filter(marking_status=="marked")%>%
  ggplot()+
  geom_path(data = captiva_coord_1, aes(x=long, y = lat), col="black")+
  geom_point(data=release_point, aes(x=long, y=lat),
             col="black",fill="black", size=0.85, shape=4)+
  geom_point(data=mrr_data_B%>%filter(gender=="male"),
             aes(x=long, y=lat, size=n_adults, color=marking_status),
             position = position_jitter(width = .00001, height = .00001, seed = 1), alpha=0.35)+
  scale_color_manual(name="", 
                     values = c("#a9a9a9","#00b300"),
                     label=c("Unmarked", "Marked"))+
  scale_size_area(limits = c(0,max(mrr_data_B%>%#filter(marking_status=="marked")%>%
                                     .$n_adults)), max_size = 10,
                  breaks = seq(0,max(mrr_data_B%>%#filter(marking_status=="marked")%>%
                                       .$n_adults),by=5),
                  name="")+
  theme_minimal(base_size = 20)+ggtitle("MRR-B Singlepoint Release")+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        plot.title=element_text(margin=margin(b=5), size = 14, hjust = 0.15),
        plot.subtitle=element_text(margin=margin(t=5, b = 5)))+
  facet_wrap(.~date)
MRR_data_B_plot3
#ggsave("./Field Activities/MRR/MRR/MRR_B_plot3.png",MRR_data_B_plot3,
#       dpi = 1000, width = 22, height = 20, units = "cm")

###
MRR_data_B_plot4<-ggplot(mrr_data_B%>%filter(marking_status=="marked"), 
                     aes(x=date, y=recap_rate*100, col=mean_geral_distan))+
  geom_point(position = position_jitter(width=0.25,seed=1))+
  ggtitle("MRR-B")+
  scale_color_gradient(low = "#00ff00", high = "#006700", na.value = NA,
                       name="Dispersion (m)")+
  scale_y_continuous(name="Recaptured males (%)")+
  scale_x_date(name = "Recapture date", date_breaks = "2 day",date_labels = "%m/%d")+
  geom_smooth(method = "glm", method.args=list(family="quasibinomial"),
              se = TRUE, col="#00b300", alpha=0.25)+
  theme_bw(base_size = 20)+theme(legend.position = c(0.7, 0.90),
                                 legend.background = element_rect(fill=NA))+
  guides(color = guide_colourbar(barwidth = 10, barheight = 1,direction = "horizontal"))
MRR_data_B_plot4
#ggsave("./Field Activities/MRR/MRR/MRR - II/MRR_B_plot4.pdf", MRR_data_B_plot4, dpi=1000, width=28, height=18, units="cm")

##### MRR-C1 Purple #####
mrr_data_C1<-mrr_data%>%filter(round=="C1")

mrr_data_C1%>%group_by(marking_status, date)%>%
  summarise(min=min(n_adults),
            max=max(n_adults),
            mean=mean(n_adults),
            sum=sum(n_adults),
            length=length(n_adults))

mrr_table_C1<-mrr_data_C1%>%filter(gender=="male", marking_status=="marked")%>%
  #group_by(date)%>%
  summarise(mean_adult=mean(n_adults),
            sem=sd(n_adults)/sqrt(length(n_adults)),
            sum_total=sum(n_adults),
            recap_rate=(sum(n_adults)/MRR_C1$released))

kable(mrr_table_C1)

mrr_table_C11<-mrr_data_C1%>%filter(gender=="male")%>%
  group_by(marking_status)%>%
  summarise(mean_recap=mean(n_adults),
            sem=sd(n_adults)/sqrt(length(n_adults)),
            sum_total=sum(n_adults),
            recap_rate=(sum(n_adults)/MRR_C1$released))

kable(mrr_table_C11)

mrr_data_C1$recap_rate<-mrr_data_C1$n_adults/MRR_C1$released

mrr_data_C1_male<-mrr_data_C1%>%
  filter(gender=="male")%>%
  group_by(date,trap_ID, long, lat, marking_status)%>%
  summarise(n_adults=sum(n_adults))

mrr_data_C1_plot<-ggplot()+ # this is not working using macOS Mojave
  geom_path(data = captiva_coord_1, aes(x=long, y = lat), col="black")+
  geom_point(data=release_multipoints, aes(x=long, y=lat),
             col="black",fill="black", size=1, shape=25, alpha=0.25)+
  geom_point(data=mrr_data_C1_male,
             aes(x=long, y=lat, col=marking_status, size=n_adults), 
             position = position_jitterdodge(jitter.width = .00001, jitter.height = .00001, dodge.width = .0001, seed = 1))+
  scale_size_continuous(name = "Recaptured", 
                        limits = c(0,max(mrr_data_C1_male$n_adults)),
                        breaks = seq(0,max(mrr_data_C1_male$n_adults),by=50))+
  scale_color_manual(name="", labels=c("Marked", "Unmarked"),values=c("#8b0000", "#008b8b"))+
  theme_minimal(base_size = 14)+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.position = "bottom", legend.box = "vertical",
        plot.title=element_text(margin=margin(b=5), size = 14, hjust = 0.15),
        plot.subtitle=element_text(margin=margin(t=5, b = 5)))+
  labs(title = 'Colletion: {closest_state}') +
  transition_states(mrr_data_C1_male$date,
                    transition_length = 0.5)+
  ease_aes("sine-in-out")

#animate(mrr_data_C1_plot, fps = 5, width = 250, height = 350)
#anim_save("./Field Activities/MRR/MRR/MRR - III/mrr_C1.gif")


### Graphs ###
MRR_data_C1_plot1<-ggplot(mrr_data_C1, aes(x=date, y=n_adults, col=marking_status))+
  geom_point(position = position_jitter(width=0.1,height=0.1,seed=1), alpha=0.5)+
  labs(title = "MRR-C1")+
  scale_color_manual(name="Marking status", values = c("#ffa500", "#A9A9A9"))+
  scale_y_continuous(trans = "log10", name = "Recaptured males (log10)")+
  scale_x_date(name = "Recapture date", date_breaks = "2 day",date_labels = "%m/%d")+
  geom_smooth(method = "glm", method.args=list(family="poisson"), se = TRUE, alpha=0.25)+
  theme_bw(base_size = 20)+theme(legend.position="bottom")
MRR_data_C1_plot1
#ggsave("./Field Activities/MRR/MRR/MRR_C1_plot1.png",MRR_data_C1_plot1, dpi = 1000,
#      width = 28, height = 18, units = "cm")

###
MRR_data_C1_plot2<-ggplot(mrr_data_C1, aes(x=date, y=n_adults, col=marking_status))+
  geom_point(position = position_jitter(width=0.08,height=0.08,seed=1))+
  labs(x="Recapture date", y="Recaptured males (log10)")+ggtitle("MRR-C1 Purple multipoints")+
  scale_color_manual(name="Marking status", values = c("#ffa500", "#A9A9A9"))+
  scale_x_date(date_breaks = "2 day",date_labels = "%m/%d")+
  scale_y_continuous(trans = "log10")+
  geom_smooth(method = "glm", method.args=list(family="poisson"),
              se = TRUE, alpha=0.25)+
  theme_bw(base_size = 20)+theme(legend.position="bottom")+facet_grid(.~gender)
MRR_data_C1_plot2
#ggsave("./Field Activities/MRR/MRR/MRR - III/MRR_C1_plot2.pdf",MRR_data_C1_plot2, dpi = 1000,
#       width = 28, height = 18, units = "cm")

###
MRR_data_C1_plot3<-mrr_data_C1%>%#filter(marking_status=="marked")%>%
  ggplot()+
  geom_path(data = captiva_coord_1, aes(x=long, y = lat), col="black")+
  geom_point(data=release_multipoints, aes(x=long, y=lat),
             col="black",fill="black", size=0.85, shape=4)+
  geom_point(data=mrr_data_C1%>%filter(gender=="male"),
             aes(x=long, y=lat, size=n_adults, color=marking_status),
             position = position_jitter(width = .00001, height = .00001, seed = 1), alpha=0.35)+
  scale_color_manual(name="", 
                     values = c("#a9a9a9","#ffa500"),
                     label=c("Unmarked", "Marked"))+
  scale_size_area(limits = c(0,max(mrr_data_C1%>%#filter(marking_status=="marked")%>%
                                     .$n_adults)), max_size = 8,
                  breaks = seq(0,max(mrr_data_C1%>%#filter(marking_status=="marked")%>%
                                       .$n_adults), by=50),
                  name="")+
  theme_minimal(base_size = 20)+ggtitle("MRR-C1")+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        plot.title=element_text(margin=margin(b=5), size = 14, hjust = 0.15),
        plot.subtitle=element_text(margin=margin(t=5, b = 5)))+
  facet_wrap(.~date)
MRR_data_C1_plot3
#ggsave("./Field Activities/MRR/MRR/MRR-C1.1_plot3.png", MRR_data_C1_plot3, dpi=1000, width=28,
#       height=24, units="cm")

#### MRR-C2 Multicolor #####
mrr_data_C2<-mrr_data%>%filter(round=="C2")

mrr_data_C2%>%group_by(marking_status, date)%>%
  summarise(min=min(n_adults),
            max=max(n_adults),
            mean=mean(n_adults),
            sum=sum(n_adults),
            length=length(n_adults))

mrr_table_C2<-mrr_data_C2%>%filter(gender=="male", marking_status!="marked")%>%
  group_by(marking_status,date)%>%
  summarise(mean_adult=mean(n_adults),
            sem=sd(n_adults)/sqrt(length(n_adults)),
            sum_total=sum(n_adults),
            reca_rate=(sum(n_adults)/MRR_C2$released))

kable(mrr_table_C2)

mrr_table_C21<-mrr_data_C2%>%filter(gender=="male")%>%
  group_by(marking_status)%>%
  summarise(mean_recap=mean(n_adults),
            sem=sd(n_adults)/sqrt(length(n_adults)),
            sum_total=sum(n_adults),
            recap_rate=(sum(n_adults)/MRR_C2$released))

kable(mrr_table_C21)

mrr_data_C2$recap_rate<-mrr_data_C2$n_adults/MRR_C2$released

mrr_data_C2_male<-mrr_data_C2%>%
  filter(gender=="male")%>%
  group_by(date,trap_ID, long, lat, marking_status)%>%
  summarise(n_adults=sum(n_adults))

mrr_data_C2_plot<-ggplot()+
  geom_path(data = captiva_coord_1, aes(x=long, y = lat), col="black")+
  geom_point(data=release_multipoints, aes(x=long, y=lat),
             col="black",fill="black", size=4, shape=25)+
  geom_point(data=mrr_data_C2_male,
             aes(x=long, y=lat, col=marking_status, size=n_adults), 
             position = position_jitter(width = .00001, height = .00001, seed = 1))+
  scale_size_continuous(name = "Recaptured", 
                        limits = c(0,max(mrr_data_C2_male$n_adults)),
                        breaks = seq(0,max(mrr_data_C2_male$n_adults),by=50))+
  scale_color_manual(name="", labels=c("Marked", "Unmarked"),values=c("#0000FF", "#A9A9A9"))+
  theme_minimal(base_size = 20)+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.position = "bottom",
        plot.title=element_text(margin=margin(b=5), size = 14, hjust = 0.15),
        plot.subtitle=element_text(margin=margin(t=5, b = 5)))+
  labs(title = 'Colletion: {closest_state}') +
  transition_states(mrr_data_C2_male$date,
                    transition_length = 0.5)+
  ease_aes("sine-in-out")

#animate(mrr_data_C2_plot, fps = 10, width = 750, height = 750)
#anim_save("./Field Activities/MRR/MRR/MRR - III/mrr_C2.gif")


### Graphs ###
MRR_data_C2_plot1<-ggplot(mrr_data_C2, aes(x=date, y=n_adults, col=marking_status))+
  geom_point(position = position_jitter(width=0.1,height=0.1,seed=1), alpha=0.5)+
  labs(title = "MRR-C2")+
  scale_color_manual(name="Marking status", values = c("#00b300", "#A9A9A9"))+
  scale_y_continuous(trans = "log10", name = "Recaptured males (log10)")+
  scale_x_date(name = "Recapture date", date_breaks = "2 day",date_labels = "%m/%d")+
  geom_smooth(method = "glm", method.args=list(family="poisson"), se = TRUE, alpha=0.25)+
  theme_bw(base_size = 20)+theme(legend.position="bottom")
MRR_data_C2_plot1
#ggsave("./Field Activities/MRR/MRR/MRR_C2_plot1.png",MRR_data_C2_plot1, dpi = 1000,
#      width = 28, height = 18, units = "cm")

###
MRR_data_C2_plot2<-ggplot(mrr_data_C2, aes(x=date, y=n_adults, col=marking_status))+
  geom_point(position = position_jitter(width=0.08,height=0.08,seed=1))+
  labs(x="Recapture date", y="Recaptured males (log10)")+ggtitle("MRR-C1 Multipoints")+
  scale_color_manual(name="Marking status", values = c("#ffa500", "#A9A9A9"))+
  scale_x_date(date_breaks = "2 day",date_labels = "%m/%d")+
  scale_y_continuous(trans = "log10")+
  geom_smooth(method = "glm", method.args=list(family="poisson"),
              se = TRUE, alpha=0.25)+
  theme_bw(base_size = 20)+theme(legend.position="bottom")+facet_grid(.~gender)
MRR_data_C2_plot2
#ggsave("./Field Activities/MRR/MRR/MRR - III/MRR_C2_plot2.pdf",MRR_data_C2_plot2, dpi = 1000,
#       width = 28, height = 18, units = "cm")

###
MRR_data_C2_plot3<-mrr_data_C2%>%#filter(marking_status=="marked")%>%
  ggplot()+
  geom_path(data = captiva_coord_1, aes(x=long, y = lat), col="black")+
  geom_point(data=release_multipoints, aes(x=long, y=lat),
             col="black",fill="black", size=0.85, shape=4)+
  geom_point(data=mrr_data_C2%>%filter(gender=="male"),
             aes(x=long, y=lat, size=n_adults, color=marking_status),
             position = position_jitter(width = .00001, height = .00001, seed = 1), alpha=0.5)+
  scale_color_manual(name="", 
                     values = c("#a9a9a9","#00b300"),
                     label=c("Unmarked", "Marked"))+
  scale_size_area(limits = c(0,max(mrr_data_C2%>%#filter(marking_status=="marked")%>%
                                     .$n_adults)),  max_size = 8,
                  breaks = seq(0,max(mrr_data_C2%>%#filter(marking_status=="marked")%>%
                                       .$n_adults),by=50),
                  name="")+
  theme_minimal(base_size = 20)+ggtitle("MRR-C2")+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        plot.title=element_text(margin=margin(b=5), size = 14, hjust = 0.15),
        plot.subtitle=element_text(margin=margin(t=5, b = 5)))+
  facet_wrap(.~date)
MRR_data_C2_plot3
#ggsave("./Field Activities/MRR/MRR/MRR-C2.1_plot3.png", MRR_data_C2_plot3, dpi=1000, width=36,
#       height=38, units="cm")

#### MRR-D #####
mrr_data_D<-mrr_data%>%filter(round=="D")

mrr_data_D%>%group_by(marking_status, date)%>%
  summarise(min=min(n_adults),
            max=max(n_adults),
            mean=mean(n_adults),
            sum=sum(n_adults),
            length=length(n_adults))

mrr_table_D<-mrr_data_D%>%filter(marking_status!="marked", gender=="male")%>%
  group_by(date)%>%
  summarise(mean_adult=mean(n_adults),
            sem=sd(n_adults)/sqrt(length(n_adults)),
            sum_total=sum(n_adults),
            reca_rate=(sum(n_adults)/MRR_D$released))

kable(mrr_table_D)

mrr_table_D1<-mrr_data_D%>%filter(gender=="male")%>%
  group_by(marking_status)%>%
  summarise(mean_recap=mean(n_adults),
            sem=sd(n_adults)/sqrt(length(n_adults)),
            sum_total=sum(n_adults),
            recap_rate=(sum(n_adults)/MRR_D$released))

kable(mrr_table_D1)

mrr_data_D$recap_rate<-mrr_data_D$n_adults/MRR_D$released

mrr_data_D_male<-mrr_data_D%>%
  filter(gender=="male")%>%
  group_by(date,trap_ID, long, lat, marking_status)%>%
  summarise(n_adults=sum(n_adults))

mrr_data_D_plot<-ggplot()+
  geom_path(data = captiva_coord_1, aes(x=long, y = lat), col="black")+
  geom_point(data=release_multipoints, aes(x=long, y=lat),
             col="black",fill="black", size=4, shape=25)+
  geom_point(data=mrr_data_D_male,
             aes(x=long, y=lat, col=marking_status, size=n_adults), 
             position = position_jitter(width = .00001, height = .00001, seed = 1))+
  scale_size_continuous(name = "Recaptured", 
                        limits = c(0,max(mrr_data_D_male$n_adults)),
                        breaks = seq(0,max(mrr_data_D_male$n_adults),by=50))+
  scale_color_manual(name="", labels=c("Marked", "Unmarked"),values=c("#0000FF", "#A9A9A9"))+
  theme_minimal(base_size = 20)+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.position = "bottom",
        plot.title=element_text(margin=margin(b=5), size = 14, hjust = 0.15),
        plot.subtitle=element_text(margin=margin(t=5, b = 5)))+
  labs(title = 'Colletion: {closest_state}') +
  transition_states(mrr_data_D_male$date,
                    transition_length = 0.5)+
  ease_aes("sine-in-out")

#animate(mrr_data_D_plot, fps = 10, width = 750, height = 750)
#anim_save("./Field Activities/MRR/MRR/MRR - III/mrr_D.gif")


### Graphs ###
MRR_data_D_plot1<-ggplot(mrr_data_D, aes(x=date, y=n_adults, col=marking_status))+
  geom_point(position = position_jitter(width=0.1,height=0.1,seed=1), alpha=0.5)+
  labs(title = "MRR-D")+
  scale_color_manual(name="Marking status", values = c("#ff00a5", "#A9A9A9"))+
  scale_y_continuous(trans = "log10", name = "Recaptured males (log10)")+
  scale_x_date(name = "Recapture date", date_breaks = "2 day",date_labels = "%m/%d")+
  geom_smooth(method = "glm", method.args=list(family="poisson"), se = TRUE, alpha=0.25)+
  theme_bw(base_size = 20)+theme(legend.position="bottom")
MRR_data_D_plot1
#ggsave("./Field Activities/MRR/MRR/MRR_D_plot1.png",MRR_data_D_plot1, dpi = 1000,
#      width = 28, height = 18, units = "cm")

###
MRR_data_D_plot2<-ggplot(mrr_data_D, aes(x=date, y=n_adults, col=marking_status))+
  geom_point(position = position_jitter(width=0.08,height=0.08,seed=1))+
  labs(x="Recapture date", y="Recaptured males (log10)")+ggtitle("MRR-C1 Multipoints")+
  scale_color_manual(name="Marking status", values = c("#ffa500", "#A9A9A9"))+
  scale_x_date(date_breaks = "2 day",date_labels = "%m/%d")+
  scale_y_continuous(trans = "log10")+
  geom_smooth(method = "glm", method.args=list(family="poisson"),
              se = TRUE, alpha=0.25)+
  theme_bw(base_size = 20)+theme(legend.position="bottom")+facet_grid(.~gender)
MRR_data_D_plot2
#ggsave("./Field Activities/MRR/MRR/MRR - III/MRR_D_plot2.pdf",MRR_data_D_plot2, dpi = 1000,
#       width = 28, height = 18, units = "cm")

###
MRR_data_D_plot3<-mrr_data_D%>%#filter(marking_status=="marked")%>%
  ggplot()+
  geom_path(data = captiva_coord_1, aes(x=long, y = lat), col="black")+
  geom_point(data=release_multipoints, aes(x=long, y=lat),
             col="black",fill="black", size=0.85, shape=4)+
  geom_point(data=mrr_data_D%>%filter(gender=="male"),
             aes(x=long, y=lat, size=n_adults, color=marking_status),
             position = position_jitter(width = .00001, height = .00001, seed = 1), alpha=0.5)+
  scale_color_manual(name="", 
                     values = c("#a9a9a9","#ff00a5"),
                     label=c("Unmarked", "Marked"))+
  scale_size_area(limits = c(0,max(mrr_data_D%>%#filter(marking_status=="marked")%>%
                                     .$n_adults)),  max_size = 8,
                  breaks = seq(0,max(mrr_data_D%>%#filter(marking_status=="marked")%>%
                                       .$n_adults),by=7),
                  name="")+
  theme_minimal(base_size = 20)+ggtitle("MRR-D")+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        plot.title=element_text(margin=margin(b=5), size = 14, hjust = 0.15),
        plot.subtitle=element_text(margin=margin(t=5, b = 5)))+
  facet_wrap(.~date)
MRR_data_D_plot3
#ggsave("./Field Activities/MRR/MRR/MRR-D.1_plot3.png", MRR_data_D_plot3, dpi=1000, width=26,
#       height=28, units="cm")

#### MRR-E #####
mrr_data_E<-mrr_data%>%filter(round=="E")

mrr_data_E%>%group_by(marking_status, date)%>%
  summarise(min=min(n_adults),
            max=max(n_adults),
            mean=mean(n_adults),
            sum=sum(n_adults),
            length=length(n_adults))

mrr_table_E<-mrr_data_E%>%filter(marking_status=="marked",gender=="male")%>%
  group_by(date)%>%
  summarise(mean_adult=mean(n_adults),
            sem=sd(n_adults)/sqrt(length(n_adults)),
            sum_total=sum(n_adults),
            recap_rate=(sum(n_adults)/MRR_E$released))

kable(mrr_table_E)

mrr_table_E1<-mrr_data_E%>%filter(gender=="male")%>%
  group_by(marking_status)%>%
  summarise(mean_recap=mean(n_adults),
            sem=sd(n_adults)/sqrt(length(n_adults)),
            sum_total=sum(n_adults),
            recap_rate=(sum(n_adults)/MRR_E$released))

kable(mrr_table_E1)

mrr_data_E$recap_rate<-mrr_data_E$n_adults/MRR_E$released

mrr_data_E_male<-mrr_data_E%>%
  filter(gender=="male")%>%
  group_by(date,trap_ID, long, lat, marking_status)%>%
  summarise(n_adults=sum(n_adults))

mrr_data_E_plot<-ggplot()+
  geom_path(data = captiva_coord_1, aes(x=long, y = lat), col="black")+
  geom_point(data=release_multipoints, aes(x=long, y=lat),
             col="black",fill="black", size=4, shape=25)+
  geom_point(data=mrr_data_E_male,
             aes(x=long, y=lat, col=marking_status, size=n_adults), 
             position = position_jitter(width = .00001, height = .00001, seed = 1))+
  scale_size_continuous(name = "Recaptured", 
                        limits = c(0,max(mrr_data_E_male$n_adults)),
                        breaks = seq(0,max(mrr_data_E_male$n_adults),by=50))+
  scale_color_manual(name="", labels=c("Marked", "Unmarked"),values=c("#0000FF", "#A9A9A9"))+
  theme_minimal(base_size = 20)+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.position = "bottom",
        plot.title=element_text(margin=margin(b=5), size = 14, hjust = 0.15),
        plot.subtitle=element_text(margin=margin(t=5, b = 5)))+
  labs(title = 'Colletion: {closest_state}') +
  transition_states(mrr_data_E_male$date,
                    transition_length = 0.5)+
  ease_aes("sine-in-out")

#animate(mrr_data_E_plot, fps = 10, width = 750, height = 750)
#anim_save("./Field Activities/MRR/MRR/MRR - III/mrr_F.gif")


### Graphs ###
MRR_data_E_plot1<-ggplot(mrr_data_E, aes(x=date, y=n_adults, col=marking_status))+
  geom_point(position = position_jitter(width=0.1,height=0.1,seed=1),alpha=0.5)+
  labs(title = "MRR-E")+
  scale_color_manual(name="Marking status", values = c("#00ffa5", "#A9A9A9"))+
  scale_y_continuous(trans = "log10", name = "Recaptured males (log10)")+
  scale_x_date(name = "Recapture date", date_breaks = "1 day",date_labels = "%m/%d")+
  geom_smooth(method = "glm", method.args=list(family="poisson"), se = TRUE, alpha=0.25)+
  theme_bw(base_size = 20)+theme(legend.position="bottom")
MRR_data_E_plot1
#ggsave("./Field Activities/MRR/MRR/MRR_E_plot1.png",MRR_data_E_plot1, dpi = 1000,
#      width = 28, height = 18, units = "cm")

###
MRR_data_E_plot2<-ggplot(mrr_data_E, aes(x=date, y=n_adults, col=marking_status))+
  geom_point(position = position_jitter(width=0.08,height=0.08,seed=1))+
  labs(x="Recapture date", y="Recaptured males (log10)")+ggtitle("MRR-E Multipoints")+
  scale_color_manual(name="Marking status", values = c("#00ffa5", "#A9A9A9"))+
  scale_x_date(date_breaks = "2 day",date_labels = "%m/%d")+
  scale_y_continuous(trans = "log10")+
  geom_smooth(method = "glm", method.args=list(family="poisson"),
              se = FALSE, alpha=0.25)+
  theme_bw(base_size = 20)+theme(legend.position="bottom")+facet_grid(.~gender)
MRR_data_E_plot2
#ggsave("./Field Activities/MRR/MRR/MRR - III/MRR_E_plot2.pdf",MRR_data_E_plot2, dpi = 1000,
#       width = 28, height = 18, units = "cm")

###
MRR_data_E_plot3<-mrr_data_E%>%#filter(marking_status=="marked")%>%
  ggplot()+
  geom_path(data = captiva_coord_1, aes(x=long, y = lat), col="black")+
  geom_point(data=release_multipoints, aes(x=long, y=lat),
             col="black",fill="black", size=0.85, shape=4)+
  geom_point(data=mrr_data_E%>%filter(gender=="male"),
             aes(x=long, y=lat, size=n_adults, color=marking_status),
             position = position_jitter(width = .00001, height = .00001, seed = 1), alpha=0.5)+
  scale_color_manual(name="", 
                     values = c("#a9a9a9","#00ffA5"),
                     label=c("Unmarked", "Marked"))+
  scale_size_area(limits = c(0,max(mrr_data_E%>%#filter(marking_status=="marked")%>%
                                     .$n_adults)),  max_size = 8,
                  breaks = seq(0,max(mrr_data_E%>%#filter(marking_status=="marked")%>%
                                       .$n_adults),by=5),
                  name="")+
  theme_minimal(base_size = 20)+ggtitle("MRR-E")+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        plot.title=element_text(margin=margin(b=5), size = 14, hjust = 0.15),
        plot.subtitle=element_text(margin=margin(t=5, b = 5)))+
  facet_wrap(.~date)
MRR_data_E_plot3
#ggsave("./Field Activities/MRR/MRR/MRR-E.1_plot3.png", MRR_data_E_plot3, dpi=1000, width=20,
#       height=18, units="cm")

#### MRR-F #####
mrr_data_F<-mrr_data%>%filter(round=="F")

mrr_data_F%>%group_by(marking_status, date)%>%
  summarise(min=min(n_adults),
            max=max(n_adults),
            mean=mean(n_adults),
            sum=sum(n_adults),
            length=length(n_adults))

mrr_table_F<-mrr_data_F%>%filter(marking_status=="marked", gender=="male")%>%
  #group_by(date)%>%
  summarise(mean_adult=mean(n_adults),
            sem=sd(n_adults)/sqrt(length(n_adults)),
            sum_total=sum(n_adults),
            reca_rate=(sum(n_adults)/MRR_F$released))

kable(mrr_table_F)

mrr_table_F1<-mrr_data_F%>%filter(gender=="male")%>%
  group_by(marking_status)%>%
  summarise(mean_recap=mean(n_adults),
            sem=sd(n_adults)/sqrt(length(n_adults)),
            sum_total=sum(n_adults),
            recap_rate=(sum(n_adults)/MRR_F$released))

kable(mrr_table_F1)

mrr_data_F$recap_rate<-mrr_data_F$n_adults/MRR_F$released

mrr_data_F_male<-mrr_data_F%>%
  filter(gender=="male")%>%
  group_by(date,trap_ID, long, lat, marking_status)%>%
  summarise(n_adults=sum(n_adults))

mrr_data_F_plot<-ggplot()+
  geom_path(data = captiva_coord_1, aes(x=long, y = lat), col="black")+
  geom_point(data=release_multipoints, aes(x=long, y=lat),
             col="black",fill="black", size=4, shape=25)+
  geom_point(data=mrr_data_F_male,
             aes(x=long, y=lat, col=marking_status, size=n_adults), 
             position = position_jitter(width = .00001, height = .00001, seed = 1))+
  scale_size_continuous(name = "Recaptured", 
                        limits = c(0,max(mrr_data_F_male$n_adults)),
                        breaks = seq(0,max(mrr_data_F_male$n_adults),by=50))+
  scale_color_manual(name="", labels=c("Marked", "Unmarked"),values=c("#0000FF", "#A9A9A9"))+
  theme_minimal(base_size = 20)+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.position = "bottom",
        plot.title=element_text(margin=margin(b=5), size = 14, hjust = 0.15),
        plot.subtitle=element_text(margin=margin(t=5, b = 5)))+
  labs(title = 'Colletion: {closest_state}') +
  transition_states(mrr_data_F_male$date,
                    transition_length = 0.5)+
  ease_aes("sine-in-out")

#animate(mrr_data_F_plot, fps = 10, width = 750, height = 750)
#anim_save("./Field Activities/MRR/MRR/MRR - III/mrr_F.gif")

### Graphs ###
MRR_data_F_plot1<-ggplot(mrr_data_F, aes(x=date, y=n_adults, col=marking_status))+
  geom_point(position = position_jitter(width=0.1,height=0.1,seed=1), alpha=0.5)+
  labs(title = "MRR-F")+
  scale_color_manual(name="Marking status", values = c("#00a5ff", "#A9A9A9"))+
  scale_y_continuous(trans = "log10", name = "Recaptured males (log10)")+
  scale_x_date(name = "Recapture date", date_breaks = "2 day",date_labels = "%m/%d")+
  geom_smooth(method = "glm", method.args=list(family="poisson"), se = TRUE, alpha=0.25)+
  theme_bw(base_size = 20)+theme(legend.position="bottom")
MRR_data_F_plot1
#ggsave("./Field Activities/MRR/MRR/MRR_F_plot1.png",MRR_data_F_plot1, dpi = 1000,
#      width = 28, height = 18, units = "cm")

###
MRR_data_F_plot2<-ggplot(mrr_data_F, aes(x=date, y=n_adults, col=marking_status))+
  geom_point(position = position_jitter(width=0.08,height=0.08,seed=1))+
  labs(x="Recapture date", y="Recaptured males (log10)")+ggtitle("MRR-C1 Multipoints")+
  scale_color_manual(name="Marking status", values = c("#00a5ff", "#A9A9A9"))+
  scale_x_date(date_breaks = "3 day",date_labels = "%m/%d")+
  scale_y_continuous(trans = "log10")+
  geom_smooth(method = "glm", method.args=list(family="poisson"),
              se = FALSE, alpha=0.25)+
  theme_bw(base_size = 20)+theme(legend.position="bottom")+facet_grid(.~gender)
MRR_data_F_plot2
#ggsave("./Field Activities/MRR/MRR/MRR - III/MRR_F_plot2.pdf",MRR_data_F_plot2, dpi = 1000,
#       width = 28, height = 18, units = "cm")

###
MRR_data_F_plot3<-mrr_data_F%>%#filter(marking_status=="marked")%>%
  ggplot()+
  geom_path(data = captiva_coord_1, aes(x=long, y = lat), col="black")+
  geom_point(data=release_multipoints, aes(x=long, y=lat),
             col="black",fill="black", size=0.85, shape=4)+
  geom_point(data=mrr_data_F%>%filter(gender=="male"),
             aes(x=long, y=lat, size=n_adults, color=marking_status),
             position = position_jitter(width = .00001, height = .00001, seed = 1), alpha=0.5)+
  scale_color_manual(name="", 
                     values = c("#a9a9a9","#00a5ff"),
                     label=c("Unmarked", "Marked"))+
  scale_size_area(limits = c(0,max(mrr_data_F%>%#filter(marking_status=="marked")%>%
                                     .$n_adults)),  max_size = 8,
                  breaks = seq(0,max(mrr_data_F%>%#filter(marking_status=="marked")%>%
                                       .$n_adults),by=3),
                  name="")+
  theme_minimal(base_size = 20)+ggtitle("MRR-F")+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        plot.title=element_text(margin=margin(b=5), size = 14, hjust = 0.15),
        plot.subtitle=element_text(margin=margin(t=5, b = 5)))+
  facet_wrap(.~date)
MRR_data_F_plot3
#ggsave("./Field Activities/MRR/MRR/MRR-F.1_plot3.png", MRR_data_F_plot3, dpi=1000, width=28,
#       height=26, units="cm")