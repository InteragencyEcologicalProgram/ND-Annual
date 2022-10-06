#Summer-Fall Habitat Action Report
#Data analysis of zooplankton, phytoplankton, and water quality from 2020 monitoring
#Laura Twardochleb
#Updated 10/6/22

rm(list = ls())

#packages
lapply(c("CDECRetrieve","tidyverse", "lubridate", "viridis", "ggpubr"), require, character.only = TRUE)

#set directory
setwd("~/NDFA/2022/SFHA_Report")

#"during" flow pulse period is Sept. 11-14

#define sampling regions
downstream<-c("RVB","RYI",  "LIB", "PRS", "BL5") 
upstream<-c("I80", "RD22", "RMB","RCS", "LIS", "STTD")

##################### Discrete WQ ###############################################################################################
#DO, sp cond., pH, turbidity, secchi, temperature

physical_wq<-read_csv("NDFS_WQ_Data_2021.csv", na = c("", "N/A", "NA")) #need to write pre-cleaning steps into this script

#data cleaning: filter to NDFA, 2020 samples, change to long format for wq measurements
physical_wq2<-physical_wq%>%filter(Measuring.Program.Name%in%c("NDFA", "NDFS", "Shared"))%>% #filter for NDFS data
  filter(`WDL_SAM_COLLECTION_DATE`>= "8/2/2021"| `WDL_SAM_COLLECTION_DATE`< "10/14/2021")%>% #efilter to NDFS monitoring season
  pivot_longer(c("secchi","water.temp","DO.probe","sp.cond","EC","pH", "turb"), names_to = "Analyte", values_to ="Result") #transform data to long format

#assign regions
physical_wq2$Region<-ifelse(physical_wq2$Station.Name %in% downstream, "Downstream",
                    ifelse(physical_wq2$Station.Name %in% upstream, "Upstream", "Middle Sac River"))
physical_wq2$Region<-factor(physical_wq2$Region, levels=c("Upstream", "Downstream", "Middle Sac River"))

#assign flow pulse periods- for automated reporting, we could define pulse periods directly from the CDEC flow data
before<-c("8/2/2021", "8/3/2021", "8/30/2021", "8/31/2021")
during<-c("9/13/2021", "9/14/2021", "9/27/2021", "9/28/2021")

physical_wq2$Pulse_period<-ifelse(physical_wq2$WDL_SAM_COLLECTION_DATE %in% before, "Before",
                           ifelse(physical_wq2$WDL_SAM_COLLECTION_DATE %in% during, "During", "After"))
physical_wq2$Pulse_period<-factor(physical_wq2$Pulse_period, levels=c("Before", "During", "After"))

#summarize values by region and flow pulse period
phys_wq_means<-physical_wq2%>%group_by(Analyte, Region, Pulse_period)%>%summarize(wq_mean=mean(Result, na.rm=TRUE), wq_sd=sd(Result, na.rm=TRUE))

#rename analytes
phys_wq_means$Analyte[phys_wq_means$Analyte=="DO.probe"]<-"DO"
phys_wq_means$Analyte[phys_wq_means$Analyte=="secchi"]<-"Secchi"
phys_wq_means$Analyte[phys_wq_means$Analyte=="sp.cond"]<-"Cond"
phys_wq_means$Analyte[phys_wq_means$Analyte=="turb"]<-"Turb"
phys_wq_means$Analyte[phys_wq_means$Analyte=="water.temp"]<-"Temp"

#plot water quality by region- change to free y-scales
wq_plot<-phys_wq_means%>%filter(Region!="Middle Sac River")%>%filter(Analyte!="EC")%>%ggplot(aes(x=Pulse_period, y=wq_mean, fill=Region))+geom_bar(stat="identity")+geom_errorbar(aes(ymin=0, ymax=wq_mean+wq_sd))+facet_grid(cols=vars(Region), rows=vars(Analyte), scales = "free_y")+
  scale_fill_viridis(discrete = TRUE, option="A")+
  xlab("Region")+ylab("Mean value")+
  theme_bw()+theme(legend.position = "none")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18), axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), strip.text=element_text(size=15))
ggsave(filename="wq.png", wq_plot, height=6, width=10, dpi=600)

#run t-tests on physical water quality data
#group data by analyte and calculate site means
wq_site_means<-physical_wq2%>%filter(Region!="Middle Sac River")%>%group_by(Region, Station.Name, Analyte)%>%summarize(mean=mean(Result, na.rm=TRUE), sd=sd(Result, na.rm=TRUE))

#subset and test by analyte
DO<-filter(wq_site_means, Analyte=="DO.probe")
hist(log(DO$mean)) #no log-transformation
DO.t.test<-t.test(mean~Region, DO, alternative="two.sided", var.equal=FALSE)
DO.t.test #downstream higher

pH<-filter(wq_site_means, Analyte=="pH")
hist(log(pH$mean))
hist(pH$mean)#no log-transformation
pH.t.test<-t.test(mean~Region, pH, alternative="two.sided", var.equal=FALSE)
pH.t.test #not significantly different

secchi<-filter(wq_site_means, Analyte=="secchi")
hist(log(secchi$mean))
hist(secchi$mean)#no log-transformation
secchi.t.test<-t.test(mean~Region, secchi, alternative="two.sided", var.equal=FALSE)
secchi.t.test #downstream higher

sp.cond<-filter(wq_site_means, Analyte=="sp.cond")
hist(log(sp.cond$mean))
hist(sp.cond$mean)#no log-transformation
sp.cond.t.test<-t.test(mean~Region, sp.cond, alternative="two.sided", var.equal=FALSE)
sp.cond.t.test #upstream higher

turb.FNU<-filter(wq_site_means, Analyte=="turb.FNU")
hist(log(turb.FNU$mean))
hist(turb.FNU$mean)#no log-transformation
turb.FNU.t.test<-t.test(mean~Region, turb.FNU, alternative="two.sided", var.equal=FALSE)
turb.FNU.t.test #upstream higher

water.temp<-filter(wq_site_means, Analyte=="water.temp")
hist(log(water.temp$mean))
hist(water.temp$mean)#no log-transformation
water.temp.t.test<-t.test(mean~Region, water.temp, alternative="two.sided", var.equal=FALSE)
water.temp.t.test #upstream higher

################### save script workspace image #######################################################################
save.image("SFHA_report.RData")
