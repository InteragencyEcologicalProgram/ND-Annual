rm(list = ls())

lapply(c("tidyverse", "lubridate", "viridis", "ggpubr", "remotes", "scales", "data.table","rstatix","gridExtra","janitor","lme4","psych","nlme","lattice","lmerTest","visreg","emmeans","multcomp","multcompView","readxl"), require, character.only = TRUE)

setwd("C:/Users/mminer/Documents/R_Data_Analysis/NDFS/WaterQuality_2021")

# "during" flow pulse period is Sept. 11-14

# read in excel file - kind of funky, but need to tell R the type of data in each column and skip first 
  # lets convert from excel to CSV before loading into R - juan
physical_wq <- read_excel("data/NDFS_WQ_Data_2021.xlsx", 
                           col_types = c("numeric", "numeric", "text", "text", "numeric", "text", "date", 
                                          "date", "numeric", "text", "text","text", "text", "text", "text", "text", 
                                          "text", "text", "numeric", "numeric", "text", "text", "text", "text", "text", 
                                          "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", 
                                          "text", "text", "text", "text", "text","text", "text", "text", "text", "text", 
                                          "text", "text", "text", "text"), skip = 1)

# remove uncessary column - artifact of excel format
physical_wq <- physical_wq[-c(1),]

# rename fore ease and R formatting compliance 
physical_wq <- physical_wq %>% 
  rename(program=3, station=`Station Name`, samcode =`Sampling Number (event)`, date = `WDL SAM_COLLECTION_DATE`,do = DO.probe) 

# format quantiative values as numeric 
physical_wq1 <- physical_wq %>% mutate_at(c("secchi", "water.temp", "do", "sp.cond", "EC", "pH", "turb"), as.numeric)

# reformat structure
physical_wq2<-physical_wq1 %>%pivot_longer(cols = c("secchi","water.temp","do","sp.cond","EC","pH", "turb"), names_to = "Analyte", values_to ="Result")

# clean & select relevant data for dataframe
wq <- dplyr::select(physical_wq2, "program", 'station', 'samcode', 'date', `Collection Time`, 'Analyte', 'Result')

# assign regions 
downstream<-c("RVB","RYI",  "LIB", "PRS", "BL5") 
upstream<-c("I80", "RD22", "RMB","RCS", "LIS", "STTD", "WWT", "DWT") 


wq$region<-ifelse(wq$station %in% downstream, "Downstream",
                            ifelse(wq$station %in% upstream, "Upstream", "Middle Sac River"))

# order regions 
wq$region<-factor(wq$region, levels=c("Upstream", "Downstream", "Middle Sac River"))


# assign flow pulse periods - pulse flow occured between 9/11 and 9/14 2021

before<-c("2021-08-02", "2021-08-03", "2021-08-16", "2021-08-17", "2021-08-30", "2021-08-31") 
during<-c("2021-09-13", "2021-09-14") 

wq <- wq %>% 
  mutate(pulse_period = case_when(
    date <= ymd_hms("2021-09-11 00:00:00") ~ "Before",
    date >= ymd_hms ("2021-09-11 01:00:00") &
      date <= ymd_hms("2021-09-15 00:00:00") ~ "During",
    date >= ymd_hms("2021-09-15 01:00:00") ~ "After"))

wq$pulse_period <- factor(wq$pulse_period, levels = c("Before", "During", "After"))

# rename analytes 
wq$Analyte[wq$Analyte=="do"]<-"DO"
wq$Analyte[wq$Analyte=="secchi"]<-"Secchi"
wq$Analyte[wq$Analyte=="sp.cond"]<-"Cond"
wq$Analyte[wq$Analyte=="turb"]<-"Turb"
wq$Analyte[wq$Analyte=="water.temp"]<-"Temp"

wq$station = as.factor(wq$station)

# summarize values by region and flow pulse period
wq_means<-wq %>%
  group_by(Analyte, region) %>%
  summarize(wq_mean=mean(Result, na.rm=TRUE), wq_sd=sd(Result, na.rm=TRUE))

# Plot WQ - removing EC and middle sac factor levels for plotting purposes only 
wq_plot <- droplevels(wq[!wq$region == "Middle Sac River", ])
wq_plot <- droplevels(wq_plot[!wq_plot$Analyte == "EC", ])

#remove inaccurate/suspicious data 
wq_plot <- wq_plot %>% dplyr::filter(!(Analyte == "Cond" & Result == 10000.00))
wq_plot <- wq_plot %>% dplyr::filter(!(Analyte == "Turb" & Result == 240.50))

# create accurate/pretty labels for wq figure below
wq.labs <- c("Conductivity \n(uS/cm)", "Dissolved \nOxygen (mg/L)","pH", "Secchi \nDepth (m)","Temperature \n(°C)","Turbidity \n(FNU)")
names(wq.labs) <- c("Cond","DO","pH","Secchi","Temp","Turb")

wq_fig <- ggplot(wq_plot, aes(x=pulse_period, y=Result, fill=pulse_period)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar")+
  scale_fill_viridis(discrete = TRUE, option="A") +
  labs(x = "Flow Pulse Period", y = "Mean Value", fill = "Flow Pulse Period") +
  theme_bw() + 
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15), 
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), strip.text=element_text(size=12))+
    facet_grid(Analyte~region, scales = "free_y", labeller = labeller(Analyte = wq.labs))
wq_fig
# ggsave(filename = "physical_wq_2021_4.png", height = 8.5, width = 6, units = "in",dpi = 600)


######### Summary stats ###############
wq_site_means<-wq %>%
  filter(region!="Middle Sac River") %>% 
  group_by(station, Analyte, region) %>%
  summarize(mean=mean(Result, na.rm=TRUE), sd=sd(Result, na.rm=TRUE))

wq_region_means<-wq %>%
  filter(region!="Middle Sac River") %>% 
  group_by(region, Analyte) %>%
  summarize(mean=mean(Result, na.rm=TRUE), sd=sd(Result, na.rm=TRUE))

wq_period_means<-wq %>%
  filter(region!="Middle Sac River") %>% 
  group_by(region, Analyte, pulse_period) %>%
  summarize(mean=mean(Result, na.rm=TRUE), sd=sd(Result, na.rm=TRUE))
# write_xlsx(wq_period_means, "output/physical_wq_means_by_period.xlsx")

##### T-TESTS #####
# T-Tests were not  used for reporting results in 2021/2022 - instead mixed anovas were used, but leaving in this code for future reference 
# subset and test by analyte
DO <- filter(wq_site_means, Analyte=="DO")
hist(log(DO$mean)) #no log-transformation
DO.t.test<-t.test(mean~region, DO, alternative="two.sided", var.equal=FALSE)
DO.t.test #downstream significantly higher

pH<-filter(wq_site_means, Analyte=="pH")
hist(log(pH$mean))
hist(pH$mean)#no log-transformation
pH.t.test<-t.test(mean~region, pH, alternative="two.sided", var.equal=FALSE)
pH.t.test #not significantly different

secchi<-filter(wq_site_means, Analyte=="Secchi")
hist(log(secchi$mean))
hist(secchi$mean)#no log-transformation
secchi.t.test<-t.test(mean~region, secchi, alternative="two.sided", var.equal=FALSE)
secchi.t.test #downstream significantly higher

sp.cond<-filter(wq_site_means, Analyte=="Cond")
hist(log(sp.cond$mean))
hist(sp.cond$mean)#no log-transformation
sp.cond.t.test<-t.test(mean~region, sp.cond, alternative="two.sided", var.equal=FALSE)
sp.cond.t.test #upstream significantly higher

turb <-filter(wq_site_means, Analyte=="Turb") #laura had this seperated into FNU. but I'm not seeing that variable in teh script. reducing to turb. 
hist(log(turb$mean))
hist(turb$mean)#no log-transformation
turb.t.test<-t.test(mean~region, turb, alternative="two.sided", var.equal=FALSE)
turb.t.test #upstream higher, but not significantly... almost, but not quite

water.temp<-filter(wq_site_means, Analyte=="Temp")
hist(log(water.temp$mean))
hist(water.temp$mean)#no log-transformation
water.temp.t.test<-t.test(mean~region, water.temp, alternative="two.sided", var.equal=FALSE)
water.temp.t.test #upstream significantly higher


# T-TEST SUMMARY: 
# Dissolved oxygen and water clarity (secchi) were significantly higher downstream while conductivity and water temperature were signficantly greater in the upstream region. Other parameters, including pH and turbidity did not differ significantly between regions. 


##### ANOVAS ##### Recode for two way interactions 11/2/22

##### Two-Way Mixed Effects ANOVAS with random effect of staion ############################
# FOCUS ON RESULTS IN TWO WAY ANOVA - POST HOC TESTS ARE NOT INTUITIVE - IF GOING BY POST HOC TEST USE COMPARISON OF CL

########## Dissolved Oxygen: 
raw_DO <- filter(wq, Analyte=="DO")
raw_DO <- droplevels(raw_DO[!raw_DO$region == "Middle Sac River", ])
naov_do <- lmer(Result ~ region*pulse_period +(1|station), raw_DO, REML = T)
  anova(naov_do) #region and period are significant 
  
# Post hoc tests   
  do_posthoc <- lsmeans(naov_do, ~region+pulse_period, adjust = "sidak") #post hoc test on region and period (signifcant from aov)
  contrast(do_posthoc) # contrast() produces t-values that are the square root of the F statistic. 
  confint(do_posthoc, adjust = "sidak")
  
# find model residuals  
  df.residual(naov_do)
  
do_sum <- raw_DO %>% 
  dplyr::group_by(pulse_period, region) %>% 
  summarise(mean = mean(Result, na.rm = T),
            sd = sd(Result, na.rm = T))

############# pH: 
raw_pH <- filter(wq, Analyte=="pH")
raw_pH <- droplevels(raw_pH[!raw_pH$region == "Middle Sac River", ])

naov_pH <- lmer(Result ~ region*pulse_period + (1|station), raw_pH, REML = T)
  anova(naov_pH) #region is significant, period is NOT
  
# post hoc test
  ph_posthoc <- lsmeans(naov_pH, ~region, adjust = "sidak") #since pulse period was not significant, should it be included here
  contrast(ph_posthoc)
  
pH_sum <- raw_pH %>% 
  dplyr::group_by(pulse_period, region) %>% 
  summarise(mean = mean(Result, na.rm = T),
            sd = sd(Result, na.rm = T)) #pH appeared to remain consistent between flow pulse periods across regions. 

############## Secchi: 
raw_secchi <- filter(wq, Analyte=="Secchi")
raw_secchi <- droplevels(raw_secchi[!raw_secchi$region == "Middle Sac River", ])

naov_secchi <- lmer(Result~region*pulse_period + (1|station), raw_secchi, REML = T)
  anova(naov_secchi) #region and pulse period are independently significant, but the interaction is also sig.

# post hoc test  
  ph_secchi <- lsmeans(naov_secchi, ~region*pulse_period, adjust="sidak")
  contrast(ph_secchi) #every comparison is significant 
  
secchi_sum <- raw_secchi %>% 
  dplyr::group_by(pulse_period, region) %>% 
  summarise(mean = mean(Result, na.rm = T),
            sd = sd(Result, na.rm = T))


##### Temp:
raw_temp <- filter(wq, Analyte=="Temp")
raw_temp <- droplevels(raw_temp[!raw_temp$region == "Middle Sac River", ])

naov_temp <- lmer(Result~region*pulse_period +(1|station), raw_temp, REML = T)
  anova(naov_temp) #region and period independently significant 
  
#post hoc test  
  ph_temp <- lsmeans(naov_temp, ~region*pulse_period, adjust="sidak")
  contrast(ph_temp)
  
temp_sum <- raw_temp %>% 
  dplyr::group_by(pulse_period, region) %>% 
  summarise(mean = mean(Result, na.rm = T),
              sd = sd(Result, na.rm = T))


########### Conductivity: 
raw_cond <- filter(wq, Analyte=="Cond")
raw_cond <- droplevels(raw_cond[!raw_cond$region == "Middle Sac River", ])


naov_cond <- lmer(Result~region*pulse_period + (1|station), raw_cond, REML = T)
  anova(naov_cond) #only region is sig
  ph_cond <- lsmeans(naov_cond, ~region, adjust="sidak")
  ph_cond
  contrast(ph_cond) #upstream after was significantly different than everything else? 
  emmip(naov_cond, region ~ pulse_period) 
  emmip(naov_cond, pulse_period ~ region) 
  
cond_sum <- raw_cond %>% 
  dplyr::group_by(pulse_period, region) %>% 
  summarise(mean = mean(Result, na.rm = T),
              sd = sd(Result, na.rm = T)) 
  
########### Turbidity: 
raw_turb <- filter(wq, Analyte=="Turb")
raw_turb <- droplevels(raw_turb[!raw_turb$region == "Middle Sac River", ])


naov_turb <- lmer(Result~region*pulse_period + (1|station), raw_turb, REML = T)
  anova(naov_turb) #nothing is sig? that doesnt seem right
  ph_turb <- lsmeans(naov_turb, ~region*pulse_period, adjust="sidak")
  ph_turb  
  contrast(ph_turb)
