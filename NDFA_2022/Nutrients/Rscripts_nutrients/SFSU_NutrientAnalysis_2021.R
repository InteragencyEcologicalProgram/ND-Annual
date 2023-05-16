rm(list = ls())
#setwd("C:/Users/mminer/Documents/R_Data_Analysis/NDFS/Nutrients_2021")

#libraries 
library(tidyverse) 
library(ggplot2)
library(lubridate)
library(data.table)
library(rstatix)
library(gridExtra)
library(janitor)
library(ggpubr)
library(viridis)
library(lme4)
library(psych)
library(nlme)
library(lattice)
library(lmerTest)
library(visreg)
library(emmeans)
library(multcomp)
library(multcompView)
library(viridis)
library(car)
library(knitr)

#load data
sfsu_nut <- read.csv("C:/Users/jtorre/Desktop/Github_Repos/ND-Annual/NDFA_2022/Nutrients/Data_nutrients/SFSU_nutrients_2021.csv")

#select relevant columns, clean df 
sfsu_nut <- dplyr::select(sfsu_nut, "Date","Time","Station","silica","nitrate_nitrite","nitrite","orthophosphate","ammonia","ammonia_conversion_mg_L")

# format date
sfsu_nut$Date = as.POSIXct(sfsu_nut$Date, format = "%m/%d/%Y", tz = "UTC")


#assign region
downstream<-c("RVB","RVB_T0", "RYI",  "LIB", "PRS", "BL5") 
upstream<-c("I80TD", "I80-TD", "RD22", "RMB", "LIS", "STTD", "WWT", "DWT")
sfsu_nut$region <- ifelse(sfsu_nut$Station %in% downstream, "Downstream",
                          ifelse(sfsu_nut$Station %in% upstream, "Upstream", "Middle Sac River"))

sfsu_nut$region <- factor(sfsu_nut$region, levels = c("Upstream","Downstream","Middle Sac River")) # Done like this to organize them in plots? -Juan

sfsu_nut$Station = as.factor(sfsu_nut$Station)

#assign period 
nut <- sfsu_nut %>% 
  mutate(period = case_when(
    Date <= ymd_hms("2021-09-11 00:00:00") ~ "Before",
    Date >= ymd_hms ("2021-09-11 01:00:00") & Date <= ymd_hms("2021-09-15 00:00:00") ~ "During",
    Date >= ymd_hms("2021-09-15 01:00:00") ~ "After"))

nut$period=as.factor(nut$period)

nut %>% count(Station)

nut <- nut %>% reorder_levels(period, order = c("Before", "During", "After"))
nut <- droplevels(nut[!nut$region=="Middle Sac River", ])
nut <- nut %>% mutate_at(c("silica","nitrate_nitrite","nitrite","orthophosphate","ammonia","ammonia_conversion_mg_L"),as.numeric)

# pivot 
nut2 <- nut %>% pivot_longer(cols = c("silica","nitrate_nitrite","nitrite","orthophosphate","ammonia","ammonia_conversion_mg_L"), names_to = "analyte", values_to = "result")


# Subset by analyte for stats tests --> all of these are in umol NOT mg/L
amm <- nut2[nut2$analyte == "ammonia", ]
amm_statsum <- amm %>% 
  group_by(period, region) %>% 
  get_summary_stats(result, type = "mean_sd") 

silica <- nut2[nut2$analyte == "silica", ]
silica_statsum <- silica %>% 
  group_by(period, region) %>% 
  get_summary_stats(result, type = "mean_sd")

nit <- nut2[nut2$analyte == "nitrate_nitrite", ]
  nit_statsum <- nit %>% 
     group_by(period, region) %>% get_summary_stats(result, type = "mean_sd")
  
phos <- nut2[nut2$analyte == "orthophosphate", ]
  phos_statsum <- phos %>% 
    group_by(period, region) %>% get_summary_stats(result, type = "mean_sd")

# summarize everything
nut_statsum <- nut %>% 
  filter(region!="Middle Sac River") %>% 
  group_by(period,region) %>% 
  get_summary_stats(type = "mean_sd")
  # Will use this method for the output into word docs


#format for plotting: 
#not plotting ammonia mg/L or nitrite
nut_plot <- nut2[nut2$analyte!="ammonia_conversion_mg_L", ]
nut_plot <- nut_plot[nut_plot$analyte!="nitrite", ]

# Make pretty lables 
nut.labs <- c("Ammonia","Nitrate/ Nitrite","ortho-\nPhosphate","Silica")
names(nut.labs) <- c("ammonia", "nitrate_nitrite", "orthophosphate","silica")

#plot 
umol <- ggplot(nut_plot, aes(x = period, y=result, fill = period))+
  geom_boxplot()+
  stat_boxplot(geom = "errorbar")+
  labs(y="Analyte Concentration (μmol)", x="Flow Pulse Period", fill = "Flow Period")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  guides(fill= "none")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15), 
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), 
        strip.text=element_text(size=12))+
  # ylim(0,30)+
  scale_fill_viridis(option = "A", discrete = TRUE) +
  facet_grid(analyte~region, scales = "free_y",labeller = labeller(analyte = nut.labs))
# ggsave(filename = "output/SFSU_nutrients_2021_umol.png", height = 8.5, width = 6, units = "in",dpi = 600)

# Nested anovas lmer(Result ~ region*period + (1|station), nutrients, REML = TRUE)
# Overall, nutrient values were greater and demonstrated increased variability in the upstream region for all measured analytes across pulse flow periods. 
# Ammonia: 
naov_amm <- lmer(result~region*period + (1|Station), amm, REML = T)
anova(naov_amm) #neither region OR period is significant?!
hist(residuals(naov_amm)) #normal, no transformation 
  emmip(naov_amm, region ~ period) 
  emmip(naov_amm, period ~ region) 
# ammonia concentrations were notably lower downstream throughout all pulse periods,  but did not differ statistically between regions or periods. in the downstream region, ammonia was highest after the pusle period. In the upstream region, ammonia was highest before the pulse period. in the downstream region, ammonia was highest after the pusle period.
  
  emmeans(naov_amm, pairwise ~ period | region) #post hoc not necessary becuase nothing in anova was sig
  with(amm, interaction.plot(period, region, result)) 

  
# Silica:   
naov_si <- lmer(result~region*period + (1|Station), silica, REML = T)
anova(naov_si) #region significant 
hist(residuals(naov_si)) # normalish 
df.residual(naov_si)
  emmip(naov_si, region ~ period) 
  emmip(naov_si, period ~ region)
  emmeans(naov_si, pairwise ~ region | period) #CI overlap, but still getting significant p value for differences between us and ds
  
# Silica concentrations were signficnatly higher upstream. Silica concentrations dipped slightly in both regions following the pulse period. 
  
  with(silica, interaction.plot(period, region, result)) 

naov_nit <- lmer(result~region*period + (1|Station), nit, REML = T)
anova(naov_nit) # period but not region is significant, also a significant interaction term 
df.residual(naov_nit)
hist(residuals(naov_nit))
  emmip(naov_nit, region ~ period)# nit mean lower in ds region during all periods, increases throughout all periods in upstream region? 
  emmip(naov_nit, period ~ region)
  emmeans(naov_nit, pairwise ~ period | region)
#nitrate/nitrite differed between pulse flow periods, increasing upstream between each period (concentrations were greater during and after the pulse flow). concentratiosn remained consistently low downstream. Significant interaction term seems to indicate the that concentrations were different between regions, but also t hat in the upstream region in particular pulse flow period had an effect on concentrations. 
  
naov_phos <- lmer(result~region*period + (1|Station), phos, REML = T)  
anova(naov_phos) # region sig 
df.residual(naov_phos)
hist(residuals(naov_phos))
  emmip(naov_phos, region ~ period) #almost parallel
  emmip(naov_phos, period ~ region)
  emmeans(naov_phos,pairwise ~period |region) # individual contrasts have sig overlapping CIs. 
  
  #ortho-phosphate concentrations were signficnatnly greater in the upstream region 
