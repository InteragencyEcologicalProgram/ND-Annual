# resources: https://www.datanovia.com/en/lessons/anova-in-r/
  # https://www.datanovia.com/en/lessons/kruskal-wallis-test-in-r/

### daily notes:
# EOD 9/23 still couldn't transform data to meet normality assumptions of ANOVA - Performed Kruskal-Wallis test and posthoc tests instead.

# Laura noted that non parametric tests are not as robust, suggested retyring ANOVA using mixed effect model. Fixed effect = region, random effect = station. 

# EOD 10/12 still struggling to make nested/mixed effect anova work... beginning to think it all needs to be done in separate dataframes by analyte type.

rm(list = ls())
setwd("C:/Users/mminer/Documents/R_Data_Analysis/NDFS/Nutrients_2021")

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

# data 
nutrients <- read_csv("C:/Users/mminer/Documents/R_Data_Analysis/NDFS/Nutrients_2021/nutrients_2021.csv")
str(nutrients)

# Reformat & Refactor
nutrients <- nutrients %>% mutate_at(c("station","Analyte","RelativetoRptLimit","period","region"), as.factor)

nutrients$Analyte <- recode_factor(nutrients$Analyte, "Dissolved Ammonia" = "Ammonia",
                                   'Dissolved Calcium' = "Calcium", 
                                   'Chlorophyll a' = "Chlorophyll", 
                                   "Dissolved Nitrate + Nitrite"= "Nitrate/Nitrite",
                                   "Dissolved Organic Carbon" = "DOC", 
                                   "Dissolved Organic Nitrogen" = "DON",
                                   "Dissolved ortho-Phosphate" = "orthoPhosphate", 
                                   "Dissolved Silica (SiO2)" = "Silica")

nutrients <- nutrients %>% reorder_levels(period, order = c("Before", "During", "After"))

# Data Transformations ------------------------------ 
nutrients$log_result <- log(nutrients$Result)

nutrients$log10_result <- log10(nutrients$Result)

# Subset by analyte and check for normality: individual analytes really left skewed
# Ammonia
nutrients_ammonia <- nutrients[nutrients$Analyte == "Ammonia", ]
#untransformed 
amm <- ggdensity(nutrients_ammonia$Result, fill = "yellow")
#log_result
ammlog <- ggdensity(nutrients_ammonia$log_result, fill = "lightgray") 
#log10_result
ammlog10 <- ggdensity(nutrients_ammonia$log10_result, fill = "steelblue")
ggarrange(amm, ammlog,ammlog10, 
            labels = c("  untransformed", "      log","      log10"), nrow = 1) 

      #in depth normality checking:
      nutrients_ammonia_model <- lm(log10_result ~ station, data = nutrients_ammonia)
      ggqqplot(residuals(nutrients_ammonia_model))
      shapiro_test(residuals(nutrients_ammonia_model))


# Phosphate 
nutrients_phosphate <- nutrients[nutrients$Analyte == "orthoPhosphate", ]  
#untransformed 
p <- ggdensity(nutrients_phosphate$Result, fill = "yellow")
#log_result
plog <- ggdensity(nutrients_phosphate$log_result, fill = "lightgray") 
#log10_result
plog10 <- ggdensity(nutrients_phosphate$log10_result, fill = "steelblue")

pall <- ggarrange(p, plog,plog10,labels = c("  untransformed", "      log","      log10"), nrow = 1)

    #more normality checking        
    nutreints_phosphate_model <- lm(log10_result ~ station, data = nutrients_phosphate)
    ggqqplot(residuals(nutreints_phosphate_model))
    shapiro_test(residuals(nutreints_phosphate_model))

    
# Nitrate/Nitrite 
nutrients_N <- nutrients[nutrients$Analyte == "Nitrate/Nitrite", ]    
# untransformed 
n <- ggdensity(nutrients_N$Result, fill = "yellow")
#log_result
nlog <- ggdensity(nutrients_N$log_result, fill = "lightgray") 
#log10_result
nlog10 <- ggdensity(nutrients_N$log10_result, fill = "steelblue")
ggarrange(n, nlog,nlog10, labels = c("  untransformed", "      log","      log10"), nrow = 1)

nutrients_Chla <- nutrients[nutrients$Analyte == "Chlorophyll", ]
ggdensity(nutrients_Chla$Result)

nutrients_DON <- nutrients[nutrients$Analyte == "DON", ]
# Summary Stats -----------------------
nutrient_sumstat <- nutrients %>% janitor::clean_names() %>% 
  group_by(analyte,station) %>% 
  get_summary_stats(result, type = "mean_sd")

nutrient_sumstat_period <- nutrients %>% janitor::clean_names() %>% 
  group_by(analyte,period) %>% 
  get_summary_stats(result, type = "mean_sd")

ammonia_sumstat <- nutrients_ammonia %>% janitor::clean_names() %>% 
  group_by(station) %>% 
  get_summary_stats(result,  type = "mean_sd")

DON_sumstat <- nutrients_DON %>% 
  group_by(region) %>% 
  get_summary_stats(Result, type = "mean_sd")

# Visualize ------------------------
bp_nuts <- ggplot(nutrients, aes(x = period, y=Result, fill = period))+
  geom_boxplot()+
  stat_boxplot(geom = "errorbar")+
  labs(y="Analyte Concentration", x="Flow Period", fill = "Fow Period")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  guides(fill= "none")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18), 
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), 
        strip.text=element_text(size=15))+
  # ylim(0,30)+
  scale_fill_viridis(option = "A", discrete = TRUE) +
  facet_grid(Analyte~region, scales = "free_y")
bp_nuts


ggplot(nutrients_DON, aes(x = period, y=log_result, fill = period))+
  geom_boxplot()+
  stat_boxplot(geom = "errorbar")+
  labs(y="Analyte Concentration", x="Flow Period", fill = "Fow Period")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  guides(fill= "none")+
  theme_bw()+

# ANOVA Assumptions: ---------------
# 1. Independence 
# 2. No significant outliers
# 3. Normality
# 4. Equal variance

# Outliers: 4 total, only appear at RMB
nutrients %>% 
  janitor::clean_names() %>% 
  group_by(period) %>% 
  identify_outliers(result)

# Normality: nothing is normal.... 
model <- lm(Result ~ station, data = nutrients_ammonia) #linear model with different predictor variables - swap variable after ~ 
model_log <- lm(log_result ~ station, data = nutrients)
model_log10 <- lm(log10_result ~ station, data = nutrients) 
#QQ plot of model residuals 
ggqqplot(residuals(model)) #looks non normal... 
#Shaprio-Wilk test of normality, pval = <0.01 --> DOES NOT MEET ASSUMPTIONS OF NORMALITY
nutrients %>% janitor::clean_names() %>% group_by(analyte,region) %>% shapiro_test(result)

ggqqplot(nutrients, "log_result", facet.by = "Analyte") #normal enough? 

# Equal Variance: 
plot(model, 1)
plot(model_log, 1)
plot(model_log10, 1)

#### ANALYSES: 
# Nested ANOVA -----------------
# fixed effect = region, random effect = station

####### Attempt #1 with full dataframe: 
naov <- lmer(Result ~ region + (1|station), nutrients, REML = TRUE)
  summary(naov)
  anova(naov) # from anova table - region is most important factor but non significant? 
  rand(naov)

#visualize
ggplot(nutrients, aes(x=factor(region), y=Result, fill=station))+geom_boxplot()

#model checking
hist(residuals(naov)) #non normal
plot(fitted(naov),residuals(naov)) #seems to be patterned.. 

#post hoc test comparison of means - trying multiple options... 
ph <- glht(naov, linfct = mcp(region = "Tukey"))
mcs <- summary(ph, test = adjusted("single-step"))
emmeans(naov, specs = pairwise ~ region, adjust="sidak")

####### Attempt #2 with dataframes subset by analyte - better approach
# Ammonia data
naov_amm <- lmer(Result ~ region + (1|station), nutrients_ammonia, REML = TRUE)
  summary(naov_amm)
  anova(naov_amm) #region non significant?
  confint(naov_amm)
  hist(residuals(naov_amm))

# Phosphate data     
naov_phos <- lmer(Result ~ region + (1|station),
                    nutrients_phosphate, REML = TRUE)
  summary(naov_phos)
  anova(naov_phos) #region non significant?
  hist(residuals(naov_phos))
  
  # Nitrate/Nitrite data     
naov_n <- lmer(Result ~ region + (1|station),
                    nutrients_N, REML = TRUE)
  summary(naov_n)
  anova(naov_n) #region non significant?
  hist(residuals(naov_n))
  
#Chlorophyll 
naov_chla <- lmer(Result ~ region + (1|station),
                  nutrients_Chla, REML = TRUE)
summary(naov_chla)
anova(naov_chla) #region non significant?
hist(residuals(naov_chla))
  
# Two Way Mixed ANOVA -----------------  
#trying two way mixed anova in rstatix
maov_amm <- anova_test(data = nutrients_ammonia, dv = Result, wid = station, between = region) 
  
  #wid = specifies sample identifier, between = between-subjects grouping variable, within = within subjects grouping variable 
  #but this isn't a repeated measures ANOVA... 
  get_anova_table(maov_amm) #region non significant predictor of ammonia concentrations 
  
  
  
# Nonparametric Tests: Kruskal-Wallis Test -------------------
# Assumptions of K-W Test: 
  # 1. Trts are independent of one another. 
  # 2. Response variable is ordinal or continuous. 
  # 3. Samples are random.   
  
nutrients %>% kruskal_test(Result~station)
nutrients %>% kruskal_effsize(Result~station)
  
  
## Ammonia
nutrients_ammonia %>% kruskal_test(Result~station) #ammonia varies significantly by station - but not by period or region
nutrients_ammonia %>% kruskal_effsize(Result~station) #large effect of ammonia by station
#Pairwise comparison (cbonferroni too conservative - what are the other options, tukey?) --> “holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”
pwc_ammonia <- nutrients_ammonia %>% dunn_test(Result~station, p.adjust.method = "bonferroni")
View(pwc_ammonia)

## Nitrate + Nitrite
nutrients_N %>% kruskal_test(Result~region) #Nitrogen varies significantly by station - but not by period or region, but only at waste water sites..
nutrients_N %>% kruskal_effsize(Result~station) #large effect of nitrogen by station
#Pairwise comparison (cbonferroni too conservative - what are the other options, tukey?)
pwc_N <- nutrients_N %>% dunn_test(Result~station, p.adjust.method = "bonferroni")  
View(pwc_N)
## Phosphate 
nutrients_phosphate %>% kruskal_test(Result~region) # varies significantly by station, 
nutrients_phosphate %>% kruskal_effsize(Result~station)

pwc_P <- nutrients_phosphate %>% dunn_test(Result~station, p.adjust.method = "bonferroni") 
View(pwc_P)