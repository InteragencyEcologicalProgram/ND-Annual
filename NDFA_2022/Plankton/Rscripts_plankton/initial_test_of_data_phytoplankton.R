#### Set working directory and load libraries ####
setwd("C:/Users/jtorre/Desktop/Github_Repos/ND-Annual/NDFA_2022/Plankton/Rscripts_plankton")

library(tidyverse)
library(patchwork)
library(viridis)
library(janitor)
library(rstatix)
library(ggpubr)
library(broom)
library(AICcmodavg)
library(lme4)
library(car)

# CONSTRUCTING THE DF TO ANALYZE

## Phytoplankton data to analyze ####

### Pulling the phytoplankton data sets ####
NDFS_phytoplankton_data_2021 <- read_csv("../Data_plankton/20211222_DWR_YOLO_BYPASS_PHYTOS_OUTPUT.csv", show_col_types = FALSE) %>% # Phytoplankton data
  clean_names()
phytoplankton_tax_data <- read_csv("C:/Users/jtorre/Desktop/Github_Repos/ND-Annual/NDFA_2022/Plankton/Data_plankton/phyto_group_classification.csv") %>% # Phytoplankton taxonomic data
  clean_names()

### Removing duplicate sites by creating a site object. Eventually this will be part of the metadata pipeline

  # this goes into the automatic filtering
NDFS_station_codes <- c("BL5", "I80", "LIB", "LIS", "RCS", "RD22", "RMB", "RVB", "RYI", "STTD")

NDFS_phytoplankton_data_2021 <- NDFS_phytoplankton_data_2021[NDFS_phytoplankton_data_2021$station_code %in% NDFS_station_codes,]

NDFS_phytoplankton_data_2021$sample_date <- as.Date(NDFS_phytoplankton_data_2021$sample_date, format="%m/%d/%Y")

### Joining the phytoplankton taxonomic data to the phytoplankton data
NDFS_phytoplankton_data_2021 <- left_join(NDFS_phytoplankton_data_2021, phytoplankton_tax_data, by = "genus") 

## Zooplankton data to analyze ####

### Appending the zooplankton files to each other
file_names <- dir("C:/Users/jtorre/Desktop/Github_Repos/ND-Annual/NDFA_2022/Plankton/Data_plankton/20211115 YB Zooplankton ID Data Sheet - SAMPLES 67-143 - Contract 4600012453_EditedCOPY", 
                  pattern = "*150um.csv", 
                  full.names = TRUE)

NDFS_zooplankton_data_2021 <- do.call(rbind, lapply(file_names, read.csv)) %>%
  clean_names()

NDFS_zooplankton_data_2021$date <- as.Date(NDFS_zooplankton_data_2021$date)#, format = "%m/%d/%y")

### subsetting and formatting the data. Will be part of metadata pipeline
NDFS_zooplankton_data_2021 <- NDFS_zooplankton_data_2021[NDFS_zooplankton_data_2021$station %in% NDFS_station_codes,]

NDFS_zooplankton_data_2021$date <- as.Date(NDFS_zooplankton_data_2021$date)

## joining the zooplankton lab and field data 
zoop_field_data <- read_csv("C:/Users/jtorre/Desktop/Github_Repos/ND-Annual/NDFA_2022/Plankton/Data_plankton/NDFS_zoop_field_data.csv", show_col_types = FALSE, skip = 1) %>%
  clean_names()

zoop_field_data <- zoop_field_data %>%
  dplyr::select("measuring_program_short_name", "sampling_event_date", "sampling_area_number", "flow_meter_start_150", "flow_meter_end_150", "net_type") #%>%
  #rename("date" = "sampling_event_date", "station" = "sampling_area_number")

zoop_field_data$sampling_event_date <- as.Date(zoop_field_data$sampling_event_date, format="%m/%d/%Y")

zoop_field_data <- zoop_field_data[zoop_field_data$measuring_program_short_name %in% c("NDFS", "Shared"),] # what are the cutoff for the dates?

zoop_field_data <- zoop_field_data[zoop_field_data$net_type == 150,]
# CONSTRUCTING THE PLANKTON METADATA 

NDFS_zooplankton_data_2021 <- left_join(NDFS_zooplankton_data_2021, zoop_field_data, by = c("date" = "sampling_event_date", "station" = "sampling_area_number")) %>%
  dplyr::select(c("project", "station", "date", "flow_meter_start_150", "flow_meter_end_150","category", "taxon", "count", "subsample", "v1_ml", "sub1_ml", "v2_ml", "sub2_ml"))#, "sampling_time", "region", "site_region"))

# CONSTRUCTING THE PLANKTON METADATA 

###<><><><><><><><><<> Here is where I left off <><><><><><

## Phytoplankton metadata ####

NDFS_phytoplankton_data_2021 <- NDFS_phytoplankton_data_2021 %>%
  rename(total_cells = number_of_cells_per_unit) %>%
  mutate(sampling_time = case_when( # CREATING "sampling_time" COLUMN
    sample_date < as.Date("2021-09-11") ~ "Before",
    sample_date >= as.Date("2021-09-11") & sample_date <= ("2021-09-14") ~ "During",
    sample_date > ("2021-09-14") ~ "After",
    TRUE ~ "Error")) %>% # Error if any dates are wrong 
  mutate(region = case_when( # CREATING "region" COLUMN
    station_code == "RMB"|station_code == "RCS" ~ "Colusa Drain/Ridge Cut",
    station_code == "RD22"|station_code == "I80" ~ "Upper Yolo Bypass",
    station_code == "LIS"|station_code == "STTD" ~ "Lower Yolo Bypass",
    station_code == "BL5"|station_code == "LIB" ~ "Cache Slough Complex",
    station_code == "RYI"|station_code == "RVB" ~ "Lower Sac River",
    station_code == "SHR" ~ "Sherwood Harbor",
    TRUE ~ "Error")) %>% # Error if any station codes are wrong 
  mutate(site_region = case_when(
    region == "Colusa Drain/Ridge Cut" | region == "Upper Yolo Bypass" | region == "Lower Yolo Bypass" ~ "Upstream region",
    region == "Cache Slough Complex" | region == "Lower Sac River" ~ "Downstream region",
    region == "Sherwood Harbor" ~ "Control region",
    TRUE ~ "Error"))%>% #Error if any regions are wrong
  #rename(NDFS_phytoplankton_data_2021, total_cells = number_of_cells_per_unit) %>% #Correction of phytoplankton data per Ted Flynn memo
  mutate(biov_avg = apply(NDFS_phytoplankton_data_2021[,c(31:40)], 1, mean, na.rm = TRUE)) %>% # Average all 10 biovolume measurements for each taxon (BioV.Avg)
  mutate(density = unit_abundance * factor) %>% # Calculate unit abundance per mL (Density)
  mutate(biov_per_mL = total_cells * biov_avg * factor) # Calculate Biovolume density (per mL) per Flynn_Brown memo

## Zooplankton metadata ####

NDFS_zooplankton_data_2021 <- NDFS_zooplankton_data_2021 %>%
  mutate(sampling_time = case_when( # CREATING "sampling_time" COLUMN
    date < as.Date("2021-09-11") ~ "Before",
    date >= as.Date("2021-09-11") & date <= ("2021-09-14") ~ "During",
    date > ("2021-09-14") ~ "After",
    TRUE ~ "Error")) %>% # Error if any dates are wrong 
  mutate(region = case_when(
    station == "RMB"|station == "RCS" ~ "Colusa Drain/Ridge Cut",
    station == "RD22"|station == "I80" ~ "Upper Yolo Bypass",
    station == "LIS"|station == "STTD" ~ "Lower Yolo Bypass",
    station == "BL5"|station == "LIB" ~ "Cache Slough Complex",
    station == "RYI"|station == "RVB" ~ "Lower Sac River",
    station == "SHR" ~ "Sherwood Harbor",
    TRUE ~ "Error")) %>%# Error if any station codes are wrong 
  mutate(site_region = case_when(
    region == "Colusa Drain/Ridge Cut" | region == "Upper Yolo Bypass" | region == "Lower Yolo Bypass" ~ "Upstream region",
    region == "Cache Slough Complex" | region == "Lower Sac River" ~ "Downstream region",
    region == "Sherwood Harbor" ~ "Control region",
    TRUE ~ "Error"))

NDFS_zooplankton_data_2021$CPUE <- #if Meso// eventually this will connect to the pipeline above for automation
  ifelse( # copy over explanation
    (NDFS_zooplankton_data_2021$category != "MICROZOOPLANKTON & NAUPLII"), 
    (NDFS_zooplankton_data_2021$count/((NDFS_zooplankton_data_2021$sub1_ml)/(NDFS_zooplankton_data_2021$v1_ml))/
       (((3.14*0.25)/4)*((abs(NDFS_zooplankton_data_2021$flow_meter_end_150-NDFS_zooplankton_data_2021$flow_meter_start_150)*57560)/999999))
    ),
    ifelse(
      (NDFS_zooplankton_data_2021$category == "MICROZOOPLANKTON & NAUPLII"), 
      NDFS_zooplankton_data_2021$count/((NDFS_zooplankton_data_2021$sub2_ml)/(NDFS_zooplankton_data_2021$v2_ml))/
        (((3.14*0.25)/4)*((abs(NDFS_zooplankton_data_2021$flow_meter_end_150-NDFS_zooplankton_data_2021$flow_meter_start_150)*57560)/999999)),
      0 # This should be error instead, so we can detect where there is an issue.
    )) 

NDFS_zooplankton_data_2021$CPUE[is.na(NDFS_zooplankton_data_2021$CPUE)] <- 0
NDFS_zooplankton_data_2021$log_CPUE <- log(NDFS_zooplankton_data_2021$CPUE)


# CONSTRUCTING THE PLANKTON PLOTS

# remaking figure 49 ####


### Creating factor levels for the plots
NDFS_phytoplankton_data_2021$sampling_time <- factor(NDFS_phytoplankton_data_2021$sampling_time, levels = c("Before", "During", "After"))
NDFS_phytoplankton_data_2021$site_region <- factor(NDFS_phytoplankton_data_2021$site_region, levels = c("Upstream region", "Downstream region"))
NDFS_phytoplankton_data_2021$group <- factor(NDFS_phytoplankton_data_2021$group, levels = c("Diatoms", "Green Algae", "Cyanobacteria", "Cryptophytes", "Golden Algae", "Other"))
NDFS_zooplankton_data_2021$category <- factor(NDFS_zooplankton_data_2021$category, levels =c("Calanoids", "Microzooplankton & Nauplii", "Cyclopoids", "Cladocera", "Macrozooplankton", "Harpacticoids"))

phytoplankton_boxplot <- ggplot(data = NDFS_phytoplankton_data_2021, aes(x = sampling_time, y = log(biov_per_mL), fill = sampling_time)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot() +
  #geom_hline(yintercept = 0) +
  scale_fill_viridis(option = "A", discrete = TRUE) +
  #scale_x_discrete(limits = sampling_time) +
  labs(x = element_blank(),
       y = expression(paste("Log Phytoplankton Biovolume (","  ",mu, m^3, "/", mL,")", sep="")),
       fill = "Sampling Period",
       tag = "A") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none") + 
  facet_grid(. ~ site_region) 

  

NDFS_zooplankton_data_2021$sampling_time <- factor(NDFS_zooplankton_data_2021$sampling_time, levels = c("Before", "During", "After"))
NDFS_zooplankton_data_2021$site_region <- factor(NDFS_zooplankton_data_2021$site_region, levels = c("Upstream region", "Downstream region"))

zooplankton_boxplot <- ggplot(data = NDFS_zooplankton_data_2021, aes(x = sampling_time, y = log(CPUE), fill = sampling_time)) + # I need to get rid of the top labels still
  stat_boxplot(geom = "errorbar")+
  geom_boxplot() +
  #geom_hline(yintercept = 0) +
  scale_fill_viridis(option = "A", discrete = TRUE) +
  #scale_x_discrete(limits = sampling_time)+
  labs(x = element_blank(),
       y = expression(paste("Log Zooplankton CPUE  (","  ",no., "/", m^2,")", sep="")),
       fill = "Sampling Period",
       tag = "B")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") + 
  facet_grid(. ~ site_region) 

combined_boxplot_49 <- ggarrange(phytoplankton_boxplot, zooplankton_boxplot, # I need to annotate each plot still, could potentially do this in powerpoint 
                              ncol = 1, nrow = 2, 
                              legend = "right", common.legend = TRUE)

# remaking figure 50 ####

NDFS_zooplankton_data_2021 <- NDFS_zooplankton_data_2021 %>%
  mutate(category = case_when(
         category == "MICROZOOPLANKTON & NAUPLII" ~ "Microzooplankton & Nauplii",
         category == "CYCLOPOIDS" ~ "Cyclopoids",
         category == "CALANOIDS" ~ "Calanoids",
         category == "CLADOCERA" ~ "Cladocera",
         category == "HARPACTICOIDS" ~ "Harpacticoids",
         category == "MACROZOOPLANKTON" ~ "Macrozooplankton",
         TRUE ~ "Error")) 

NDFS_zooplankton_data_2021 <- NDFS_zooplankton_data_2021 %>%
  filter(category != "Macrozooplankton")

zooplankton_tax_group_plot <- ggplot(data = NDFS_zooplankton_data_2021, aes(x = sampling_time, y = log(CPUE), fill = sampling_time)) +
  stat_boxplot(geom = "errorbar")+
  geom_boxplot() +
  scale_fill_viridis(option = "A", discrete = TRUE) +
  #scale_x_discrete(limits = sampling_time)+
  labs(x = "Sampling Period",
       y = expression(paste("Log Zooplankton CPUE  ("," ",no., "/", m^2,")", sep="")),
       fill = "Sampling Period",
       tag = "B")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_wrap(. ~ category)

NDFS_phytoplankton_data_2021$group <- replace(NDFS_phytoplankton_data_2021$group, NDFS_phytoplankton_data_2021$group == "Dinoflagellates", "Other")
NDFS_phytoplankton_data_2021$group <- replace(NDFS_phytoplankton_data_2021$group, is.na(NDFS_phytoplankton_data_2021$group), "Other")
NDFS_phytoplankton_data_2021 <- NDFS_phytoplankton_data_2021 %>%
  filter(group != 'Other') %>%
  filter(group != "Golden Algae")

phytoplankton_tax_group_plot <- ggplot(data = NDFS_phytoplankton_data_2021, aes(x = sampling_time, y = log(biov_per_mL), fill = sampling_time)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot() +
  scale_fill_viridis(option = "A", discrete = TRUE) +
  #scale_x_discrete(limits = sampling_time) +
  labs(x = "Sampling Period",
       y = expression(paste("Log Phytoplankton Biovolume ("," " ,mu, m^3, "/", mL,")", sep="")),
       fill = "Sampling Period",
       tag = "A") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+#,
        #axis.text.x = element_blank()) +
  facet_wrap(. ~ group)

combined_group_plot_50 <- ggarrange(phytoplankton_tax_group_plot, zooplankton_tax_group_plot, # I need to annotate each plot still, could potentially do this in powerpoint 
                                 ncol = 2, nrow = 1, 
                                 legend = "right", common.legend = TRUE)





# ALTERNATE VERSIONS OF FIGURE 50 ####
# version 2
phytoplankton_tax_group_plot2 <- ggplot(data = NDFS_phytoplankton_data_2021, aes(x = sampling_time, y = log(biov_per_mL), fill = sampling_time)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot() +
  scale_fill_viridis(option = "A", discrete = TRUE) +
  #scale_x_discrete(limits = sampling_time) +
  labs(x = "Sampling Period",
       y = expression(paste("Log Phytoplankton Biovolume ("," " ,mu, m^3, "/", mL,")", sep="")),
       fill = "Sampling Period") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_grid(group ~ site_region)

zooplankton_tax_group_plot2 <- ggplot(data = NDFS_zooplankton_data_2021, aes(x = sampling_time, y = log(CPUE), fill = sampling_time)) +
  stat_boxplot(geom = "errorbar")+
  geom_boxplot() +
  scale_fill_viridis(option = "A", discrete = TRUE) +
  #scale_x_discrete(limits = sampling_time)+
  labs(x = "Sampling Period",
       y = expression(paste("Log Zooplankton Density CPUE  ("," ",no., "/", m^3,")", sep="")),
       fill = "Sampling Period")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_grid(category ~ site_region)

combined_group_plot_502 <- ggarrange(phytoplankton_tax_group_plot2, zooplankton_tax_group_plot2, # I need to annotate each plot still, could potentially do this in powerpoint 
                                    ncol = 2, nrow = 1, 
                                    legend = "right", common.legend = TRUE)

# version 3
phytoplankton_tax_group_plot3 <- ggplot(data = NDFS_phytoplankton_data_2021, aes(x = site_region, y = log(biov_per_mL), fill = site_region)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot() +
  scale_fill_viridis(option = "A", discrete = TRUE) +
  #scale_x_discrete(limits = sampling_time) +
  labs(x = "Sampling Period",
       y = expression(paste("Log Phytoplankton Biovolume ("," " ,mu, m^3, "/", mL,")", sep="")),
       fill = "Sampling Period") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_wrap(. ~ group)

zooplankton_tax_group_plot3 <- ggplot(data = NDFS_zooplankton_data_2021, aes(x = site_region, y = log(CPUE), fill = site_region)) +
  stat_boxplot(geom = "errorbar")+
  geom_boxplot() +
  scale_fill_viridis(option = "A", discrete = TRUE) +
  #scale_x_discrete(limits = sampling_time)+
  labs(x = "Sampling Period",
       y = expression(paste("Log Zooplankton Density CPUE  ("," ",no., "/", m^3,")", sep="")),
       fill = "Sampling Period")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_wrap(. ~ category)

combined_group_plot_503 <- ggarrange(phytoplankton_tax_group_plot3, phytoplankton_tax_group_plot, zooplankton_tax_group_plot3, zooplankton_tax_group_plot, # I need to annotate each plot still, could potentially do this in powerpoint 
                                     ncol = 2, nrow = 2, 
                                     legend = "right")


# STATISTICAL ANALYSIS ####
# ANOVA analysis for figure 49####
# Make the data frames for each "group" of phytoplankton, to automate this I will make a for loop function
phyto_groups_ANOVA_df <- NDFS_phytoplankton_data_2021 %>%
  group_by(site_region, station_code, sampling_time) %>%
  summarize(total_biovolume = mean(biov_per_mL))

interaction_two_way_phyto_groups <- lmer(log(total_biovolume) ~ sampling_time*site_region + (1|station_code), phyto_groups_ANOVA_df, REML = TRUE)
#summary(interaction_two_way_phyto_groups)
Anova(interaction_two_way_phyto_groups, 
      type = "III")

# interaction_two_way_phyto_groups <- aov(log(total_biovolume) ~ sampling_time*site_region + station_code, data = phyto_groups_ANOVA_df)  
# summary(interaction_two_way_phyto_groups)

zoop_categories_ANOVA_df <- NDFS_zooplankton_data_2021 %>%
  group_by(site_region, station, sampling_time) %>%
  summarize(regional_total_cpue = mean(CPUE)) 

interaction_two_way_zoop_categories <- lmer(log(regional_total_cpue) ~ sampling_time*site_region + (1|station), zoop_categories_ANOVA_df, REML = TRUE)  
summary(interaction_two_way_zoop_categories)
Anova(interaction_two_way_zoop_categories, 
      type = "III")
  
# for figure 50: ####
  #one-way ANOVA analysis with sampling period as independent variable for each of the subplots
  # Phytoplankton
  # Diatoms one-way ANOVA
  diatoms <- NDFS_phytoplankton_data_2021 %>%
    filter(group == "Diatoms") %>%
    group_by(sampling_time, site_region, station_code) %>%
    summarize(total_biovolume = mean(biov_per_mL))
  
  one_way_diatoms = lm(log(total_biovolume) ~ sampling_time, data = diatoms)
  summary(one_way_diatoms)  
  Anova(one_way_diatoms,
        type = "II") 
  
  # Cryptophyte one-way ANOVA
  cryptophytes <- NDFS_phytoplankton_data_2021 %>%
    filter(group == "Cryptophytes") %>%
    group_by(sampling_time,site_region, station_code) %>%
    summarize(total_biovolume = mean(biov_per_mL))
  
  one_way_cryptophytes = lm(log(total_biovolume) ~ sampling_time, data = cryptophytes)
  summary(one_way_cryptophytes)  
  Anova(one_way_cryptophytes,
        type = "II") 
  
  # Cyanobacteria one-way ANOVA
  cyanobacteria <- NDFS_phytoplankton_data_2021 %>%
    filter(group == "Cyanobacteria") %>%
    group_by(sampling_time,site_region, station_code) %>%
    summarize(total_biovolume = mean(biov_per_mL))
  
  one_way_cyanobacteria = lm(log(total_biovolume) ~ sampling_time, data = cyanobacteria)
  summary(one_way_cyanobacteria)  
  Anova(one_way_cyanobacteria,
        type = "II") 
  
  # Green Algae one-way ANOVA
  green_algae <- NDFS_phytoplankton_data_2021 %>%
    filter(group == "Green Algae") %>%
    group_by(sampling_time,site_region, station_code) %>%
    summarize(total_biovolume = mean(biov_per_mL))
  
  one_way_green_algae = lm(log(total_biovolume) ~ sampling_time, data = green_algae)
  summary(one_way_green_algae)  
  Anova(one_way_green_algae,
        type = "II")
  
  
  #Zooplankton
  # Calanoids one-way ANOVA
  calanoids <- NDFS_zooplankton_data_2021 %>%
    filter(category == "Calanoids") %>%
    group_by(site_region, sampling_time, station) %>%
    summarize(regional_total_cpue = mean(CPUE) + 1) 
  
  one_way_calanoids = lm(log(regional_total_cpue) ~ sampling_time, data = calanoids)
  summary(one_way_calanoids)  
  Anova(one_way_calanoids,
        type = "II")
  
  # Cladocera one-way ANOVA
  cladocera <- NDFS_zooplankton_data_2021 %>%
    filter(category == "Cladocera") %>%
    group_by(site_region, sampling_time, station) %>%
    summarize(regional_total_cpue = mean(CPUE)) 
  
  one_way_cladocera = lm(log(regional_total_cpue) ~ sampling_time, data = cladocera)
  summary(one_way_cladocera)  
  Anova(one_way_cladocera,
        type = "II")
  
  # Cyclopoids one-way ANOVA
  cyclopoids <- NDFS_zooplankton_data_2021 %>%
    filter(category == "Cyclopoids") %>%
    group_by(site_region, sampling_time, station) %>%
    summarize(regional_total_cpue = mean(CPUE)) 
  
  one_way_cyclopoids = lm(log(regional_total_cpue) ~ sampling_time, data = cyclopoids)
  summary(one_way_cyclopoids)  
  Anova(one_way_cyclopoids,
        type = "II")
  
  # Harpacticoids one-way ANOVA
  harpacticoids <- NDFS_zooplankton_data_2021 %>%
    filter(category == "Harpacticoids") %>%
    group_by(site_region, sampling_time, station) %>%
    summarize(regional_total_cpue = mean(CPUE) + 1) # quick and dirty fix to the problem I had before.
  
  one_way_harpacticoids = lm(log(regional_total_cpue) ~ sampling_time, data = harpacticoids)
  summary(one_way_harpacticoids)  
  Anova(one_way_harpacticoids,
        type = "II")
  
  # Microzooplanton & Nauplii one-way ANOVA
  microzoop_nauplii <- NDFS_zooplankton_data_2021 %>%
    filter(category == "Microzooplankton & Nauplii") %>%
    group_by(site_region, sampling_time, station) %>%
    summarize(regional_total_cpue = mean(CPUE)) # quick and dirty fix to the problem I had before.
  
  one_way_microzoop_nauplii = lm(log(regional_total_cpue) ~ sampling_time, data = microzoop_nauplii)
  summary(one_way_microzoop_nauplii)  
  Anova(one_way_microzoop_nauplii,
        type = "II")
  

  