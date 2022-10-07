#### Set working directory and load libraries ####
setwd("C:/Users/jtorre/Desktop/Github_Repos/ND-Annual/Plankton/Rscripts_plankton")

library(tidyverse)
library(patchwork)
library(viridis)

#### PHYTOPLANKTON DATA ####

#### PHYTOPLANKTON DATA ####
phytoplankton_data <- read_csv("../Data_plankton/20211222_DWR_YOLO_BYPASS_PHYTOS_OUTPUT.csv", show_col_types = FALSE)
NDFS_phytoplankton_data_2021 <- phytoplankton_data

# Removing duplicate sites by creating a site object.
NDFS_station_codes <- c("BL5", "I80", "LIB", "LIS", "RCS", "RD22", "RMB", "RVB", "RYI", "STTD")

NDFS_phytoplankton_data_2021 <- NDFS_phytoplankton_data_2021[NDFS_phytoplankton_data_2021$StationCode %in% NDFS_station_codes,]

# Adding column Sampling.Date
# Merge boat and land dates into single day
## Convert Date column into date format
NDFS_phytoplankton_data_2021$SampleDate <- as.Date(NDFS_phytoplankton_data_2021$SampleDate, format="%m/%d/%Y")

rename(NDFS_phytoplankton_data_2021, TotalCells = "Number of Cells per unit")

# CONSTRUCTING THE PHYTO PLANKTON METADATA 
NDFS_phytoplankton_data_2021 <- NDFS_phytoplankton_data_2021 %>%
  rename(TotalCells = `Number of cells per unit`) %>%
  mutate(Sampling.Time = case_when( # CREATING "Sampling.Time" COLUMN
    SampleDate < as.Date("2021-09-11") ~ "Before",
    SampleDate >= as.Date("2021-09-11") & SampleDate <= ("2021-09-14") ~ "During",
    SampleDate > ("2021-09-14") ~ "After",
    TRUE ~ "Error")) %>% # Error if any dates are wrong 
  mutate(Region = case_when( # CREATING "Region" COLUMN
    StationCode == "RMB"|StationCode == "RCS" ~ "Colusa Drain/Ridge Cut",
    StationCode == "RD22"|StationCode == "I80" ~ "Upper Yolo Bypass",
    StationCode == "LIS"|StationCode == "STTD" ~ "Lower Yolo Bypass",
    StationCode == "BL5"|StationCode == "LIB" ~ "Cache Slough Complex",
    StationCode == "RYI"|StationCode == "RVB" ~ "Lower Sac River",
    StationCode == "SHR" ~ "Sherwood Harbor",
    TRUE ~ "Error")) %>% # Error if any station codes are wrong 
  mutate(Site.Region = case_when(
    Region == "Colusa Drain/Ridge Cut" | Region == "Upper Yolo Bypass" | Region == "Lower Yolo Bypass" ~ "Upstream Region",
    Region == "Cache Slough Complex" | Region == "Lower Sac River" ~ "Downstream Region",
    Region == "Sherwood Harbor" ~ "Control Region",
    TRUE ~ "Error"))%>%
  mutate(BioV.Avg = apply(NDFS_phytoplankton_data_2021[,c(31:40)], 1, mean, na.rm = TRUE)) %>% # Average all 10 biovolume measurements for each taxon (BioV.Avg)
  mutate(Density = `Unit Abundance` * Factor) %>% # Calculate unit abundance per mL (Density)
  mutate(BioV.per.mL = TotalCells * BioV.Avg * Factor) # Calculate Biovolume density (per mL) per Flynn_Brown memo

## making the phytoplankton plots ####
regions = c("Cache Slough Complex", "Lower Sac River", "Colusa Drain/Ridge Cut", "Upper Yolo Bypass", "Lower Yolo Bypass", "Sherwood Harbor")
sampling.times = c("Before", "During", "After")

# Still need to fix this plot: Need to fix the order of the boxplots and the legend. Also need to add letters to each one and use patchwork to combine them.
magma(3)

phytoplankton_boxplot <- ggplot(data = NDFS_phytoplankton_data_2021, aes(x = Sampling.Time, y = log(BioV.per.mL), fill = Sampling.Time)) +#, fill = Sampling.Time)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot() +
  scale_fill_viridis(option = "A", discrete = TRUE) +
  scale_x_discrete(limits = sampling.times) +
  labs(x = "Sampling Period",
       y = expression(paste("Log Phytoplankton Biovolume (",mu, m^3, "/", mL,")", sep="")),
       fill = "Sampling Period") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  # theme(axis.text.x = element_text(angle = 0, hjust = 1))
  facet_grid(. ~ Site.Region)

#### Making the zooplankton metadata ####

# Pull in the desired CSVs and append them to each other
file_names <- dir("C:/Users/jtorre/Desktop/Github_Repos/ND-Annual/Plankton/Data_plankton/20211115 YB Zooplankton ID Data Sheet - SAMPLES 67-143 - Contract 4600012453_EditedCOPY",  pattern = "*150um.csv", full.names = TRUE)

NDFS_zooplankton_data_2021 <- do.call(rbind, lapply(file_names, read.csv))

NDFS_zooplankton_data_2021 <- NDFS_zooplankton_data_2021[NDFS_zooplankton_data_2021$station %in% NDFS_station_codes,]

NDFS_zooplankton_data_2021$date <- as.Date(NDFS_zooplankton_data_2021$date)

NDFS_zooplankton_data_2021 <- NDFS_zooplankton_data_2021 %>%
  mutate(Sampling.Time = case_when( # CREATING "Sampling.Time" COLUMN
    date < as.Date("2021-09-11") ~ "Before",
    date >= as.Date("2021-09-11") & date <= ("2021-09-14") ~ "During",
    date > ("2021-09-14") ~ "After",
    TRUE ~ "Error")) %>% # Error if any dates are wrong 
  mutate(Region = case_when(
    station == "RMB"|station == "RCS" ~ "Colusa Drain/Ridge Cut",
    station == "RD22"|station == "I80" ~ "Upper Yolo Bypass",
    station == "LIS"|station == "STTD" ~ "Lower Yolo Bypass",
    station == "BL5"|station == "LIB" ~ "Cache Slough Complex",
    station == "RYI"|station == "RVB" ~ "Lower Sac River",
    station == "SHR" ~ "Sherwood Harbor",
    TRUE ~ "Error")) %>%# Error if any station codes are wrong 
  mutate(Site.Region = case_when(
    Region == "Colusa Drain/Ridge Cut" | Region == "Upper Yolo Bypass" | Region == "Lower Yolo Bypass" ~ "Upstream Region",
    Region == "Cache Slough Complex" | Region == "Lower Sac River" ~ "Downstream Region",
    Region == "Sherwood Harbor" ~ "Control Region",
    TRUE ~ "Error"))
# making the zooplankton plots ####

## Merging lab and field data 
zoop_field_data <- read_csv("../Data_plankton/NDFS_zoop_field_data.csv", show_col_types = FALSE, skip = 1)

zoop_field_data$`Sampling Event Date` <- as.Date(zoop_field_data$`Sampling Event Date`, format="%m/%d/%Y")

zoop_field_data <- zoop_field_data[zoop_field_data$`Measuring program short name` %in% c("NDFS", "Shared"),] 
  
zoop_field_data <- zoop_field_data[zoop_field_data$`Net Type` == 150,]

zoop_field_data <- zoop_field_data %>%
  select("Sampling Event Date", "Sampling Area Number", "Flow Meter Start (150)", "Flow Meter End (150)") %>%
  rename("date" = "Sampling Event Date", "station" = "Sampling Area Number")


NDFS_zooplankton_data_2021 <- left_join(NDFS_zooplankton_data_2021, zoop_field_data, by = c("date" = "date", "station" = "station")) %>%
  select(c("project", "station", "date", "Flow Meter Start (150)", "Flow Meter End (150)","category", "taxon", "count", "subsample", "v1_ml", "sub1_ml", "v2_ml", "sub2_ml", "Sampling.Time", "Region", "Site.Region"))

NDFS_zooplankton_data_2021$CPUE <- #if Meso
  ifelse(
    (NDFS_zooplankton_data_2021$category != "MICROZOOPLANKTON & NAUPLII"), 
    (NDFS_zooplankton_data_2021$count/((NDFS_zooplankton_data_2021$sub1_ml)/(NDFS_zooplankton_data_2021$v1_ml))/
       (((3.14*0.25)/4)*((abs(NDFS_zooplankton_data_2021$`Flow Meter End (150)`-NDFS_zooplankton_data_2021$`Flow Meter Start (150)`)*57560)/999999))
    ),
    ifelse(
      (NDFS_zooplankton_data_2021$category == "MICROZOOPLANKTON & NAUPLII"), 
      NDFS_zooplankton_data_2021$count/((NDFS_zooplankton_data_2021$sub2_ml)/(NDFS_zooplankton_data_2021$v2_ml))/
        (((3.14*0.25)/4)*((abs(NDFS_zooplankton_data_2021$`Flow Meter End (150)`-NDFS_zooplankton_data_2021$`Flow Meter Start (150)`)*57560)/999999)),
      0
    )) 

NDFS_zooplankton_data_2021$CPUE[is.na(NDFS_zooplankton_data_2021$CPUE)] <- 0

#NDFS_zooplankton_data_2021$Sampling.Time <- factor(NDFS_zooplankton_data_2021$Sampling.Time, levels = c("Before", "During", "After"))

zooplankton_boxplot <- ggplot(data = NDFS_zooplankton_data_2021, aes(x = Sampling.Time, y = log(CPUE), fill = Sampling.Time)) +#, fill = Sampling.Time)) +
  stat_boxplot(geom = "errorbar")+
  geom_boxplot() +
  scale_fill_viridis(option = "A", discrete = TRUE) +
  scale_x_discrete(limits = sampling.times)+
  labs(x = "Sampling Period",
       y = expression(paste("Log Zooplankton Density CPUE  ("," ",no., "/", m^3,")", sep="")),
       fill = "Sampling Period")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  # theme(axis.text.x = element_text(angle = 0, hjust = 1))
  facet_grid(. ~ Site.Region)

# remaking figure 50 ####

zooplankton_tax_group_plot <- ggplot(data = NDFS_zooplankton_data_2021, aes(x = Sampling.Time, y = log(CPUE), fill = Sampling.Time)) +#, fill = Sampling.Time)) +
  stat_boxplot(geom = "errorbar")+
  geom_boxplot() +
  scale_fill_viridis(option = "A", discrete = TRUE) +
  scale_x_discrete(limits = sampling.times)+
  labs(x = "Sampling Period",
       y = expression(paste("Log Zooplankton Density CPUE  ("," ",no., "/", m^3,")", sep="")),
       fill = "Sampling Period")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  # theme(axis.text.x = element_text(angle = 0, hjust = 1))
  facet_wrap(. ~ category)

# Remaking figure 50, but for phytoplankton ####
phytoplankton_tax_data <- read_csv("C:/Users/jtorre/Desktop/Github_Repos/ND-Annual/Plankton/Data_plankton/phyto_taxonomy_WoRMS_complete.csv")

NDFS_phytoplankton_data_2021 <- left_join(NDFS_phytoplankton_data_2021, phytoplankton_tax_data, by = "Genus") 

NDFS_phytoplankton_data_2021$Sampling.Time <- factor(NDFS_phytoplankton_data_2021$Sampling.Time, levels = c("Before", "During", "After"))

phytoplankton_tax_group_plot <- ggplot(data = NDFS_phytoplankton_data_2021, aes(x = Sampling.Time, y = log(BioV.per.mL), fill = Sampling.Time)) +#, fill = Sampling.Time)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot() +
  scale_fill_viridis(option = "A", discrete = TRUE) +
  scale_x_discrete(limits = sampling.times) +
  labs(x = "Sampling Period",
       y = expression(paste("Log Phytoplankton Biovolume ("," " ,mu, m^3, "/", mL,")", sep="")),
       fill = "Sampling Period") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
    #axis.text.x = element_text(angle = 0, hjust = 1))
  facet_wrap(. ~ Group)
