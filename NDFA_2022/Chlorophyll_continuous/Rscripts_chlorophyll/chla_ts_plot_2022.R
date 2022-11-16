#Project: NDFS
#Script to import data from NWIS, prepare and plot average chlorophyll fluorescence data.
#LIS, STTD, I80 and RD22 data were downloaded from Hydstra and were formatted & cleaned in separate R file called Clean_Hydstra
#9/21/2022
#Author: Elena Huynh, adapted script written by Dave Bosworth

# Load packages
library(tidyverse)
library(lubridate)
library(dataRetrieval)
library(scales)
library(patchwork)


rm(list=ls()) #clear all objects in environment

setwd("Z:/Water Quality Evaluation Section/Projects and Programs/North Delta FLow Action/2022/NDFS 2022 report")
getwd()

##### 1. Import and download data ---------------------------------------------

# Download LIB and RYI data from USGS using dataRetrieval package
# Create vectors for start and end dates
start_date <- "2022-07-27" 
end_date <- "2022-10-03" 

lib <- readNWISuv("11455315", parameterCd = "32316", start_date, end_date, tz = "Etc/GMT+8")
ryi <- readNWISuv("11455385", parameterCd = "32316", start_date, end_date, tz = "Etc/GMT+8")


#import cleaned data saved from Hydstra
lis <- read.csv("LIS_chl_clean.csv")
i80 <- read.csv("I80_chl_clean.csv")
sttd <- read.csv("STTD_chl_clean.csv") 
rd22 <- read.csv("RD22_chl_clean.csv") 
  
#import data from CEMP
rvb <- read.csv("RVB_raw_chl.csv", skip = 4)


##### 2. Prepare data for plot ------------------------------------------------

# Prepare each data frame with the same column names so they can be joined together
# USGS data:
usgs_clean <- list("LIB" = lib, "RYI" = ryi) %>% 
  map(
    ~select(.x, dateTime, contains("32316")) %>% 
      rename(
        DateTime = dateTime,
        Chla_Qual = ends_with("_cd")
      ) %>% 
      rename(Chla = contains("32316"))
  ) %>% 
  bind_rows(.id = "StationCode") %>% 
  as_tibble()

#WQES convert DateTime to POSIXCT and convert Chla_Qual to character
lis_clean <- lis %>% 
  mutate(DateTime = ymd_hms(DateTime, tz = "Etc/GMT+8"))
lis_clean$Chla_Qual <- as.character(lis_clean$Chla_Qual)

i80_clean <- i80 %>% 
  mutate(DateTime = ymd_hms(DateTime, tz = "Etc/GMT+8"))
i80_clean$Chla_Qual <- as.character(i80_clean$Chla_Qual)

sttd_clean <- sttd %>% 
  mutate(DateTime = ymd_hms(DateTime, tz = "Etc/GMT+8"))
sttd_clean$Chla_Qual <- as.character(sttd_clean$Chla_Qual)

rd22_clean <- rd22 %>% 
  mutate(
    DateTime = ymd_hms(DateTime, tz = "Etc/GMT+8"))
rd22_clean$Chla_Qual <- as.character(rd22_clean$Chla_Qual)

#CEMP 
rvb_clean <- rvb %>% 
  select(
    DateTime = DATE,
    Chla = VALUE,
    Chla_Qual = `QAQC.Flag`
  ) %>% 
  mutate(
    DateTime = mdy_hm(DateTime, tz = "Etc/GMT+8"),
    StationCode = "RVB"
  ) 
  


# Bind data together
df_chla <- bind_rows(usgs_clean, lis_clean, rd22_clean, sttd_clean, i80_clean, rvb_clean)

# Create a vector for the factor order of StationCode
sta_order <- c(
  "RD22",
  "I80",
  "LIS",
  "STTD",
  "LIB",
  "RYI",
  "RVB"
)

# Calculate daily averages and prepare for plot
df_chla_daily_avg <- df_chla %>% 
  mutate(Date = date(DateTime)) %>% 
  group_by(StationCode, Date) %>%
  drop_na() %>%
  summarize(Daily_avg = mean(Chla)) %>% 
  # Fill in missing dates with NA values for geom_line to not interpolate data gaps
  complete(Date = seq.Date(min(Date), max(Date), by = "1 day")) %>%
  ungroup() %>% 
  mutate(
    # Add a grouping variable for Region
    Region = if_else(
      StationCode %in% c("RD22", "I80", "LIS", "STTD"),
      "Upstream",
      "Downstream"
    ),
    # Apply factor order to StationCode
    StationCode = factor(StationCode, levels = sta_order)
  )

# 3. Create Plot ----------------------------------------------------------

# Create a list for custom formatting of time-series plots
ts_cust_format <- list(
  theme_light(),
  theme(
    strip.text = element_text(color = "black"),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  ),
  scale_x_date(
    name = "Date",
    breaks = breaks_pretty(10),
    labels = label_date_short(),
    # define limits of the x-axes to make them consistent
    limits = c(ymd(start_date), ymd(end_date)),
    expand = expansion(mult = 0.01)
  ),
  scale_color_viridis_d(
    name = "Station:", 
    option = "plasma", 
    end = 0.95
  )
)

# Define max and min values for the y-axes of the time-series plots
y_min <- floor(min(df_chla_daily_avg$Daily_avg, na.rm = TRUE))
y_max <- ceiling(max(df_chla_daily_avg$Daily_avg, na.rm = TRUE))

# Create time-series plots of continuous WQ Data by Region
plt_upstream <- df_chla_daily_avg %>% 
  filter(Region == "Upstream") %>% 
  ggplot(aes(x = Date, y = Daily_avg, color = StationCode)) +
  geom_line() +
  ts_cust_format +
  scale_x_date(labels = label_date_short(c(NA, "%b", "%d", NA)))+
  scale_y_continuous(
    name = "Chlorophyll (µg/L)", 
    limits = c(y_min, y_max), 
    labels = label_comma()
  ) +
  ggtitle("Upstream Region") +
  annotate(
    "rect",
    xmin = ymd("2022-09-21"), 
    xmax = ymd("2022-09-22"),
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.2,
    fill = "grey50"
  )
plt_upstream

plt_downstream <- df_chla_daily_avg %>% 
  filter(Region == "Downstream") %>% 
  ggplot(aes(x = Date, y = Daily_avg, color = StationCode)) +
  geom_line() +
  ts_cust_format +
  scale_x_date(labels = label_date_short(c(NA, "%b", "%d", NA)))+
  scale_y_continuous(
    name = NULL, 
    limits = c(y_min, y_max), #customized y max based on this year's data
    labels = label_comma()
  ) +
  ggtitle("Downstream Region") + 
  annotate(
    "rect",
    xmin = ymd("2022-09-21"), 
    xmax = ymd("2022-09-22"),
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.2,
    fill = "grey50"
  ) 

plt_downstream

# Combine plots together
plt_upstream + plt_downstream

# Export plot
ggsave(
  "chla_plot_2022.jpg",
  width = 8.25, 
  height = 4, 
  units = "in", 
  dpi = 300
)


