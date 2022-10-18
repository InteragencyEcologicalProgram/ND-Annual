# Example Script used to import, prepare and plot daily avg chlorophyll fluorescence data
# Author: Dave Bosworth
# Contact: David.Bosworth@water.ca.gov

# Load packages
library(tidyverse)
library(lubridate)
library(dataRetrieval)
library(scales)
library(patchwork)


# 1. Import and download data ---------------------------------------------

# Download LIB and RYI data from USGS using dataRetrieval package
# Create vectors for start and end dates
start_date <- "2021-07-20" 
end_date <- "2021-10-04" 

lib <- readNWISuv("11455315", parameterCd = "32316", start_date, end_date, tz = "Etc/GMT+8")
ryi <- readNWISuv("11455385", parameterCd = "32316", start_date, end_date, tz = "Etc/GMT+8")

# Import data from WDL
lis <- read_csv("Chlorophyll/B9156000_Chlorophyll_Raw_LIS.csv", skip = 8)

# Import data from CEMP
rvb <- read_csv("Chlorophyll/RVB.csv", skip = 4)


# 2. Prepare data for plot ------------------------------------------------

# Prepare each data frame so that they can be bound together
# USGS:
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

# WDL:
lis_clean <- lis %>% 
  rename(
    DateTime = `Date Time`,
    Chla = `Chlorophyll Raw  (UT)`,
    Chla_Qual = `Quality Code`
  ) %>% 
  mutate(
    DateTime = mdy_hm(DateTime, tz = "Etc/GMT+8"),
    StationCode = "LIS"
  )

# CEMP
rvb_clean <- rvb %>% 
  select(
    DateTime = DATE,
    Chla = VALUE,
    Chla_Qual = `QAQC Flag`
  ) %>% 
  mutate(
    DateTime = mdy_hm(DateTime, tz = "Etc/GMT+8"),
    StationCode = "RVB"
  )

# Bind data together
df_chla <- bind_rows(usgs_clean, lis_clean, rvb_clean)

# Create a vector for the factor order of StationCode
sta_order <- c(
  "RCS",
  "RD22",
  "I80",
  "LIS",
  "TOE",
  "STTD",
  "LIBCUT",
  "LIB",
  "RYI",
  "RVB"
)

# Calculate daily averages and prepare for plot
df_chla_daily_avg <- df_chla %>% 
  mutate(Date = date(DateTime)) %>% 
  group_by(StationCode, Date) %>% 
  summarize(Daily_avg = mean(Chla)) %>% 
  # Fill in missing dates with NA values for geom_line to not interpolate data gaps
  complete(Date = seq.Date(min(Date), max(Date), by = "1 day")) %>%
  ungroup() %>% 
  mutate(
    # Add a grouping variable for Region
    Region = if_else(
      StationCode %in% c("RCS", "RD22", "I80", "LIS", "TOE", "STTD"),
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
  ),
  annotate(
    "rect",
    xmin = ymd("2021-09-11"), 
    xmax = ymd("2021-09-15"),
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.2,
    fill = "grey50"
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
  scale_y_continuous(
    name = "Chlorophyll (µg/L)", 
    limits = c(y_min, y_max), 
    labels = label_comma()
  ) +
  ggtitle("Upstream Region")
  
plt_downstream <- df_chla_daily_avg %>% 
  filter(Region == "Downstream") %>% 
  ggplot(aes(x = Date, y = Daily_avg, color = StationCode)) +
  geom_line() +
  ts_cust_format +
  scale_y_continuous(
    name = NULL, 
    limits = c(y_min, y_max), 
    labels = label_comma()
  ) +
  ggtitle("Downstream Region")

# Combine plots together
plt_upstream + plt_downstream

# Export plot
ggsave(
  "Chlorophyll/chla_plot_2021.jpg",
  width = 7, 
  height = 4, 
  units = "in", 
  dpi = 300
)

# >>>> REMAINING QUESTION:
# 1) I made the y-axis scales the same for the two regions. Do we want to change that?

