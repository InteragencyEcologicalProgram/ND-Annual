#Project: NDFS 
#Purpose: Import and clean real time monitoring data.
#Data were imported from Hydstra and converted to .xlsx for cleaning
#Author: Elena Huynh
#9/21/2022

library(tidyverse)
library(dplyr)
library(readxl)

rm(list=ls()) #clear all objects in environment


#import data, skip first 2 rows with station and parameter identification info
lis <- read_excel("Y:/Projects and Programs/North Delta Flow Action/2022/NDFS 2022 report/LIS_Hyds_chl.xlsx", skip = 2)
i80 <- read_excel("Y:/Projects and Programs/North Delta Flow Action/2022/NDFS 2022 report/I80_Hyds_chl.xlsx", skip = 2)
sttd <- read_excel("Y:/Projects and Programs/North Delta Flow Action/2022/NDFS 2022 report/STTD_Hyds_chl.xlsx", skip = 2) 
rd22 <- read_excel("Y:/Projects and Programs/North Delta Flow Action/2022/NDFS 2022 report/RD22_Hyds_chl.xlsx",  skip = 2) 


#remove comments column, rename columns, standardize the first reading and add StationCode column (and move it to 1st column), keep data that are good quality (=1)
lis_clean <- lis %>% rename("DateTime" = "Date",  "Chla" = "Point",  "Chla_Qual" = "Qual") %>%
  filter(DateTime >= "2022-07-27 06:15:00") %>%
  mutate(StationCode = "LIS", .before = 1) %>%
  filter(Chla_Qual == "1")

i80_clean <- i80 %>% rename("DateTime" = "Date",  "Chla" = "Point",  "Chla_Qual" = "Qual") %>%
  filter(DateTime >= "2022-07-27 06:15:00") %>%
  mutate(StationCode = "I80", .before = 1) %>%
  filter(Chla_Qual == "1")

sttd_clean <- sttd %>% rename("DateTime" = "Date",  "Chla" = "Point",  "Chla_Qual" = "Qual") %>%
  filter(DateTime >= "2022-07-27 06:15:00") %>%
  mutate(StationCode = "STTD", .before = 1) %>%
  filter(Chla_Qual == "1")

rd22_clean <- rd22 %>% rename("DateTime" = "Date",  "Chla" = "Point",  "Chla_Qual" = "Qual") %>%
  filter(DateTime >= "2022-07-27 06:15:00") %>%
  mutate(StationCode = "RD22", .before = 1) %>%
  filter(Chla_Qual == "1")


#save clean data frames
write.csv(lis_clean, file = "Y:/Projects and Programs/North Delta Flow Action/2022/NDFS 2022 report/LIS_chl_clean.csv", row.names=F)
write.csv(i80_clean, file = "Y:/Projects and Programs/North Delta Flow Action/2022/NDFS 2022 report/I80_chl_clean.csv", row.names=F)
write.csv(sttd_clean, file = "Y:/Projects and Programs/North Delta Flow Action/2022/NDFS 2022 report/STTD_chl_clean.csv", row.names=F)
write.csv(rd22_clean, file = "Y:/Projects and Programs/North Delta Flow Action/2022/NDFS 2022 report/RD22_chl_clean.csv", row.names=F)
