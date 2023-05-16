rm(list = ls()) #removes objects in workspace 

#packages
# remotes::install_github("flowwest/CDECRetrieve") #initially required to load package since update. download directly from github? https://github.com/FlowWest/CDECRetrieve
lapply(c("CDECRetrieve","tidyverse", "lubridate", "viridis", "ggpubr", "remotes", "scales", "data.table"), require, character.only = TRUE)


setwd("C:/Users/mminer/Documents/R_Data_Analysis/NDFS/Flow_2022")

### 2022 LIS Flow Data ########
#call LIS metadata
cdec_datasets("lis") 

#cfs data is sensor #20. 
lis_flows22 <- CDECRetrieve::cdec_query("LIS", "20", "E", "2022-06-01", "2022-10-31")
lis_flows21 <- CDECRetrieve::cdec_query("LIS", "20", "E", "2021-06-01", "2021-10-31")
lis_flows20 <- CDECRetrieve::cdec_query("LIS", "20", "E", "2020-06-01", "2020-10-31")
lis_flows19 <- CDECRetrieve::cdec_query("LIS", "20", "E", "2019-06-01", "2019-10-31") #last wet year for reference
lis_flows16 <- CDECRetrieve::cdec_query("LIS", "20", "E", "2016-06-01", "2016-10-31")

#LIS 2022 format and select####
lis_flows22$location_id=as.factor(lis_flows22$location_id)
lis_flows22=lis_flows22[,c("datetime" ,"location_id","parameter_value")]

lis_flows22$datetime2<-lis_flows22$datetime
lis_flows22_2<-lis_flows22%>%mutate(Month = month(datetime), #create a month and year variable
                                   Year = year(datetime), 
                                   Day = day(datetime), 
                                   Hour= hour(datetime))%>%
  separate(datetime2, into=c("date", "time"), sep=" ")%>%
  group_by(Month, Day)%>%
  mutate(daily_mean_flow=mean(parameter_value, na.rm=TRUE))%>%
  distinct(date, daily_mean_flow, .keep_all=TRUE)

#LIS 2021 format and select ####
lis_flows21$location_id=as.factor(lis_flows21$location_id)
lis_flows21=lis_flows21[,c("datetime" ,"location_id","parameter_value")]

lis_flows21$datetime2<-lis_flows21$datetime
lis_flows21_2<-lis_flows21%>%mutate(Month = month(datetime), #create a month and year variable
                                    Year = year(datetime), 
                                    Day = day(datetime), 
                                    Hour= hour(datetime))%>%
  separate(datetime2, into=c("date", "time"), sep=" ")%>%
  group_by(Month, Day)%>%
  mutate(daily_mean_flow=mean(parameter_value, na.rm=TRUE))%>%
  distinct(date, daily_mean_flow, .keep_all=TRUE)

#LIS 2020 format and select ###################
lis_flows20$location_id=as.factor(lis_flows20$location_id)
lis_flows20=lis_flows20[,c("datetime" ,"location_id","parameter_value")]

lis_flows20$datetime2<-lis_flows20$datetime
lis_flows20_2<-lis_flows20%>%mutate(Month = month(datetime), #create a month and year variable
                                    Year = year(datetime), 
                                    Day = day(datetime), 
                                    Hour= hour(datetime))%>%
  separate(datetime2, into=c("date", "time"), sep=" ")%>%
  group_by(Month, Day)%>%
  mutate(daily_mean_flow=mean(parameter_value, na.rm=TRUE))%>%
  distinct(date, daily_mean_flow, .keep_all=TRUE)


#LIS 2019 foramt and select ####
lis_flows19$location_id=as.factor(lis_flows19$location_id)
lis_flows19=lis_flows19[,c("datetime" ,"location_id","parameter_value")]

lis_flows19$datetime2<-lis_flows19$datetime
lis_flows19_2<-lis_flows19%>%mutate(Month = month(datetime), #create a month and year variable
                                    Year = year(datetime), 
                                    Day = day(datetime), 
                                    Hour= hour(datetime))%>%
  separate(datetime2, into=c("date", "time"), sep=" ")%>%
  group_by(Month, Day)%>%
  mutate(daily_mean_flow=mean(parameter_value, na.rm=TRUE))%>%
  distinct(date, daily_mean_flow, .keep_all=TRUE)

#calculating pulse flow period ####
#consecutive days of net positive flow: kinda works... flow period should just be Sept 21 and 22, 2022
lis_flow_period22 <- lis_flows22_2 %>% 
  group_by(return_rleid = {return_rleid = rle(daily_mean_flow > 0); rep(seq_along(return_rleid$lengths), return_rleid$lengths)}) %>%
  mutate(consec_net_pos = ifelse(daily_mean_flow <= 0, NA, seq_along(return_rleid))) %>% 
  ungroup() #%>% 
  # select(-return_rleid)

#2021 flow period September 11-14. also net positive flow at the end of october 
lis_flow_period21 <- lis_flows21_2 %>% 
  group_by(return_rleid = {return_rleid = rle(daily_mean_flow > 0); rep(seq_along(return_rleid$lengths), return_rleid$lengths)}) %>%
  mutate(consec_net_pos = ifelse(daily_mean_flow <= 0, NA, seq_along(return_rleid))) %>% 
  ungroup() #%>% 
# select(-return_rleid)

lis_flow_period20 <- lis_flows20_2 %>% 
  group_by(return_rleid = {return_rleid = rle(daily_mean_flow > 0); rep(seq_along(return_rleid$lengths), return_rleid$lengths)}) %>%
  mutate(consec_net_pos = ifelse(daily_mean_flow <= 0, NA, seq_along(return_rleid))) %>% 
  ungroup()

lis_flow_period19 <- lis_flows19_2 %>% 
  group_by(return_rleid = {return_rleid = rle(daily_mean_flow > 0); rep(seq_along(return_rleid$lengths), return_rleid$lengths)}) %>%
  mutate(consec_net_pos = ifelse(daily_mean_flow <= 0, NA, seq_along(return_rleid))) %>% 
  ungroup() 

# PLOTS ####
#ORIGINAL: plot all flow data within NDFS monitoring period: ####
a <- lis_flows22_2%>%
  ggplot(aes(x=datetime,y=daily_mean_flow))+
  geom_point(pch=20, size=1, color="black")+
  geom_smooth(size=1, color="blue", method="loess", span=0.1)+
  geom_hline(yintercept = 0, linetype = 3)+
  scale_x_datetime("Date (2022)",breaks=date_breaks("1 week"), labels=date_format("%b %d"))+
  ylab("Flow (CFS)")+
  theme(legend.position = c(0.8, 0.9))+
  annotate("rect", xmin = as.POSIXct("2022-09-21"), xmax = as.POSIXct("2022-09-22"),
           ymin = -Inf, ymax = Inf,fill="dimgrey",alpha = .2)+
  # annotate("text", x=lis_flows2$datetime[1], y=30, label="a", size=8)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size=20))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x=element_text(size=rel(1.5)), axis.text.y=element_text(size=rel(1.5)), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))+
  theme(axis.text.x=element_text(angle=45, hjust = 1))
a
# ggsave("output/lis_flow2022.png", width = 7, height = 4, units = "in")
#plot possible flow window  
b <- lis_flows22_2%>%filter(datetime>"2022-09-01"&datetime<"2022-09-30")%>%
  ggplot(aes(x=datetime,y=daily_mean_flow))+
  geom_point(pch=20,size=1, color="black")+
  geom_smooth(size=1, color="blue",method="loess", span=0.3)+
  scale_x_datetime("Date (2022)",breaks=date_breaks("72 hour"), labels=date_format("%m/%d")) +
  ylab("Flow (CFS)")+
  theme(legend.position = c(0.8, 0.9))+
  annotate("rect", xmin = as.POSIXct("2022-09-21"), xmax = as.POSIXct("2022-09-22"),
             ymin = -Inf, ymax = Inf,fill="darkgrey",alpha = .2)+
  # annotate("text", x=lis_flows2$datetime[94], y=47, label="b", size=8)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size=20))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x=element_text(size=rel(1.5)), axis.text.y=element_text(size=rel(1.5)), axis.title.x = element_text(size=20), axis.title.y = element_text(size=15))+
  theme(axis.text.x=element_text(angle=45, hjust = 1))
b 
  
#arrange pretty
ab <- ggarrange(a,b, labels = c("A","B"), nrow = 2)
ab
    
# save 
# ggsave("output/flow_comparison.png", width = 7, height = 4, units="in")
  
# Multi-year Plot #####
#combine month-day#
lis_flows19_2$md <- format(as.Date(lis_flows19_2$date), "%m/%d")
lis_flows19_2$md = as.POSIXct(lis_flows19_2$md, format = "%m/%d")

lis_flows20_2$md <- format(as.Date(lis_flows20_2$date), "%m/%d")
lis_flows20_2$md = as.POSIXct(lis_flows20_2$md, format = "%m/%d")

lis_flows21_2$md <- format(as.Date(lis_flows21_2$date), "%m/%d")
  lis_flows21_2$md = as.POSIXct(lis_flows21_2$md, format = "%m/%d") 

lis_flows22_2$md <- format(as.Date(lis_flows22_2$date), "%m/%d")
  lis_flows22_2$md = as.POSIXct(lis_flows22_2$md, format = "%m/%d")
  

lis_flows19_2$Year=as.factor(lis_flows19_2$Year)
lis_flows20_2$Year=as.factor(lis_flows20_2$Year)
lis_flows21_2$Year=as.factor(lis_flows21_2$Year)
lis_flows22_2$Year=as.factor(lis_flows22_2$Year)

lims <- as.POSIXct(strptime(c("2022-06-01", "2022-10-31"),format = "%Y-%m-%d"))

comparison_plot <- ggplot()+
  geom_line(data= lis_flows22_2, aes(x=md, y= daily_mean_flow, color = Year), size=1)+
  geom_line(data= lis_flows21_2, aes(x=md, y= daily_mean_flow, color = Year), size=1)+
  geom_line(data = lis_flows20_2, aes(x=md, y=daily_mean_flow, color = Year), size = 1)+
  geom_line(data= lis_flows19_2, aes(x=md, y= daily_mean_flow, color = Year), size=1)+
  geom_hline(yintercept = 0, linetype = 3)+
  scale_x_datetime("",breaks=date_breaks("1 week"), labels=date_format("%b %d"), limits = lims)+
  # scale_color_viridis(option = "A") +
  scale_color_viridis(option = "D", discrete = TRUE)+
  ylab("Flow (CFS)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x=element_text(size=rel(1.5)), axis.text.y=element_text(size=rel(1.5)), 
        axis.title.x = element_text(size=20), axis.title.y = element_text(size=15))+
  theme(axis.text.x=element_text(angle=45, hjust = 1))+
  theme(legend.position = "bottom", legend.title = element_blank(), legend.direction = "horizontal", legend.box.spacing = unit(-0.75, "cm"), legend.text = element_text(size = rel(1.25)))
# ggsave("output/flow_comparison.png", width = 7, height = 4, units="in")
# ggplot()+
#   geom_smooth(data= lis_flows22_2, aes(x=md, y= daily_mean_flow, color = Year), se = FALSE)+
#   geom_smooth(data = lis_flows21_2, aes(x=md, y= daily_mean_flow, color = Year), se = FALSE)+
#   geom_smooth(data = lis_flows20_2, aes(x= md, y= daily_mean_flow, color = Year), se = FALSE)+
#   geom_smooth(data = lis_flows19_2, aes(x=md, y= daily_mean_flow, color = Year), se = FALSE)+
# #add last year there was a sac river action - 2016?
#   scale_x_datetime("",breaks=date_breaks("1 week"), labels=date_format("%b%d"), limits = lims)+
#   # scale_color_viridis(option = "A") +
#   scale_color_viridis(option = "D", discrete = TRUE)+
#   ylab("Flow (CFS)")+
#   theme_linedraw()+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   theme(axis.text.x=element_text(size=rel(1.5)), axis.text.y=element_text(size=rel(1.5)), 
#         axis.title.x = element_text(size=20), axis.title.y = element_text(size=15))+
#   theme(axis.text.x=element_text(angle=45, hjust = 1))+
#   theme(legend.position = "bottom", legend.title = element_blank(), legend.direction = "horizontal", legend.box.spacing = unit(-0.75, "cm"), legend.text = element_text(size = rel(1.25)))

# Identify pulse periods 
lis21_nosum<- lis_flows21 %>% 
  mutate(pulse_period = case_when(
    datetime2 <= ymd_hms("2021-09-11 00:00:00") ~ "Before",
    datetime2 >= ymd_hms ("2021-09-11 00:00:30") &
      datetime2 <= ymd_hms("2021-09-15 00:00:00") ~ "During",
    datetime2 >= ymd_hms("2021-09-15 00:00:30") ~ "After"))

lis21_means<-lis21_nosum %>%
  group_by(pulse_period) %>%
  summarize(cfs_mean=mean(parameter_value, na.rm=TRUE), cfs_sd=sd(parameter_value, na.rm=TRUE))

lis22_nosum <- lis_flows22 %>% 
  mutate(pulse_period = case_when(
    datetime2 <= ymd_hms("2022-09-21 00:00:00") ~ "Before",
    datetime2 >= ymd_hms ("2022-09-21 00:00:00") &
      datetime2 <= ymd_hms("2022-09-23 00:00:00") ~ "During",
    datetime2 >= ymd_hms("2022-09-23 00:00:30") ~ "After"))

lis22_means<-lis22_nosum %>%
  group_by(pulse_period) %>%
  summarize(cfs_mean=mean(parameter_value, na.rm=TRUE), cfs_sd=sd(parameter_value, na.rm=TRUE))

# write_csv(lis_flows19, "output/lis_flows19.csv")
