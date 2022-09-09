#### Set working directory and load libraries ####
wd_path <- file.path("C:", "Users", "jtorre", "Desktop", "Github_Repos", "ND-Annual", "Plankton", "Rscripts_plankton", fsep = "/")
setwd(wd_path)

library(tidyverse)

#### PHYTOPLANKTON DATA ####
#phyto_data_path <- file.path("..", "Data_plantkon", "20211222_DWR_YOLO_BYPASS_PHYTOS_OUTPUT.csv", fsep = "/")
phytoplankton_data <- read_csv("../Data_plankton/20211222_DWR_YOLO_BYPASS_PHYTOS_OUTPUT.csv", show_col_types = FALSE)



