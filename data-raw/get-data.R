library(rgdal)  # Geospatial Data Abstraction Library 
library(dplyr)

shape <- readOGR("data-raw/", layer="newDistribuidoras", encoding="utf-8")

d1 <- read.csv("./data-raw/ANEEL407_Bias_Corrected.csv")  # Efficiency variables.
d2 <- read.csv("./data-raw/ANEEL_Ambientais_NEW.csv")  # Environment variables.

dados <- full_join(d1, d2,  by = c("Codigo" = "d_Code"))

shape <- merge(shape, dados, by.x = "Codig", by.y="Codigo", all.x=FALSE)
aneelshape <- subset(shape, Codig != "D05") ## Retira "Roraima"

devtools::use_data(aneelshape, overwrite = FALSE)
