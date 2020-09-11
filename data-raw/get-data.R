shape <- rgdal::readOGR("data-raw/", layer = "distribuidoras", encoding = "latin-1")

df_dea <- read.csv("./data-raw/aneel-407-bias-corrected.csv", encoding = "latin-1")
df_env <- read.csv("./data-raw/aneel-environmental-data.csv", encoding = "latin-1")

df_full <- dplyr::full_join(df_dea, df_env, by = c("Codigo" = "d_Code"))

shape <- merge(shape, df_full, by.x = "Codig", by.y = "Codigo", all.x = FALSE)
aneelshape <- subset(shape, Codig != "D05")

usethis::use_data(aneelshape, overwrite = FALSE, version = 3)
