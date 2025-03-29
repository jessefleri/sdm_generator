defs <- c(
  "Annual Mean Temperature",
  "Mean Diurnal Range",
  "Isothermality",
  "Temperature Seasonality",
  "Max Temperature of Warmest Month",
  "Min Temperature of Coldest Month",
  "Temperature Annual Range",
  "Mean Temperature of Wettest Quarter",
  "Mean Temperature of Driest Quarter",
  "Mean Temperature of Warmest Quarter",
  "Mean Temperature of Coldest Quarter",
  "Annual Precipitation",
  "Precipitation of Wettest Month",
  "Precipitation of Driest Month",
  "Precipitation Seasonality",
  "Precipitation of Wettest Quarter",
  "Precipitation of Driest Quarter",
  "Precipitation of Warmest Quarter",
  "Precipitation of Coldest Quarter"
)

bios <- str_c("bio_", 1:19)

df <- data.frame(bios, defs)

write.csv(df, "bioclim_lookup.csv", row.names = FALSE)
