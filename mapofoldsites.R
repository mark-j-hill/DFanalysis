library(tidyverse)
hpculled <- read_xlsx("HP_ALL_Culled.xlsx", sheet = "Combined Storm")
excel_sheets("HP_ALL_Culled.xlsx")

locations <- read_xlsx("HPLatLong.xlsx")


# convert dms to dd -------------------------------------------------------

locations <- separate(locations, Lat, c("degrees","minutes","seconds"), sep= " ")
locations$degrees <-  as.numeric(locations$degrees)
locations$minutes <- as.numeric(locations$minutes)
locations$seconds <- as.numeric(locations$seconds)
locations <- mutate(locations, Lat = degrees + minutes/60 + seconds/3600)
locations <- locations[,-c(2:4)]

locations <- separate(locations, Long, c("degrees","minutes","seconds"), sep= " ")
locations$degrees <-  as.numeric(locations$degrees)
locations$minutes <- as.numeric(locations$minutes)
locations$seconds <- as.numeric(locations$seconds)
locations <- mutate(locations, Long = degrees + minutes/60 + seconds/3600)
locations <- locations[,-c(2:5)]
locations$Long <- locations$Long * -1
#join location data with sample data

names(locations)[1] <- "LOCCODE"
hpjoined <- left_join(hpculled, locations)


hpjoined <- subset(hpjoined, (!is.na(hpjoined[,23])))
hpjoined <- hpjoined[,c(23,22,1:21)]

#hpjoined <- usmap_transform(hpjoined)
#plot_usmap("counties", include = c("MS"))+
#  geom_point(data = hpjoined,
#            aes(x= LongDD.1, y= LatDD.1, fill= Turbidity ))



hpjoined <- st_as_sf(hpjoined, coords = c("Long", "Lat"))
mapview(hpjoined)


m <- leaflet(hpjoined) %>% 
  addProviderTiles(providers$OpenTopoMap) %>% 
  addMarkers()

m
