
## A. Lucore 06/30/2020
## Written to look at some EPA F2F data

if (!require(pacman)) {
  # do i have this package?
  install.packages("pacman")  # install it if i don't
  library(pacman)  # call this packge, cause i'm gonna need it
}
## use p_load from pacman package to do the above with much less coding
p_load(
  dataRetrieval, # For USGS water data
  lubridate,     # For Dates
  tidyverse      # For cleaning data
)

# library(dataRetrieval)# For USGS water data
# library(lubridate)    # For Dates
# library(beepr)
# library(tidyverse)    # For cleaning data

# clear envrionment
remove(list = ls())


dat.sonde <- read.csv("data/Sonde data.csv",
                       stringsAsFactors = FALSE,
                       na.strings = c("","NA"), # replaces blanks with NAs
                       skip = 2
  ) %>% 
  select(
    c(1:4, 6, 7, 9, 10, 17, 18, 23)
      # drops unneeded columns
  ) %>% 
  rename_all(
    tolower
      # sets column names to only use lower case letters
  ) %>% 
  separate(
    col = date.time,
    into = c("date", "time"),
    sep = " ",
    remove = FALSE,
      # seperates date and time into their own columns
  ) %>% 
    separate(
    col = location.name,
    into = c("code", "site"),
    sep = " ",
    extra = "merge" # since there is more than 1 space, this merges the other two levels
  ) %>% 
  select(
    -c(3)
  ) %>%
  mutate(
    date.time = mdy_hm(date.time),
    date = mdy(date),
    site = factor(site)
  ) %>% 
  mutate(
    date.time = round_date(date.time, "15 mins")
  )

# Adding USGS codes
dat.sonde$usgs <- ifelse(
  dat.sonde$code == "113A30", "341550090391300",
  ifelse(
    dat.sonde$code == "111B40", "07288521",
    ifelse(
      dat.sonde$code == "111A14", "07288048", NA
    )
  )
)

# Fix Lat
dat.sonde$latitude <- ifelse(
  dat.sonde$code == "113A30", "34.264099",
  ifelse(
    dat.sonde$code == "111B40", "33.548207",
    ifelse(
      dat.sonde$code == "111A14", "34.209083", 
      ifelse(
        dat.sonde$code == "202d61", "34.14084", NA
      )
    )
  )
)

# Fix Long
dat.sonde$longitude <- ifelse(
  dat.sonde$code == "113A30", "-90.653717",
  ifelse(
    dat.sonde$code == "111B40", "-90.672393",
    ifelse(
      dat.sonde$code == "111A14", "-90.70025", 
      ifelse(
        dat.sonde$code == "202d61", "-90.653074", NA
      )
    )
  )
)

last <- max(dat.sonde$date)

flowRB <- readNWISuv(
  siteNumbers = "07288048",  
  parameterCd = "00060",     # code for instant height
  startDate = "2020-05-28",
  endDate = last,            # This can also be left blank and it will default to the latest date available
  tz = "America/Chicago"
) %>% 
  renameNWISColumns() %>%    # Gives the columns names instead of just X and X.2
  select(c(2:4)) %>%         # Keeps the data we actually care about
  rename(
    usgs = site_no,
    date.time = dateTime,
    RB = Flow_Inst
    )


flowPB <- readNWISuv(
  siteNumbers = "07288521",
  parameterCd = "00060",     # code for instant height
  startDate = "2020-05-28",
  endDate = last,            # This can also be left blank and it will default to the latest date available
  tz = "America/Chicago"
) %>% 
  renameNWISColumns() %>%    # Gives the columns names instead of just X and X.2
  select(c(2:4)) %>%         # Keeps the data we actually care about
  rename(
    usgs = site_no,
    date.time = dateTime,
    PB = Flow_Inst
    )

flowOS <- readNWISuv(
  siteNumbers = "341550090391300",
  parameterCd = "00060",     # code for instant height
  startDate = "2020-05-28",
  endDate = last,            # This can also be left blank and it will default to the latest date available
  tz = "America/Chicago"
) %>% 
  renameNWISColumns() %>%    # Gives the columns names instead of just X and X.2
  select(c(2:4)) %>%         # Keeps the data we actually care about
  rename(
    usgs = site_no,
    date.time = dateTime,
    OS = Flow_Inst
    )


dat.sonde <- merge(
  x = dat.sonde,
  y = flowOS,
  by = c("date.time", "usgs"),
  all.x = TRUE
) %>%
  merge(y = flowPB,
        by = c("date.time", "usgs"),
        all.x = TRUE
  ) %>%
  merge(y = flowRB,
        by = c("date.time", "usgs"),
        all.x = TRUE)

dat.sonde <- dat.sonde %>% 
  mutate(flow = coalesce(OS, PB, RB)) %>% 
  select(-c(OS, PB, RB)); remove(flowOS, flowPB, flowRB, last)

dat.sonde <- dat.sonde %>% 
  mutate(code = as.factor(code))
