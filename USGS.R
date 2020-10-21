 

# ..Library ---------------------------------------------------------------

library(dataRetrieval)
library(tidyverse)

# Parametere codes
pcod <- parameterCdFile
colnames(pcod)[1] <- "parm_cd"



# PORTER BAYOU 07288521 ---------------------------------------------------
# ..What's available ------------------------------------------------------------------------

# Shows what data the USGS collects at the site.
  #the service arguement differentiates bewteen instantaneous measurments (uv) and periodic meausurments (qw)
  #leave it blank to show all types of data
what <- whatNWISdata(
  #service  = "uv",
  #service = "qw",
  siteNumbers = "07288068"
)

# This combines the list of parameter codes and the list of measurments taken at the site so they are labled
what <- merge(
  x = pcod,
  y = what,
  by = "parm_cd",
  all = FALSE
) %>% 
  select(c(1, 3, 5:6, 8:9, 11:12, 18:19, 27:29))


view(what)


# ..Discharge -------------------------------------------------------------

# REMEMBER
# readNWISuv() is used for instantaneous data, and 
# readNWISqw() is used for periodic data.

dis <- readNWISuv(
  siteNumbers = "07288521",
  parameterCd = "00060",     # code for instant discharge
  startDate = "2020-03-01",
  endDate = "2020-03-31",    # This can also be left blank and it will default to the latest date available
  tz = "America/Chicago"
) %>% 
  renameNWISColumns() %>%    # Gives the columns names instead of just X and X.2
  select(c(2:4)) %>%         # Keeps the data we actually care about
  rename(date.time = dateTime)

# Check the "what" data frame to see the units the function returns

ggplot()+
  geom_point(
    data = dis,
    aes(
      x = dis$dateTime,
      y = dis$Flow_Inst
    )
  )



# ..Total nitrogen --------------------------------------------------------

tn <- readNWISqw(
  siteNumbers = "07288521",
  parameterCd = "71887",     # Code for total nitrogen
  startDate = "2009-11-17",  # This is the earliest date for this parameter
  endDate = Sys.Date(),      # Sys.Date() just returns the current date and time
  tz = "America/Chicago"
)

# This returns a lot more data than the readNWISuv(), so we get rid of a lot more

tn <- tn %>% 
  select(-c(1, 3:18, 23:33))

# The actual values are in the "result_va" column

ggplot()+
  geom_point(
    data = tn,
    aes(
      x = tn$startDateTime,
      y = tn$result_va
    )
  )



# OVERCUP SLOUGH 341550090391300  -----------------------------------------


# RICHES BAYOU 07288048 ---------------------------------------------------

what <- whatNWISdata(
  #service  = "uv",
  #service = "qw",
  siteNumbers = "07288048"
)

# This combines the list of parameter codes and the list of measurments taken at the site so they are labled
what <- merge(
  x = pcod,
  y = what,
  by = "parm_cd",
  all = FALSE
) %>% 
  select(c(1, 3, 5:6, 8:9, 11:12, 18:19, 27:29))

view(what)

# ..Gage height (ft) ------------------------------------------------------

depth <- readNWISuv(
  siteNumbers = "07288048",
  parameterCd = "00065",     # code for instant discharge
  startDate = "2020-05-28",
  # endDate = "2020-03-31",    # This can also be left blank and it will default to the latest date available
  tz = "America/Chicago"
) %>% 
  renameNWISColumns() %>%    # Gives the columns names instead of just X and X.2
  select(c(2:4)) %>%         # Keeps the data we actually care about
  rename(date.time = dateTime)
