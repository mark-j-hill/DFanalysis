# HEADER ------------------------------------------------------------------
# A. Lucore, M. Hill 08/17/2020
# Written as the conclusion of NIFA data cleaning and manipulating for the 
# water quality lab's NIFA project.

# This script has 4 major parts. They are:
# Read and clean
# Logic Check
# Pair data
# Reduction potential

# More detailed documentation will be given below in each section but briefly:
# Data is read in and manipulated so that R will be able to work with it. Our
# logic columns are checked and any data that falls outside of our acceptable 
# range is removed (the logic columns are then removed). A new data frame 
# containing only events where paired samplers collected samples is created, 
# allowing for a closer comparison to be made. Finally, from the paired
# dataframe, a dataframe showing "reduction potential" is created showing the 
# difference and percent difference between experimental and control plots.


# ..Library ---------------------------------------------------------------
rm(list= ls())
library(lubridate)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(ggsignif)
library(stats)
library(zoo)
library(scales)
library(gridExtra)
library(grid)
library(mosaic)
library(mosaicData)
library(patchwork)
library(mice)
# READ AND CLEAN ----------------------------------------------------------

# This is where the data is initially read in. Nothing special here, just
# some basic tidyverse stuff.
setwd("/Users/mark/Water Quality/DeltaFARManalysis")
dat <- read.csv("DF_final.csv",
                stringsAsFactors = FALSE,
                na.strings = c("", "NA") # replaces blanks with NAs
               # skip = 1 # skips to the first row for data
) %>% 
mutate( 
    # changes the formatting of some of the columns to make R happier
    Date = lubridate::mdy(Date),
    # makes the dates as dates R will recognize
    sampleid = as.factor(sampleid),
    # Sets the sampleid column as a factor instead of text strings
    farm = as.factor(farm),
    # Sets the Farm column as a factor instead of text strings
    treatment = as.factor(treatment),
    # Sets the Treatment column as a factor instead of text strings
    TUR_ntu = as.numeric(TUR_ntu),      
    # Sets turbidity as a number instead of text strings 
    #(I have no idea why it wasn't already)
    sample_type = as.factor(sample_type)
  )

# replace underscores in the column headers with periods
names(dat) <- gsub(
  x = names(dat), 
  pattern = "\\_", 
  replacement = "."
)

dat$season <- ifelse(
  month(dat$Date) == 11 | month(dat$Date) ==  12 | month(dat$Date) == 1| 
    month(dat$Date) == 2 |  month(dat$Date) == 3 | month(dat$Date) == 4,
  "cover",
  "cash")

dat <- dat %>% 
  mutate(season = as.factor(season))


# add location data -------------------------------------------------------
 
df_samplesites <- read.csv("df_samplesites.csv")
dat <- left_join(dat, df_samplesites)

# change locations for fields that moved
dat[dat$sampleid == "ARR1",]$lat <-
  ifelse(
    dat[dat$sampleid == "ARR1",]$Date >="2019-10-18",
    33.7475234, # New location 
    33.743923   # old location 
  )
dat[dat$sampleid == "ARR1",]$long <-
  ifelse(
    dat[dat$sampleid == "ARR1",]$Date >="2019-10-18",
    -90.439047, # New location 
    -90.434174   # old location 
  )

dat[dat$sampleid == "ARR2",]$lat <-
  ifelse(
    dat[dat$sampleid == "ARR2",]$Date >="2019-10-18",
    33.747566, # New location 
    33.7439064   # old location 
  )
dat[dat$sampleid == "ARR2",]$long <-
  ifelse(
    dat[dat$sampleid == "ARR2",]$Date >="2019-10-18",
    -90.441001, # New location 
    -90.431302   # old location 
  )

dat[dat$sampleid == "STU1",]$lat <-
  ifelse(
    dat[dat$sampleid == "STU1",]$Date >="2019-12-5",
    33.8177149, # New location 
    33.7921411   # old location 
  )
dat[dat$sampleid == "STU1",]$long <-
  ifelse(
    dat[dat$sampleid == "STU1",]$Date >="2019-12-5",
    -90.252741, # New location 
    -90.267364   # old location 
  )

dat[dat$sampleid == "STU2",]$lat <-
  ifelse(
    dat[dat$sampleid == "STU2",]$Date >="2019-12-5",
    33.8194787, # New location 
    33.7955973   # old location 
  )
dat[dat$sampleid == "STU2",]$long <-
  ifelse(
    dat[dat$sampleid == "STU2",]$Date >="2019-12-5",
    -90.258794, # New location 
    -90.272795   # old location 
  )


specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
dat$lat <- specify_decimal(dat$lat, 7)
dat$long <- specify_decimal(dat$long, 7)

# LOGIC CHECK -------------------------------------------------------------

# Here the standard columns are checked and any data from dates where the 
# standards didn't meet QAQC objectives are removed.

# First, replace concentration data on bad dates
# list of columns that have concentration data
conc <- c(10,11,13:15,17,19,21,23,25,27,31:39)
  
# This basically says, "For all the columns with concentration data, if the 
# QAQC for the day was bad, set the concentration data to 'NA', if not, then
# do not change it". All of the loops follow the same basic pattern.
for (i in conc) {
  dat[,i] <- ifelse(
    dat[,2] == FALSE,
    NA,
    dat[,i]
  )
}

# Then, replace concentration on individual dates
# list of logical column numbers besides qaqc.check
logic <- c(12, 16, 18, 20, 22, 24, 26)

# Since the check columns are behind the columns they check, we can just change
# the column before the logical one.
for (i in logic) {
  dat[,i - 1] <- ifelse(
    dat[,i] == FALSE,
    NA,
    dat[,i - 1]
  )
}


# replace NO3 values b/c they are tied to NOx check
# we do this seperately b/c it doesn't fit the pattern of the for() loop above
dat[,27] <- ifelse(
  dat[,24] == FALSE,
  NA,
  dat[,27]
)

# Next we replace data in the loads columns. Since they are spaced differently,
# it has to be done in 3 different loops.

# replace tss loads
for (i in c(21, 24, 27)) { 
  dat[,10 + i] <- ifelse(
    is.na(dat[,10]) == TRUE,
    NA,
    dat[,10 + i]
  )
}

# replace tn loads
for (i in c(17, 20, 23)) { 
  dat[,15 + i] <- ifelse(
    is.na(dat[,15]) == TRUE,
    NA,
    dat[,15 + i]
  )
}

# replace tip loads
for (i in c(16, 19, 22)) { 
  dat[,17 + i] <- ifelse(
    is.na(dat[,17]) == TRUE,
    NA,
    dat[,17 + i]
  )
}

# Finally, we remove the columns that have the logical checks (because I
# the rest of the code relies on absolute places in the dataframe).
# remove logical columns
dat <- dat %>% 
  select(-c(2, 12, 16, 18, 20, 22, 24, 26))

rm(conc, i, logic)


# Impute missing values ---------------------------------------------------
#dat <- dat %>%
#  mutate(
#    farm = as.character(farm),
#    treatment = as.character(treatment)
#  )
#
#dat$treatment <-
#  ifelse(
#    is.na(dat$treatment),
#    "INST",
#    dat$treatment
#  )
#dat$farm <-
#  ifelse(
#    is.na(dat$farm),
#    "INST",
#    dat$farm
#  )
#
#dat <- mutate(
#  dat,
#  farm = as.factor(farm),
#  treatment = as.factor(treatment)
#)
#
# ..Sediment ---------------------------

impdat <-
  mice(
    data = dat[,c(9, 10, 14)],
    m = 5,
    maxit = 50,
    method = 'pmm',
    seed = 1234,
    printFlag = FALSE
  )

set.seed(1234)  # set seed for random imputed dataset selection
compdat <- complete(impdat, round(runif(1, 1, 5), 0))
# randomly select imputed dataset and combine with complete dataset

# inspect imputed data
# imputed datasets in magenta; observed in blue. are values plausible?
#stripplot(impdat)
# do imputed points fall in likely area of distrib.?

dat[,c(9, 10, 14)] <- compdat[,1:ncol(compdat)]

# ..Ortho P -------------------------


impdat <-
  mice(
    data = dat[,c(14, 16)],
    m = 5,
    maxit = 50,
    method = 'pmm',
    seed = 1234,
    printFlag = FALSE
  )

set.seed(1234)  # set seed for random imputed dataset selection
compdat <- complete(impdat, round(runif(1, 1, 5), 0))
# randomly select imputed dataset and combine with complete dataset

# inspect imputed data
# imputed datasets in magenta; observed in blue. are values plausible?
stripplot(impdat)
# do imputed points fall in likely area of distrib.?

dat[,c(14, 16)] <- compdat[,1:ncol(compdat)]

# ..Nitrogen ------------------------------------

impdat <-
  mice(
    data = dat[,c(11:13, 15, 17)],
    m = 5,
    maxit = 50,
    method = 'pmm',
    seed = 1234,
    printFlag = FALSE
  )

set.seed(1234)  # set seed for random imputed dataset selection
compdat <- complete(impdat, round(runif(1, 1, 5), 0))
# randomly select imputed dataset and combine with complete dataset

# inspect imputed data
# imputed datasets in magenta; observed in blue. are values plausible?
stripplot(impdat)
# do imputed points fall in likely area of distrib.?

dat[,c(11:13, 15, 17)] <- compdat[,1:ncol(compdat)]

# ..No3 -----------------------------------

impdat <-
  mice(
    data = dat[,c(11:13, 15, 19)],
    m = 5,
    maxit = 50,
    method = 'pmm',
    seed = 1234,
    printFlag = FALSE
  )

set.seed(1234)  # set seed for random imputed dataset selection
compdat <- complete(impdat, round(runif(1, 1, 5), 0))
# randomly select imputed dataset and combine with complete dataset

# inspect imputed data
# imputed datasets in magenta; observed in blue. are values plausible?
stripplot(impdat)
# do imputed points fall in likely area of distrib.?

dat[,c(11:13, 15, 19)] <- compdat[,1:ncol(compdat)]

rm(impdat)
rm(compdat)

# repair MDLs -------------------------------------------------------------
if (sum(dat$TSS.mgl < 5, na.rm = TRUE) / sum(!is.na(dat$TSS.mgl)) < 0.05) {
  print("Less than 5% of TSS values are below MDL")
  print(sum(dat$TSS.mgl < 5, na.rm = TRUE) / sum(!is.na(dat$TSS.mgl)) *100)
  dat <- dat %>% mutate(TSS.mgl = replace(TSS.mgl, TSS.mgl < 5, 5/2))
} else {
  print("More than 5% of TSS values are below MDL and")
  print("the values have remained unchanged")
  print(dat %>% count(TSS.mgl < 5))
  print(sum(dat$TSS.mgl < 5, na.rm = TRUE) / sum(!is.na(dat$TSS.mgl)) *100)}


if (sum(dat$TUR.ntu < .02, na.rm = TRUE) / sum(!is.na(dat$TUR.ntu)) < 0.05) {
  print("Less than 5% of TUR values are below MDL")
  print(sum(dat$TUR.ntu < .02, na.rm = TRUE) / sum(!is.na(dat$TUR.ntu)) *100)
  dat <- dat %>% mutate(TUR.ntu = replace(TUR.ntu, TUR.ntu < .02, .02/2))
} else {
  print("More than 5% of TUR values are below MDL and")
  print("the values have remained unchanged")
  print(dat %>% count(TUR.ntu < .02))
  print(sum(dat$TUR.ntu < .02, na.rm = TRUE) / sum(!is.na(dat$TUR.ntu)) *100)}


if (sum(dat$NO3.NO2.mgl < .23, na.rm = TRUE) / sum(!is.na(dat$NO3.NO2.mgl)) < 0.05) {
  print("Less than 5% of NO3.NO2 values are below MDL")
  print(sum(dat$NO3.NO2.mgl < .23, na.rm = TRUE) / sum(!is.na(dat$NO3.NO2.mgl)) *100)
  dat <- dat %>% mutate(NO3.NO2.mgl = replace(NO3.NO2.mgl, NO3.NO2.mgl < .23, .23/2))
} else {
  print("More than 5% of NO3.NO2 values are below MDL and")
  print("the values have remained unchanged")
  print(dat %>% count(NO3.NO2.mgl < .23))
  print(sum(dat$NO3.NO2.mgl < .23, na.rm = TRUE) / sum(!is.na(dat$NO3.NO2.mgl)) *100)}


if (sum(dat$TKN.mgl < .77, na.rm = TRUE) / sum(!is.na(dat$TKN.mgl)) < 0.05) {
  print("Less than 5% of TKN values are below MDL")
  print(sum(dat$TKN.mgl < .77, na.rm = TRUE) / sum(!is.na(dat$TKN.mgl)) *100)
  dat <- dat %>% mutate(TKN.mgl = replace(TKN.mgl, TKN.mgl < .77, .77/2))
} else {
  print("More than 5% of TKN values are below MDL and")
  print("the values have remained unchanged")
  print(dat %>% count(TKN.mgl < .77))
  print(sum(dat$TKN.mgl < .77, na.rm = TRUE) / sum(!is.na(dat$TKN.mgl)) * 100)
  }


if (sum(dat$TN.mgl < 1, na.rm = TRUE) / sum(!is.na(dat$TN.mgl)) < 0.05) {
  print("Less than 5% of TN values are below MDL")
  print(sum(dat$TN.mgl < 1, na.rm = TRUE) / sum(!is.na(dat$TN.mgl)) *100)
  dat <- dat %>% mutate(TN.mgl = replace(TN.mgl, TN.mgl < 1, 1/2))
} else {
  print("More than 5% of TN values are below MDL and")
  print("the values have remained unchanged")
  print(dat %>% count(TN.mgl < 1))
  print(sum(dat$TN.mgl < 1, na.rm = TRUE) / sum(!is.na(dat$TN.mgl)) *100)
  }


if (sum(dat$TIP.mgl < .05, na.rm = TRUE) / sum(!is.na(dat$TIP.mgl)) < 0.05) {
  print("Less than 5% of TIP values are below MDL")
  print(sum(dat$TIP.mgl < .05, na.rm = TRUE) / sum(!is.na(dat$TIP.mgl)) *100)
  dat <- dat %>% mutate(TIP.mgl = replace(TIP.mgl, TIP.mgl < .05, .05/2))
} else {
  print("More than 5% of TIP values are below MDL and")
  print("the values have remained unchanged")
  print(dat %>% count(TIP.mgl < .05))
  print(sum(dat$TIP.mgl < .05, na.rm = TRUE) / sum(!is.na(dat$TIP.mgl)) *100)}


if (sum(dat$NH3.mgl < .01, na.rm = TRUE) / sum(!is.na(dat$NH3.mgl)) < 0.05) {
  print("Less than 5% of NH3 values are below MDL")
  print(sum(dat$NH3.mgl < .01, na.rm = TRUE) / sum(!is.na(dat$NH3.mgl)) *100)
  dat <- dat %>% mutate(NH3.mgl = replace(NH3.mgl, NH3.mgl < .01, .01/2))
} else {
  print("More than 5% of NH3 values are below MDL and")
  print("the values have remained unchanged")
  print(dat %>% count(NH3.mgl < .01))
  print(sum(dat$NH3.mgl < .01, na.rm = TRUE) / sum(!is.na(dat$NH3.mgl)) *100)}


if (sum(dat$OrthoP.mgl < .01, na.rm = TRUE) / sum(!is.na(dat$OrthoP.mgl)) < 0.05) {
  print("Less than 5% of Ortho Phosphorus values are below MDL")
  print(sum(dat$OrthoP.mgl < .01, na.rm = TRUE) / sum(!is.na(dat$OrthoP.mgl)) *100)
  dat <- dat %>% mutate(OrthoP.mgl = replace(OrthoP.mgl, OrthoP.mgl < .01, .01/2))
} else {
  print("More than 5% of Ortho Phosphorus values are below MDL and")
  print("the values have remained unchanged")
  print(dat %>% count(OrthoP.mgl < .01))
  print(sum(dat$OrthoP.mgl < .01, na.rm = TRUE) / sum(!is.na(dat$OrthoP.mgl)) *100)}


if (sum(dat$NOx.mgl < .007, na.rm = TRUE) / sum(!is.na(dat$NOx.mgl)) < 0.05) {
  print("Less than 5% of NOx values are below MDL")
  print(sum(dat$NOx.mgl < .007, na.rm = TRUE) / sum(!is.na(dat$NOx.mgl)) *100)
  dat <- dat %>% mutate(NOx.mgl = replace(NOx.mgl, NOx.mgl < .007, .007/2))
} else {
  print("More than 5% of NOx values are below MDL and")
  print("the values have remained unchanged")
  print(dat %>% count(NOx.mgl < .007))
  print(sum(dat$NOx.mgl < .007, na.rm = TRUE) / sum(!is.na(dat$NOx.mgl)) *100)}


if (sum(dat$NO2.mgl < .05, na.rm = TRUE) / sum(!is.na(dat$NO2.mgl)) < 0.05) {
  print("Less than 5% of NO2 values are below MDL")
  print(sum(dat$NO2.mgl < .05, na.rm = TRUE) / sum(!is.na(dat$NO2.mgl)) *100)
  dat <- dat %>% mutate(NO2.mgl = replace(NO2.mgl, NO2.mgl < .05, .05/2))
} else {
  print("More than 5% of NO2 values are below MDL and")
  print("the values have remained unchanged")
  print(dat %>% count(NO2.mgl < .05))
  print(sum(dat$NO2.mgl < .05, na.rm = TRUE) / sum(!is.na(dat$NO2.mgl)) *100)}


if (sum(dat$NO3.mgl < .0007, na.rm = TRUE) / sum(!is.na(dat$NO3.mgl)) < 0.05) {
  print("Less than 5% of NO3 values are below MDL")
  print(sum(dat$NO3.mgl < .0007, na.rm = TRUE) / sum(!is.na(dat$NO3.mgl)) *100)
  dat <- dat %>% mutate(NO3.mgl = replace(NO3.mgl, NO3.mgl < .0007, .0007/2))
} else {
  print("More than 5% of NO3 values are below MDL and")
  print("the values have remained unchanged")
  print(dat %>% count(NO3.mgl < .0007))
  print(sum(dat$NO3.mgl < .0007, na.rm = TRUE) / sum(!is.na(dat$NO3.mgl)) *100)}

#another round of imput for ammonia and tkn
#NH3 ---------------

dat$NH3.mgl <-
  ifelse(
    dat$NH3.mgl < .01,
    NA,
    dat$NH3.mgl
  )

impdat <-
  mice(
    data = dat[,c(11:13, 15, 17)],
    m = 5,
    maxit = 50,
    method = 'pmm',
    seed = 1234,
    printFlag = FALSE
  )

set.seed(1234)  # set seed for random imputed dataset selection
compdat <- complete(impdat, round(runif(1, 1, 5), 0))
# randomly select imputed dataset and combine with complete dataset

# inspect imputed data
# imputed datasets in magenta; observed in blue. are values plausible?
stripplot(impdat)
# do imputed points fall in likely area of distrib.?

dat[,c(11:13, 15, 17)] <- compdat[,1:ncol(compdat)]

# TKN-------------------
print(
  paste(
    sum(dat$TKN.mgl < 0, na.rm = TRUE) / sum(!is.na(dat$TKN.mgl)) *100,
    "% of TKN values are below 0"
  )
)


dat$TKN.mgl <-
  ifelse(
    dat$TKN.mgl < 0,
    NA,
    dat$TKN.mgl
  )

summary(dat$TKN.mgl)

impdat <-
  mice(
    data = dat[,c(11:13, 15, 17)],
    m = 5,
    maxit = 50,
    method = 'pmm',
    seed = 1234,
    printFlag = FALSE
  )

set.seed(1234)  # set seed for random imputed dataset selection
compdat <- complete(impdat, round(runif(1, 1, 5), 0))
# randomly select imputed dataset and combine with complete dataset

# inspect imputed data
# imputed datasets in magenta; observed in blue. are values plausible?
stripplot(impdat)
# do imputed points fall in likely area of distrib.?

dat[,c(11:13, 15, 17)] <- compdat[,1:ncol(compdat)]



# Calculate Loads ---------------------------------------------------------

#first we'll get rid of loads in the csv that don't reflect our adjusted mdl/2 values

dat <- dat[-c(22:31)]

# introduce a list with the field acreages and merge it with our dat so we can calculate loads, housekeeping

a <- as.data.frame(cbind(c("ARR1", "ARR2",  "PBR1", "PBR2", "STU1", "STU2", "MOS1", "MOS2", "DCDC1", "DCDC2",
                           "SCH1", "SCH2", "PRE1", "PRE2", "MUZ1", "MUZ2", "SIM1", "SIM2", "MUR1", "MUR2", "CAR1", "CAR2"),
          c(14.55, 15.48, 52.39, 23.89, 73.04, 34.24, 19.23, 19.26, 20.33, 20.73,
             10.38, 9.12, 36.62, 37.37, 19.29, 17.83, 14.29, 17.59, 19.04, 16.13, 22.36, 23.4)))


a <- as.data.frame(a) %>% 
  transform(V2= as.numeric(as.character(V2)))

dat <- merge(
  dat,
  a,
  by.x = "sampleid",
  by.y = "V1",
  all.x = TRUE,
  all.y = FALSE
)
names(dat)[names(dat) == "V2"] <- "acres"

# change acres where we moved fields
dat[dat$sampleid == "ARR1",]$acres <-
  ifelse(
    dat[dat$sampleid == "ARR1",]$Date >="2019-10-18",
    7.84, # New field size
    14.55   # old field size
  )

dat[dat$sampleid == "ARR2",]$acres <-
  ifelse(
    dat[dat$sampleid == "ARR2",]$Date >="2019-10-18",
    6.94, # New field size
    15.48   # old field size
  )

dat[dat$sampleid == "STU1",]$acres <-
  ifelse(
    dat[dat$sampleid == "STU1",]$Date >="2019-12-5",
    29.97, # New field size
    73.04   # old field size
  )

dat[dat$sampleid == "STU2",]$acres <-
  ifelse(
    dat[dat$sampleid == "STU2",]$Date >="2019-12-5",
    15.85, # New field size
    34.24   # old field size
  )

dat <- dat[c(2,1,3:29)]
#calculate loads
dat$in.runoff.ac <- dat$acft.discharge * 12 / dat$acres
dat$TSS.kg <- dat$TSS.mgl * dat$event.disch.l / 1000000
dat$TN.kg <- dat$TN.mgl * dat$event.disch.l / 1000000
dat$TIP.kg <- dat$TIP.mgl * dat$event.disch.l / 1000000
dat$TSS.lb <- dat$TSS.kg * 2.2046
dat$TN.lb <- dat$TN.kg * 2.2046
dat$TIP.lb <- dat$TIP.kg * 2.2046
dat$TSS.lbac <- dat$TSS.lb / dat$acres
dat$TN.lbac <- dat$TN.lb / dat$acres
dat$TIP.lbac <- dat$TIP.lb / dat$acres

#rearrange columns
dat <- dat[c(1:21,30:39,22:29)]

# PAIR DATA ---------------------------------------------------------------

# This is function I wrote to make things a little easier to read.
# You can open it to read a little bit more. 
dupedates <- function(data) {
  filter(data,
         duplicated(Date)|duplicated(Date, fromLast = TRUE) == TRUE
  )
  # filters and keeps where dates are duplicated
  # "== TRUE" means that there are two identical dates in the Date column
}                                                               

# ..Subsetting and keeping the pairs
  # We have to subset each farm and store them as a seperate data frames.                         
  # Then we join them afterwards and clean up the enverionment                                  
  # We need to use these 3 functions I found on the internet so I'm not 
  # sure how they work      

#....Generic form
'%=%' = function(l, r, ...) UseMethod('%=%')                                                    

#....Binary Operator
'%=%.lbunch' = function(l, r, ...) {
  Envir = as.environment(-1)
  
  if (length(r) > length(l))
    warning("RHS has more args than LHS. Only first", length(l), "used.")
  
  if (length(l) > length(r))  {
    warning("LHS has more args than RHS. RHS will be repeated.")
    r <- extendToMatch(r, l)
  }
  
  for (II in 1:length(l)) {
    do.call('<-', list(l[[II]], r[[II]]), envir=Envir)
  }
} 

#....Grouping the left hand side
g = function(...) {
  List = as.list(substitute(list(...)))[-1L]
  class(List) = 'lbunch'
  return(List)
}

#....Assigning and finding duplicate dates
g(arr, car, dcdc, mos, mur, muz, sch, sim, stu, pbr, pre) %=% list(
  subset(dat, farm == "ARR") %>% 
    dupedates(),
  subset(dat, farm == "CAR") %>% 
    dupedates(),
  subset(dat, farm == "DCDC") %>% 
    dupedates(),
  subset(dat, farm == "MOS") %>% 
    dupedates(),
  subset(dat, farm == "MUR") %>% 
    dupedates(),
  subset(dat, farm == "MUZ") %>% 
    dupedates(),
  subset(dat, farm == "SCH") %>% 
    dupedates(),
  subset(dat, farm == "SIM") %>% 
    dupedates(),
  subset(dat, farm == "STU") %>% 
    dupedates(),
  subset(dat, farm == "PBR") %>% 
    dupedates(),
  subset(dat, farm == "PRE") %>% 
    dupedates()
)


#..Combining all of the these into one data frame and cleaning the environment

dat.pair <- rbind(arr, car, dcdc, mos, mur, muz, sch, sim, stu, pbr, pre)
  # Combines all the individual farm data frames into
rm(arr, car, dcdc, mos, mur, muz, sch, sim, stu, pbr, pre)
  # removes the farm data frames


dat.pair <- dat.pair[ order(dat.pair$Date, dat.pair$sampleid),]
head(dat.pair[,c(1:5,9:11)])


# REDUCTION POTENTIAL calculations -----------------------------------------------------

# Now we use the paired dataframe to create this reduction potential dataframe.

# Here we create the new data frame and change the format of the paired data
# from a "long" format into a "wide" format so instead of every sample having
# its own row, the sample pairs each have their own row. This currently only 
# looks at concentration data but that can easily be changed.
dat.reduct <- dat.pair %>% 
  pivot_wider(
    # this is how we want to identify how to split up our rows
    id_cols = c(Date, farm, season),
    # this is the factor we split our columns by
    names_from = treatment,
    names_sep = ".",
    # this is the data we want to see in our new dataframe
    # this is also where you can add or remove columns you want to see.
    values_from = names(dat.pair[,c(9:31)])
  ) %>% 
  as.data.frame()

# .. calculating reductions
# Here is where we actually claculate the reduction values and percentages.
# It goes column by column and first calculates the difference between the
# FBM and CCMT (FBM minus CCMT) and then calculates the percent the CCMT reduced 
# the nutrient (difference divided by FBM). 

# A NEGATIVE DIFFERENCE INDICATES AN INCREASE FROM FBM TO CCMT. 
# A NEGATIVE PERCENTAGE INDICATES AN INCREASE FROM FBM TO CCMT. 
# You can change this below if you'd like.

# The difference columns are named nutrient_dif, and the percentages columns 
# are named nutrient_per. The layout of the reduction columns will always be
# nutrient1_dif, nutrient1_per, nutrient2_dif, nutrient2_per... 

for (i in 4:48) {
  if (i %% 2 == 0) {  
    # This is the difference code.
    
    # This creates a new column at the end of the dataframe and gives it the value 
    # of the FBM - CCMT. To switch the order as noted above use the commented out 
    # line instead. This will make the percentage follow the same pattern.
    dat.reduct[,ncol(dat.reduct) + 1] <- dat.reduct[,i + 1] - dat.reduct[,i]
#   dat.reduct[,ncol(dat.reduct) + 1] <- dat.reduct[,i] - dat.reduct[,i + 1]
      # the way I make it always a new column at the end is the function
      # ncol() returns the number of columns in the dataframe, so adding 1 to this
      # number adds a column to the end. But afterwards, I have to use ncol() again
    # This names the column from where the data was taken. 
    # For names(dat.reduct)[i] the [i] has to be outside of the parenthesies
    names(dat.reduct)[ncol(dat.reduct)] <- names(dat.reduct)[i]
    # the name is changed by replaceing "CCMT" with "dif" in the column title
    names(dat.reduct)[ncol(dat.reduct)] <- gsub(
      pattern = "CCMT",
      replacement = "dif",
      x = names(dat.reduct)[ncol(dat.reduct)]
    )
    
    # Then it moves to find the percent reduction in the same way as above. 
    # If you did change the difference code so that a negative difference means a 
    # reduction, you don't need to change this line.
    dat.reduct[,ncol(dat.reduct) + 1] <- dat.reduct[,ncol(dat.reduct)] / dat.reduct[,i + 1] * 100
    names(dat.reduct)[ncol(dat.reduct)] <- names(dat.reduct)[i]
    names(dat.reduct)[ncol(dat.reduct)] <- gsub(
      pattern = "CCMT",
      replacement = "per",
      x = names(dat.reduct)[ncol(dat.reduct)]
    )
  }
}


# summary stuff -----------------------------------------------------------

# total number of samples
length(dat$sampleid)

# number of samples by type
dat$sample.type[dat$sample.type == "EOF"] %>% length()
dat$sample.type[dat$sample.type == "INST"] %>% length()

# number of samples by treatment
length(dat$treatment[!is.na(dat$treatment) & dat$treatment == "CCMT"])
      #length(dat$treatment[dat$sample.type =="EOF" & dat$treatment == "CCMT"])
length(dat$treatment[!is.na(dat$treatment) & dat$treatment == "FBM"])

# number of EOF samples by season
length(dat$season[dat$sample.type == "EOF" & dat$season == "cover"])
length(dat$season[dat$sample.type == "EOF" & dat$season == "cash"])

# number of paired sample events
length(dat.pair$sampleid)/2
length(dat.pair$sampleid[dat.pair$season == "cover"])/2
length(dat.pair$sampleid[dat.pair$season == "cash"])/2

### statistical summary of analytes in dat
favstats(dat[,9]~ treatment, data = dat[dat$sample.type == "EOF",])

#favstats for the ten variables below

for(i in c(9,10,11,13,14,16,22,29,30,31)){
  
   print(favstats(dat[i] ~ treatment, data = dat[dat$sample.type== "EOF",]))
            
  } 

b <- for(i in c(9,10,11,13,14,16,22,29,30,31)){
  
  print(favstats(dat[,i][dat$sample.type== "EOF" & dat$treatment== "CCMT"]))
  
} 

### makes two data frame then combines to make a table of the favstats for all analytes (needs names somehow though)
exampleDFccmt <- data.frame()
for(i in c(9,10,11,13,14,16,22,29,30,31)){
  
  x = favstats(dat[,i][dat$sample.type== "EOF" & dat$treatment== "CCMT"])
  exampleDFccmt <- rbind(exampleDFccmt, x)
}

exampleDFfbm <- data.frame()
for(i in c(9,10,11,13,14,16,22,29,30,31)){
  
  x = favstats(dat[,i][dat$sample.type== "EOF" & dat$treatment== "FBM"])
  exampleDFfbm <- rbind(exampleDFfbm, x)
}

ex <- do.call("rbind", Map("rbind", split(exampleDFccmt, 1:nrow(exampleDFccmt)), split(exampleDFfbm, 1:nrow(exampleDFfbm))))

write.csv(ex, "summary.csv")

rm(exampleDF)


  
### makes the table into a graphical object which you can then plot
c <- (favstats(TSS.mgl ~ treatment, data = dat[dat$sample.type== "EOF",]))
d <- tableGrob(c)

(p1 + p2) / d


for(i in c(9,10,11,13,14,16,22,29,30,31)){
  print(names(dat[i]))
   print( favstats(dat[,i] ~ treatment, data = dat[dat$sample.type== "EOF",])
  )
}

c <- do.call(rbind, dfapply(dat, favstats, select= c(9:35), groups = "treatment")) %>% write.table


### builds a list of all the returns of the for loop
df <- list()

for(i in 9:19) {
  df <- mosaic::favstats(dat[,i])
  df1$i <- names(dat[i])
  df[[i -8]] <- df1
}
df <- do.call(
  rbind,
  df
)
df <- df[,c(10,1:9)]
as_tibble(df)

#and treatment as a factor
df <- list()

for(i in 9:19) {
  df1 <- mosaic::favstats(dat[,i] ~ dat$treatment)
  df1$i <- names(dat[i])
  df[[i -8]] <- df1
}
df <- do.call(
  rbind,
  df
)
df <- df[,c(11,1:10)]

### visualizing variable distributions --------------------------------------

#TSS.mgl
p1 <- ggdensity(dat[dat$sample.type == "EOF", "TSS.mgl"],
                xlab = names(dat["TSS.mgl"]))

p2 <- ggplot(dat[dat$sample.type == "EOF",], 
             aes(y= TSS.mgl, x= "All year", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 3600))+
  stat_compare_means()


p4 <- ggplot(dat[dat$sample.type == "EOF" & dat$season == "cover",], 
             aes(y= TSS.mgl, x= "Cover season", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 3600))+
  stat_compare_means()

p3 <- ggplot(dat.pair[dat$sample.type == "EOF",], 
             aes(y= TSS.mgl, x= " paired All year", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 3600))+
  stat_compare_means()

p5 <- ggplot(dat[dat$sample.type == "EOF",], aes(Date, TSS.mgl))+
  geom_point()+
  stat_smooth(span = 0.3)+
  geom_hline(data = dat, yintercept = median(dat$TSS.mgl, na.rm= T))+
  scale_x_date(date_labels = "%m-%Y", date_breaks= "2 months")

grid.arrange(p1,p2, p3, p4, p5, layout_matrix= rbind(c(1,2,3,4), c(5,5,5,5)),
             top= textGrob("TSS.mgl", gp= gpar(fontsize= 20, font=4)))

(p1+p2+p3+p4)/ p5 + plot_layout(design = layout)

layout <- "
AABBCCDD
AABBCCDD
EEEEEEEE
"



#TUR.ntu 
p1 <- ggdensity(dat[dat$sample.type == "EOF", "TUR.ntu"],
                xlab = names(dat["TUR.ntu"]))

p2 <- ggplot(dat[dat$sample.type == "EOF",], 
             aes(y= TUR.ntu, x= "All year", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 3600))+
  stat_compare_means()


p4 <- ggplot(dat[dat$sample.type == "EOF" & dat$season == "cover",], 
             aes(y= TUR.ntu, x= "Cover season", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 3600))+
  stat_compare_means()

p3 <- ggplot(dat.pair[dat$sample.type == "EOF",], 
             aes(y= TUR.ntu, x= " paired All year", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 3600))+
  stat_compare_means()

p5 <- ggplot(dat[dat$sample.type == "EOF",], aes(Date, TUR.ntu))+
  geom_point()+
  stat_smooth(span = 0.3)+
  geom_hline(data = dat, yintercept = median(dat$TUR.ntu, na.rm= T))+
  scale_x_date(date_labels = "%m-%Y", date_breaks= "2 months")

grid.arrange(p1,p2, p3, p4, p5, layout_matrix= rbind(c(1,2,3,4), c(5,5,5,5)),
             top= textGrob("TUR.ntu", gp= gpar(fontsize= 20, font=4)))

#NO3.NO2.mgl

p1 <- ggdensity(dat[dat$sample.type == "EOF", "NO3.NO2.mgl"],
                xlab = names(dat["NO3.NO2.mgl"]))

p2 <- ggplot(dat[dat$sample.type == "EOF",], 
             aes(y= NO3.NO2.mgl, x= "All year", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 7.5))+
  stat_compare_means()


p4 <- ggplot(dat[dat$sample.type == "EOF" & dat$season == "cover",], 
             aes(y= NO3.NO2.mgl, x= "Cover season", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 7.5))+
  stat_compare_means()

p3 <- ggplot(dat.pair[dat$sample.type == "EOF",], 
             aes(y= NO3.NO2.mgl, x= " paired All year", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 7.5))+
  stat_compare_means()

p5 <- ggplot(dat[dat$sample.type == "EOF",], aes(Date, NO3.NO2.mgl))+
  geom_point()+
  stat_smooth(span = 0.3)+
  geom_hline(data = dat, yintercept = median(dat$NO3.NO2.mgl, na.rm= T))+
  scale_x_date(date_labels = "%m-%Y", date_breaks= "2 months")

grid.arrange(p1,p2, p3, p4, p5, layout_matrix= rbind(c(1,2,3,4), c(5,5,5,5)),
             top= textGrob("NO3.NO2.mgl", gp= gpar(fontsize= 20, font=4)))
#TN.mgl
p1 <- ggdensity(dat[dat$sample.type == "EOF", "TN.mgl"],
                xlab = names(dat["TN.mgl"]))

p2 <- ggplot(dat[dat$sample.type == "EOF",], 
             aes(y= TN.mgl, x= "All year", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 7.5))+
  stat_compare_means()


p4 <- ggplot(dat[dat$sample.type == "EOF" & dat$season == "cover",], 
             aes(y= TN.mgl, x= "Cover season", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 7.5))+
  stat_compare_means()

p3 <- ggplot(dat.pair[dat$sample.type == "EOF",], 
             aes(y= TN.mgl, x= " paired All year", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 7.5))+
  stat_compare_means()

p5 <- ggplot(dat[dat$sample.type == "EOF",], aes(Date, TN.mgl))+
  geom_point()+
  stat_smooth(span = 0.3)+
  geom_hline(data = dat, yintercept = median(dat$TN.mgl, na.rm= T))+
  scale_x_date(date_labels = "%m-%Y", date_breaks= "2 months")+
  

grid.arrange(p1,p2, p3, p4, p5, layout_matrix= rbind(c(1,2,3,4), c(5,5,5,5)),
             top= textGrob("TN.mgl", gp= gpar(fontsize= 20, font=4)))

#TIP.mgl
p1 <- ggdensity(dat[dat$sample.type == "EOF", "TIP.mgl"],
                xlab = names(dat["TIP.mgl"]))

p2 <- ggplot(dat[dat$sample.type == "EOF",], 
             aes(y= TIP.mgl, x= "All year", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 7.5))+
  stat_compare_means()


p4 <- ggplot(dat[dat$sample.type == "EOF" & dat$season == "cover",], 
             aes(y= TIP.mgl, x= "Cover season", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 7.5))+
  stat_compare_means()

p3 <- ggplot(dat.pair[dat$sample.type == "EOF",], 
             aes(y= TIP.mgl, x= " paired All year", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 7.5))+
  stat_compare_means()

p5 <- ggplot(dat[dat$sample.type == "EOF",], aes(Date, TIP.mgl))+
  geom_point()+
  stat_smooth(span = 0.3)+
  geom_hline(data = dat, yintercept = median(dat$TIP.mgl, na.rm= T))+
  scale_x_date(date_labels = "%m-%Y", date_breaks= "2 months")

grid.arrange(p1,p2, p3, p4, p5, layout_matrix= rbind(c(1,2,3,4), c(5,5,5,5)),
             top= textGrob("TIP.mgl", gp= gpar(fontsize= 20, font=4)))
#OrthoP
p1 <- ggdensity(dat[dat$sample.type == "EOF", "OrthoP.mgl"],
                xlab = names(dat["OrthoP.mgl"]))

p2 <- ggplot(dat[dat$sample.type == "EOF",], 
             aes(y= OrthoP.mgl, x= "All year", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0,.75))+
  stat_compare_means()


p4 <- ggplot(dat[dat$sample.type == "EOF" & dat$season == "cover",], 
             aes(y= OrthoP.mgl, x= "Cover season", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, .75))+
  stat_compare_means()

p3 <- ggplot(dat.pair[dat$sample.type == "EOF",], 
             aes(y= OrthoP.mgl, x= " paired All year", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, .75))+
  stat_compare_means()

p5 <- ggplot(dat[dat$sample.type == "EOF",], aes(Date, OrthoP.mgl))+
  geom_point()+
  stat_smooth(span = 0.3)+
  geom_hline(data = dat, yintercept = median(dat$OrthoP.mgl, na.rm= T))+
  scale_x_date(date_labels = "%m-%Y", date_breaks= "2 months")

grid.arrange(p1,p2, p3, p4, p5, layout_matrix= rbind(c(1,2,3,4), c(5,5,5,5)),
             top= textGrob("OrthoP.mgl", gp= gpar(fontsize= 20, font=4)))

#in.runoff.ac
p1 <- ggdensity(dat[dat$sample.type == "EOF", "in.runoff.ac"],
                xlab = names(dat["in.runoff.ac"]))

p2 <- ggplot(dat[dat$sample.type == "EOF",], 
             aes(y= in.runoff.ac, x= "All year", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 2))+
  stat_compare_means()


p4 <- ggplot(dat[dat$sample.type == "EOF" & dat$season == "cover",], 
             aes(y= in.runoff.ac, x= "Cover season", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 2))+
  stat_compare_means()

p3 <- ggplot(dat.pair[dat$sample.type == "EOF",], 
             aes(y= in.runoff.ac, x= " paired All year", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 2))+
  stat_compare_means()

p5 <- ggplot(dat[dat$sample.type == "EOF",], aes(Date, in.runoff.ac))+
  geom_point()+
  stat_smooth(span = 0.3)+
  geom_hline(data = dat, yintercept = median(dat$in.runoff.ac, na.rm= T))+
  scale_x_date(date_labels = "%m-%Y", date_breaks= "2 months")

grid.arrange(p1,p2, p3, p4, p5, layout_matrix= rbind(c(1,2,3,4), c(5,5,5,5)),
             top= textGrob("in.runoff.ac", gp= gpar(fontsize= 20, font=4)))

#TSS.lbac
p1 <- ggdensity(dat[dat$sample.type == "EOF", "TSS.lbac"],
                xlab = names(dat["TSS.lbac"]))

p2 <- ggplot(dat[dat$sample.type == "EOF",], 
             aes(y= TSS.lbac, x= "All year", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 7.5))+
  stat_compare_means()


p4 <- ggplot(dat[dat$sample.type == "EOF" & dat$season == "cover",], 
             aes(y= TSS.lbac, x= "Cover season", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 7.5))+
  stat_compare_means()

p3 <- ggplot(dat.pair[dat$sample.type == "EOF",], 
             aes(y= TSS.lbac, x= " paired All year", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 7.5))+
  stat_compare_means()

p5 <- ggplot(dat[dat$sample.type == "EOF",], aes(Date, TSS.lbac))+
  geom_point()+
  stat_smooth(span = 0.3)+
  geom_hline(data = dat, yintercept = median(dat$TSS.lbac, na.rm= T))+
  scale_x_date(date_labels = "%m-%Y", date_breaks= "2 months")

grid.arrange(p1,p2, p3, p4, p5, layout_matrix= rbind(c(1,2,3,4), c(5,5,5,5)),
             top= textGrob("TSS.lbac", gp= gpar(fontsize= 20, font=4)))

#TN.lbac
p1 <- ggdensity(dat[dat$sample.type == "EOF", "TN.lbac"],
                xlab = names(dat["TN.lbac"]))

p2 <- ggplot(dat[dat$sample.type == "EOF",], 
             aes(y= TN.lbac, x= "All year", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 1))+
  stat_compare_means()


p4 <- ggplot(dat[dat$sample.type == "EOF" & dat$season == "cover",], 
             aes(y= TN.lbac, x= "Cover season", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 1))+
  stat_compare_means()

p3 <- ggplot(dat.pair[dat$sample.type == "EOF",], 
             aes(y= TN.lbac, x= " paired All year", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 1))+
  stat_compare_means()

p5 <- ggplot(dat[dat$sample.type == "EOF",], aes(Date, TN.lbac))+
  geom_point()+
  stat_smooth(span = 0.3)+
  geom_hline(data = dat, yintercept = median(dat$TN.lbac, na.rm= T))+
  scale_x_date(date_labels = "%m-%Y", date_breaks= "2 months")

grid.arrange(p1,p2, p3, p4, p5, layout_matrix= rbind(c(1,2,3,4), c(5,5,5,5)),
             top= textGrob("TN.lbac", gp= gpar(fontsize= 20, font=4)))

#TIP.lb.ac
p1 <- ggdensity(dat[dat$sample.type == "EOF", "TIP.lbac"],
                xlab = names(dat["TIP.lbac"]))

p2 <- ggplot(dat[dat$sample.type == "EOF",], 
             aes(y= TIP.lbac, x= "All year", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 1.5))+
  stat_compare_means()


p4 <- ggplot(dat[dat$sample.type == "EOF" & dat$season == "cover",], 
             aes(y= TIP.lbac, x= "Cover season", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 1.5))+
  stat_compare_means()

p3 <- ggplot(dat.pair[dat$sample.type == "EOF",], 
             aes(y= TIP.lbac, x= " paired All year", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 1.5))+
  stat_compare_means()

p5 <- ggplot(dat[dat$sample.type == "EOF",], aes(Date, TIP.lbac))+
  geom_point()+
  stat_smooth(span = 0.3)+
  geom_hline(data = dat, yintercept = median(dat$TIP.lbac, na.rm= T))+
  scale_x_date(date_labels = "%m-%Y", date_breaks= "2 months")

grid.arrange(p1,p2, p3, p4, p5, layout_matrix= rbind(c(1,2,3,4), c(5,5,5,5)),
             top= textGrob("TIP.lbac", gp= gpar(fontsize= 20, font=4)))












#boxplot
boxplot(dat$TSS.mgl~ dat$treatment, outline= F)

#density plots for analytes

for (a in c(9, 10, 11, 13, 14, 16, 22, 30)) {

     print(ggdensity(dat[dat$sample.type == "EOF", a ], 
                     xlab = names(dat[a])))
}



#quantile plots

x <- dat$TSS.mgl[dat$treatment == "CCMT"]
y <- dat$TSS.mgl[dat$treatment == "FBM"]
qqplot(x,y,
       xlab = "CCMT",
       ylab = "FBM"
      
      )
  qqline(y)







  
  
#time vs mean plot



ggplot(dat[dat$sample.type == "EOF",], aes(Date, TSS.mgl))+
  geom_point()+
  stat_smooth(span = 0.3)+
  geom_hline(data = dat, yintercept = median(dat$TSS.mgl, na.rm= T))+
  scale_x_date(date_labels = "%m-%Y", date_breaks= "2 months")

ggplot(dat[dat$sample.type == "EOF",], aes(Date, TUR.ntu))+
  geom_point()+
  stat_smooth(span = 0.3)+
  geom_hline(data = dat, yintercept = median(dat$TUR.ntu, na.rm= T))+
  scale_x_date(date_labels = "%m-%Y", date_breaks= "2 months")

ggplot(dat[dat$sample.type == "EOF",], aes(Date, NO3.NO2.mgl))+
  geom_point()+
  stat_smooth(span = 0.3)+
  geom_hline(data = dat, yintercept = median(dat$NO3.NO2.mgl, na.rm= T))+
  scale_x_date(date_labels = "%m-%Y", date_breaks= "2 months")
multiplot(a, b, c, cols = 1)

ggplot(dat[dat$sample.type == "EOF",], aes(Date, TN.mgl))+
  geom_point()+
  stat_smooth(span = 0.3)+
  geom_hline(data = dat, yintercept = median(dat$TN.mgl, na.rm= T))+
  scale_x_date(date_labels = "%m-%Y", date_breaks= "2 months")

ggplot(dat[dat$sample.type == "EOF",], aes(Date, TIP.mgl))+
  geom_point()+
  stat_smooth(span = 0.3)+
  geom_hline(data = dat, yintercept = median(dat$TIP.mgl, na.rm= T))+
  scale_x_date(date_labels = "%m-%Y", date_breaks= "2 months")

ggplot(dat[dat$sample.type == "EOF",], aes(Date, OrthoP.mgl))+
  geom_point()+
  stat_smooth(span = 0.3)+
  geom_hline(data = dat, yintercept = median(dat$OrthoP.mgl, na.rm= T))+
  scale_x_date(date_labels = "%m-%Y", date_breaks= "2 months")



#for (i in c(9,10,11,13,14,16)) {
#  
#print(
#  ggplot(dat[dat$sample.type == "EOF",], aes(Date, dat[,i]))+
#  geom_point()+
#  stat_smooth(span = 0.3)+
#  geom_hline(data = dat, yintercept = median(dat[,i], na.rm= T))+
#  scale_x_date(date_labels = "%m-%Y", date_breaks= "2 months")
#   )
#}



### Analyte plots across all sites

analyteplot <- function(analyte,Title,Ylab = "mg/L",Xlab,Ylimt){
  ggplot()+ 
    geom_boxplot(data = dat, 
                 aes(x = sampleid, 
                     y = analyte, 
                     fill= farm)
    )+
    ggtitle(Title)+
    ylab(Ylab)+
    xlab(Xlab)+
    # coord_cartesian(ylim = c())+
    coord_cartesian(ylim= c(0,Ylimt))+
    theme_bw()
  
}


analyteplot(dat$TSS.mgl,"Total suspended solids","mg/L","site",10000)

analyteplot(dat$TUR.ntu,"Turbidity","NTU","site",12000)

analyteplot(dat$NO3.NO2.mgl,"NO3.NO2.mgl","mg/L","site", 15)

analyteplot(dat$TN.mgl,"Total Nitrogen","mg/L","site",20)

analyteplot(dat$TIP.mgl,"Total inorganic Phosphorus","mg/L","site",25)

analyteplot(dat$OrthoP.mgl,"Ortho Phosphorus","mg/L","site",2)

analyteplot(dat$in.runoff.ac,"Inches of runoff per acre","in","site",3)

analyteplot(dat$TSS.lbac,"Soil losses","pounds/ac","site",1000)

analyteplot(dat$TN.lbac,"Nitrogen discharge","pounds/acre","site",3)

analyteplot(dat$TIP.lbac,"Phosphorus discharge","pounds/acre","site",4.5)







# Average reduction potentials  -------------------------------------------

# percent reduction of TSS exploration
fivenum(dat.reduct$TSS.mgl.per, na.rm=T)
mean(dat.reduct$TSS.mgl.per, trim = 0, na.rm = T)
sd(dat.reduct$TSS.mgl.per, na.rm = T)
median(dat.reduct$TSS.mgl.per, na.rm= T)
boxplot(dat.reduct$TSS.mgl.per)
  ## a few outliers can completely skew this data. below plots the data w/o the outliers
boxplot(dat.reduct$TSS.mgl.per, outline= FALSE)



##during the entire season

par(mfrow = c(3,3))

for ( i in c(51,53,55,59,61,65,77, 93, 95)){
  cat(c(names(dat.reduct[i]), "five number summary \n"))
  print(fivenum(dat.reduct[, i], na.rm=T))
  print(paste("Mean  =", mean(dat.reduct[, i], trim = 0, na.rm = T)))
  print(paste("standard deviation =", sd(dat.reduct[, i], na.rm = T)))
  
  #boxplot(dat.reduct[,i], main = names(dat.reduct[i]))
  ## a few outliers completely skew this data. below plot the data w/o the outliers
  boxplot(dat.reduct[, i], 
          outline= FALSE,
          na.rm= T, 
          main = names(dat.reduct[i]), 
          ylab = "(+) represents reduction from treatment"
          )
}
title("Overall Reduction potential", outer = T)


## during the cover season
  for ( i in c(51,53,55,59,61,65,77, 93, 95)){
cat(c(names(dat.reduct[i]), "five number summary \n"))
 print(fivenum(dat.reduct[dat.reduct$season == "cover", i], na.rm=T))
 print(paste("Mean  =", mean(dat.reduct[dat.reduct$season == "cover", i], trim = 0, na.rm = T)))
 print(paste("standard deviation =", sd(dat.reduct[dat.reduct$season == "cover", i], na.rm = T)))

 #boxplot(dat.reduct[,i], main = names(dat.reduct[i]))
 ## a few outliers completely skew this data. below plot the data w/o the outliers
 boxplot(dat.reduct[dat.reduct$season == "cover", i],
         outline= FALSE, 
         na.rm= T, 
         main = names(dat.reduct[i]),
         ylab = "(+) represents reduction from treatment"
         )
 
}
title( "Reduction potential cover season", outer = T)

par(mfrow = c(1,1))

    #dev.off()

### makes 3x3 plot of reduction potentials
par(mfrow = c(3,3))
boxplot(dat.reduct[,51], outline= FALSE, na.rm= T, main = names(dat.reduct[51]) , ylab = "(+) represents reduction from treatment")
boxplot(dat.reduct[,53], outline= FALSE, na.rm= T, main = names(dat.reduct[53]) , ylab = "(+) represents reduction from treatment")
boxplot(dat.reduct[,55], outline= FALSE, na.rm= T, main = names(dat.reduct[55]) , ylab = "(+) represents reduction from treatment")
boxplot(dat.reduct[,59], outline= FALSE, na.rm= T, main = names(dat.reduct[59]) , ylab = "(+) represents reduction from treatment")
boxplot(dat.reduct[,61], outline= FALSE, na.rm= T, main = names(dat.reduct[61]) , ylab = "(+) represents reduction from treatment")
boxplot(dat.reduct[,65], outline= FALSE, na.rm= T, main = names(dat.reduct[65]) , ylab = "(+) represents reduction from treatment")
boxplot(dat.reduct[,77], outline= FALSE, na.rm= T, main = names(dat.reduct[77]) , ylab = "(+) represents reduction from treatment")
boxplot(dat.reduct[,93], outline= FALSE, na.rm= T, main = names(dat.reduct[93]) , ylab = "(+) represents reduction from treatment")
boxplot(dat.reduct[,95], outline= FALSE, na.rm= T, main = names(dat.reduct[95]) , ylab = "(+) represents reduction from treatment")
title("Percent reduction", outer = T)




par(mfrow = c(3,3))
for (i in c(51,53,55,57,59,61,63,65,67)){

rd <- boxplot(dat.reduct[,i], outline= FALSE, na.rm= T, main = names(dat.reduct[i]) , ylab = "(+) represents reduction from treatment")
rdmean <- with(rd, stats[3,], names, pos= 0)
text(1, rdmean, paste(round(rdmean,0)),  adj = 10)

}


dev.off()

# scratch -----------------------------------------------------------------
data <- dat$NO3.NO2.mgl

logccmt <- log10(dat$NO3.NO2.mgl[dat$treatment == "CCMT"])
logfbm <- log10(dat$NO3.NO2.mgl[dat$treatment == "FBM"])

t.test(logccmt, logfbm)

plot(density(logccmt, na.rm= T))
lines(density(logfbm, na.rm= T))

ggdensity(logccmt)+ stat_overlay_normal_density(color = "red")



a <- cbind(c("ARR1", "ARR2",  "PBR1", "PBR2", "STU1", "STU2", "MOS1", "MOS2", "DCDC1",
             "DCDC2", "SCH1", "SCH2", "PRE1", "PRE2", "MUZ1", "MUZ2", "SIM1", "SIM2", "MUR1", "MUR2", "CAR1", "CAR2"),
           c(21.24, 22.35, 76.2, 34.72, 105.35, 49.64, 27.83, 27.63, 28.78, 31.52,
             10.38, 9.12, 36.62, 37.37, 19.29, 17.83, 14.29, 17.59, 19.04, 16.13, 22.36, 23.4))
a <- as.factor(a[,1])

b <- cbind(c("ARR1","ARR2","STU1","STU2"), c(11.41, 10.12, 30.2, 23.08))

datex <- merge(
  dat,
  a,
  by.x = "sampleid",
  by.y = "V1",
  all.x = TRUE,
  all.y = FALSE
)
names(datex)[names(datex) == "V2"] <- "acres"

datex <- datex[,c(2,1,3:37)]


datex <- datex %>% 
  mutate( acres = as.numeric(acres)) %>% 
  mutate(acres = replace(acres, sampleid== "ARR1", 21.24))


datex <- datex[datex$sampleid == "ARR1"] %>% 
  mutate( acres = as.numeric(acres)) %>% 
  mutate(acres = replace(acres, Date > "2019-10-18", 11.41))


dat <- mutate(subset(dat,
              subset= dat$Date > "2019-10-18" & dat$sampleid == "ARR1", 
              !is.null(dat$acres), 
              11.41))



dat[dat$sampleid == "ARR1",]$acres <-
  ifelse(
    dat[dat$sampleid == "ARR1",]$Date >="2019-10-18",
    100, # New field size
    1    # old field size
  )



#qqplots
ggplot(dat[dat$sample.type== "EOF",])+
  stat_qq(aes(sample= dat$TSS.mgl))+
  stat_qq_line()
##### composite

p1 <- ggplot(dat[dat$sample.type == "EOF",], 
             aes(y= TSS.mgl, x= "All year", colour= treatment))+
  geom_boxplot(outlier.shape = NA )+
  coord_cartesian(ylim= c(0,4500))+
  stat_compare_means(label.y = 4000, lab)

p2 <- ggplot(dat.pair[dat$sample.type == "EOF",], 
             aes(y= TSS.mgl, x= " paired All year", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  coord_cartesian(ylim= c(0,4500))+
  stat_compare_means(label.y = 4000)


p3 <- ggplot(dat[dat$sample.type == "EOF" & dat$season == "cover",], 
             aes(y= TSS.mgl, x= "Cover season", colour= treatment))+
  geom_boxplot(outlier.shape = NA)+
  coord_cartesian(ylim= c(0,4500))+
  stat_compare_means(label.y = 4000)


p4 <- ggplot(dat[dat$sample.type == "EOF",], aes(Date, TSS.mgl))+
  geom_point()+
  stat_smooth(span = 0.3)+
  geom_hline(data = dat, yintercept = median(dat$TSS.mgl, na.rm= T))+
  scale_x_date(date_labels = "%m-%Y", date_breaks= "2 months")+
  rotate_x_text(angle = 45)


(p1 + p2 + p3) / p4

grid.arrange(p1,p2, p3, p4, layout_matrix= rbind(c(1,2,3), c(4,4,4)),
             top= textGrob("TSS.mgl", gp= gpar(fontsize= 20, font=4)))



#### try to use map to run ggplot boxplot in a loop
expl <- dat.reduct %>% 
  select(TSS.mgl.per, TUR.ntu.per, NO3.NO2.mgl.per,TKN.mgl.per,TN.mgl.per,TIP.mgl.per,NH3.mgl.per,OrthoP.mgl.per,NOx.mgl.per,in.runoff.ac.per) %>% names()
  


#names(dat.pair)[,c(51,53,55,57,59,61,63,65,67,77)]

expl = set_names(expl)

box_fun = function(y) {
  #ylimt = boxplot.stats(dat.reduct$y)$stats[c(1,5)]
  ylimt <- boxplot.stats(.data[[y]])$stats[c(1,5)]
  ggplot(dat.reduct, aes(y = .data[[y]] )) +
    geom_boxplot() +
    theme_bw() +
    labs(y = y)+
    coord_cartesian(ylim = c(ylimt))
  #boxplot(dat.reduct[[y]],
   #       outline = FALSE,
    #      main= names(.data[[y]])
     #     )
  
  }


test_plots <- map(expl, ~box_fun(.x))

test_plots


### summarize sample count by month
dat$month <- month(dat$Date)
dat$year <- year(dat$Date)
table(dat$month, dat$year)

### compare sonde TUR to lab TUR

ggplot(dat)+
  geom_line(aes(dat$Date, dat$TUR.ntu, alpha= 0.2))+
  geom_line(aes(dat$Date, dat$turbidity.sonde.ntu), color= "red")

p1 <- ggplot(dat[dat$sample.type == "INST",])+
  geom_boxplot(mapping = aes(y= TUR.ntu))
  
p2 <- ggplot(dat)+
  geom_boxplot(mapping = aes(y= dat$turbidity.sonde.ntu))

p1 | p2
p1
p2


## scatterplot of nutrients
ggplot(dat[dat$sample.type == "EOF",])+
  geom_jitter(aes(x= treatment, y= NO3.NO2.mgl, color= treatment))

ggplot(dat[dat$sample.type == "EOF",],
       aes(y= TN.mgl, x= treatment))+
  geom_boxplot()+
  stat_compare_means()+
  theme_pubr()+
  coord_cartesian(ylim = c(0,25))

ggplot(dat.pair[dat.pair$sample.type== "EOF",],aes(log(TN.mgl), fill= treatment))+
  geom_density(position = "stack")+
  coord_cartesian(xlim = c(0,5))

##lm
dat$logtur <- log10(dat$TUR.ntu)
dat$logtss <- log10(dat$TSS.mgl)

ggdensity(dat$logtss, data = dat[dat$sample.type== "EOF",])

ggplot(dat[dat$sample.type== "EOF",])+
  geom_density(aes(x=logtss))
 
ggplot(dat[dat$sample.type== "EOF",],
       aes(x= logtss, y= logtur))+  
  geom_point()+
  geom_smooth(method = "lm")

tnlm <- lm(logtur ~ logtss,
           data = dat[dat$sample.type == "EOF",])

autoplot(tnlm)
anova(tnlm)
summary(tnlm) 