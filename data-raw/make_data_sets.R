

##################################### load packages ######################################

library(magrittr)
library(dplyr)
library(lowliner)
library(stringr)
library(devtools)

oldwd <- setwd("data-raw/")


################################ Read and transform data #################################


files <- dir()
files %<>% keep(~ substring(., str_length(.) - 3, ) == ".txt")

dfs <- map(files, . %>% read.csv(check.names = FALSE, row.names = 1)) 
   

####################### Save datasets as required for lazy loading #######################

indata <- dfs[[1]]
inpop2 <- dfs[[2]]
inpop1 <- dfs[[3]]

setwd(oldwd)
use_data(indata, inpop1, inpop2)


