# Data Prep

# clear environment
rm(list = ls())
gc()

#######################
# Load libraries
#######################
library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)
library(RPostgreSQL)
source("Functions/helper_functions.R")

dir_out <- "Output"

#######################
# Load data
#######################

family <- read.csv( "Data/W_Susquehanna_Peterson_3pass_Network.csv", stringsAsFactors = FALSE)
colnames(family)[1] = "child_name"
family = cbind( family, "child_b"=1:nrow(family) )

df_fish <- read.csv("Data/West_Susquehanna_3pass.csv", stringsAsFactors = FALSE)
df_fish_MR <- read.csv("Data/Lincoln_Petersen_MR.csv", stringsAsFactors = FALSE)


# combine with 3-pass data
df_fish$method <- "removal"
df_fish_MR$method <- "MR"
df_fish <- dplyr::bind_rows(df_fish, df_fish_MR)


#---------- check that equal effort for 1st pass of both methods ---------



#######################
# Prep Trout Data
#######################

obs_sites_network <- family %>%
  dplyr::select(child_name) %>%
  distinct() %>%
  dplyr::filter(!(grepl("N_", child_name)))

# organize
df_fish <- df_fish %>%
  dplyr::mutate(site = GIS_Key,
                site_visit = paste0(site, "_", date),
                site_year = paste0(site, "_", year),
                pass = EffortNumber)

# assign stage class
df_bkt <- df_fish %>%
  dplyr::filter(species == "Brook Trout")

# Assign Stage
df_bkt_adult <- df_fish %>%
  dplyr::mutate(stage = ifelse(sizebin >= 100, "adult", "yoy")) %>%
  dplyr::filter(species == "Brook Trout" & stage == "adult") %>%
  dplyr::group_by(site_visit, site, year, date, pass, method) %>%
  dplyr::summarise(count = sum(catch))

df_bkt_yoy <- df_fish %>%
  dplyr::mutate(stage = ifelse(sizebin < 100, "yoy", "adult")) %>%
  dplyr::filter(species == "Brook Trout" & stage == "yoy") 
df_bkt_yoy <- df_bkt_yoy %>%
  dplyr::group_by(site_visit, site, year, date, pass, method) %>%
  dplyr::summarise(count = sum(catch))

# expand to 3-pass for all visits
foo <- df_bkt_adult %>%
  dplyr::ungroup() %>%
  dplyr::select(site_visit, site, date, year, pass, method, count)
df_bkt_3_pass <- tidyr::complete(ungroup(foo), nesting(site_visit, site, date, year, method), pass, fill = list(count = 0))

bar <- df_bkt_yoy %>%
  dplyr::ungroup() %>%
  dplyr::select(site_visit, site, date, year, pass, method, count)
df_bkt_3_pass_yoy <- tidyr::complete(ungroup(bar), nesting(site_visit, site, date, year, method), pass, fill = list(count = 0))

# test that expands to correct size
dim(df_bkt_adult)
dim(df_bkt_yoy)
dim(df_bkt_3_pass)[1] == length(unique(df_bkt_adult$site_visit)) * length(unique(df_bkt_adult$pass))
str(df_bkt_3_pass)
head(df_bkt_3_pass, 20)
dim(df_bkt_3_pass)[1] == dim(df_bkt_3_pass_yoy)[1]

# any sites visited multiple days in a year
length(unique(df_bkt_adult$site_visit))
length(unique(df_fish$site_year)) # yes

# interesting site 103-1 was sampled 2 days in a row in 1996 and have about half as many on the second day.

# convert long pass format to wide format
df_bkt_3_pass <- df_bkt_3_pass %>%
  dplyr::mutate(pass = paste0("pass_", pass)) %>%
  tidyr::spread(key = pass, value = count) 
dim(df_bkt_3_pass)

df_bkt_3_pass_yoy <- df_bkt_3_pass_yoy %>%
  dplyr::mutate(pass = paste0("pass_", pass)) %>%
  tidyr::spread(key = pass, value = count) 
dim(df_bkt_3_pass_yoy)


# compare captures in pass 1 to pass 2 to see if equal effort
df_bkt_3_pass$diff <- (df_bkt_3_pass$pass_1 - df_bkt_3_pass$pass_2)/df_bkt_3_pass$pass_1
ggplot(df_bkt_3_pass, aes(diff, fill = method)) + geom_histogram(binwidth=.5, alpha=.5, position="identity")
summarise(group_by(df_bkt_3_pass, method), min = min(diff), mean = mean(diff), max=max(diff))


df_bkt_3_pass_yoy$diff <- df_bkt_3_pass_yoy$pass_1 - df_bkt_3_pass_yoy$pass_2
ggplot(df_bkt_3_pass_yoy, aes(diff, fill = method)) + geom_histogram(binwidth=.5, alpha=.5, position="identity") #+ xlim(-20,300)
summarise(group_by(df_bkt_3_pass_yoy, method), min = min(diff), mean = mean(diff), max=max(diff))


######## Effort seems fine of first pass using Lincoln-Peterson #######
