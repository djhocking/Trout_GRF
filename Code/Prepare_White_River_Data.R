# Data Prep

#######################
# Load libraries
#######################
library(dplyr)
library(lubridate)
library(tidyr)

#######################
# Load data
#######################

load( "Data/White_River_Network.RData")
colnames(family)[1] = "child_name"
family = cbind( family, "child_b"=1:nrow(family) )

df_bkt <- readRDS("Data/White_River_Trout.Rdata")

df_locations <- read.csv("Data/upperWhitePoints.csv", header = TRUE, stringsAsFactors = FALSE)
df_locations$featureid <- as.numeric(gsub(",", "", df_locations$featureid))

#######################
# Prep Trout Data
#######################

# First try use three passes for adults only for years with block nets
df_bkt_reduced <- df_bkt %>%
  dplyr::filter(net == TRUE & species == "BKT" & stage == "adult") %>%
  dplyr::group_by(visit, location_id, date, year, month, pass) %>%
  dplyr::summarise(count = sum(count),
                   length_sample = mean(length_sample),
                   width = mean(width))

# just use the most recent date a site was sampled to avoid temporal trends at sites
# recent_site_visits <- df_bkt_reduced %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(location_id) %>%
#   dplyr::summarise(date = max(date)) %>%
#   dplyr::mutate(visit = paste0(location_id, "_", date))
# 
 df_visit <- df_bkt_reduced %>%
   dplyr::select(visit, location_id, date, year, month, length_sample, width)
# 
# df_recent_visits <- left_join(recent_site_visits, df_visit) %>%
#   distinct() # i don't understand why this was needed but only way not duplicated -- maybe due to grouping of one or both of the dataframes - groupings seems to cause problems, maybe should always ungroup at end of pipes.

# expand to 3-pass for all visits
foo <- df_bkt_reduced %>%
  dplyr::ungroup() %>%
  dplyr::select(visit, location_id, year, pass, count)

df_bkt_3_pass <- tidyr::complete(ungroup(foo), c(visit, location_id, year), pass, fill = list(count = 0))

# test that expands to correct size
dim(df_bkt_reduced)
dim(df_bkt_3_pass)[1] == length(unique(df_bkt_reduced$visit)) * length(unique(df_bkt_reduced$pass))
str(df_bkt_3_pass)

# convert long pass format to wide format
df_bkt_3_pass <- df_bkt_3_pass %>%
  dplyr::mutate(pass = paste0("pass_", pass)) %>%
  tidyr::spread(key = pass, value = count) 

# reduce to most recent site-visit
#df_trout <- left_join(df_recent_visits, df_bkt_3_pass)


############ Pull covariates from DB ############
library(RPostgreSQL)
library(tidyr)

# get featureid for database
df <- left_join(family, df_locations[ , c("featureid", "location_id")], by = c("child_name" = "location_id")) %>%
  dplyr::mutate(featureid = ifelse(is.na(featureid), gsub("^N_", "", child_name), featureid))

featureids <- unique(df$featureid)
featureid_string <- paste(shQuote(featureids), collapse = ", ")

# connect to database
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, dbname='sheds', host='felek.cns.umass.edu', user=options('SHEDS_USERNAME'), password=options('SHEDS_PASSWORD'))

param_list <- c("forest", 
                "herbaceous", 
                "agriculture", 
                "devel_hi", 
                "developed",
                "AreaSqKM",  
                "allonnet",
                "alloffnet",
                "surfcoarse", 
                "srad", 
                "dayl", 
                "swe")

cov_list_string <- paste(shQuote(param_list), collapse = ", ")

qry_covariates <- paste0("SELECT * FROM covariates WHERE zone='upstream' AND variable IN (", cov_list_string, ") AND featureid IN (", featureid_string, ");")
rs <- dbSendQuery(con, qry_covariates)
df_covariates_long <- fetch(rs, n=-1)

# transform from long to wide format
df_covariates_upstream <- tidyr::spread(df_covariates_long, variable, value) %>%
  dplyr::mutate(featureid = as.character(featureid),
                impound = AreaSqKM * allonnet)

# reconnect to database if lost
if(isPostgresqlIdCurrent(con) == FALSE) {
  drv <- dbDriver("PostgreSQL")
  con <- dbConnect(drv, dbname='sheds', host='felek.cns.umass.edu', user=options('SHEDS_USERNAME'), password=options('SHEDS_PASSWORD'))
}

########## pull daymet records ##########
years <- unique(df_bkt$year)
# need to use previous fall and summer weather so need to get previous year if available (daymet only goes back to 1980)
if(min(years) > 1980) {
  years <- c(min(years)-1, years)
}
year_string <- paste(shQuote(years), collapse = ", ")

qry_daymet <- paste0("SELECT featureid, date, tmax, tmin, prcp, swe, (tmax + tmin) / 2.0 AS airTemp, date_part('year', date) AS year FROM daymet WHERE featureid IN (", featureid_string, ") AND date_part('year', date) IN (", year_string, ") ;")
rs <- dbSendQuery(con, statement = qry_daymet)
climateData <- fetch(rs, n=-1)


df_climate <- climateData %>%
  group_by(featureid, year) %>%
  arrange(featureid, year, date) %>%
  dplyr::mutate(month = month(date),
                season = ifelse(month %in% c(1,2,3), "winter",
                                ifelse(month %in% c(4,5,6), "spring", 
                                       ifelse(month %in% c(7,8,9), "summer",
                                              ifelse(month %in% c(10,11,12), "fall", NA))))) %>%
  group_by(featureid, year, season) %>%
  dplyr::summarise(precip_mean = mean(prcp, na.rm = TRUE),
                   precip_sd = sd(prcp, na.rm = TRUE),
                   prceip_max = max(prcp, na.rm = TRUE),
                   temp_mean = mean(airtemp, na.rm = TRUE),
                   temp_sd = sd(airtemp, na.rm = TRUE),
                   temp_max = max(airtemp, na.rm = TRUE))

library(tidyr)
precip_mean <- df_climate %>%
  ungroup() %>%
  group_by(featureid, year) %>%
  dplyr::select(featureid, year, season, precip_mean) %>%
  tidyr::spread(season, precip_mean) %>%
  dplyr::rename(prcp_mean_fall = fall,
                prcp_mean_spring = spring,
                prcp_mean_summer = summer,
                prcp_mean_winter = winter) %>%
  dplyr::arrange(featureid, year)
precip_mean$prcp_mean_summer_1 <- c(NA_real_, precip_mean$prcp_mean_summer[1:(nrow(precip_mean)-1)])
precip_mean <- precip_mean %>%
  dplyr::mutate(prcp_mean_summer_1 = ifelse(year == min(precip_mean$year), NA_real_, prcp_mean_summer_1))

precip_sd <- df_climate %>%
  ungroup() %>%
  group_by(featureid, year) %>%
  dplyr::select(featureid, year, season, precip_sd) %>%
  tidyr::spread(season, precip_sd)
temp_mean <- df_climate %>%
  ungroup() %>%
  group_by(featureid, year) %>%
  dplyr::select(featureid, year, season, temp_mean) %>%
  tidyr::spread(season, temp_mean) %>%
  dplyr::rename(temp_mean_fall = fall,
                temp_mean_spring = spring,
                temp_mean_summer = summer,
                temp_mean_winter = winter) %>%
  dplyr::arrange(featureid, year)
temp_mean$temp_mean_summer_1 <- c(NA_real_, temp_mean$temp_mean_summer[1:(nrow(temp_mean)-1)])
temp_mean <- temp_mean %>%
  dplyr::mutate(temp_mean_summer_1 = ifelse(year == min(temp_mean$year), NA_real_, temp_mean_summer_1))
temp_sd <- df_climate %>%
  ungroup() %>%
  group_by(featureid, year) %>%
  dplyr::select(featureid, year, season, temp_sd) %>%
  tidyr::spread(season, temp_sd)


# combine covariate data back in
df_trout <- left_join(df_bkt_3_pass, df_visit) %>% 
  distinct() %>%
  dplyr::left_join(dplyr::select(df_locations, location_id, featureid)) %>%
  left_join(temp_mean) %>%
  left_join(precip_mean)

df <- left_join(family, df_trout, by = c("child_name" = "location_id"))
str(df)

# need to add dates and other variables for nodes without measurements for covariates
df[is.na(df$date), "date"] <- "2010-08-01"

df <- df %>%
  dplyr::mutate(year = year(date),
                month = month(date),
                length_sample = ifelse(is.na(length_sample), 100, length_sample),
                width = ifelse(is.na(width), 10, width)) # width should really be spatially interpolated but it won't affect the model

str(df)
summary(df)
head(df, 20)
tail(df, 20)

c_ip <- dplyr::select(df, starts_with("pass"))
year <- as.factor(df$year)
dummies <- model.matrix(~year)
std <- function(data, var) {
  var_std <- (data[ , c(var)] - mean(data[ , c(var)], na.rm = TRUE)) / sd(data[ , c(var)], na.rm = TRUE)
}
length_std <- (df$length_sample - mean(df$length_sample)) / sd(df$length_sample)
width_std <- (df$width - mean(df$width)) / sd(df$width)
temp_summer_std <- std(df, "temp_mean_summer_1")
prcp_winter_std <- std(df, "prcp_mean_winter")
X_ij = data.frame(dummies[ , 2:ncol(dummies)], length_std, width_std, temp_summer_std, prcp_winter_std) #[ , 2:ncol(dummies)]

df_stds <- data.frame(parameter = c("length", "width"), means = c(mean(df$length_sample), mean(df$width)), sds = c(sd(df$length_sample), sd(df$width)))

t_i <- df$year

save(c_ip, X_ij, df_stds, t_i, df, family, file = "Data/Prepared_Data_White_River.RData")
