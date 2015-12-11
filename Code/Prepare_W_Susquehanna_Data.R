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

# just use first pass from lincoln-peterson MR data
df_fish_MR <- df_fish_MR %>%
  dplyr::filter(EffortNumber == 1)

# combine with 3-pass data
df_fish$method <- "removal"
df_fish_MR$method <- "MR"
df_fish <- dplyr::bind_rows(df_fish, df_fish_MR)

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

# make paths
if(!exists(file.path(getwd(), dir_out, "/Figures/Stage/BKT"))) dir.create(file.path(getwd(), dir_out, "/Figures/Stage/BKT"), recursive = T)

# DON'T OVERWRITE
if(!exists(file.path(getwd(), dir_out, "/Figures/Stage/BKT/yoy_assignment_bkt.csv"))) {
  write.csv(unique(df_bkt$site_visit), file = file.path(getwd(), dir_out, "/Figures/Stage/BKT/yoy_assignment_bkt.csv"), row.names = F)
}

# export histograms
if(FALSE) {
  j <- 0
  for(i in 1:length(unique(df_bkt$site_visit))) {
    bar <- df_bkt %>% dplyr::filter(site_visit == unique(df_bkt$site_visit)[i])
    if(sum(bar$catch) >= 100) {
      j <- j + 1
      g <- ggplot(bar, aes(sizebin)) + geom_histogram(binwidth = 5, aes(weight = catch)) + ggtitle(paste0("BKT ", unique(df_bkt$site_visit)[i])) + xlim(50, 175) #+ geom_density(aes(y = 5*..count..)) 
      ggsave(filename = paste0(dir_out, "/Figures/Stage/BKT/", unique(df_bkt$site_visit)[i], ".png"))
    }
    print(j)
  }
}

# large size bins make it impossible to separate out YOY and 1+ very well. I will use 100 mm and larger as adults and below 100m as YOY

# Covariates each site-visit
df_covs_visit <- df_fish %>%
  group_by(site_visit, site, year, date) %>%
  dplyr::summarise(length_sample = mean(length, na.rm = T),
                   width = mean(width, na.rm = T),
                   effort = mean(SiteEffortHours, na.rm = T))

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
df_bkt_3_pass <- tidyr::complete(ungroup(foo), c(site_visit, site, date, year, method), pass, fill = list(count = 0))

bar <- df_bkt_yoy %>%
  dplyr::ungroup() %>%
  dplyr::select(site_visit, site, date, year, pass, method, count)
df_bkt_3_pass_yoy <- tidyr::complete(ungroup(bar), c(site_visit, site, date, year, method), pass, fill = list(count = 0))

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

# combine covariate data back in
df_trout <- left_join(df_covs_visit, df_bkt_3_pass) %>% distinct()
df_trout_yoy <- left_join(df_covs_visit, df_bkt_3_pass_yoy) %>% distinct()

# Fill zeros for sites where no fish were caught but it was visited
df_trout[is.na(df_trout)] <- 0
df_trout_yoy[is.na(df_trout_yoy)] <- 0

# replace 0 in second and third pass MR with NA
df_trout[which(df_trout$method == "MR"), c("pass_2", "pass_3")] <- NA
df_trout_yoy[which(df_trout_yoy$method == "MR"), c("pass_2", "pass_3")] <- NA

# Don't use sites over 1km in length as untrustworthy
df_trout[which(df_trout$length_sample >= 1000), c("pass_1", "pass_2", "pass_3")] <- NA
df_trout_yoy[which(df_trout_yoy$length_sample >= 1000), c("pass_1", "pass_2", "pass_3")] <- NA

# only use first visit each year? (possible that previous visit if close in time could change the abundance through emigration) - plus not sure how to handle it 
df_trout <- dplyr::arrange(df_trout, site, year, date)
df_trout$shift_site <- c(NA, df_trout$site[1:nrow(df_trout)-1])
df_trout$shift_year <- c(NA, df_trout$year[1:nrow(df_trout)-1])
df_trout <- df_trout %>%
  dplyr::filter(!(site == shift_site & year == shift_year))

df_trout_yoy <- dplyr::arrange(df_trout_yoy, site, year, date)
df_trout_yoy$shift_site <- c(NA, df_trout_yoy$site[1:nrow(df_trout_yoy)-1])
df_trout_yoy$shift_year <- c(NA, df_trout_yoy$year[1:nrow(df_trout_yoy)-1])
df_trout_yoy <- df_trout_yoy %>%
  dplyr::filter(!(site == shift_site & year == shift_year))

df <- left_join(family, df_trout, by = c("child_name" = "site"))
str(df)
df_yoy <- left_join(family, df_trout_yoy, by = c("child_name" = "site"))
str(df_yoy)

foo <- unique(family$child_name)
bar <- unique(c(df_trout$site, df_trout_yoy$site))

sna <- bar %in% foo

# need to add dates and other variables for nodes without measurements for covariates - do for other covariates after standardizing
df[is.na(df$date), "date"] <- "2010-08-01"

df <- df %>%
  dplyr::mutate(year = year(date),
                month = month(date)) # width should really be spatially interpolated but it won't affect the model

str(df)
summary(df)
head(df, 20)
tail(df, 20)


df_yoy[is.na(df_yoy$date), "date"] <- "2010-08-01"

df_yoy <- df_yoy %>%
  dplyr::mutate(year = year(date),
                month = month(date),
                length_sample = ifelse(is.na(length_sample), 300, length_sample),
                width = ifelse(is.na(width), 5, width),
                effort = ifelse(is.na(effort), mean(df_yoy$effort, na.rm = TRUE), effort)) # width should really be spatially interpolated but it won't affect the model

# check for weird distances (should all be positive)
nearby <- dplyr::filter(family, dist_b < 0.01) # nodes less than 10 m apart
df <- df %>%
  dplyr::mutate(dist_b = ifelse(dist_b < 0.001, 0.001, dist_b))
df_yoy <- df_yoy %>%
  dplyr::mutate(dist_b = ifelse(dist_b < 0.001, 0.001, dist_b))

############ Pull covariates from DB ############
# get featureid for database
df_loc <- read.csv("Data/occupancySitesSusqWest_Threepass.csv", header = TRUE, stringsAsFactors = FALSE)
df_loc$featureid <- as.numeric(gsub(",", "", df_loc$FEATUREID))

df_loc2 <- read.csv("Data/occupancySitesSusqWest_Petersen.csv", header = TRUE, stringsAsFactors = FALSE)
df_loc2$featureid <- as.numeric(gsub(",", "", df_loc2$FEATUREID))

dim(df_loc)
length(unique(df_loc$GIS_Key))
dim(df_loc2)
length(unique(df_loc2$GIS_Key))

df_loc2 <- df_loc2 %>%
  dplyr::filter(!(GIS_Key %in% unique(df_loc$GIS_Key)))

df_locs <- bind_rows(df_loc, df_loc2) %>%
  distinct()
dim(df_locs)
length(unique(df_locs$GIS_Key))

df <- left_join(df, df_locs[ , c("featureid", "GIS_Key")], by = c("child_name" = "GIS_Key")) %>%
  dplyr::mutate(featureid = ifelse(is.na(featureid), gsub("^N_", "", child_name), featureid))

df_yoy <- left_join(df_yoy, df_locs[ , c("featureid", "GIS_Key")], by = c("child_name" = "GIS_Key")) %>%
  dplyr::mutate(featureid = ifelse(is.na(featureid), gsub("^N_", "", child_name), featureid))

featureids <- unique(df$featureid)
featureid_string <- paste(shQuote(featureids), collapse = ", ")

# connect to database
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, dbname='sheds', host='osensei.cns.umass.edu', user=options('SHEDS_USERNAME'), password=options('SHEDS_PASSWORD'))

param_list <- c("forest", 
                "herbaceous", 
                "agriculture", 
                "devel_hi", 
                "developed",
                "AreaSqKM",  
                "allonnet",
                "alloffnet",
                "surfcoarse"
)

cov_list_string <- paste(shQuote(param_list), collapse = ", ")

qry_covariates <- paste0("SELECT * FROM covariates WHERE zone='upstream' AND variable IN (", cov_list_string, ") AND featureid IN (", featureid_string, ");")
rs <- dbSendQuery(con, qry_covariates)
df_covariates_long <- fetch(rs, n=-1)

# transform from long to wide format
df_covariates_upstream <- tidyr::spread(df_covariates_long, variable, value) %>%
  dplyr::mutate(featureid = as.character(featureid),
                impound = AreaSqKM * allonnet)

##########################

df <- left_join(df, df_covariates_upstream)
df_yoy <- left_join(df_yoy, df_covariates_upstream)


dbClearResult(res = rs)
dbDisconnect(con)
dbUnloadDriver(drv)
#dbListConnections(drv)

gc()

########## pull daymet records ##########
years <- unique(df_bkt$year)
# need to use previous fall and summer weather so need to get previous year if available (daymet only goes back to 1980)
if(min(years) > 1980) {
  years <- c(min(years)-1, years)
}
year_string <- paste(shQuote(years), collapse = ", ")

########## Set Up Parallel Processing ##########
library(foreach)
library(doParallel)
library(RPostgreSQL)

# size of chunks
catchmentid <- featureids
n.catches <- length(catchmentid)
chunk.size <- 50
n.loops <- ceiling(n.catches / chunk.size)

# set up parallel backend & make database connection available to all workers
nc <- min(c(detectCores(), 15)) # use the maximum number of cores minus 1 or up to 15 because max 16 database connections
cl <- makePSOCKcluster(nc)
registerDoParallel(cl)

# setup to write out to monitor progress
logFile = paste0("Diagnostics/log_file.txt")
cat("Monitoring progress of prediction loop in parallel", file=logFile, append=FALSE, sep = "\n")

########## Run Parallel Loop ########## 
# start loop
climate <- foreach(i = 1:n.loops, 
                   .inorder=FALSE, 
                   .combine = rbind,
                   .packages=c("DBI", 
                               "RPostgreSQL",
                               "dplyr",
                               "tidyr",
                               "lubridate")#,
                   #.export=ls(envir=globalenv())
) %dopar% {
  
  #try({
  ########## Set up database connection ##########
  drv <- dbDriver("PostgreSQL")
  con <- dbConnect(drv, dbname='NHDHRDV1_daymet_SusqWest', host='osensei.cns.umass.edu', user=options('OSENSEI_USERNAME'), password=options('SHEDS_PASSWORD'))
  
  # write start of each iteration
  cat(paste("Starting iteration", i, " of ", n.loops, "at", Sys.time(), "\n"), file = logFile, append = TRUE)
  
  ########## Set group of featureid to calculate for ##########
  k <- i*chunk.size
  if(k <= n.catches) {
    catches <- catchmentid[(1+(i-1)*chunk.size):k]
  } else {
    catches <- catchmentid[(1+(i-1)*chunk.size):n.catches]
  }
  catches_string <- paste(catches, collapse = ', ')
  
  # reconnect to database if lost
  #   if(isPostgresqlIdCurrent(con) == FALSE) {
  #     drv <- dbDriver("PostgreSQL")
  #     con <- dbConnect(drv, dbname='sheds', host='ecosheds.org', user=options('SHEDS_USERNAME'), password=options('SHEDS_PASSWORD'))
  #   }
  
  ########## pull daymet records ##########
  qry_daymet <- paste0("SELECT featureid, date, tmax, tmin, prcp, (tmax + tmin) / 2.0 AS airTemp, date_part('year', date) AS year FROM daymet WHERE featureid IN (", catches_string, ") ;")
  rs <- dbSendQuery(con, statement = qry_daymet)
  climateData <- fetch(rs, n=-1)
  dbDisconnect(con)
  dbUnloadDriver(drv)
  
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
  # get summer and fall precip from the previous year to predict current year abudance
  precip_mean$prcp_mean_summer_1 <- c(NA_real_, precip_mean$prcp_mean_summer[1:(nrow(precip_mean)-1)])
  precip_mean <- precip_mean %>%
    dplyr::mutate(prcp_mean_summer_1 = ifelse(year == min(precip_mean$year), NA_real_, prcp_mean_summer_1))
  precip_mean$prcp_mean_fall_1 <- c(NA_real_, precip_mean$prcp_mean_fall[1:(nrow(precip_mean)-1)])
  precip_mean <- precip_mean %>%
    dplyr::mutate(prcp_mean_fall_1 = ifelse(year == min(precip_mean$year), NA_real_, prcp_mean_fall_1))
  
  precip_sd <- df_climate %>%
    ungroup() %>%
    group_by(featureid, year) %>%
    dplyr::select(featureid, year, season, precip_sd) %>%
    tidyr::spread(season, precip_sd) %>%
    dplyr::rename(prcp_sd_fall = fall,
                  prcp_sd_spring = spring,
                  prcp_sd_summer = summer,
                  prcp_sd_winter = winter) %>%
    dplyr::arrange(featureid, year)
  # get summer and fall precip from the previous year to predict current year abudance
  precip_sd$prcp_sd_summer_1 <- c(NA_real_, precip_sd$prcp_sd_summer[1:(nrow(precip_sd)-1)])
  precip_sd <- precip_sd %>%
    dplyr::mutate(prcp_sd_summer_1 = ifelse(year == min(precip_sd$year), NA_real_, prcp_sd_summer_1))
  precip_sd$prcp_sd_fall_1 <- c(NA_real_, precip_sd$prcp_sd_fall[1:(nrow(precip_sd)-1)])
  precip_sd <- precip_sd %>%
    dplyr::mutate(prcp_sd_fall_1 = ifelse(year == min(precip_sd$year), NA_real_, prcp_sd_fall_1))
  
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
  # get the previous summer
  temp_mean$temp_mean_summer_1 <- c(NA_real_, temp_mean$temp_mean_summer[1:(nrow(temp_mean)-1)])
  temp_mean <- temp_mean %>%
    dplyr::mutate(temp_mean_summer_1 = ifelse(year == min(temp_mean$year), NA_real_, temp_mean_summer_1))
  temp_mean$temp_mean_fall_1 <- c(NA_real_, temp_mean$temp_mean_fall[1:(nrow(temp_mean)-1)])
  temp_mean <- temp_mean %>%
    dplyr::mutate(temp_mean_fall_1 = ifelse(year == min(temp_mean$year), NA_real_, temp_mean_fall_1))
  
  temp_sd <- df_climate %>%
    ungroup() %>%
    group_by(featureid, year) %>%
    dplyr::select(featureid, year, season, temp_sd) %>%
    tidyr::spread(season, temp_sd) %>%
    dplyr::rename(temp_sd_fall = fall,
                  temp_sd_spring = spring,
                  temp_sd_summer = summer,
                  temp_sd_winter = winter) %>%
    dplyr::arrange(featureid, year)
  # get the previous summer
  temp_sd$temp_sd_summer_1 <- c(NA_real_, temp_sd$temp_sd_summer[1:(nrow(temp_sd)-1)])
  temp_sd <- temp_sd %>%
    dplyr::mutate(temp_sd_summer_1 = ifelse(year == min(temp_sd$year), NA_real_, temp_sd_summer_1))
  temp_sd$temp_sd_fall_1 <- c(NA_real_, temp_sd$temp_sd_fall[1:(nrow(temp_sd)-1)])
  temp_sd <- temp_sd %>%
    dplyr::mutate(temp_sd_fall_1 = ifelse(year == min(temp_sd$year), NA_real_, temp_sd_fall_1))
  
  
  climate_summary <- temp_mean %>%
    left_join(temp_sd) %>%
    left_join(precip_mean) %>%
    left_join(precip_sd)
  
  return(climate_summary)
  #})
  
  #saveRDS(metrics.lat.lon, file = paste0(data_dir, "/derived_site_metrics.RData"))
  #write.table(metrics.lat.lon, file = paste0(data_dir, "/derived_site_metrics.csv"), sep = ',', row.names = F)
  
} # end dopar

stopCluster(cl)

dbClearResult(res = rs)
dbDisconnect(con)
dbUnloadDriver(drv)
dbListConnections(drv)

gc()

saveRDS(climate, file = file.path(dir_out, "climate.RData"))

climate$featureid <- as.character(climate$featureid)
df <- left_join(df, climate)
df_yoy <- left_join(df_yoy, climate)

# remove data from pre-1981 because no daymet records for previous summer and fall

df <- df %>%
  dplyr::filter(year > 1980)
df_yoy <- df_yoy %>%
  dplyr::filter(year >1980)

################################################

# Separate count data
c_ip <- dplyr::select(df, starts_with("pass"))
c_ip_yoy <- dplyr::select(df_yoy, starts_with("pass"))

# YOY not captured well until maybe June so make earlier samples NA for YOY.
df_yoy$date <- as.Date(df_yoy$date)
df_yoy$month <- month(df_yoy$date)
c_ip_yoy[which(df_yoy$month < 6), c("pass_1", "pass_2", "pass_3")] <- NA

# pull out years if want to use as fixed effects
year <- as.factor(df$year)
dummies <- model.matrix(~year)

# make covariate list
cov_list <- c("length_sample", "width", "effort", "surfcoarse", "forest", "AreaSqKM", "impound", names(climate))
cov_list <- cov_list[which(!(cov_list %in% c("featureid", "year")))]

# extract standardized covariates
X_ij <- df %>%
  dplyr::select(one_of(cov_list)) %>%
  mutate_each_(funs(std_vec), vars = cov_list)

# assign all missing values the mean
X_ij[is.na(X_ij)] <- 0

# save means and sds used for standardization
means <- df %>%
  dplyr::select(one_of(cov_list)) %>%
  apply(MARGIN = 2, FUN = mean, na.rm = TRUE)
sds <- df %>%
  dplyr::select(one_of(cov_list)) %>%
  apply(MARGIN = 2, FUN = sd, na.rm = TRUE)
df_stds <- data.frame(var = names(means), mean = as.numeric(means), sd = as.numeric(sds))

# year vector for autoregressive
t_i <- df$year

# save
save(c_ip, c_ip_yoy, X_ij, df_stds, t_i, df, df_yoy, family, file = "Data/Prepared_Data_W_Susquehanna_Combined.RData")
