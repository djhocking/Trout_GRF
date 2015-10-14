# Data Prep

#######################
# Load libraries
#######################
library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)

# clear environment
rm(list = ls())
gc()

dir_out <- "Output"

#######################
# Load data
#######################

family <- read.csv( "Data/West_Susquehanna_3pass_Network.csv", stringsAsFactors = FALSE)
colnames(family)[1] = "child_name"
family = cbind( family, "child_b"=1:nrow(family) )

df_fish <- read.csv("Data/West_Susquehanna_3pass.csv", stringsAsFactors = FALSE)

#######################
# Prep Trout Data
#######################

obs_sites_network <- family %>%
  dplyr::select(child_name) %>%
  distinct() %>%
  dplyr::filter(!(grepl("N_", child_name)))

# assign stage class
df_fish <- df_fish %>%
  dplyr::mutate(site = GIS_Key,
                site_visit = paste0(site, "_", date),
                site_year = paste0(site, "_", year),
                pass = EffortNumber)

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

# First try use three passes for adults only for years with block nets
df_fish <- df_fish %>%
  dplyr::group_by(site_visit, site, year, date, pass) %>%
  dplyr::mutate(length_sample = mean(length, na.rm = T),
                   width = mean(width, na.rm = T),
                   effort = mean(SiteEffortHours, na.rm = T)) # effort the total over all 3 passes

df_covs_visit <- df_fish %>%
  group_by(site_visit, site, year, date) %>%
  dplyr::summarise(length_sample = mean(length_sample, na.rm = T),
                   width = mean(width, na.rm = T),
                   effort = mean(effort, na.rm = T))

df_bkt_adult <- df_fish %>%
  dplyr::mutate(stage = ifelse(sizebin >= 100, "adult", "yoy")) %>%
  dplyr::filter(species == "Brook Trout" & stage == "adult") %>%
  dplyr::group_by(site_visit, site, year, date, pass) %>%
  dplyr::summarise(count = sum(catch))

df_bkt_yoy <- df_fish %>%
  dplyr::mutate(stage = ifelse(sizebin < 100, "adult", "yoy")) %>%
  dplyr::filter(species == "Brook Trout" & stage == "yoy") %>%
  dplyr::group_by(site_visit, site, year, date, pass) %>%
  dplyr::summarise(count = sum(catch))

# expand to 3-pass for all visits
foo <- df_bkt_adult %>%
  dplyr::ungroup() %>%
  dplyr::select(site_visit, site, date, year, pass, count)

bar <- df_bkt_yoy %>%
  dplyr::ungroup() %>%
  dplyr::select(site_visit, site, date, year, pass, count)

df_bkt_3_pass <- tidyr::complete(ungroup(foo), c(site_visit, site, date, year), pass, fill = list(count = 0))

df_bkt_3_pass_yoy <- tidyr::complete(ungroup(bar), c(site_visit, site, date, year), pass, fill = list(count = 0))

# test that expands to correct size
dim(df_bkt_adult)
dim(df_bkt_3_pass)[1] == length(unique(df_bkt_adult$site_visit)) * length(unique(df_bkt_adult$pass))
str(df_bkt_3_pass)
head(df_bkt_3_pass, 20)
dim(df_bkt_3_pass)[1] == dim(df_bkt_3_pass_yoy)[1]

# any sites visited multiple days in a year
length(unique(df_bkt_adult$site_visit))
length(unique(df_bkt$site_year)) # yes

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
df_trout <- left_join(df_bkt_3_pass, df_covs_visit) %>% distinct()
df_trout_yoy <- left_join(df_bkt_3_pass_yoy, df_covs_visit) %>% distinct()

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

# need to add dates and other variables for nodes without measurements for covariates
df[is.na(df$date), "date"] <- "2010-08-01"

df <- df %>%
  dplyr::mutate(year = year(date),
                month = month(date),
                length_sample = ifelse(is.na(length_sample), 300, length_sample),
                width = ifelse(is.na(width), 5, width),
                effort = ifelse(is.na(effort), mean(df$effort, na.rm = TRUE), effort)) # width should really be spatially interpolated but it won't affect the model

str(df)
summary(df)
head(df, 20)
tail(df, 20)

df[is.na(df$date), "date"] <- "2010-08-01"
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
library(RPostgreSQL)
library(tidyr)

# get featureid for database
df_loc <- read.csv("Data/occupancySitesSusqWest_Threepass.csv", header = TRUE, stringsAsFactors = FALSE)
df_loc$featureid <- as.numeric(gsub(",", "", df_loc$FEATUREID))

df <- left_join(df, df_loc[ , c("featureid", "GIS_Key")], by = c("child_name" = "GIS_Key")) %>%
  dplyr::mutate(featureid = ifelse(is.na(featureid), gsub("^N_", "", child_name), featureid))

df_yoy <- left_join(df_yoy, df_loc[ , c("featureid", "GIS_Key")], by = c("child_name" = "GIS_Key")) %>%
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
qry_daymet <- paste0("SELECT featureid, date, tmax, tmin, prcp, dayl, srad, vp, swe, (tmax + tmin) / 2.0 AS airTemp FROM daymet WHERE featureid IN (", featureid_string, ") ;")
rs <- dbSendQuery(con, statement = qry_daymet)
climateData <- fetch(rs, n=-1)

dbClearResult(res = rs)
dbDisconnect(con)
dbUnloadDriver(drv)
#dbListConnections(drv)

climateData <- climateData %>%
  group_by(featureid, year) %>%
  arrange(featureid, year, date) %>%
  dplyr::mutate(month = month(date),
                season = ifelse(month %in% c(1,2,3), "winter",
                                ifelse(month %in% c(4,5,6), "spring", 
                                       ifelse(month %in% c(7,8,9), "summer",
                                              ifelse(month %in% c(10,11,12), "fall", NA))))) %>%
  group_by(featureid, year, season)
  
  # fall precip (previous fall)
  dplyr::mutate(precip)
  
  # winter precip (include december of previous year?)
  
  # summer precip
  
  # fall temp (previous year)
  
  # summer temp
  
# summer temp (previous year)
  
  mutate(impoundArea = AreaSqKM * allonnet,
         airTempLagged1 = lag(airTemp, n = 1, fill = NA),
         temp5p = rollapply(data = airTempLagged1, 
                            width = 5, 
                            FUN = mean, 
                            align = "right", 
                            fill = NA, 
                            na.rm = T),
         temp7p = rollapply(data = airTempLagged1, 
                            width = 7, 
                            FUN = mean, 
                            align = "right", 
                            fill = NA, 
                            na.rm = T),
         prcp2 = rollsum(x = prcp, 2, align = "right", fill = NA),
         prcp7 = rollsum(x = prcp, 7, align = "right", fill = NA),
         prcp30 = rollsum(x = prcp, 30, align = "right", fill = NA))


df <- left_join(df, df_covariates_upstream)
df_yoy <- left_join(df_yoy, df_covariates_upstream)

gc()
################################################

c_ip <- dplyr::select(df, starts_with("pass"))
c_ip_yoy <- dplyr::select(df_yoy, starts_with("pass"))
year <- as.factor(df$year)
dummies <- model.matrix(~year)
length_std <- (df$length_sample - mean(df_covs_visit$length_sample)) / sd(df_covs_visit$length_sample)
width_std <- (df$width - mean(df_covs_visit$width)) / sd(df_covs_visit$width)
effort_std <- (df$effort - mean(df_covs_visit$effort, na.rm = T)) / sd(df_covs_visit$effort, na.rm = T)
surfcoarse_std <- (df$surfcoarse - mean(df$surfcoarse, na.rm = T)) / sd(df$surfcoarse, na.rm = T)
forest_std <- (df$forest - mean(df$forest, na.rm = T)) / sd(df$forest, na.rm = T)
area_std <- (df$AreaSqKM - mean(df$AreaSqKM, na.rm = T)) / sd(df$AreaSqKM, na.rm = T)
impound_std <- (df$impound - mean(df$impound, na.rm = T)) / sd(df$impound, na.rm = T)
X_ij = data.frame(dummies[ , 2:ncol(dummies)], 
                  length_std, 
                  width_std, 
                  effort_std, 
                  surfcoarse_std, 
                  forest_std, 
                  area_std, 
                  impound_std) #

df_stds <- data.frame(parameter = c("length", "width", "effort", "surfcoarse", "forest", "area", "impound"), 
                      means = c(mean(df_covs_visit$length_sample), 
                                mean(df_covs_visit$width), 
                                mean(df_covs_visit$effort, na.rm = T),
                                mean(df$surfcoarse, na.rm = T),
                                mean(df$forest, na.rm = T),
                                mean(df$AreaSqKM, na.rm = T),
                                mean(df$impound, na.rm = T)), 
                      sds = c(sd(df_covs_visit$length_sample), 
                              sd(df_covs_visit$width), 
                              sd(df_covs_visit$effort),
                              sd(df$surfcoarse, na.rm = T),
                              sd(df$forest, na.rm = T),
                              sd(df$AreaSqKM, na.rm = T),
                              sd(df$impound, na.rm = T)))

t_i <- df$year

save(c_ip, c_ip_yoy, X_ij, df_stds, t_i, df, df_yoy, family, file = "Data/Prepared_Data_W_Susquehanna.RData")
