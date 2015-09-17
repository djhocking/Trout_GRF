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

df_bkt <- read.csv("Data/West_Susquehanna_3pass.csv", stringsAsFactors = FALSE)

#######################
# Prep Trout Data
#######################

obs_sites_network <- family %>%
  dplyr::select(child_name) %>%
  distinct() %>%
  dplyr::filter(!(grepl("N_", child_name)))

# assign stage class
df_bkt <- df_bkt %>%
  dplyr::filter(species == "Brook Trout") %>%
  dplyr::mutate(site = GIS_Key,
                site_visit = paste0(site, "_", date),
                site_year = paste0(site, "_", year),
                pass = EffortNumber)

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


# large size bins make it impossible to separate out YOY and 1+ very well. I will use 125 mm and large as adults (stable part of population with potential for reproduction and 100m and below as YOY)




# First try use three passes for adults only for years with block nets
df_bkt_adult <- df_bkt %>%
  dplyr::mutate(stage = ifelse(sizebin >= 125, "adult", "yoy")) %>%
  dplyr::filter(species == "Brook Trout" & stage == "adult") %>%
  dplyr::group_by(site_visit, site, year, date, pass) %>%
  dplyr::summarise(count = sum(catch),
                   length_sample = mean(length, na.rm = T),
                   width = mean(width, na.rm = T),
                   effort = mean(SiteEffortHours, na.rm = T)) # effort the total over all 3 passes

df_covs_visit <- df_bkt_adult %>%
  group_by(site_visit, site, year, date) %>%
  dplyr::summarise(length_sample = mean(length_sample, na.rm = T),
                   width = mean(width, na.rm = T),
                   effort = mean(effort, na.rm = T))

# expand to 3-pass for all visits
foo <- df_bkt_adult %>%
  dplyr::ungroup() %>%
  dplyr::select(site_visit, site, date, year, pass, count)

df_bkt_3_pass <- tidyr::complete(ungroup(foo), c(site_visit, site, date, year), pass, fill = list(count = 0))

# test that expands to correct size
dim(df_bkt_adult)
dim(df_bkt_3_pass)[1] == length(unique(df_bkt_adult$site_visit)) * length(unique(df_bkt_adult$pass))
str(df_bkt_3_pass)
head(df_bkt_3_pass, 20)

# any sites visited multiple days in a year
length(unique(df_bkt_adult$site_visit))
length(unique(df_bkt$site_year)) # yes

# interesting site 103-1 was sampled 2 days in a row in 1996 and have about half as many on the second day.

# convert long pass format to wide format
df_bkt_3_pass <- df_bkt_3_pass %>%
  dplyr::mutate(pass = paste0("pass_", pass)) %>%
  tidyr::spread(key = pass, value = count) 
dim(df_bkt_3_pass)

# combine covariate data back in
df_trout <- left_join(df_bkt_3_pass, df_covs_visit) %>% distinct()

# only use first visit each year? (possible that previous visit if close in time could change the abundance through emigration) - plus not sure how to handle it 
df_trout <- dplyr::arrange(df_trout, site, year, date)
df_trout$shift_site <- c(NA, df_trout$site[1:nrow(df_trout)-1])
df_trout$shift_year <- c(NA, df_trout$year[1:nrow(df_trout)-1])

df_trout <- df_trout %>%
  dplyr::filter(!(site == shift_site & year == shift_year))

df <- left_join(family, df_trout, by = c("child_name" = "site"))
str(df)

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

# check for weird distances (should all be positive)
nearby <- dplyr::filter(family, dist_b < 0.01) # nodes less than 10 m apart
df <- df %>%
  dplyr::mutate(dist_b = ifelse(dist_b < 0.01, 0.01, dist_b))


c_ip <- dplyr::select(df, starts_with("pass"))
year <- as.factor(df$year)
dummies <- model.matrix(~year)
length_std <- (df$length_sample - mean(df_covs_visit$length_sample)) / sd(df_covs_visit$length_sample)
width_std <- (df$width - mean(df_covs_visit$width)) / sd(df_covs_visit$width)
effort_std <- (df$effort - mean(df_covs_visit$effort)) / sd(df_covs_visit$effort)
X_ij = cbind(dummies, length_std, width_std, effort_std) #[ , 2:ncol(dummies)]

df_stds <- data.frame(parameter = c("length", "width", "effort"), means = c(mean(df$length_sample), mean(df$width), mean(df$effort)), sds = c(sd(df$length_sample), sd(df$width), sd(df$effort)))

t_i <- df$year

save(c_ip, X_ij, df_stds, t_i, df, family, file = "Data/Prepared_Data_W_Susquehanna.RData")
