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

# combine covariate data back in
df_trout <- left_join(df_bkt_3_pass, df_visit) %>% distinct()

df <- left_join(family, df_trout, by = c("child_name" = "location_id")) %>%
  dplyr::select(child_name, parent_b, dist_b, child_b, NodeLat, NodeLon, pass_1, pass_2, pass_3, date, length_sample, width)
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
length_std <- (df$length_sample - mean(df$length_sample)) / sd(df$length_sample)
width_std <- (df$width - mean(df$width)) / sd(df$width)
X_ij = data.frame(dummies[ , 2:ncol(dummies)], length_std, width_std) #[ , 2:ncol(dummies)]

df_stds <- data.frame(parameter = c("length", "width"), means = c(mean(df$length_sample), mean(df$width)), sds = c(sd(df$length_sample), sd(df$width)))

t_i <- df$year

save(c_ip, X_ij, df_stds, t_i, df, family, file = "Data/Prepared_Data_White_River.RData")
