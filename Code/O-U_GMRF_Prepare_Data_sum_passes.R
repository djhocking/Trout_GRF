# Data Prep

#######################
# Load libraries
#######################
library(dplyr)
library(lubridate)

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

# First try use sum of three passes for adults only
df_bkt_reduced <- df_bkt %>%
  dplyr::filter(net == TRUE & species == "BKT" & stage == "adult") %>%
  dplyr::group_by(visit, location_id, date, year, month) %>%
  dplyr::summarise(count = sum(count),
                   length_sample = mean(length_sample),
                   width = mean(width))

# just use the most recent year a site was sampled to avoid temporal trends at sites
first_site_visits <- df_bkt_reduced %>%
  dplyr::ungroup() %>%
  dplyr::group_by(location_id) %>%
  dplyr::summarise(year = max(year))

df_trout <- left_join(first_site_visits, df_bkt_reduced)

df <- left_join(family, df_trout, by = c("child_name" = "location_id")) %>%
  dplyr::select(child_name, parent_b, dist_b, child_b, NodeLat, NodeLon, count, date, length_sample, width)
str(df)

df[is.na(df$date), "date"] <- "2010-08-01"

df <- df %>%
  dplyr::mutate(year = year(date),
                month = month(date),
                length_sample = ifelse(is.na(length_sample), 100, length_sample),
                width = ifelse(is.na(width), 10, width)) # width should really be spatially interpolated but it won't affect the model


c_i <- df$count
year <- as.factor(df$year)
dummies <- model.matrix(~year)
length_std <- (df$length_sample - mean(df$length_sample)) / sd(df$length_sample)
width_std <- (df$width - mean(df$width)) / sd(df$width)
X_ij = cbind(dummies[ , 2:ncol(dummies)], length_std, width_std)

save(c_i, X_ij, df, file = "Data/Prepared_Data.RData")