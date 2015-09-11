

#----------- Generate Data -------------
n_sites <- 30
sites <- 1:n_sites

# abundance at each site

N <- rpois(n_sites, lambda = 10)

# probability of detection
p_mean <- 0.5

# generate 1st captures (assume all marked and returned)
captures1 <- rep(NA, times = n_sites)
for(i in 1:n_sites) {
  captures1[i] <- rbinom(1, N[i], p_mean)
}

# generate captures on pass 2, assume equal capturability and mixing
captures2 <- rep(NA, times = n_sites)
for(i in 1:n_sites) {
  captures2[i] <- rbinom(1, N[i], p_mean)
}

# generate recaptures, probability that an individual captured in time 2 was also catpured in time 1 given equal catchability
recaptures2 <- rep(NA, times = n_sites)
for(i in 1:n_sites) {
  recaptures2[i] <- min(rbinom(1, captures2[i], captures1[i] / N[i]), captures1[i])
}

# loop at the data
dat <- data.frame(site = 1:length(captures1), captures1, captures2, recaptures2)


#-------------- Traditional Estimates ------------
# Lincoln-Petersen Estimator
lp_est <- (captures1 * captures2) / recaptures2 # will fail if 0 recaptures

# estimate detection
lp_detect <- recaptures2 / captures2

# Chapman Adjustment to the L-P Estimator
chapman_est <- ((captures1 + 1) * (captures2 + 1)) / (recaptures2 + 1)

(out_mat <- cbind(N, captures1, captures2, recaptures2, lp_est, lp_detect, chapman_est))

apply(out_mat, MARGIN = 2, mean)


#--------- Convert to individual capture histories ------------

df_hist <- NULL
for(i in 1:nrow(dat)) {
  uni_caps <- dat$captures1[i] + dat$captures2[i] - dat$recaptures2[i]
  hist1 <- c(rep(1, times = dat$captures1[i]), rep(0, times = uni_caps - dat$captures1[i]))
  hist2 <- c(rep(1, times = dat$recaptures2[i]), rep(0, times = uni_caps - dat$recaptures2[i]))
  site <- rep(dat$site[i], times = length(hist1))
  foo <- cbind(site, hist1, hist2)
  df_hist <- rbind(df_hist, foo)
}
df_hist <- as.data.frame(df_hist)
library(dplyr)
df_hist <- df_hist %>%
  mutate(hist2 = ifelse(hist1 == 0, 1, hist2))

df_hist$hist <- paste0(df_hist$hist1, df_hist$hist2)

df_hist_sum <- df_hist %>%
  group_by(site, hist) %>%
  dplyr::summarise(count = n())

sites <- unique(df_hist_sum$site)
histories <- c("11", "10", "01")
combos <- expand.grid(site = sites, hist = histories)
df_hist_expand <- left_join(combos, df_hist_sum) %>%
  dplyr::arrange(site, hist)
df_hist_expand[is.na(df_hist_expand)] <- 0

table(df_hist_expand$site, df_hist_expand$hist)
library(tidyr)
spread(df_hist_expand, key = hist, value = count)

#--------- Function to convert capture histories to multinomial cell probabilities ---------
# Function based on unmarked cap-recap vignette

crPiFun <- function

#----------estimate
library(marked)



# library(Rcapture)
# The closedp.t command with method ”Mt” gives exactly the same results and standard errors as the Lincoln-Petersen estimator; while the Schnabel estimator and standard error is given by closedp.bc command. I confirm that all the calculated results are correct.



