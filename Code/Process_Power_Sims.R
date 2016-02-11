mean_N <- 50
n_years_vec <- c(2, 4, 8, 10, 15, 20)
n_years <- max(n_years_vec)
sample_sites_vec <- c(10, 25, 50, 100, nrow(family))
p <- c(0.75, 0.75, 0.75)
theta <- 0.1
SD <- 0.05
rhot <- 0.8
SD_t <- 0.1
theta_st <- 0.1
SD_st <- 0.05
rho <- 0.5

n_sim <- 100
counter <- 0
for(i in 1:n_sim) {
  for(b in 1:length(sample_sites_vec)) {
    for(t in 1:length(n_years_vec)) {
      for(s in 1:2) {
        counter <- counter + 1
        if(file.exists(paste0("Output/Power_Sim/Data/sim_", i, "_st_", s-1, "_sites_", sample_sites_vec[b], "_years_", n_years_vec[t], ".RData"))) {
          load(file = paste0("Output/Power_Sim/Data/sim_", i, "_st_", s-1, "_sites_", sample_sites_vec[b], "_years_", n_years_vec[t], ".RData"))
          df_N <- data.frame(N_i = network$N_i, N_hat = mod$Report$N_ip[,1])
          if(class(mod) == "try-error") {
            dat[counter, "iter"] <- i
            dat[counter, "n_sites"] <- sample_sites_vec[b]
            dat[counter, "n_years"] <- n_years_vec[t]
            dat[counter, "spatialTF"] <- s - 1
            dat[counter, "mean_N"] <- mean(network$N_i)
            dat[counter, "mean_N_est"] <- NA_real_
            dat[counter, "N_se"] <- NA_real_
            dat[counter, "RMSE"] <- NA_real_
            dat[counter, "theta"] <- theta
            dat[counter, "theta_hat"] <- NA_real_
            dat[counter, "rhot"] <- rhot
            dat[counter, "rhot_hat"] <- NA_real_
            dat[counter, "theta_st"] <- theta_st
            dat[counter, "theta_st_hat"] <- NA_real_
          } else {
            N_se <- mod$SD$sd[which(names(mod$SD$value) == "mean_N")]
            
            dat[counter, "iter"] <- i
            dat[counter, "n_sites"] <- sample_sites_vec[b]
            dat[counter, "n_years"] <- n_years_vec[t]
            dat[counter, "spatialTF"] <- s - 1
            dat[counter, "mean_N"] <- mean(network$N_i)
            dat[counter, "mean_N_est"] <- mod$Report$mean_N
            if(N_se == "NaN") {
              dat[counter, "N_se"] <- NA_real_
            } else {
              dat[counter, "N_se"] <- N_se
            }
            dat[counter, "RMSE"] <- rmse(df_N$N_i - df_N$N_hat)
            dat[counter, "theta"] <- theta
            dat[counter, "theta_hat"] <- mod$Report$theta
            dat[counter, "rhot"] <- rhot
            dat[counter, "rhot_hat"] <- mod$Report$rhot
            dat[counter, "theta_st"] <- theta_st
            dat[counter, "theta_st_hat"] <- mod$Report$theta_st
          }
        } else {
          dat[counter, "iter"] <- i
          dat[counter, "n_sites"] <- sample_sites_vec[b]
          dat[counter, "n_years"] <- n_years_vec[t]
          dat[counter, "spatialTF"] <- s - 1
          dat[counter, "mean_N"] <- NA_real_
          dat[counter, "mean_N_est"] <- NA_real_
          dat[counter, "N_se"] <- NA_real_
          dat[counter, "RMSE"] <- NA_real_
          dat[counter, "theta"] <- theta
          dat[counter, "theta_hat"] <- NA_real_
          dat[counter, "rhot"] <- rhot
          dat[counter, "rhot_hat"] <- NA_real_
          dat[counter, "theta_st"] <- theta_st
          dat[counter, "theta_st_hat"] <- NA_real_
        }# end if file exists statement
      }
    }
  }
}

df_sims <- dat %>%
  dplyr::mutate(converged = ifelse(is.na(N_se) | N_se > mean_N_est | is.na(mean_N_est), 0, 1))

saveRDS(df_sims, file = "Output/Power_Sim/STsim_Results.RData")

str(df_sims)
df_sim_summary <- df_sims %>%
  dplyr::group_by(n_sites, n_years, spatialTF) %>%
  dplyr::select(-iter) %>%
  dplyr::summarise_each(funs(mean(., na.rm = T)))
as.data.frame(df_sim_summary %>% dplyr::filter(spatialTF == 1)) 

# what predicts prob of convergence
library(lme4)
df_sims_s <- scale(df_sims[ ,2:ncol(df_sims)])
df_sims_s <- data.frame(iter = df_sims$iter, df_sims_s, converge = df_sims$converged)
converge_lm <- glmer(converge ~ n_sites + n_years + spatialTF + n_years * spatialTF + (1|iter), data = df_sims_s, family = "binomial")
summary(converge_lm)

df_sims %>%
  dplyr::group_by(n_sites, spatialTF) %>%
  dplyr::select(-iter) %>%
  dplyr::summarise_each(funs(mean(., na.rm = T)))

df_sims %>%
  dplyr::group_by(n_years, spatialTF) %>%
  dplyr::select(-iter) %>%
  dplyr::summarise_each(funs(mean(., na.rm = T)))

SD_array <- array(NA, c(6, 3, 100))


SD_means <- data.frame(param = names(mod$SD$value), 
                       est = as.numeric(mod$SD$value), 
                       sd = as.numeric(mod$SD$sd), stringsAsFactors = F)

# N_est <- SD_means %>%
#   dplyr::filter(grepl("N", param))

N_est <- mod$Report$N_ip[ , 1]

plot(network$N_i, N_est)
abline(0,1, col = "blue")


N_hat[[i]] <- data.frame(N_hat = mod$Report$N_ip[ , 1], N_SD = N_est$sd, lambda_hat = mod$Report$lambda_ip[ , 1], pass_1 = network$c_ip[ , 1], stringsAsFactors = F)
N_hat[[i]] <- N_hat[[i]] %>%
  dplyr::mutate(N_hat = ifelse(is.na(pass_1), lambda_hat, N_hat))