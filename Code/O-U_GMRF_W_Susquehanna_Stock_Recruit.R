# Run yoy model with previous year adult (stock-recruit) and vice versa

# limitation - no good time series so only estimates of abundance the previous year not accounting for overdispersion (random site-year effects which can be large for YOY)

rm(list = ls())
gc()

library(tidyr)
library(dplyr)
library(TMB)
library(minqa)
library(ggplot2)
source("Functions/Input_Functions.R")
source("Functions/simOUGMRF.R")
source("Functions/runOUGMRF.R")
source("Functions/summary_functions.R")

data_dir <- "Output"

yoy <- readRDS(file.path(data_dir, "W_Susquehanna_YOY_Summary.Rds"))
adult <- readRDS(file.path(data_dir, "W_Susquehanna_Summary.Rds"))

# problem of how to covert to match year and site
lambda_dt <- data.frame(child_name = adult$family$child_name, adult$Report5$lambda_dt, stringsAsFactors = FALSE)
names(lambda_dt) <- c("child_name", seq(from = 1982, to = 1982+ncol(lambda_dt)-2, by = 1))
lambda_adult <- tidyr::gather(lambda_dt, year, lambda_adult, -child_name) %>%
  distinct()
lambda_adult$child_name <- as.character(lambda_adult$child_name)
lambda_adult$year <- as.integer(as.character(lambda_adult$year))
lambda_adult$year_next <- lambda_adult$year + 1
adult$df <- left_join(adult$df, dplyr::select(lambda_adult, -year), by = c("child_name"="child_name", "year"="year_next"))
adult$df <- adult$df %>%
  dplyr::rename(adult_prev = lambda_adult)
adult$df <- adult$df[which(!is.na(adult$df$adult_prev)), ]

lambda_dt <- data.frame(child_name = yoy$family$child_name, yoy$Report5$lambda_dt, stringsAsFactors = FALSE)
names(lambda_dt) <- c("child_name", seq(from = 1982, to = 1982+ncol(lambda_dt)-2, by = 1))
lambda_yoy <- tidyr::gather(lambda_dt, year, lambda_yoy, -child_name) %>%
  distinct()
lambda_yoy$child_name <- as.character(lambda_yoy$child_name)
lambda_yoy$year <- as.integer(as.character(lambda_yoy$year))
lambda_yoy$year_next <- lambda_yoy$year + 1
yoy$df_yoy <- left_join(yoy$df_yoy, dplyr::select(lambda_yoy, -year), by = c("child_name"="child_name", "year"="year_next"))
yoy$df_yoy <- yoy$df_yoy %>%
  dplyr::rename(yoy_prev = lambda_yoy)
yoy$df_yoy <- yoy$df_yoy[which(!is.na(yoy$df_yoy$yoy_prev)), ]

covs <- dplyr::select(adult$df, child_name, year, forest, surfcoarse, temp_mean_fall_1, temp_mean_winter, temp_mean_spring, prcp_mean_fall_1, prcp_mean_winter, prcp_mean_spring, adult_prev)
covs <- covs %>%
  dplyr::left_join(dplyr::select(yoy$df_yoy, child_name, year, yoy_prev))


offset <-  as.numeric(yoy$df_yoy$length_sample)

Version = "OU_GMRF_v1i"

# Compile
if(FALSE) {
  dyn.unload(dynlib(paste0("Code/", Version)))
  file.remove( paste0("Code/", Version,c(".o",".dll")) )
}
compile( paste0("Code/", Version,".cpp") )

# initial site-years to get SD for lambda (cycle through just for best model)
Calc_lambda_ip <- rep(0, length.out = nrow(adult$df))


#######################
# YOY -> Adults (recruitment)
#######################

#----------------- Temporal + Spatiotemporal ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=0, "TemporalTF"=1, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1, "abundTF"=0)

# covariates
X_ij <- as.matrix(dplyr::select(covs, -child_name, -year, -adult_prev, -yoy_prev))

# Make inputs
Inputs <- makeInput(family = adult$family, df = adult$df, c_ip = as.matrix(dplyr::select(adult$df, starts_with("pass"))), options = Options_vec, X = X_ij, t_i = adult$df$year, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

adult_t_st <- runOUGMRF(inputs = Inputs)

#----------------- Spatiotemporal Only ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=0, "TemporalTF"=0, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1, "abundTF"=0)

# covariates
X_ij <- as.matrix(dplyr::select(covs, -child_name, -year, -adult_prev, -yoy_prev))

# Make inputs
Inputs <- makeInput(family = adult$family, df = adult$df, c_ip = as.matrix(dplyr::select(adult$df, starts_with("pass"))), options = Options_vec, X = X_ij, t_i = adult$df$year, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

adult_st <- runOUGMRF(inputs = Inputs)

#----------------- Spatiotemporal Only ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=0, "TemporalTF"=0, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1, "abundTF"=0)

# covariates
X_ij <- as.matrix(dplyr::select(covs, -child_name, -year, -adult_prev))

# Make inputs
Inputs <- makeInput(family = adult$family, df = adult$df, c_ip = as.matrix(dplyr::select(adult$df, starts_with("pass"))), options = Options_vec, X = X_ij, t_i = adult$df$year, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

recruit <- runOUGMRF(inputs = Inputs)

str(recruit) # fails

# try building up
Options_vec = c("SpatialTF"=0, "TemporalTF"=0, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1, "abundTF"=0)

# covariates
X_ij <- as.matrix(dplyr::select(covs, forest, yoy_prev))

# Make inputs
Inputs <- makeInput(family = adult$family, df = adult$df, c_ip = as.matrix(dplyr::select(adult$df, starts_with("pass"))), options = Options_vec, X = X_ij, t_i = adult$df$year, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

recruit_st <- runOUGMRF(inputs = Inputs)

str(recruit_st)

#----------------- Temporal Spatiotemporal Recruit ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=0, "TemporalTF"=1, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1, "abundTF"=0)

# covariates
X_ij <- as.matrix(dplyr::select(covs, -child_name, -year, -adult_prev))

# Make inputs
Inputs <- makeInput(family = adult$family, df = adult$df, c_ip = as.matrix(dplyr::select(adult$df, starts_with("pass"))), options = Options_vec, X = X_ij, t_i = adult$df$year, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

recruit_t_st <- runOUGMRF(inputs = Inputs)

str(recruit_t_st) # fails

# try building up
Options_vec = c("SpatialTF"=0, "TemporalTF"=1, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1, "abundTF"=0)

# covariates
X_ij <- as.matrix(dplyr::select(covs, forest, yoy_prev))

# Make inputs
Inputs <- makeInput(family = adult$family, df = adult$df, c_ip = as.matrix(dplyr::select(adult$df, starts_with("pass"))), options = Options_vec, X = X_ij, t_i = adult$df$year, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

recruit_simple_t_st <- runOUGMRF(inputs = Inputs)


# check for convergence
adult_aic <- data.frame(model = c("t-st",
                     "st",
                     "st-recruit",
                     "t-st-recruit",
                     "t-st-recuit-simple"),
           converge = c(adult_t_st$opt$convergence, 
             adult_st$opt$convergence,
             recruit$opt$convergence,
             #recruit_simple$opt$convergence,
             recruit_t_st$opt$convergence,
             recruit_simple_t_st$opt$convergence),
           sd_max = c(max(adult_t_st$SD$sd),
             max(adult_st$SD$sd),
             max(recruit$SD$sd),
            # max(recruit_simple$SD$sd),
             max(recruit_t_st$SD$sd),
             max(recruit_simple_t_st$SD$sd)),
           AIC = c(adult_t_st$opt$AIC,
                   adult_st$opt$AIC,
                   recruit$opt$AIC,
                 #  recruit_simple$opt$AIC,
                   recruit_t_st$opt$AIC,
                   recruit_simple_t_st$opt$AIC)
) %>%
  dplyr::arrange(AIC)



# AIC comparison

data.frame(model = c("climate",
                     "recruit"),
           AIC = c(adult_st$opt$AIC,
                   recruit_simple$opt$AIC))

#######################
# Adults -> YOY (Stock-Recruit)
#######################

#----------------- Spatiotemporal Only ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=0, "TemporalTF"=0, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1, "abundTF"=0)

# covariates
X_ij <- as.matrix(dplyr::select(covs, forest, temp_mean_fall_1, temp_mean_winter, temp_mean_spring, prcp_mean_fall_1, prcp_mean_winter))

# Make inputs
Inputs <- makeInput(family = yoy$family, df = yoy$df_yoy, c_ip = as.matrix(dplyr::select(yoy$df_yoy, starts_with("pass"))), options = Options_vec, X = X_ij, t_i = yoy$df_yoy$year, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

yoy_st <- runOUGMRF(inputs = Inputs)

#----------------- Spatiotemporal Only ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=0, "TemporalTF"=0, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1, "abundTF"=0)

# covariates
X_ij <- as.matrix(dplyr::select(covs, forest, temp_mean_fall_1, temp_mean_winter, temp_mean_spring, prcp_mean_fall_1, prcp_mean_winter, adult_prev))

# Make inputs
Inputs <- makeInput(family = yoy$family, df = yoy$df_yoy, c_ip = as.matrix(dplyr::select(yoy$df_yoy, starts_with("pass"))), options = Options_vec, X = X_ij, t_i = yoy$df_yoy$year, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

SR <- runOUGMRF(inputs = Inputs)


#----------------- Temporal + Spatiotemporal Only ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=0, "TemporalTF"=1, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1, "abundTF"=0)

# covariates
X_ij <- as.matrix(dplyr::select(covs, forest, temp_mean_fall_1, temp_mean_winter, temp_mean_spring, prcp_mean_fall_1, prcp_mean_winter))

# Make inputs
Inputs <- makeInput(family = yoy$family, df = yoy$df_yoy, c_ip = as.matrix(dplyr::select(yoy$df_yoy, starts_with("pass"))), options = Options_vec, X = X_ij, t_i = yoy$df_yoy$year, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

yoy_t_st <- runOUGMRF(inputs = Inputs)

#----------------- Temporal + Spatiotemporal Only ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=0, "TemporalTF"=1, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1, "abundTF"=0)

# covariates
X_ij <- as.matrix(dplyr::select(covs, forest, temp_mean_fall_1, temp_mean_winter, temp_mean_spring, prcp_mean_fall_1, prcp_mean_winter, adult_prev))

# Make inputs
Inputs <- makeInput(family = yoy$family, df = yoy$df_yoy, c_ip = as.matrix(dplyr::select(yoy$df_yoy, starts_with("pass"))), options = Options_vec, X = X_ij, t_i = yoy$df_yoy$year, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

SR_t_st <- runOUGMRF(inputs = Inputs)



# check for convergence
yoy_aic <- data.frame(model = c("yoy_st",
                     "yoy_t_st",
                     "stock_recruit_st",
                     "stock_recruit_t_st"),
           converge = c(yoy_st$opt$convergence,
                        yoy_t_st$opt$convergence,
             SR$opt$convergence,
             SR_t_st$opt$convergence),
           SD = c(max(yoy_st$SD$sd),
                  max(yoy_t_st$SD$sd),
             max(SR$SD$sd),
             max(SR_t_st$SD$sd)),
           AIC = c(yoy_st$opt$AIC,
                   yoy_t_st$opt$AIC,
                   SR$opt$AIC,
                   SR_t_st$opt$AIC)
) %>%
  dplyr::arrange(AIC)


yoy_aic
adult_aic

save.image("Output/W_Susquehanna_Stock_Recruit.RData")

# compare predictions from stock-recruit and base model

