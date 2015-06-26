

#setwd("C:/Users/James.Thorson/Desktop/UW Hideaway/Collaborations/2015 -- river network GMRF/exploratory data/")
#######################
# Load libraries
#######################
library(dplyr)
library(TMB)
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


#######################
# Fit in TMB
#######################

Version = "OU_GMRF_v1b"
# v1a -- Original version
# v1b -- added covariates matrix X_ij
#setwd( TmbFile )
if(FALSE){
  dyn.unload(dynlib(paste0("Code/", Version)))
  file.remove( paste0("Code/", Version,c(".o",".dll")) )
}
compile( paste0("Code/", Version,".cpp") )

# Make inputs
if(Version=="OU_GMRF_v1a") Data = list( "n_i"=length(c_i), "n_b"=nrow(family), "c_i"=c_i, "d_i"=family[,'child_b']-1, "parent_b"=family[,'parent_b']-1, "child_b"=family[,'child_b']-1, "dist_b"=family[,'dist_b'])
if(Version=="OU_GMRF_v1b") Data = list( "n_i"=length(c_i), "n_b"=nrow(family), "c_i"=c_i, "d_i"=family[,'child_b']-1, "X_ij"=X_ij, "parent_b"=family[,'parent_b']-1, "child_b"=family[,'child_b']-1, "dist_b"=family[,'dist_b'])
if(Version=="OU_GMRF_v1a") Params = list( "log_theta"=log(1), "log_SD"=log(1), "log_mean"=log(1), "Epsiloninput_d"=rnorm(Data$n_b))
if(Version=="OU_GMRF_v1b") Params = list( "log_theta"=log(1), "log_SD"=log(1), "log_mean"=log(1), "gamma_j"=rep(0,ncol(Data$X_ij)), "Epsiloninput_d"=rnorm(Data$n_b))
Random = c( "Epsiloninput_d" )
Map = NULL

# Make object
dyn.load( dynlib(paste0("Code/", Version) ))
obj <- MakeADFun(data=Data, parameters=Params, random=Random, map=Map, hessian=FALSE, inner.control=list(maxit=1000) )

# First run
obj$fn( obj$par )
# Check for parameters that don't do anything
which( obj$gr( obj$par )==0 )

# Run model
opt = nlminb(start=obj$env$last.par.best[-c(obj$env$random)], objective=obj$fn, gradient=obj$gr, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14) )
opt[["final_gradient"]] = obj$gr( opt$par )

# Get standard errors
Report = obj$report()
Sdreport = sdreport( obj )

# compare predicted vs. observed on original scale
c_est <- exp(Report[["Epsiloninput_d"]]+Report[["log_mean"]]+Report[["eta_i"]])
plot(c_i, c_est)
abline( a=0, b=1, lty="dotted")

