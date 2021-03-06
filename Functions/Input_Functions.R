rmatrix <- function( nrow=1, ncol=1, mean=0, sd=1, ... ){
  Return = matrix( rnorm(nrow*ncol,mean=mean,sd=sd), nrow=nrow, ncol=ncol, ...)
  return( Return )
}


# Make inputs
makeInput <- function(family, df = NULL, c_i = NULL, c_ip = NULL, options, X, t_i, version, CalcSD_lambda_ip, offset_i = NULL, spatial_equal = FALSE) {
  
  # included so don't have to redo for simulations
  if(is.null(df)) {
    df <- family
    if(nrow(df) != nrow(c_ip)) {
      stop("family must have the same number of rows as c_ip unless df specified")
    }
  }
  
  # convert 3-pass counts to abundance using Carle & Strub 1978 methods
  if(Version=="OU_GMRF_v1b") {
    #source("http://www.rforge.net/FSA/InstallFSA.R")
    if(!"FSA" %in% installed.packages()) {
      if (!require('devtools')) install.packages('devtools')
      require('devtools')
      devtools::install_github('droglenc/FSA')
    }
    library(FSA)
    c_obs <- dplyr::filter(c_ip, !is.na(pass_1))
    N_cs <- as.data.frame(matrix(NA, nrow(c_obs), 8))
    for(i in 1:dim(c_obs)[1]){
      foo <- FSA::removal(catch = as.vector(c_obs[i, ]), method = c("CarleStrub"), just.ests = TRUE)
      N_cs[i, ] <- foo
    } # end carle-strub for loop
    names(N_cs) <- names(foo)
    data.frame(c_obs, c_sum = rowSums(c_obs), N_cs)
    foo <- c_ip %>%
      dplyr::select(pass_1) %>%
      dplyr::filter(is.na(pass_1))
    c_i <- as.vector(dplyr::bind_rows(foo, dplyr::select(N_cs, No))$No)
  } # end if version b statement
  
  if(Version=="OU_GMRF_v1a") Data = list( "n_i"=length(c_i), "n_b"=nrow(family), "c_i"=c_i, "d_i"=family[,'child_b']-1, "parent_b"=family[,'parent_b']-1, "child_b"=family[ , 'child_b']-1, "dist_b"=family[,'dist_b'])
  if(Version=="OU_GMRF_v1b") Data = list( "n_i"=length(c_i), "n_b"=nrow(family), "c_i"=c_i, "d_i"=family[,'child_b']-1, "X_ij"=X_ij, "parent_b"=family[,'parent_b']-1, "child_b"=child_b-1, "dist_b"=family[,'dist_b'])
  if(Version%in%c("OU_GMRF_v1c","OU_GMRF_v1d")) Data = list( "n_i"=dim(c_ip)[1], "n_b"=nrow(family), "c_ip"=as.matrix(c_ip), "d_i"=df[,'child_b']-1, "X_ij"=X_ij, "parent_b"=family[ ,'parent_b']-1, "child_b"=family[ ,'child_b']-1, "dist_b"=family[,'dist_b'])
  if(Version%in%c("OU_GMRF_v1e")) Data = list( "n_i"=dim(c_ip)[1], "n_b"=nrow(family), "n_t"=length(YearSet), "c_ip"=as.matrix(c_ip), "d_i"=df[,'child_b']-1, "X_ij"=X_ij, "t_i"=t_i-min(t_i), "parent_b"=family[ ,'parent_b']-1, "child_b"=family[ ,'child_b']-1, "dist_b"=family[,'dist_b'])
  if(Version%in%c("OU_GMRF_v1g","OU_GMRF_v1f", "OU_GMRF_v1h", "OU_GMRF_v1i")) {
    YearSet = min(t_i):max(t_i)
    n_sd = length(Calc_lambda_ip[which(Calc_lambda_ip != 0)])
    if(is.null(offset_i)) {
      offset_i <- rep(1, length.out = nrow(c_ip))
      warning("No offsets included")
    }
    Data = list( "Options_vec"=options, "n_sd"=n_sd, "CalcSD_lambda_ip"=CalcSD_lambda_ip, "n_i"=dim(c_ip)[1], "n_b"=nrow(family), "n_t"=length(YearSet), "c_ip"=as.matrix(c_ip), "d_i"=df[,'child_b']-1, "X_ij"=X, "t_i"=t_i-min(t_i), "parent_b"=family[ ,'parent_b']-1, "child_b"=family[ ,'child_b']-1, "dist_b"=family[,'dist_b'], "offset_i"=offset_i) # d_i and child_b redundant?
  }
  
  ############### Sanity checks on inputs ##############
  # Please add more here
  
  # If any distances are negative or below threshold, fix at lower bound
  if( any(Data$dist_b < (max(Data$dist_b,na.rm=TRUE)/1e4),na.rm=TRUE) ) {
    Data$dist_b = ifelse( Data$dist_b<(max(Data$dist_b,na.rm=TRUE)/1e4), (max(Data$dist_b,na.rm=TRUE)/1e4), Data$dist_b)
    warning("Negative distances fixed to lower bound")
  }
  
  # check for correct dimensions
  if(nrow(Data$X_ij) != nrow(as.data.frame(Data$c_ip))) {
    stop("covariates X_ij must be of the length as the observed data c_ip")
  }
  
  # If any covariates are missing
  if( any(is.na(Data$X_ij)) ) {
    Data$X_ij = ifelse(is.na(Data$X_ij), 0, Data$X_ij)
    warning("Missing covariate values: replaced with mean")
  }
  
  # offset must be a vector not a dataframe and of length n_i
  if(class(Data$offset_i) != "numeric") {
    stop("offset must be class numeric")
  }
  
  if(length(Data$offset_i) != nrow(c_ip)) {
    stop("offset must be numeric vector of length == nrow(c_ip)")
  }
  
  # replace any missing offset values with the median
  if( any(is.na(Data$offset_i)) ) {
    Data$offset_i[is.na(Data$offset_i)] <- median(Data$offset_i, na.rm = TRUE)
    warning("Missing offsets replaced with median")
  }
  # check that all offsets are positive
  if(any(Data$offset_i <= 0)) {
    stop("All offsets need to be positive values")
  }

  
  
  ###### make dist min = 10 m
  #Data$dist_b = ifelse(Data$dist_b < 0.01, 0.01, Data$dist_b)
  ######
  # Check for Options_vec combos that don't make sense
  
  ######################################################
  
  if(Version=="OU_GMRF_v1a") Params = list( "log_theta"=log(1), "log_SD"=log(1), "log_mean"=log(1), "Epsiloninput_d"=rnorm(Data$n_b))
  if(Version=="OU_GMRF_v1b") Params = list( "log_theta"=log(1), "log_SD"=log(1), "log_mean"=log(1), "gamma_j"=rep(0,ncol(Data$X_ij)), "Epsiloninput_d"=rnorm(Data$n_b))
  if(Version=="OU_GMRF_v1c") Params = list( "log_theta"=log(1), "log_SD"=log(1), "log_mean"=log(1), "gamma_j"=rep(0,ncol(Data$X_ij)), "log_detectrate"=log(0.2), "log_extradetectrate_i"=log(rep(1,Data$n_i)), "Epsiloninput_d"=rnorm(Data$n_b))
  if(Version=="OU_GMRF_v1d") Params = list( "log_theta"=log(1), "log_SD"=log(1), "log_mean"=log(1), "log_extradetectionSD"=log(1), "gamma_j"=rep(0,ncol(Data$X_ij)), "log_detectrate"=log(0.2), "log_extradetectrate_i"=log(rep(1,Data$n_i)), "Epsiloninput_d"=rnorm(Data$n_b))
  if(Version=="OU_GMRF_v1e") Params = list( "log_theta"=log(1), "log_SD"=log(1), "log_mean"=log(1), "log_extradetectionSD"=log(1), "rhot"=0, "log_sigmat"=log(1), "gamma_j"=rep(0,ncol(Data$X_ij)), "log_detectrate"=log(0.2), "log_extradetectrate_i"=log(rep(1,Data$n_i)), "Epsiloninput_d"=rnorm(Data$n_b), "Deltainput_t"=rnorm(Data$n_t))
  if(Version=="OU_GMRF_v1f") Params = list( "log_theta"=log(1), "log_SD"=log(1), "log_theta_sp"=log(1), "log_SD_st"=log(1), "rho_sp"=0, "log_mean"=log(1), "log_extradetectionSD"=log(1), "rhot"=0, "log_sigmat"=log(1), "gamma_j"=rep(0,ncol(Data$X_ij)), "log_detectrate"=log(0.2), "log_extradetectrate_i"=rnorm(Data$n_i,sd=0.01), "Epsiloninput_d"=rnorm(Data$n_b,sd=0.01), "Deltainput_t"=rnorm(Data$n_t,sd=0.01), "Nu_dt"=rmatrix(Data$n_b,Data$n_t,sd=0.01))
  
  if(Version%in%c("OU_GMRF_v1h","OU_GMRF_v1g")) Params = list( "log_theta"=log(1), "log_SD"=log(1), "log_theta_st"=log(1), "log_SD_st"=log(1), "rho_st"=0, "log_sigmaIID"=log(1), "log_mean"=log(1), "log_extradetectionSD"=log(1), "rhot"=0, "log_sigmat"=log(1), "gamma_j"=rep(0,ncol(Data$X_ij)), "log_detectrate"=log(0.2), "log_extradetectrate_i"=rnorm(Data$n_i,sd=0.01), "lognormal_overdispersed_i"=rnorm(Data$n_i,sd=0.01), "Epsiloninput_d"=rnorm(Data$n_b,sd=0.01), "Deltainput_t"=rnorm(Data$n_t,sd=0.01), "Nu_dt"=rmatrix(Data$n_b,Data$n_t,sd=0.01))
  
  if(Version%in%c("OU_GMRF_v1i")) Params = list( "log_theta_vec"=rep(log(1),2), "log_SD"=log(1), "log_SD_st"=log(1), "rho_st"=0, "log_sigmaIID"=log(1), "log_mean"=log(1), "log_extradetectionSD"=log(1), "rhot"=0, "log_sigmat"=log(1), "gamma_j"=rep(0,ncol(Data$X_ij)), "log_detectrate"=log(0.2), "log_extradetectrate_i"=rnorm(Data$n_i,sd=0.01), "lognormal_overdispersed_i"=rnorm(Data$n_i,sd=0.01), "Epsiloninput_d"=rnorm(Data$n_b,sd=0.01), "Deltainput_t"=rnorm(Data$n_t,sd=0.01), "Nu_dt"=rmatrix(Data$n_b,Data$n_t,sd=0.01))
  
  if(Version%in%c("OU_GMRF_v1a","OU_GMRF_v1b")) Random = c( "Epsiloninput_d" )
  if(Version%in%c("OU_GMRF_v1c","OU_GMRF_v1d")) Random = c( "Epsiloninput_d", "log_extradetectrate_i" )
  if(Version%in%c("OU_GMRF_v1e")) Random = c( "Epsiloninput_d", "log_extradetectrate_i", "Deltainput_t" )
  if(Version%in%c("OU_GMRF_v1f")) Random = c( "Epsiloninput_d", "log_extradetectrate_i", "Deltainput_t", "Nu_dt" )
  if(Version%in%c("OU_GMRF_v1g", "OU_GMRF_v1h", "OU_GMRF_v1i")) Random = c( "Epsiloninput_d", "log_extradetectrate_i", "Deltainput_t", "Nu_dt", "lognormal_overdispersed_i" )
  if(Version%in%c("OU_GMRF_v1h", "OU_GMRF_v1i")){
    if( sum(unlist(Options_vec[c("SpatialTF", "TemporalTF", "SpatiotemporalTF", "DetectabilityTF", "OverdispersedTF")]))==0 ) Random = NULL
  }
  
  # Turn off random effects if desired
  Map = list()
  if( Version%in%c("OU_GMRF_v1h","OU_GMRF_v1g","OU_GMRF_v1f","OU_GMRF_v1e","OU_GMRF_v1d") & Options_vec[["SpatialTF"]]==FALSE ){
    Map[["log_theta"]] = factor(NA)
    Map[["log_SD"]] = factor(NA)
    Params[["Epsiloninput_d"]] = rep(0,length(Params[["Epsiloninput_d"]]))
    Map[["Epsiloninput_d"]] = factor( rep(NA,length(Params[["Epsiloninput_d"]])) )
  }
  if( Version%in%c("OU_GMRF_v1i") & Options_vec[["SpatialTF"]]==FALSE ){
    Map[["log_theta_vec"]] = factor(c(NA,1))
    Map[["log_SD"]] = factor(NA)
    Params[["Epsiloninput_d"]] = rep(0,length(Params[["Epsiloninput_d"]]))
    Map[["Epsiloninput_d"]] = factor( rep(NA,length(Params[["Epsiloninput_d"]])) )
  }
  if( Version%in%c("OU_GMRF_v1h", "OU_GMRF_v1i","OU_GMRF_v1g","OU_GMRF_v1f","OU_GMRF_v1e","OU_GMRF_v1d") & Options_vec[["TemporalTF"]]==FALSE ){
    Map[["rhot"]] = factor(NA)
    Map[["log_sigmat"]] = factor(NA)
    Params[["Deltainput_t"]] = rep(0,length(Params[["Deltainput_t"]]))
    Map[["Deltainput_t"]] = factor( rep(NA,length(Params[["Deltainput_t"]])) )
  }
  if( Version%in%c("OU_GMRF_v1h","OU_GMRF_v1g","OU_GMRF_v1f","OU_GMRF_v1e","OU_GMRF_v1d") & Options_vec[["SpatiotemporalTF"]]==FALSE ){
    Map[["log_theta_st"]] = factor(NA)
    Map[["log_SD_st"]] = factor(NA)
    Map[["rho_st"]] = factor(NA)
    Params[["Nu_dt"]] = array(0,dim(Params[["Nu_dt"]]))
    Map[["Nu_dt"]] = factor( array(NA,dim(Params[["Nu_dt"]])) )
  }
  if( Version%in%c("OU_GMRF_v1i") & Options_vec[["SpatiotemporalTF"]]==FALSE ){
    Map[["log_theta_vec"]] = factor(c(1,NA))
#    Map[["log_theta_st"]] = factor(NA)
    Map[["log_SD_st"]] = factor(NA)
    Map[["rho_st"]] = factor(NA)
    Params[["Nu_dt"]] = array(0,dim(Params[["Nu_dt"]]))
    Map[["Nu_dt"]] = factor( array(NA,dim(Params[["Nu_dt"]])) )
  }
  if( Version%in%c("OU_GMRF_v1h", "OU_GMRF_v1i","OU_GMRF_v1g","OU_GMRF_v1f","OU_GMRF_v1e","OU_GMRF_v1d") & Options_vec[["DetectabilityTF"]]==FALSE ){
    Map[["log_extradetectionSD"]] = factor(NA)
    Map[["log_extradetectrate_i"]] = factor( rep(NA,Data$n_i) )
    Params[["log_extradetectrate_i"]] = rep(0,Data$n_i)
  }
  if( Version%in%c("OU_GMRF_v1h", "OU_GMRF_v1i") & Options_vec[["OverdispersedTF"]]==FALSE ){
    Map[["log_sigmaIID"]] = factor(NA)
    Map[["lognormal_overdispersed_i"]] = factor( rep(NA,Data$n_i) )
    Params[["lognormal_overdispersed_i"]] = rep(0,Data$n_i)
  }
  
  if( Version%in%c("OU_GMRF_v1i") & Options_vec[["SpatialTF"]]==FALSE  & Options_vec[["SpatiotemporalTF"]]==FALSE ){
    Map[["log_theta_vec"]] = factor(c(NA,NA))
  }
  
  if( Version%in%c("OU_GMRF_v1i") & Options_vec[["SpatialTF"]]==TRUE  & Options_vec[["SpatiotemporalTF"]]==TRUE & spatial_equal == TRUE){
    Map[["log_theta_vec"]] = factor(c(1,1))
  }
  
  input_list <- list(Data = Data, Params = Params, Map = Map, Random = Random)
  
  return(input_list)
} # end input function
