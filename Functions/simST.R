# Generate Data
simST <- function(family, theta = 0.25, SD = 0.1, rhot = 0.5, SD_t = 1, theta_st = 0.5, SD_st = 0.05, mean_N = 10, n_years = 10, rho = 0.8, gamma_j, X_ij, p = c(0.75, 075, 0.75), spatial = TRUE, temporal = TRUE, spatiotemporal = TRUE) {#, sample_pct = NULL, sample_n = NULL, sample_years = NULL) {
  
  # load( "sample.RData")
  # colnames(family)[1] = "child_name"
  # family = cbind( family, "child_b"=1:nrow(family) )
  
  library( mvtnorm )
  
  #   # checks
  #   if(is.null(sample_years)) {
  #     sample_years <- n_years
  #   }
  #   if(sample_years > n_years) {
  #     stop("sample years must be less than or equal to n_years")
  #   }
  
  # ------------- Deterministic --------------------
  eta_i <- gamma_j * X_ij
  eta_i <- rowSums(eta_i)
  log_mean <- log(mean_N)
  
  #------------ Spatial -------------
  if(spatial) {
  # object
  rho_b = condSD_b = x_b = rep(NA, nrow(family))
  
  # seed at top of network
  WhichRoot = which( is.na(family[,'parent_b']) )
  condSD_b[WhichRoot] = sqrt( SD^2 / 2*theta )
  x_b[WhichRoot] = rnorm(1, mean=0, sd=condSD_b[WhichRoot])
  rho_b[WhichRoot] = 0

  # Loop through network
  while( TRUE ){
    for(i in 1:nrow(family)){
      if( is.na(x_b[i]) ){
        SimulatedNodes = which(!is.na(x_b))
        Match = match( family[i,'parent_b'], SimulatedNodes ) # Which
        if(length(Match)==1){
          condSD_b[i] = sqrt( SD^2/(2*theta) * (1-exp(-2*theta*family[i,'dist_b'])) )
          rho_b[i] = exp(-theta * family[i,'dist_b']);
          x_b[i] = rho_b[i]*x_b[SimulatedNodes[Match]] + rnorm(1, mean=0, sd=condSD_b[i])
        }
      }
    }
    # Stopping condition
    if( all(!is.na(x_b)) ) break()
  }
  } else {
    x_b = rep(0, nrow(family))
  }
  
  #---------------- Temporal Autocorrelation AR1 ----------
  if(temporal) {
  x_t <- c(0, rep(NA, times = n_years - 1))
  for(t in 2:n_years) {  
    x_t[t] <- rhot * x_t[t-1] + rnorm(1, 0, SD_t)
  }
  } else {
    x_t <- c(0, rep(0, times = n_years - 1))
  }
  
  #-------------- Spatiotemporal ------------------
  if(spatiotemporal) {
  # Covariance for a given site among years
  Corr_tt = rho ^ outer( 1:n_years, 1:n_years, FUN=function(a,b){abs(a-b)} )
  
  # object
  condSD_b = rep(NA, nrow(family))
  rho_b = rep(NA, nrow(family))
  x_bt = matrix(NA, nrow=nrow(family), ncol=n_years)
  
  # seed at top of network
  WhichRoot = which( is.na(family[,'parent_b']) )
  condSD_b[WhichRoot] = sqrt( SD_st^2 / 2*theta_st )
  x_bt[WhichRoot,] = rmvnorm(1, mean=rep(0,n_years), sigma=condSD_b[WhichRoot]^2 * Corr_tt )[1,]
  rho_b[WhichRoot] = 0

  # Loop through network
  while( TRUE ){
    for(i in 1:nrow(family)){
      if( is.na(x_bt[i,1]) ){
        SimulatedNodes = which(!is.na(x_bt[,1]))
        Match = match( family[i,'parent_b'], SimulatedNodes ) # Which
        if(length(Match)==1){
          condSD_b[i] = sqrt( SD_st^2/(2*theta_st) * (1-exp(-2*theta_st*family[i,'dist_b'])) )
          rho_b[i] = exp(-theta_st * family[i,'dist_b']);
          x_bt[i,] = rho_b[i]*x_bt[SimulatedNodes[Match],] + rmvnorm(1, mean=rep(0,n_years), sigma=condSD_b[i]^2 * Corr_tt )[1,]   # Cov_matrix = pointwise_variance * Corr_matrix
        }
      }
    }
    # Stopping condition
    if( all(!is.na(x_bt)) ) break()
  }
  } else {
    x_bt = matrix(0, nrow=nrow(family), ncol=n_years)
  }
  # Abundance
  log_Npred_bt <- log_mean + eta_i
  if(spatial) {
    log_Npred_bt <- log_Npred_bt + outer(x_b,rep(1,n_years))
  }
  if(temporal) {
    log_Npred_bt <- log_Npred_bt + outer(rep(1,nrow(family)),x_t)
  }
  if(spatiotemporal) {
    log_Npred_bt = log_Npred_bt + x_bt
  }
  N_i = rpois( prod(dim(x_bt)), lambda=exp(log_Npred_bt)) # organizes so i is ordered by site then year(site 1, 2, 3 for year 1, then site 1,2,3 for year 2)
  
  t_i <- rep(1:n_years, each = nrow(x_bt))
  
  #-------------- Observation Process --------------
  # Simulate binomial observation (count) process
  c_ip_full <- as.data.frame(matrix(NA, length(N_i), length(p)))
  for(a in 1:length(N_i)) {
    removals <- 0
    c_ip_full[a,1] <- rbinom(1, N_i[a], p[1])
    for(b in 2:length(p)) {
      removals <- removals + c_ip_full[a,b-1]
      c_ip_full[a,b] <- rbinom(1, N_i[a] - removals, p[b])
    }
  }
  
  return(list(x_bt = x_bt, x_t=x_t, x_b=x_b, N_i = N_i, c_ip = c_ip_full, t_i = t_i, log_Npred_bt = log_Npred_bt))
}
