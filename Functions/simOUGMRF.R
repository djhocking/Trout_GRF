
# Generate Data
simOUGMRF <- function(family, theta, SD, mean_N, gamma, X_ij, p, sample_pct, spatial = TRUE) {

  #------------ Spatial -------------
  # object
  condSD_b = x_b = rep(NA, nrow(family))
  
  # seed at top of network
  WhichRoot = which( is.na(family[,'parent_b']) )
  condSD_b[WhichRoot] = sqrt( SD^2 / 2*theta )
  x_b[WhichRoot] = rnorm(1, mean=0, sd=condSD_b[WhichRoot])
  
  # Loop through network
  while( TRUE ){
    for(i in 1:nrow(family)){
      if( is.na(x_b[i]) ){
        SimulatedNodes = which(!is.na(x_b))
        Match = match( family[i,'parent_b'], SimulatedNodes ) # Which
        if(length(Match)==1){
          condSD_b[i] = sqrt( SD^2/(2*theta) * (1-exp(-2*theta*family[i,'dist_b'])) )
          x_b[i] = x_b[SimulatedNodes[Match]] + rnorm(1, mean=0, sd=condSD_b[i])
        }
      }
    }
    # Stopping condition
    if( all(!is.na(x_b)) ) break()
  }
  
  #-------------- Spatiotemporal -------------------
  
  # ------------- Deterministic --------------------
  # Covariates
  eta_i <- gamma_j * X_ij
  eta_i <- rowSums(eta_i)
  log_mean <- log(mean_N)
  if(spatial) {
    N_i = rpois( length(x_b), lambda=exp(x_b+log_mean+eta_i))
  } else {
    N_i = rpois( length(x_b), lambda=exp(log_mean+eta_i))
  }
  
  #-------------- Observation Process --------------
  # Simulate binomial observation (count) process
  c_ip_full <- matrix(NA, length(N_i), length(p))
  for(a in 1:length(N_i)) {
    for(b in 1:length(p)) {
      c_ip_full[a,b] <- rbinom(1, N_i[a], p[b])
    }
  }
  
  #---------------- Temporal Autocorrelation AR1 ----------
  # temporal - no temporal variability currently
  t_i <- rep(2000, times = length(N_i))
  
  # thin observations at random (to answer how does sampling density within network affect ability to estimate theta)
  remove_pct <- 1 - sample_pct
  rows <- 1:nrow(c_ip_full)
  remove_rows <- sample(rows, size = trunc(length(rows)*remove_pct), replace = FALSE)
  c_ip <- c_ip_full
  c_ip[remove_rows, ] <- NA
  
  return(list(x_b = x_b, N_i = N_i, c_ip = c_ip, t_i = t_i))
}