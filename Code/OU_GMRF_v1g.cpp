// Space time 
#include <TMB.hpp>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// Objective function
template<class Type>
Type objective_function<Type>::operator() ()
{
  
  // Load GMRF functions
  using namespace density;
  
  // Options
  DATA_FACTOR( Options_vec );
  // Slot 0 : Include spatial variation (0=No, 1=Yes)
  // Slot 1 : Include temporal variation (0=No, 1=Yes)
  // Slot 2 : Include spatiotemporal variation (0=No, 1=Yes)
  
  // Sizes
  DATA_INTEGER(n_i);   // Number of data points
  DATA_INTEGER(n_b);	// Number of branches in acyclic graph (b) (presumably n_b = n_d, i.e., including infinitely-long branch for root)
  DATA_INTEGER(n_t);  // Number of years

  // Data
  DATA_MATRIX(c_ip);       	// Count data
  DATA_FACTOR(d_i);         // Branch number b for each observationp
  DATA_MATRIX(X_ij);      // Covariate matrix for observation i and covariate j
  DATA_FACTOR(t_i);      // Covariate matrix for observation i and covariate j

  // Inputs regarding branched network
  DATA_FACTOR(parent_b); // index of parent
  DATA_FACTOR(child_b); // index of child for branch b
  DATA_VECTOR(dist_b);  // distance to parent
  
  // Fixed effects
  PARAMETER(log_theta);             // Autocorrelation (i.e. density dependence)
  PARAMETER(log_SD);
  PARAMETER(log_theta_st);         // Spatial correlation of spatiotemporal error
  PARAMETER(log_SD_st);         // Magnitude of Spatiotemporal error
  PARAMETER(rho_st);              // Temporal correlation of spatiomteporal error
  PARAMETER(log_sigmaIID);
  PARAMETER(log_mean);      
  PARAMETER(log_extradetectionSD);
  PARAMETER(rhot);                  // Purely temporal autocorrelation
  PARAMETER(log_sigmat);
  PARAMETER_VECTOR(gamma_j);
  PARAMETER(log_detectrate);

  // Random effects
  PARAMETER_VECTOR(log_extradetectrate_i);
  PARAMETER_VECTOR(lognormal_overdispersed_i);
  PARAMETER_VECTOR(Epsiloninput_d);  // Spatial process variation
  PARAMETER_VECTOR(Deltainput_t);
  PARAMETER_MATRIX(Nu_dt);

  // objective function -- joint negative log-likelihood 
  Type jnll = 0; 
  vector<Type> jnll_comp(6);
  jnll_comp.setZero();
  int n_j = X_ij.col(0).size();
  
  // Derived parameters
  Type SDinput = exp(log_SD);
  Type theta = exp(log_theta);
  Type SDinput_st = exp(log_SD_st);
  Type theta_st = exp(log_theta_st);
  Type detectrate = exp(log_detectrate);
  Type sigmaIID = exp(log_sigmaIID);
  vector<Type> extradetectrate_i(n_i);
  extradetectrate_i = exp(log_extradetectrate_i);
  Type extradetectionSD = exp(log_extradetectionSD);
  Type sigmat = exp(log_sigmat);
  vector<Type> Delta_t(n_t);
  Delta_t = Deltainput_t * sigmat;
  
  // Detection probability
  matrix<Type> detectprob_ip(n_i,3);
  for (int i=0; i<n_i; i++){
    detectprob_ip(i,0) = 1.0 - exp(-1 * (detectrate * extradetectrate_i(i)));
    detectprob_ip(i,1) = (1-detectprob_ip(i,0)) * (1.0 - exp(-1 * (detectrate * extradetectrate_i(i))));
    detectprob_ip(i,2) = (1-detectprob_ip(i,0)) * (1-detectprob_ip(i,1)) * (1.0 - exp(-1 * (detectrate * extradetectrate_i(i))));
  }  
  
  // Probability of GRF on network -- SPATIAL
  vector<Type> rho_b(n_b); 
  vector<Type> SD_b(n_b); 
  vector<Type> SDinput_b(n_b); 
  for(int b=0; b<n_b; b++){
    if( isNA(dist_b(b)) ){
      // Correlation between i and parent(i) as distance -> INF
      rho_b(b) = 0; 
      // SD of Ornstein-Uhlenbeck process as distance -> INF
      SDinput_b(b) = SDinput / pow(2*theta, 0.5);
      // conditional probability
      if(Options_vec(0)==1) jnll_comp(0) -= dnorm(Epsiloninput_d(child_b(b)), Type(0.0), SDinput_b(b), true);
    }
    if( !isNA(dist_b(b)) ){
      // Correlation between i and parent(i)
      rho_b(b) = exp(-theta * dist_b(b)); 
      // SD of O-U process
      SDinput_b(b) = pow( pow(SDinput,2)/(2*theta) * (1-exp(-2*theta*dist_b(b))), 0.5 );
      // conditional probability
      if(Options_vec(0)==1) jnll_comp(0) -= dnorm(Epsiloninput_d(child_b(b)), Epsiloninput_d(parent_b(b)), SDinput_b(b), true);
    }
  }
  
  // Probability of GRF on network -- SPATIOTEMPORAL
  vector<Type> rho_t_b(n_b); 
  vector<Type> SD_t_b(n_b); 
  vector<Type> SDinput_t_b(n_b); 
  vector<Type> temp_b(n_b);
  for(int b=0; b<n_b; b++){
    if( isNA(dist_b(b)) ){
      // Correlation between i and parent(i) as distance -> INF
      rho_t_b(b) = 0; 
      // SD of O-U process as distance -> INF
      SDinput_t_b(b) = SDinput_st / pow(2*theta_st, 0.5);
      // conditional probability
      temp_b = Nu_dt.row(child_b(b));
      if(Options_vec(2)==1) jnll_comp(2) += SCALE( AR1(rho_st), SDinput_t_b(b))(temp_b);
    }
    if( !isNA(dist_b(b)) ){
      // Correlation between i and parent(i)
      rho_t_b(b) = exp(-theta_st * dist_b(b)); 
      // SD of O-U process
      SDinput_t_b(b) = pow( pow(SDinput_st,2)/(2*theta_st) * (1-exp(-2*theta_st*dist_b(b))), 0.5 );
      // conditional probability
      temp_b = Nu_dt.row(child_b(b)) - rho_t_b(b)*Nu_dt.row(parent_b(b));
      if(Options_vec(2)==1) jnll_comp(2) += SCALE( AR1(rho_st), SDinput_t_b(b))(temp_b);
    }
  }
  
  // Probability of temporal variation
  if(Options_vec(1)==1) jnll_comp(1) += AR1(rhot)(Deltainput_t);
  
  // Probability of IID lognormal variation
  for (int i=0; i<n_i; i++){
    jnll_comp(5) -= dnorm( lognormal_overdispersed_i(i), Type(0.0), sigmaIID, true );
  }
  
  // Random variation in detection probability
  for (int i=0; i<n_i; i++){
    if(Options_vec(3)==1) jnll_comp(3) -= dnorm( log_extradetectrate_i(i), Type(0.0), extradetectionSD, true );
  }
  
  // Covariates
  vector<Type> eta_i(n_i);
  eta_i = X_ij * gamma_j.matrix();
  
  // Likelihood contribution from observations
  matrix<Type> lambda_ip(n_i,3);
  for (int i=0; i<n_i; i++){
  for (int p=0; p<3; p++){
    lambda_ip(i,p) = exp( log_mean + Epsiloninput_d(d_i(i)) + eta_i(i) + Delta_t(t_i(i)) + Nu_dt(d_i(i),t_i(i)) );
    if( !isNA(c_ip(i,p)) ){                
      if(Options_vec(4)==1) jnll_comp(4) -= dpois(c_ip(i,p), lambda_ip(i,p)*detectprob_ip(i,p)*exp(lognormal_overdispersed_i(i)), true);
    }
  }}

  // Add up components
  jnll = jnll_comp.sum();

  // Spatial field summaries
  REPORT( rho_b );
  REPORT( rho_st );
  REPORT( rho_t_b );
  REPORT( SD_b );
  REPORT( SDinput_b );
  REPORT( theta );
  REPORT( theta_st );
  REPORT( Epsiloninput_d );
  REPORT( log_mean );
  REPORT( detectprob_ip );
  REPORT( extradetectrate_i );
  REPORT( detectrate );
  REPORT( gamma_j );
  REPORT( eta_i);
  REPORT( lambda_ip );
  REPORT( Delta_t );
  REPORT( rhot );
  REPORT( sigmat );
  REPORT( jnll_comp );
  REPORT( jnll );
  REPORT( sigmaIID );
  REPORT( lognormal_overdispersed_i );
  REPORT( SDinput_t_b );
  REPORT( temp_b );
  REPORT( Nu_dt );
    
  ADREPORT( lambda_ip);
  
  return jnll;
}
