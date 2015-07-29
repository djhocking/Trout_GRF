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
  
  // Sizes
  DATA_INTEGER(n_i);   // Number of data points
  DATA_INTEGER(n_b);	// Number of branches in acyclic graph (b) (presumably n_b = n_d, i.e., including infinitely-long branch for root)

  // Data
  DATA_MATRIX(c_ip);       	// Count data
  DATA_FACTOR(d_i);         // Branch number b for each observationp
  DATA_MATRIX(X_ij);      // Covariate matrix for observation i and covariate j

  // Inputs regarding branched network
  DATA_FACTOR(parent_b); // index of parent
  DATA_FACTOR(child_b); // index of child for branch b
  DATA_VECTOR(dist_b);  // distance to parent
  
  // Fixed effects
  PARAMETER(log_theta);             // Autocorrelation (i.e. density dependence)
  PARAMETER(log_SD);
  PARAMETER(log_mean);      
  PARAMETER_VECTOR(gamma_j);
  PARAMETER(log_detectrate);
  PARAMETER_VECTOR(log_extradetectrate_i);

  // Random effects
  PARAMETER_VECTOR(Epsiloninput_d);  // Spatial process variation

  // objective function -- joint negative log-likelihood 
  Type jnll = 0; 
  int n_j = X_ij.col(0).size();
  
  // Derived parameters
  Type SDinput = exp(log_SD);
  Type theta = exp(log_theta);
  Type detectrate = exp(log_detectrate);
  vector<Type> extradetectrate_i(n_i);
  extradetectrate_i = exp(log_extradetectrate_i);
  
  // Detection probability
  matrix<Type> detectprob_ip(n_i,3);
  for (int i=0; i<n_i; i++){
    detectprob_ip(i,0) = 1.0 - exp(-1 * (detectrate * extradetectrate_i(i)));
    detectprob_ip(i,1) = (1-detectprob_ip(i,0)) * (1.0 - exp(-1 * (detectrate * extradetectrate_i(i))));
    detectprob_ip(i,2) = (1-detectprob_ip(i,0))*(1-detectprob_ip(i,1)) * (1.0 - exp(-1 * (detectrate * extradetectrate_i(i))));
  }  
  
  // Probability of GRF on network
  vector<Type> rho_b(n_b); 
  vector<Type> SD_b(n_b); 
  vector<Type> SDinput_b(n_b); 
  for(int b=0; b<n_b; b++){
    if( isNA(dist_b(b)) ){
      // Correlation between i and parent(i) as distance -> INF
      rho_b(b) = 0; 
      // SD of O-U process as distance -> INF
      SDinput_b(b) = SDinput / pow(2*theta, 0.5);
      // conditional probability
      jnll -= dnorm(Epsiloninput_d(child_b(b)), Type(0.0), SDinput_b(b), true);
    }
    if( !isNA(dist_b(b)) ){
      // Correlation between i and parent(i)
      rho_b(b) = exp(-theta * dist_b(b)); 
      // SD of O-U process
      SDinput_b(b) = pow( pow(SDinput,2)/(2*theta) * (1-exp(-2*theta*dist_b(b))), 0.5 );
      // conditional probability
      jnll -= dnorm(Epsiloninput_d(child_b(b)), Epsiloninput_d(parent_b(b)), SDinput_b(b), true);
    }
  }
  
  // Covariates
  vector<Type> eta_i(n_i);
  eta_i = X_ij * gamma_j.matrix();
  
  // Likelihood contribution from observations
  matrix<Type> lambda_ip(n_i,3);
  for (int i=0; i<n_i; i++){
  for (int p=0; p<3; p++){
    lambda_ip(i,p) = exp(log_mean + Epsiloninput_d(d_i(i)) + eta_i(i));
    if( !isNA(c_ip(i,p)) ){                
      jnll -= dpois(c_ip(i,p), lambda_ip(i,p)*detectprob_ip(i,p), true);
    }
  }}

  // Spatial field summaries
  REPORT( rho_b );
  REPORT( SD_b );
  REPORT( SDinput_b );
  REPORT( theta );
  REPORT( Epsiloninput_d );
  REPORT( log_mean );
  REPORT( detectprob_ip );
  REPORT( extradetectrate_i );
  REPORT( detectrate );
  REPORT( gamma_j );
  REPORT( eta_i);
  
  return jnll;
}
