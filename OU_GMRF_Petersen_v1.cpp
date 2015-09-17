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

  // Data
  DATA_VECTOR(c1_i);       	// Count data for first capture pass
  DATA_VECTOR(c2_i);         // Count data for second capture pass
  DATA_VECTOR(r2_i);      // Count data of recaptures in second pass
  
  // Fixed effects
  PARAMETER(log_theta);             // Autocorrelation (i.e. density dependence)
  PARAMETER(log_SD);
  PARAMETER(log_mean);      
  PARAMETER(log_extradetectionSD);
  PARAMETER_VECTOR(gamma_j);
  PARAMETER(log_detectrate);
  PARAMETER_VECTOR(log_extradetectrate_i);

  // Random effects
  PARAMETER_VECTOR(Epsiloninput_d);  // Spatial process variation

  // objective function -- joint negative log-likelihood 
  Type jnll = 0; 
  
  // Derived parameters
  Type SDinput = exp(log_SD);
  Type detectrate = exp(log_detectrate);
  vector<Type> extradetectrate_i(n_i);
  extradetectrate_i = exp(log_extradetectrate_i);
  Type extradetectionSD = exp(log_extradetectionSD);
  
  // Detection probability
  vector<Type> detectprob_i(n_i);
  for (int i=0; i<n_i; i++){
    detectprob_ip(i,0) = 1.0 - exp(-1 * (detectrate * extradetectrate_i(i)));
  }  
  
    // Detection probability
  matrix<Type> detectprob_ip(n_i,3);
  for (int i=0; i<n_i; i++){
    detectprob_ip(i,0) = 1.0 - exp(-1 * (detectrate * extradetectrate_i(i)));
    detectprob_ip(i,1) = (1-detectprob_ip(i,0)) * (1.0 - exp(-1 * (detectrate * extradetectrate_i(i))));
    detectprob_ip(i,2) = (1-detectprob_ip(i,0)-detectprob_ip(i,1)) * (1.0 - exp(-1 * (detectrate * extradetectrate_i(i))));
  }  
  
    // Detection probability
  matrix<Type> detectprob_ip(n_i,3);
  for (int i=0; i<n_i; i++){
    detectprob_ip(i,0) = p2 * (1-p1); // history = 01
    detectprob_ip(i,1) = p1;          // history = 10
    detectprob_ip(i,2) = p1 * p2;     // history = 11
  }  
  
  // Random variation in detection probability
  for (int i=0; i<n_i; i++){
    jnll -= dnorm( log_extradetectrate_i(i), Type(0.0), extradetectionSD, true );
  }
  
  // Covariates
  
  // Likelihood contribution from observations
  vector<Type> lambda_i(n_i);
  for (int i=0; i<n_i; i++){
    lambda_ip(i) = exp(log_mean + Epsiloninput_d(d_i(i)) + eta_i(i));
    if( !isNA(r2_i(i)) ){                
      jnll -= dpois(c_i(i), lambda_ip(i)*detectprob_ip(i), true);
    }
  }

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
  REPORT( lambda_ip );
  ADREPORT( lambda_ip);
  
  return jnll;
}
