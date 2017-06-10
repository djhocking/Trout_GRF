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
  DATA_VECTOR(CalcSD_lambda_ip);
  DATA_FACTOR(Options_vec);
  // Slot 0 : Include spatial variation (0=No, 1=Yes)
  // Slot 1 : Include temporal variation (0=No, 1=Yes)
  // Slot 2 : Include spatiotemporal variation (0=No, 1=Yes)
  // Slot 3 : Include log_extradetectrate_i variation (0=No, 1=Yes)
  // Slot 4 : Observation model (1=Poisson)
  // Slot 5 : Include lognormal_overdispersed_i (0=No, 1=Yes)
  // Slot 6 : Output SD of abundancem N_i (0=No, 1=Yes) - only use Yes for small datasets
  
  // Sizes
  DATA_INTEGER(n_i);   // Number of data points
  DATA_INTEGER(n_b);  // Number of branches in acyclic graph (b) (presumably n_b = n_d, i.e., including infinitely-long branch for root)
  DATA_INTEGER(n_t);  // Number of years
  DATA_INTEGER(n_sd); // Number of points to estimate SD
  
  // Data
  DATA_MATRIX(c_ip);         // Count data
  DATA_FACTOR(d_i);         // Branch number b for each observationp
  DATA_MATRIX(X_ij);      // Covariate matrix for observation i and covariate j
  DATA_FACTOR(t_i);      // Covariate matrix for observation i and covariate j
  DATA_VECTOR(offset_i); // Length or area of survey i to offset abundance and make comparable among site-visits
  
  // Inputs regarding branched network
  DATA_FACTOR(parent_b); // index of parent
  DATA_FACTOR(child_b); // index of child for branch b
  DATA_VECTOR(dist_b);  // distance to parent
  
  // Fixed effects
 // PARAMETER(log_theta);             // Autocorrelation (i.e. density dependence)
  PARAMETER(log_SD);
//  PARAMETER(log_theta_st);         // Spatial correlation of spatiotemporal error
  PARAMETER(log_SD_st);         // Magnitude of Spatiotemporal error
  PARAMETER(rho_st);              // Temporal correlation of spatiomteporal error
  PARAMETER(log_sigmaIID);
  PARAMETER(log_mean);      
  PARAMETER(log_extradetectionSD);
  PARAMETER(rhot);                  // Purely temporal autocorrelation
  PARAMETER(log_sigmat);
  PARAMETER_VECTOR(gamma_j);
  PARAMETER(log_detectrate);
  PARAMETER_VECTOR(log_theta_vec);
  
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
  Type theta = exp(log_theta_vec(0));
  Type SDinput_st = exp(log_SD_st);
  Type theta_st = exp(log_theta_vec(1));
  Type detectrate = exp(log_detectrate);
  Type sigmaIID = exp(log_sigmaIID);
  Type mu = exp(log_mean);
  vector<Type> extradetectrate_i(n_i);
  extradetectrate_i = exp(log_extradetectrate_i);
  Type extradetectionSD = exp(log_extradetectionSD);
  Type sigmat = exp(log_sigmat);
  vector<Type> Delta_t(n_t);
  Delta_t = Deltainput_t * sigmat;
  vector<Type> log_offset_i(n_i);
  log_offset_i = log(offset_i);
  Type SD_inf = SDinput / pow(2*theta, 0.5);
  Type SD_st_inf = SDinput_st / pow(2*theta_st, 0.5);
  
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
    if(Options_vec(5)==1) {
      jnll_comp(5) -= dnorm( lognormal_overdispersed_i(i), Type(0.0), sigmaIID, true );
    }
  }
  
  // Detection probability
  matrix<Type> detectprob_ip(n_i,3);
  for (int i=0; i<n_i; i++){
    detectprob_ip(i,0) = 1.0 - exp(-1 * (detectrate * extradetectrate_i(i)));
    detectprob_ip(i,1) = (1-detectprob_ip(i,0)) * (1.0 - exp(-1 * (detectrate * extradetectrate_i(i))));
    detectprob_ip(i,2) = (1-detectprob_ip(i,0)-detectprob_ip(i,1)) * (1.0 - exp(-1 * (detectrate * extradetectrate_i(i))));
  }  
  
  // Random variation in detection probability
  for (int i=0; i<n_i; i++){
    if(Options_vec(3)==1) jnll_comp(3) -= dnorm( log_extradetectrate_i(i), Type(0.0), extradetectionSD, true );
  }
  
  // Covariates
  vector<Type> eta_i(n_i);
  eta_i = X_ij * gamma_j.matrix(); // + log(offset_i) where offset is unstandardized length or area or it could be converted to units that make sense fish per meter vs. fish per 100 m
  
  // log-predictions
  matrix<Type> lambda_dt(n_b,n_t);
  //    matrix<Type> density_dt(n_b,n_t);
  //    matrix<Type> log_N100_dt(5,n_t); // first 5 sites with data
  int counter_d = 0;
  for(int d=0; d<n_b; d++){
    for(int t=0; t<n_t; t++){  
      lambda_dt(d,t) = exp( log_mean + Epsiloninput_d(d) + Delta_t(t) + Nu_dt(d,t) ) * 100; // * 100 for offset but should bring in as set value based on whatever the relative offset is set at (offset_denom)
      //      if( !isNA(c_ip(i,1)) ){  
      //        density_dt(d,t) = lambda_dt(d,t)*exp(lognormal_overdispersed_i(i));
      //        if( counter_d < 5 ){
      //          log_N100_dt(counter_d,t) = log(density_dt(d,t) * 100); // total hack - get SD for density of first 5 sites with observed data estimated in all years
      //          counter_d = counter_d + 1;
      //        }
      //      }
      //      else {
      //        density_dt(d,t) = lambda_dt(d,t);
      //      }
      }
    }    
 
 vector<Type> lambda_t(n_t);
  lambda_t(n_t) = 0;
 for(int t=0; t<n_t; t++){  
   for(int d=0; d<n_b; d++){
     //lambda_d_sum += lambda_dt(d, t);
     lambda_t(t) += lambda_dt(d,t) / n_b;
   }
  //lambda_t(t) = lambda_d_sum / n_b;
 }
  
  // Likelihood contribution from observations
  matrix<Type> lambda_ip(n_i,3);
  matrix<Type> N_ip(n_i,3);
  matrix<Type> chat_ip(n_i,3);
  for (int i=0; i<n_i; i++){
    for (int p=0; p<3; p++){
      lambda_ip(i,p) = exp( log_mean + Epsiloninput_d(d_i(i)) + eta_i(i) + log_offset_i(i) + Delta_t(t_i(i)) + Nu_dt(d_i(i),t_i(i)) );
      if( !isNA(c_ip(i,p)) ){                
        if(Options_vec(4)==1) {
          jnll_comp(4) -= dpois(c_ip(i,p), lambda_ip(i,p)*exp(lognormal_overdispersed_i(i))*detectprob_ip(i,p), true);
        }
      }
      N_ip(i,p) = lambda_ip(i,p)*exp(lognormal_overdispersed_i(i));
      chat_ip(i,p) = N_ip(i,p)*detectprob_ip(i,p);
    }}
  
  // predictions with SD - consider changing this for density in each year (density_dt - converted to fish per 100m)
  vector<Type> SD_report(n_sd);
  int counter = -1;
  for(int i=0; i<n_i; i++){
    if(CalcSD_lambda_ip(i)==1){
      counter = counter + 1;
      SD_report(counter) = log(N_ip(i,1));
    }
  }
  
  vector<Type> N_i(n_i);
  for(int i=0; i<n_i; i++){
    N_i(i) = N_ip(i,1);
  }
  Type mean_N = N_i.sum() / N_i.size();
  Type total_N = N_i.sum() ;
  
  Type mean_p = detectprob_ip.mean();
  
  // Add up components
  jnll = jnll_comp.sum();
  
  // Spatial field summaries
  REPORT( rho_b );
  REPORT( rho_st );
  // REPORT( rho_t_b );
  REPORT( SD_b );
  REPORT( SDinput );
  REPORT( SDinput_st );
  REPORT( SD_inf );
  REPORT( SD_st_inf );
  //  REPORT( SDinput_b );
  //  REPORT( SDinput_t_b );
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
  REPORT( Deltainput_t );
  REPORT( rhot );
  REPORT( sigmat );
  REPORT( jnll_comp );
  REPORT( jnll );
  REPORT( sigmaIID );
  REPORT( lognormal_overdispersed_i );
  //REPORT( SDinput_t_b );
  //REPORT( SDinput_st );
  REPORT( Nu_dt );
  REPORT( lambda_dt );
  REPORT( N_ip );
  REPORT( chat_ip );
  REPORT( extradetectionSD );
  REPORT( mean_N );
  REPORT( rho_t_b );
  REPORT( SDinput_t_b );
  REPORT( temp_b );
  ADREPORT( lambda_t );

  // ADREPORT( lambda_ip);
  ADREPORT( mean_N );
  ADREPORT( total_N );
  ADREPORT( log_mean );
  ADREPORT( mu );
  ADREPORT( gamma_j );
  // ADREPORT( mean_p );
  ADREPORT( detectrate );
  ADREPORT( extradetectionSD );
  ADREPORT( sigmaIID );
  ADREPORT( log_theta_vec );
  ADREPORT( theta );
  ADREPORT( SDinput );
  ADREPORT( rhot );
  ADREPORT( sigmat );
  ADREPORT( rho_st );
  ADREPORT( theta_st );
  ADREPORT( SDinput_st );
  ADREPORT( sigmaIID );
  ADREPORT( lambda_t );
  //ADREPORT( SDinput );
  //  ADREPORT( SDinput_t_b );
  //  ADREPORT( SD_report );
  // ADREPORT( log_N100_dt );
  //   if( Options_vec(6)==1 ) {
  ////    vector<Type> log_N_i(n_i);
  //    vector<Type> N_i(n_i);
  //    N_i = N_ip.col(1);
  ////    log_N_i = log(N_i);
  ////     ADREPORT( log_N_i );
  //    ADREPORT( N_i );
  //   }
  
  return jnll;
}
