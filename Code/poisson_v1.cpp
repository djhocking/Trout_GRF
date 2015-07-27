#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(Data_ij);
  
  PARAMETER_VECTOR(Alpha_j);
  PARAMETER_VECTOR(Beta_i);
  
  int n_i = Data_ij.col(0).size();
  int n_j = Data_ij.row(0).size();  
  matrix<Type> Lambda_ij(n_i,n_j);
  Type jnll = 0;
  
  for(int i=0; i<n_i; i++){
  for(int j=0; j<n_j; j++){
    Lambda_ij(i,j) = exp( Alpha_j(j) + Beta_i(i) );
    jnll -= dpois( Data_ij(i,j), Lambda_ij(i,j), true );
  }}
  
  vector<Type> Prob_j(n_j);
  Prob_j = exp(Alpha_j) / exp(Alpha_j).sum();
  
  REPORT( Lambda_ij );
  REPORT( Prob_j );
  ADREPORT( Prob_j );
  return jnll;
}

