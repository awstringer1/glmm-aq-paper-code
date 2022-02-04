#include <TMB.hpp>                                // Links in the TMB libraries
//#include <fenv.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Taken from glmmTMB: check if openmp configured, set parallelization accordingly.
#ifdef _OPENMP
  this -> max_parallel_regions = omp_get_max_threads();
  // std::cout << "OpenMP max_parallel_regions=" << this -> max_parallel_regions << "\n";
#else
  this -> max_parallel_regions = 1;
  // std::cout << "no OpenMP (max_parallel_regions=1)\n";
#endif

  // Weibull survival model with Gamma frailty.
  // This template implements the joint log likelihood for a SINGLE patient
  // (one index i).
  // Data
  DATA_VECTOR(t); // Vector of survival times for patient i, dimension m_i
  int mi = t.size();
  DATA_VECTOR(delta); // Vector of censoring indicators, delta_j = 1 if t_j is a failure time, 0 if censoring time
  DATA_MATRIX(x); // Covariate matrix for subject i. Dimension mi x p where there are p covariates

  // Parameters
  // Note the only TMB "PARAMETER" is z, because we are implementing the joint log likelihood
  // for the purpose of integrating out z for fixed values of the other params. They are therefore "DATA"
  // in TMB
  DATA_SCALAR(logmu);
  Type mu = exp(logmu);
  DATA_SCALAR(logalpha);
  Type alpha = exp(logalpha);
  DATA_VECTOR(beta);
  vector<Type> xTbeta = x*beta;
  DATA_SCALAR(logeta);
  Type eta = exp(logeta);
  PARAMETER(z); // Log-frailty term for subject i; this is the variable of integration

  // Log likelihood
  // Term 1: log density of z
  Type loglik = eta*logeta - lgamma(eta) + eta*(z - exp(z));
  // Term 2: data terms
  for (int j=0;j<mi;j++) {
    loglik += delta(j) * (logmu + logalpha + (alpha-1.0)*log(t(j)) + xTbeta(j) + z);
    loglik -= mu*pow(t(j),alpha) * exp(xTbeta(j) + z);
  }
  return loglik;
}
