
// include libraries
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
#include <string>
using namespace density;
using Eigen::SparseMatrix;
// Can also choose which likelihood to use.
// Lognormal density
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}
// Inverse gamma
template<class Type>
Type dinvgauss(Type x, Type mean, Type shape, int give_log=0){
  Type logres = 0.5*log(shape) - 0.5*log(2*M_PI*pow(x,3)) - (shape * pow(x-mean,2) / (2*pow(mean,2)*x));
  if(give_log) return logres; else return exp(logres);
}
// dcauchy for hyperparameters
template<class Type>
Type dcauchy(Type x, Type mean, Type shape, int give_log=0){
  Type logres = 0.0;
  logres-= log(M_PI);
  logres-= log(shape);
  // Note, this is unstable and should switch to log1p formulation
  logres-= log(1 + pow( (x-mean)/shape ,2));
  if(give_log) return logres; else return exp(logres);
}
// helper function to make sparse SPDE precision matrix
// Inputs:
// logkappa: log(kappa) parameter value
// logtau: log(tau) parameter value
//  M0, M1, M2: these sparse matrices are output from R::INLA::inla.spde2.matern()$param.inla$M*
template<class Type>
  SparseMatrix<Type> spde_Q(Type logkappa, Type logtauO, SparseMatrix<Type> M0, SparseMatrix<Type> M1, SparseMatrix<Type> M2) {
  SparseMatrix<Type> Q;
  Type kappa2 = exp(2. * logkappa);
  Type kappa4 = kappa2*kappa2;
  Q = pow(exp(logtauO), 2.)  * (kappa4*M0 + Type(2.0)*kappa2*M1 + M2);
  return Q;
}
template<class Type>
Type objective_function<Type>::operator() ()
{
  //=============================================================================================================
  //                                              DATA SECTION
  //=============================================================================================================
  DATA_INTEGER(likelihood);      // the likelihood function to use
  DATA_INTEGER(space);	          // form of the spatial component; 1=S; 2=ST
  DATA_INTEGER(n_t);             // number of years
  DATA_INTEGER(n_s);             // number of sites (grids)
   
  // Indices for factors
  DATA_FACTOR(year);             // year as a factor
  DATA_FACTOR(site);             // Random effect index for observation i
  DATA_FACTOR(trim);             // trimester
  DATA_FACTOR(destine);          // destine
   
  // Vectors of real data
  DATA_VECTOR(cpue);             // cpue (response variable)
  DATA_VECTOR(depth);		        // depth
  
  // SPDE objects
  DATA_SPARSE_MATRIX(M0);    // used to make gmrf precision
  DATA_SPARSE_MATRIX(M1);    // used to make gmrf precision
  DATA_SPARSE_MATRIX(M2);    // used to make gmrf precision
  
  // Prior means
  DATA_SCALAR(logtaumeanO);  // mean of prior for logtauO
  DATA_SCALAR(logkappamean); // mean of prior for logkappa
  
  
  //=============================================================================================================
  //                                              PARAMETERS SECTION
  //=============================================================================================================
  
  // Only for Gamma distribution
  //PARAMETER(shape);
  //PARAMETER(scale);

  // Fixed effects
  PARAMETER(intercept);		           // global mean cpue
  PARAMETER(beta_depth);		           // beta for depth
  
  // Parameters for factors  
  PARAMETER_VECTOR(beta_year);	        // beta for year       (factor)
  PARAMETER_VECTOR(beta_trim);	        // beta for trimester  (factor)
  PARAMETER_VECTOR(beta_destine);       // beta for destine    (factor)
    
  // Variances
  PARAMETER(logsigma);		         
  PARAMETER(logtauO);		         // spatial process
  PARAMETER(logkappa);		         // decorrelation distance (kind of)
  
  PARAMETER_VECTOR(omega_s);	     // spatial effects; n_s length /parameters for sites)
  SparseMatrix<Type> Q   = spde_Q(logkappa, logtauO, M0, M1, M2);
  
  // Priors
  Type nlp=0.0;  // negative log prior  (priors)
  nlp-= dnorm(intercept,     Type(0.0), Type(5.0), true);
  nlp-= dnorm(beta_year,     Type(0.0), Type(5.0), true).sum();
  nlp-= dnorm(beta_depth,    Type(0.0), Type(5.0), true);
  nlp-= dnorm(beta_trim,     Type(0.0), Type(5.0), true).sum();
  nlp-= dnorm(beta_destine,  Type(0.0), Type(5.0), true).sum();
   
  // Variance component
  Type sigma = exp(logsigma);
   
  nlp -= dcauchy(sigma,   Type(0.0),       Type(2.0));
   
  // Hyperpriors
  Type taumeanO = exp(logtaumeanO);
  Type kappamean = exp(logkappamean);
   
  Type tauO = exp(logtauO);
  Type kappa = exp(logkappa);
   
  nlp -= dnorm(tauO,   taumeanO,   Type(2.0), true);
  nlp -= dnorm(kappa,  kappamean , Type(2.0), true);
   
  //=============================================================================================================
  // Objective function is sum of negative log likelihood components
  using namespace density;
  int n_i = cpue.size();	             // number of observations (cpue)
  Type nll_omega=0;		               // spatial effects
   
   
  // The model predicted cpue for each observation, in natural space:
  // cpue
  vector<Type> mu(n_i);
  for( int i=0; i<n_i; i++){
    mu(i) = exp(intercept + beta_year(year(i)) + beta_depth*depth(i)  + beta_trim(trim(i)) + beta_destine(destine(i)) + omega_s(site(i))); 
  }


 
  
    // Probability of random effects
  // Space 
  if(space>0) 
    nll_omega += GMRF(Q)(omega_s);
  
    // Probability of the data, given random effects (likelihood)
  vector<Type> log_lik(n_i);
  for( int i = 0; i<n_i; i++){
    if(likelihood==1) // lognormal case
      log_lik(i) = dlognorm(cpue(i), log(mu(i)), sigma, true);
    if(likelihood==2) // gamma case
      log_lik(i) = dgamma(cpue(i), 1/pow(sigma,2), mu(i)*pow(sigma,2), true);
  }
  Type nll = -log_lik.sum(); // total NLL
   
  // cpue by year
  vector<Type> cpue_t(n_t); // one for each year
  for( int t=0; t<n_t; t++){
      cpue_t(t)=0;
      for(int s=0; s<n_s; s++){
      cpue_t(t) = exp(intercept + beta_year(t) + omega_s(s));
    }
  }


  // Simule data from the mu
  vector<Type> cpue_sim(n_i);
  for( int i=0; i<n_i; i++){
  if(likelihood==1)
      SIMULATE {
      cpue_sim(i) = rnorm(log(mu(i)), sigma);
      REPORT(cpue_sim)};
  if(likelihood==2)
      SIMULATE {
      cpue_sim(i) = rgamma(1/pow(sigma,2), mu(i)*pow(sigma,2));
      REPORT(cpue_sim);
   }
  }


  

  
  
  // Jacobian adjustment for variances
  nll -= logsigma + logtauO + logkappa;
   
  // Calculate joint negative log likelihood
  Type jnll = nll + nll_omega + nlp;
   
   
  vector<Type> preds = mu;
    
    
  // Derived quantities, given parameters
  // Geospatial
  Type Range = sqrt(8) / exp(logkappa);
  Type SigmaO = 1 / sqrt(4 * M_PI * exp(2*logtauO) * exp(2*logkappa));
   
  
  //=====================================================================================================
  //  Reporting
  REPORT(nll);
  REPORT(nll_omega);
    
  REPORT(intercept);
  REPORT(beta_year);
  REPORT(beta_depth);
  REPORT(beta_trim);
  REPORT(beta_destine);
  REPORT(omega_s);
  REPORT(preds);
  REPORT(Range);
  REPORT(SigmaO);
  REPORT(logsigma);
  REPORT(logtauO);
  REPORT(logkappa);
  REPORT(cpue_t);
  REPORT(log_lik);

  ADREPORT(intercept);		            // mean cpue
  ADREPORT(beta_year);	              // year 
  ADREPORT(beta_depth);	              // depth
  ADREPORT(beta_trim);	              // trimester
  ADREPORT(beta_destine);	            // destine
  ADREPORT(logsigma);	                // observations  
  ADREPORT(logtauO);
  ADREPORT(cpue_t);
  ADREPORT(logkappa);
    
  // Derived geospatial components
  ADREPORT(Range);		         // geostatistical range
  ADREPORT(SigmaO);		         // spatiotemporal
  return jnll;
}

