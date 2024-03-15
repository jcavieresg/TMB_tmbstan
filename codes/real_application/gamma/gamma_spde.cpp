
// include libraries
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
#include <string>
using namespace density;
using Eigen::SparseMatrix;
using namespace R_inla; //includes SPDE-spesific functions, e.g. Q_spde()


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
  logres-= log(1 + pow( (x-mean)/shape ,2));
  if(give_log) return logres; else return exp(logres);
}





template<class Type>
Type objective_function<Type>::operator() ()
{
  //=============================================================================================================
  //                                              DATA SECTION
  //=============================================================================================================
  DATA_INTEGER(likelihood);      // the likelihood function to use

  // Numerics
  DATA_VECTOR(cpue);             // cpue (response variable)
  DATA_VECTOR(depth);		        // depth
  
  // Factors
  DATA_FACTOR(year);             // year as a factor
  DATA_FACTOR(trim);             // trimester
  DATA_FACTOR(destine);          // destine
   
  // SPDE objects
  DATA_STRUCT(spde_mat, spde_t);  // see p. 8 in Lindgren et al. (2011)
  DATA_SPARSE_MATRIX(A);          // used to project from spatial mesh to data locations
  
  // Prior means
  DATA_SCALAR(tau0);  // mean of prior for logtauO
  DATA_SCALAR(kappa0); // mean of prior for logkappa
  
  
  //=============================================================================================================
  //                                              PARAMETERS SECTION
  //=============================================================================================================
  // Fixed effects
  PARAMETER(beta0);		           // global mean cpue
  PARAMETER(beta_depth);		           // beta for depth
  
  // Parameters for factors  
  PARAMETER_VECTOR(beta_year);	        // beta for year       (factor)
  PARAMETER_VECTOR(beta_trim);	        // beta for trimester  (factor)
  PARAMETER_VECTOR(beta_destine);       // beta for destine    (factor)
    
  // Variances
  PARAMETER(logsigma);		         
  PARAMETER(logtau);		         // spatial process
  PARAMETER(logkappa);		       // decorrelation distance
  
  PARAMETER_VECTOR(u);	     // spatial effects
  
  // Priors
  Type nlp=0.0;  // negative log prior  (priors)
  nlp-= dnorm(beta0,         Type(0.0), Type(5.0), true);
  nlp-= dnorm(beta_year,     Type(0.0), Type(1.0), true).sum();
  nlp-= dnorm(beta_trim,     Type(0.0), Type(1.0), true).sum();
  nlp-= dnorm(beta_destine,  Type(0.0), Type(5.0), true).sum();
   
   nlp-= dnorm(beta_depth,    Type(0.0), Type(2.0), true);
  
  // Variance component
  Type sigma = exp(logsigma);
  nlp -= dcauchy(sigma,   Type(0.0),       Type(2.0));
   
  // Hyperpriors
  Type tau = exp(logtau);
  Type kappa = exp(logkappa);
   
  nlp -= dnorm(tau,   tau0,    Type(2.0), true);
  nlp -= dnorm(kappa, kappa0 , Type(2.0), true);
  
  SparseMatrix<Type> Q = Q_spde(spde_mat, kappa);
   
  //=============================================================================================================
  int n = cpue.size();	         // number of observations (cpue)
  Type nll_u=0;		               // spatial effects
  
  
  // The model predicted cpue for each observation, in natural space:
  // cpue
  vector<Type> mu(n);
  mu  = beta0 + beta_year(year) + beta_trim(trim) + beta_destine(destine) + beta_depth*depth + A*u;

  nll_u += GMRF(Q)(u);
  
  // Likelihood
  vector<Type> log_lik(n);
  for( int i = 0; i<n; i++){
    if(likelihood==1) // lognormal case
      log_lik(i) = dlognorm(cpue(i), log(mu(i)), sigma, true);
    if(likelihood==2) // gamma case
      log_lik(i) = dgamma(cpue(i), 1/pow(sigma,2), exp(mu(i))*pow(sigma,2), true);
  }
  Type nll = -log_lik.sum(); // total NLL
   
  // Simule data from the mu
  vector<Type> cpue_sim(n);
  for( int i=0; i<n; i++){
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
  nll -= logsigma + logtau + logkappa;
   
  // Calculate joint negative log likelihood
  Type jnll = nll + nll_u + nlp;
   
   
  vector<Type> pred = mu;
    
    
  // Derived quantities
  Type rho = sqrt(8) /kappa;
  Type sigma_u = 1 / sqrt(4 * M_PI * 2*tau * 2*kappa);
   
  
  //=====================================================================================================
  //  Reporting
  REPORT(beta0);
  REPORT(beta_year);
  REPORT(beta_depth);
  REPORT(beta_trim);
  REPORT(beta_destine);
  REPORT(u);
  REPORT(pred);
  REPORT(rho);
  REPORT(sigma_u);
  REPORT(logsigma);
  REPORT(logtau);
  REPORT(logkappa);
  REPORT(log_lik);

  ADREPORT(beta0);		            // mean cpue
  ADREPORT(beta_year);	              // year 
  ADREPORT(beta_depth);	              // depth
  ADREPORT(beta_trim);	              // trimester
  ADREPORT(beta_destine);	            // destine
  ADREPORT(logsigma);	                // observations  
  ADREPORT(logtau);
  ADREPORT(logkappa);
    
  ADREPORT(u);
  return jnll;
}

