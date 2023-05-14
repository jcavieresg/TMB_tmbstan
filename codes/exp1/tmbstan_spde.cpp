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
  // Note, this is unstable and should switch to log1p formulation
  logres-= log(1 + pow( (x-mean)/shape ,2));
  if(give_log) return logres; else return exp(logres);
}


//=====================================================
// Initialize the TMB model
//=====================================================
template<class Type>
Type objective_function<Type>::operator() ()
{


//=========================
//      DATA SECTION
//=========================
  DATA_INTEGER(model); 	          // Likelihood flag
  DATA_VECTOR(y);		              // Response
  DATA_MATRIX(X);                 // Fixed effects
  DATA_FACTOR(site);		          // Random effect index for observation i
  DATA_STRUCT(spde_mat, spde_t);  // Three matrices needed for representing the GMRF, see p. 8 in Lindgren et al. (2011)
  DATA_SPARSE_MATRIX(A);          // used to project from spatial mesh to data locations
  
  // DATA_SCALAR(tau0);
  // DATA_SCALAR(kappa0);

//=========================
//   PARAMETER SECTION
//=========================
  
  // Fixed effects
  PARAMETER_VECTOR(beta);

  PARAMETER(logsigma);
  PARAMETER(logtau);
  PARAMETER(logkappa);
  
  PARAMETER_VECTOR(u);	          // spatial effects
  
  
  // Transformed parameters
  Type sigma  = exp(logsigma);
  Type tau = exp(logtau);
  Type kappa = exp(logkappa);
  

  SparseMatrix<Type> Q = Q_spde(spde_mat, kappa);

// ===================================
//               Priors
// ===================================
   Type nlp = 0.0;                  

   nlp -= dnorm(beta,     Type(0.0),     Type(10.0), true).sum();   //

// Prior on sigma
   nlp -= dnorm(sigma,   Type(0.0), Type(2.0), true);


// Prior in Hyperparameters
   nlp -= dnorm(tau,      Type(0.0),       Type(1.0), true);   //
   nlp -= dnorm(kappa,    Type(5.0),        Type(0.5), true);   //
   
   // nlp -= dnorm(tau,      tau0,       Type(5.0), true);   //
   // nlp -= dnorm(kappa,    kappa0,     Type(5.0), true);   //
   

  //=============================================================================================================
  // Objective function is sum of negative log likelihood components
  int n = y.size();	                 // number of observations 
  vector<Type> projS(n);             // value of gmrf at data points
  Type nll_u = 0.0;		               // spatial effects
  
  
  Type nll=0;	 
  
  // Linear predictor
  vector<Type> pred(n);
  projS = A * u;                     // Project S at mesh pts to data pts
  //for(int i=0; i<n; i++){
    if(model==1)
      //pred(i) = X*beta + u(site(i));
      pred  = X*beta + u(site);
    if(model==2) 
      //pred(i) = X*beta + projS(i);
      pred = X*beta + projS;
    
  //}
  
  
  //nll_u += SCALE(GMRF(Q), 1/tau)(u); // returns negative already
  nll_u += GMRF(Q)(u); // returns negative already
  
  // Probability of data conditional on fixed effect values
  for(int i=0; i<n; i++){
    // And the contribution of the likelihood
    nll -= dnorm(y(i), pred(i), sigma, true);
  }
  
  
  // Simule data from mu
  vector<Type> y_sim(n);
  for( int i=0; i<n; i++){
    SIMULATE {
      y_sim(i) = rnorm(pred(i), sigma);
      REPORT(y_sim)};
  }
  
  
  Type rho = sqrt(8.0) / kappa;
  Type sigma_s = 1.0 / sqrt (4.0 * M_PI* 2.0*tau * 2.0 * kappa );
  
  // Jacobian adjustment 
  nll -= logsigma + logtau + logkappa;
  
  // Calculate joint negative log likelihood
  Type jnll = nll + nll_u + nlp;
  

  //=====================================================================================================
  // Reporting
  REPORT(jnll);
  REPORT(beta);
  REPORT(u);
  REPORT(pred);
  REPORT(sigma);
  REPORT(tau);
  REPORT(kappa);
  REPORT(rho);		         
  REPORT(sigma_s);		         
  
  //=======================================================================================================
  // AD report (standard devations)
  ADREPORT(beta)
  ADREPORT(sigma);
  ADREPORT(tau);
  ADREPORT(kappa);
  
  // Derived geospatial components
  ADREPORT(rho);		               // geostatistical range
  ADREPORT(sigma_s);		         
  return jnll;
}

