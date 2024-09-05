// Simple linear regression with a smoothing spline

#include <TMB.hpp>

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


template<class Type>
Type objective_function<Type>::operator() ()
{

  using namespace R_inla;
  using namespace density;
  using namespace Eigen;

  
  //============================================
  //                 DATA
  //============================================
  DATA_SCALAR(model);
  
  DATA_VECTOR(y);
  DATA_VECTOR(x1);
  
  //============================================
  //                 PARAMETERS
  //============================================
  PARAMETER(beta0);            // intercept
  PARAMETER(beta1);            

  PARAMETER(logsigma_e);     // measurement noise 

  // spline things
  PARAMETER_VECTOR(x);     // spline params
  DATA_MATRIX(X);          // spline design matrix
  DATA_SPARSE_MATRIX(S);   // smoothing penalty matrix
  PARAMETER(loglambda);       // smoothing parameter


  //==========================================
  // Transformed parameters
  Type sigma_e  = exp(logsigma_e);
  Type lambda = exp(loglambda);
  
  SparseMatrix<Type> Q = lambda*S;   // precision for spline
  
  
  //==========================================  
  // Priors
  Type nlp = Type(0.0);                                 
  
  nlp-= dnorm(beta0,    Type(1.0), Type(2.0), true);
  nlp-= dnorm(beta1,    Type(1.0), Type(2.0), true);

  // Variance component
  nlp-= dcauchy(sigma_e,   Type(0.0),   Type(5.0));

  // Penalty parameter
  nlp-= dnorm(lambda, Type(0.0),   Type(1.0), true);
  //nlp-= dcauchy(lambda,    Type(5.0), Type(5.0));
  
  
  nlp-= dnorm(x, Type(0.0), Type(1.0), true).sum(); 
  
  
  // We create a vector of means
  vector<Type> mu(y.size());
  if(model == 1)
  mu = beta0 + beta1*x1 + X*x;
  if(model == 2)
  mu  = X*x;;
  
  
  // Probability of the data, given random effects (likelihood)
  vector<Type> log_lik(y.size());
  for( int i = 0; i<y.size(); i++){
      log_lik(i) = dnorm(y(i), mu(i), sigma_e, true);
  }
  Type nll = -log_lik.sum(); // total NLL

  
  nll -= Type(0.5)*1.0*loglambda - 0.5*lambda*GMRF(S).Quadform(x);

  
  // Simule data from the mu 
  vector<Type> y_sim(y.size());
  vector<Type> tps(y.size());
  tps = X*x;

if(model == 1)
   for( int i=0; i<y.size(); i++){
   SIMULATE {y_sim(i) = rnorm(beta0 + beta1*x1(i) + tps(i), sigma_e);
   REPORT(y_sim)};
}
if(model == 2)
   for( int i=0; i<y.size(); i++){
   SIMULATE {y_sim(i) = rnorm(tps(i), sigma_e);
   REPORT(y_sim)};
}

  
  
  //==========================================
  // Derived quantities
  vector<Type> pred = mu;

  // Jacobian adjustment for transformed parameters
  nll -= logsigma_e + loglambda;   // 
  
  // Calculate joint negative log likelihood
  Type jnll = nll + nlp;
  //----------------------------------
  return jnll;
  
  REPORT(log_lik);
}
