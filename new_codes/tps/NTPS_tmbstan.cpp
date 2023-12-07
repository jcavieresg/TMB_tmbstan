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


template<class Type>
Type objective_function<Type>::operator() ()
{
  
  
  //============================================
  //                 DATA
  //============================================
  //DATA_INTEGER(space);                         // Type of model
  DATA_INTEGER(likelihood);                    // the likelihood function to use
  DATA_VECTOR(y);                              // The response
  DATA_VECTOR(x1);                              
  
  DATA_MATRIX(TPS);                              // Design matrix for splines
  
  DATA_SPARSE_MATRIX(S);                       // Penalization matrix diag(S1,S2,S3,S4,S5) without storing off-diagonal zeros.
  DATA_IVECTOR(Sdims);                         // Dimensions of S1,S2,S3,S4 and S5
  DATA_SPARSE_MATRIX(tpsReport);               // Design matrix for report of splines
  
  
  //----------------------------------
  
  //============================================
  //                 PARAMETERS
  //============================================
  
  // Parameters for factors  
  PARAMETER(beta0);
  PARAMETER(beta1);
  
  PARAMETER(logsigma);                        // Penalizatiion parameter
  PARAMETER(loglambda);                       // skew normal parameter (slant)
  
  PARAMETER_VECTOR(smoothCoefs);              // Spline regression parameters
  PARAMETER(logalpha);                  // Penalization parameters
  
  
  
  //==========================================
  // Transformed parameters
  Type sigma  = exp(logsigma);
  Type lambda = exp(loglambda);
  //vector<Type> alpha = exp(logalpha);
  Type alpha = exp(logalpha);

  
  //==========================================  
  // Priors
  Type nlp = Type(0.0);                                 // negative log prior  (priors)
  
  nlp-= dnorm(beta0,    Type(0.0), Type(2.0), true);
  nlp-= dnorm(beta1,    Type(0.0), Type(2.0), true);

  // Variance component
  nlp-= dcauchy(sigma,   Type(0.0),   Type(2.0));

  // Penalty parameter
  nlp-= dnorm(lambda, Type(0.0),   Type(1.0), true);
  nlp-= dnorm(alpha, Type(0.0),    Type(1.0), true);
  
  nlp-= dnorm(smoothCoefs, Type(0.0), Type(1.0), true).sum(); 
  
  //==========================================
  // OBJETIVE FUNCTION
  Type nll = Type(0.0); 
  
  int k=0;  // Counter
  //if(space > 0)
  for(int i=0;i<Sdims.size();i++){
    int m_i = Sdims(i);
    vector<Type> smoothCoefs_i = smoothCoefs.segment(k,m_i);       // Recover betai
    SparseMatrix<Type> S_i = S.block(k,k,m_i,m_i);  // Recover Si
    //nll -= Type(0.5)*m_i*logalpha(i) - Type(0.5)*alpha(i)*GMRF(S_i).Quadform(smoothCoefs_i);
    nll -= Type(0.5)*m_i*logalpha - Type(0.5)*alpha*GMRF(S_i).Quadform(smoothCoefs_i);
    k += m_i;
  }
  
  
  
  
  
  // We create a vector of means
  vector<Type> mu(y.size());
  mu = beta0 + beta1*x1 + TPS*smoothCoefs;
  //mu = exp(beta0 + beta_year(year) + beta_depth*depth + TPS*smoothCoefs);
  
  //===============================================================
  // log-likelihood for the response
  int n = y.size();
  vector<Type> log_lik(n);
  vector<Type> z = y - mu;  
  
  for( int i = 0; i< n; i++){
    if(likelihood==1) // normal case
      log_lik(i) = dnorm(y(i), mu(i), sigma, true);
    if(likelihood==2) // skew-normal case
      log_lik(i) = dsn(z(i)/sigma, lambda, true) - log(sigma);
  }
  nll = -log_lik.sum(); // total NLL  
  
  //---------------------------------------------------------------
  
  
  // Simule data from the mu 
  vector<Type> y_sim(n);
  for( int i=0; i<n; i++){
    SIMULATE {
      y_sim(i) = rnorm(mu(i), sigma);
      REPORT(y_sim);
    }
  }
  //----------------------------------
  
  
  
  //==========================================
  // Derived quantities
  vector<Type> preds = mu;
  vector<Type> splines2D = tpsReport*smoothCoefs;
  
  // Jacobian adjustment for transformed parameters
  nll -= logsigma + logalpha + loglambda;   // add logalpha? how?
  
  // Calculate joint negative log likelihood
  Type jnll = nll + nlp;
  
  
  //=====================================================================================================
  // REPORT
  REPORT(nll);
  REPORT(jnll);
  
  REPORT(beta0);
  REPORT(beta1);
  REPORT(logalpha);
  REPORT(smoothCoefs);
  REPORT(mu);
  REPORT(log_lik);
  REPORT(splines2D);
  
  //=====================================================================================================
  // AD REPORT
  ADREPORT(beta0);
  ADREPORT(beta1);
  ADREPORT(logsigma);
  ADREPORT(loglambda);
  ADREPORT(logalpha);
  ADREPORT(smoothCoefs);
  ADREPORT(splines2D);
  return jnll;
}