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
  DATA_INTEGER(likelihood);                    // the likelihood function to use
  DATA_VECTOR(y);                              // The response
  DATA_MATRIX(X);                             // Continuous covariate
  
  // 
  DATA_MATRIX(TPS);                              // Design matrix for splines
  DATA_SPARSE_MATRIX(S);                       // Penalization matrix diag(S1,S2,S3,S4,S5) without storing off-diagonal zeros.
  DATA_IVECTOR(Sdims);                         // Dimensions of S1,S2,S3,S4 and S5
  DATA_SPARSE_MATRIX(tpsReport);               // Design matrix for report of splines
  
  

  //============================================
  //                 PARAMETERS
  //============================================
  PARAMETER_VECTOR(beta);	        

  PARAMETER_VECTOR(smoothCoefs);               // Spline regression parameters
  //PARAMETER_VECTOR(logalpha);                  // Penalization parameters
  PARAMETER(logalpha);                  // Penalization parameters
  
  
  PARAMETER(logsigma);                        
  PARAMETER(loglambda);                       
  
  

  //==========================================
  // Transformed parameters
  Type sigma  = exp(logsigma);
  Type lambda = exp(loglambda);
  //vector<Type> alpha = exp(logalpha);
  Type alpha = exp(logalpha);
  //----------------------------------
  
  
  
  //==========================================  
  // Priors
  Type nlp = Type(0.0);                                 // negative log prior  (priors)
  
  nlp-= dnorm(beta,    Type(0.0), Type(10.0), true).sum();

 // Variance component
  nlp-= dcauchy(sigma,   Type(0.0),   Type(2.0));
  
 // Penalty parameter
  nlp-= dnorm(alpha, Type(1.0),   Type(2.0), true);
  
// Slant parameter of Skew Normal
  nlp-= dcauchy(lambda,   Type(0.0),  Type(2.0));
  
  
//==========================================
  // OBJETIVE FUNCTION
  Type nll=0;
  
  int k=0;  // Counter
  for(int i=0;i<Sdims.size();i++){
    int m_i = Sdims(i);
    vector<Type> smoothCoefs_i = smoothCoefs.segment(k,m_i);       // Recover betai
    SparseMatrix<Type> S_i = S.block(k,k,m_i,m_i);  // Recover Si
    nll -= Type(0.5)*m_i*logalpha - Type(0.5)*alpha*GMRF(S_i).Quadform(smoothCoefs_i);
    k += m_i;
  }
  
  
  
// We create a vector of means
  vector<Type> mu(y.size());
  mu = X*beta + TPS*smoothCoefs;


  //===============================================================
  // log-likelihood for the response
  int n_i = y.size();
  vector<Type> z = y - mu;
  for( int i = 0; i<n_i; i++){

    // Skew-Normal case
    if(likelihood==1)
      nll -= dsn(z(i)/sigma, lambda, true) - log(sigma);
    // Normal case
    if(likelihood==2)
      nll -= dnorm(y(i), mu(i), sigma, true);
  }
  //---------------------------------------------------------------
  
  

  // Simule data from the mu (only an experimental analysis because TMB not support rdsn)
  vector<Type> y_sim(n_i);
  for( int i=0; i<n_i; i++)
    SIMULATE {
      y_sim(i) = rnorm(mu(i), sigma);
      REPORT(y_sim)
    }
  //----------------------------------
  
  
  
  //==========================================
  // Derived quantities
  vector<Type> preds = mu;
  vector<Type> splines2D = tpsReport*smoothCoefs;
  
  // Jacobian adjustment for transformed parameters
  nll -= logsigma + loglambda + logalpha; 
  
  // Calculate joint negative log likelihood
  Type jnll = nll + nlp;
  
  
  //=====================================================================================================
  // REPORT
  REPORT(nll);
  
  REPORT(beta);
  REPORT(logsigma);
  REPORT(loglambda);
  REPORT(logalpha);
  REPORT(smoothCoefs);
  REPORT(preds);
  REPORT(splines2D);
  
  //=====================================================================================================
  // AD REPORT
  ADREPORT(beta);		            
  ADREPORT(logsigma);
  ADREPORT(loglambda);	
  ADREPORT(logalpha); 
  ADREPORT(smoothCoefs);
  ADREPORT(splines2D);
  
  
  return jnll;
}