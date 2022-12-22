// Estimating parameters in a Tweedie distribution.
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
#include <string>
using namespace density;
using Eigen::SparseMatrix;
using namespace R_inla; //includes SPDE-spesific functions, e.g. Q_spde()


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
  

//==========================================
// Data section
//  DATA_INTEGER(likelihood);                    // the likelihood function to use
  DATA_VECTOR(y);
  DATA_MATRIX(X);                              // Design matrix for splines
  
  DATA_SPARSE_MATRIX(S);                       // Penalization matrix diag(S1,S2,S3,S4,S5) without storing off-diagonal zeros.
  DATA_IVECTOR(Sdims);                         // Dimensions of S1,S2,S3,S4 and S5
  DATA_SPARSE_MATRIX(tpsReport);               // Design matrix for report of splines
  

//==========================================
// Parameter section
  PARAMETER(beta0);
  PARAMETER(phi);
  PARAMETER(p);
  

  PARAMETER_VECTOR(smoothCoefs);               // Spline regression parameters
  PARAMETER_VECTOR(logalpha);                  // Penalization parameters
  
//==========================================
// Transformed parameters
  vector<Type> alpha = exp(logalpha);  
  
  
//==========================================  
// Priors
  Type nlp = Type(0.0);                                 // negative log prior  (priors)
  
  //nlp-= dnorm(beta0,    Type(1.0), alpha, true);
  nlp-= dcauchy(beta0,    Type(0.0), Type(5.0), true);
  
  // Parameter Tweddie distribution
  nlp-= dcauchy(phi,   Type(0.0),   Type(5.0));
  
  
  // Penalty parameter
  nlp-= dnorm(alpha, Type(1.0),   Type(2.0), true).sum();
  //nlp-= dnorm(alpha, Type(0.0),   Type(2.0), true).sum();
  
  
  // Prior for smoothCoefs
  nlp-= dnorm(smoothCoefs, Type(0.0), Type(2.0), true).sum(); 
  
  
//==========================================
// OBJETIVE FUNCTION
// Initialize the likelihood
Type nll = 0; 
  

int k=0;  // Counter
  for(int i=0;i<Sdims.size();i++){
    int m_i = Sdims(i);
    vector<Type> smoothCoefs_i = smoothCoefs.segment(k,m_i);       // Recover betai
    SparseMatrix<Type> S_i = S.block(k,k,m_i,m_i);  // Recover Si
    nll -= Type(0.5)*m_i*logalpha(i) - Type(0.5)*alpha(i)*GMRF(S_i).Quadform(smoothCoefs_i);

    k += m_i;
}
  
  

// We create a vector of means
  vector<Type> mu(y.size());
  mu = beta0 + X*smoothCoefs;
  
  

//===============================================================
// log-likelihood for the response
  
int n_i = y.size();
  for( int i = 0; i<n_i; i++){
  nll -= dtweedie(y(i), mu(i), phi, p, true);
  
}
//---------------------------------------------------------------
  
  
  
  
//Simule data from the mu (only an experimental analysis because TMB not support rdsn)
vector<Type> y_sim(n_i);
for( int i=0; i<n_i; i++)
  SIMULATE {
    y_sim(i) = rtweedie(mu(i), phi, p);
    REPORT(y_sim)
  }
//----------------------------------
  
  
  
  //==========================================
  // Derived quantities
  vector<Type> preds = mu;
  vector<Type> splines2D = tpsReport*smoothCoefs;
  
  // Jacobian adjustment for transformed parameters
  //nll -= logphi;   // add logalpha? how?
  
  // Calculate joint negative log likelihood
  Type jnll = nll + nlp;
  
  
  //=====================================================================================================
  // REPORT
  REPORT(jnll);
  
  REPORT(beta0);
  REPORT(phi);
  REPORT(p);
  REPORT(logalpha);
  REPORT(smoothCoefs);
  REPORT(preds);
  REPORT(splines2D);
  
  //=====================================================================================================
  // AD REPORT
  ADREPORT(beta0);		            
  ADREPORT(phi);
  ADREPORT(p);	
  ADREPORT(logalpha);
  ADREPORT(smoothCoefs);
  ADREPORT(splines2D);
  return jnll;
}
