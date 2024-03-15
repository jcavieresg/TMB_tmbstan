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
  //DATA_INTEGER(model);                    // the likelihood function to use
  DATA_VECTOR(cpue);                              // The response
  
  DATA_FACTOR(year);             // year as a factor
  DATA_FACTOR(trim);             // trimester
  DATA_FACTOR(destine);          // destine
  
  DATA_VECTOR(depth);
  
  
   
  DATA_MATRIX(TPS);                              // Design matrix for splines
  
  DATA_SPARSE_MATRIX(S);                       // Penalization matrix diag(S1,S2,S3,S4,S5) without storing off-diagonal zeros.
  DATA_IVECTOR(Sdims);                         // Dimensions of S1,S2,S3,S4 and S5
  DATA_SPARSE_MATRIX(tpsReport);               // Design matrix for report of splines
  
  
  //----------------------------------
  
  //============================================
  //                 PARAMETERS
  //============================================
  
  // Parameters for factors  
  //PARAMETER_VECTOR(beta);	        // beta for trimester  (factor)
  PARAMETER(beta0);
  PARAMETER_VECTOR(beta_year);
  PARAMETER_VECTOR(beta_trim);
  PARAMETER_VECTOR(beta_destine);
  
  PARAMETER(beta_depth);
  
  PARAMETER_VECTOR(smoothCoefs);               // Spline regression parameters
  //PARAMETER_VECTOR(logalpha);                  // Penalization parameters
  PARAMETER(loglambda);                  // Penalization parameters
  PARAMETER(logsigma);                        // Penalizatiion parameter
  
  //PARAMETER(logomega);                        // skew normal parameter (slant)
  
  

  //==========================================
  // Transformed parameters
  Type sigma  = exp(logsigma);
  Type lambda = exp(loglambda);
  //vector<Type> alpha = exp(logalpha);
  //Type lambda = exp(loglambda);
  //Type omega = exp(logomega);
  
  //==========================================  
  // Priors
  Type nlp = Type(0.0);                                 // negative log prior  (priors)
  
  nlp-= dnorm(beta0,        Type(2.0), Type(2.0), true);
  nlp-= dnorm(beta_year,    Type(0.0), Type(1.0), true).sum();
  nlp-= dnorm(beta_trim,    Type(0.0), Type(1.0), true).sum();
  nlp-= dnorm(beta_destine, Type(0.0), Type(2.0), true).sum();
  
  nlp-= dnorm(beta_depth,   Type(1.0), Type(2.0), true);
  
  // Variance component
  nlp-= dcauchy(sigma,   Type(0.0),   Type(2.0));
  
  // Penalty parameter
  //nlp-= dgamma(alpha,   Type(1.0),   Type(2.0), true).sum(); ---> WORKS!!
  nlp-= dnorm(lambda, Type(1.0),   Type(1.0), true);
  
  nlp-= dnorm(smoothCoefs, Type(0.0), Type(1.0), true).sum(); 
  

  // We create a vector of means
  int n = cpue.size();
  vector<Type> mu(cpue.size());
  mu = beta0 + beta_year(year) + beta_trim(trim) + beta_destine(destine) + beta_depth*depth + TPS*smoothCoefs;  
  
  
  
  

//==========================================
//      OBJETIVE FUNCTION
Type nll = 0.0;
//==========================================
    
//---------------------------------------------------------------
vector<Type> log_lik(n);
//vector<Type> z = cpue - mu;  
  
for( int i = 0; i< n; i++){
    if(likelihood==1) // lognormal case
      log_lik(i) = dlognorm(cpue(i), mu(i), sigma, true);
    if(likelihood==2) // Gamma case
      log_lik(i) = dgamma(cpue(i), 1/pow(sigma,2), exp(mu(i))*pow(sigma,2), true);
    //if(likelihood==3) // Skew normal case
    //  log_lik(i) = dsn(z(i)/sigma, omega, true) - log(sigma);
}
nll = -log_lik.sum(); // total NLL  
//---------------------------------------------------------------

  
int k=0;  // Counter
for(int i=0;i<Sdims.size();i++){
    int m_i = Sdims(i);
    vector<Type> smoothCoefs_i = smoothCoefs.segment(k,m_i);       // Recover betai
    SparseMatrix<Type> S_i = S.block(k,k,m_i,m_i);  // Recover Si
    //nll -= Type(0.5)*m_i*logalpha(i) - Type(0.5)*alpha(i)*GMRF(S_i).Quadform(smoothCoefs_i);
    nll -= Type(0.5)*m_i*loglambda - Type(0.5)*lambda*GMRF(S_i).Quadform(smoothCoefs_i);
    k += m_i;
}
  
    
  
// Simule data from the mu 
  vector<Type> cpue_sim(n);
  for( int i=0; i<n; i++){
    if(likelihood==1)
      SIMULATE {
        cpue_sim(i) = rnorm(log(mu(i)), sigma);
        REPORT(cpue_sim)};
    if(likelihood==2)
      SIMULATE {
        cpue_sim(i) = rgamma(1/pow(sigma,2), exp(mu(i))*pow(sigma,2));
        REPORT(cpue_sim)};
    // if(likelihood==3)
    //   SIMULATE {
    //     cpue_sim(i) = rnorm(mu(i), sigma);
    //     REPORT(cpue_sim);
    //   }
}
//----------------------------------
  
  
  
//==========================================
// Derived quantities
vector<Type> preds = mu;
vector<Type> splines2D = tpsReport*smoothCoefs;

//==========================================
// Jacobian adjustment for transformed parameters
//nll -= logsigma + logomega + loglambda;   
nll -= logsigma + loglambda;   


//==========================================
// Calculate joint negative log likelihood
Type jnll = nll + nlp;
  
  
//=====================================================================================================
// REPORT
REPORT(jnll);
  
REPORT(beta0);
REPORT(beta_year);
REPORT(beta_depth);
REPORT(logsigma);
REPORT(loglambda);
//REPORT(logomega);
REPORT(smoothCoefs);
REPORT(preds);
REPORT(log_lik);
REPORT(splines2D);
  
//=====================================================================================================
// AD REPORT
ADREPORT(beta0);		            
ADREPORT(logsigma);
ADREPORT(loglambda);	
//ADREPORT(logomega); 
ADREPORT(smoothCoefs);
ADREPORT(splines2D);
return jnll;
}