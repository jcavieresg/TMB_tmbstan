// Estimating parameters in a Tweedie distribution.
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
#include <string>

using namespace density;
using Eigen::SparseMatrix;
using namespace R_inla; //includes SPDE-spesific functions, e.g. Q_spde()


// LOGNORMAL
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
  DATA_INTEGER(likelihood); 	                 // Likelihood flag
  DATA_VECTOR(y);
  DATA_MATRIX(X);                              // Design matrix for splines
  
  DATA_SPARSE_MATRIX(S);                       // Penalization matrix diag(S1,S2,S3,S4,S5) without storing off-diagonal zeros.
  DATA_IVECTOR(Sdims);                         // Dimensions of S1,S2,S3,S4 and S5
  DATA_SPARSE_MATRIX(tpsReport);               // Design matrix for report of splines
  
  
  //==========================================
  // Parameter section
  PARAMETER(intercept); 		    // lintercept
  PARAMETER(theta); 		        // logit of zero_prob
  PARAMETER(logsigma);          // measurement noise 
  
  
  
  PARAMETER_VECTOR(smoothCoefs);               // Spline regression parameters
  //PARAMETER_VECTOR(logalpha);                  // Penalization parameters
  PARAMETER(logalpha);                  // Penalization parameters
  
  
  //======================================
  // Transformed quantities and parameters
  //======================================
  Type zero_prob = 1 / (1 + exp(-theta));
  Type sigma  = exp(logsigma);
  
  //vector<Type> alpha = exp(logalpha);
  Type alpha = exp(logalpha);
  //----------------------------------
  
  
//==========================================
// Priors
  Type nlp = Type(0.0);                                 // negative log prior  (priors)

  nlp-= dnorm(intercept,    Type(1.0), Type(2.0), true);

  nlp-= dcauchy(theta,   Type(0.0),   Type(1.0));
  nlp-= dcauchy(sigma,   Type(0.0),   Type(1.0));


  // Penalty parameter
  //nlp-= dnorm(alpha, Type(1.0),   Type(2.0), true).sum();
  nlp-= dnorm(alpha, Type(1.0),   Type(2.0), true);
  

  // Prior for smoothCoefs
  nlp-= dnorm(smoothCoefs, Type(0.0), Type(2.0), true).sum();

  
  //==========================================
  // OBJETIVE FUNCTION
  // Initialize the likelihood
  Type nll = Type(0.0); 
  
  
  int k=0;  // Counter
  for(int i=0;i<Sdims.size();i++){
    int m_i = Sdims(i);
    vector<Type> smoothCoefs_i = smoothCoefs.segment(k,m_i);       // Recover betai
    SparseMatrix<Type> S_i = S.block(k,k,m_i,m_i);  // Recover Si
    //nll -= Type(0.5)*m_i*logalpha(i) - Type(0.5)*alpha(i)*GMRF(S_i).Quadform(smoothCoefs_i);
    nll -= Type(0.5)*m_i*logalpha - Type(0.5)*alpha*GMRF(S_i).Quadform(smoothCoefs_i);
    
    k += m_i;
  }
  
  
  
  // We create a vector of means
  vector<Type> lin_pred(y.size());
  lin_pred = intercept + X*smoothCoefs;
  
  
  
// Probability of data conditional on fixed effect values
  int n = y.size();
  for(int i=0; i<n; i++){
    // If data is zero
    if(y(i)==0){
      nll -= log( zero_prob );
    } else {
      // Positive data contribution
      nll -= log( 1-zero_prob );
      // And the contribution of the likelihood
      if(likelihood==1)
        nll -= dinvgauss(y(i), exp(lin_pred(i)), sigma, true);
      else if(likelihood==2)
        nll -= dlognorm(y(i), lin_pred(i), sigma, true);
      else {
        std::cout << "Invalid likelihood specified" << std::endl;
        return 0;
      }
    }
  }
  

  

  
  
  //==========================================
  // Derived quantities
  vector<Type> pred = lin_pred;
  vector<Type> splines2D = tpsReport*smoothCoefs;
  
  // Jacobian adjustment for transformed parameters
  nll -= logsigma + logalpha;   // add logalpha? how?
  
  // Calculate joint negative log likelihood
  Type jnll = nll + nlp;
  
  
  //=====================================================================================================
  // REPORT
  REPORT(jnll);
  
  REPORT(lin_pred);
  REPORT(theta);
  REPORT(logsigma);
  REPORT(logalpha);
  REPORT(smoothCoefs);
  REPORT(pred);
  REPORT(splines2D);
  
  //=====================================================================================================
  // AD REPORT
  ADREPORT(lin_pred);		            
  ADREPORT(theta);
  ADREPORT(logsigma);	
  ADREPORT(logalpha);
  ADREPORT(smoothCoefs);
  ADREPORT(splines2D);
  return jnll;
}
