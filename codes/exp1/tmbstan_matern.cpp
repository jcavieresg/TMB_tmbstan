// 2D Gaussian process - taken from https://kaskr.github.io/adcomp/spatial_8cpp-example.html
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#include <vector>
#include <string>
using namespace density;
using Eigen::SparseMatrix;


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

// Inverse gamma
template<class Type>
Type dinvgauss(Type x, Type mean, Type shape, int give_log=0){
  Type logres = 0.5*log(shape) - 0.5*log(2*M_PI*pow(x,3)) - (shape * pow(x-mean,2) / (2*pow(mean,2)*x));
  if(give_log) return logres; else return exp(logres);
}



template<class Type>
Type objective_function<Type>::operator() ()
{
  
  
  //=========================
  //      DATA SECTION
  //=========================
  DATA_VECTOR(y);                   // Response
  DATA_MATRIX(X);                 // Model matrix for fixed effects
  DATA_MATRIX(dist_mat);            // Matrix of squared distance b/w unique locations
  

//=========================
//   PARAMETER SECTION
//=========================
// Intercept
  PARAMETER_VECTOR(beta);              

  
// Variance (nugget)
   PARAMETER(logsigma);              // variance errors
   PARAMETER(logphi);                // Spatial var
   PARAMETER(logkappa);              // Spatial decay
   
  PARAMETER_VECTOR(u);             //Spatial random effect
  
// Transformed parameters
   Type sigma = exp(logsigma);
   Type phi = exp(logphi);
   Type kappa = exp(logkappa);


//============
//  Priors
//===========
  Type nlp = 0.0;                   
 
  nlp -= dnorm(beta,     Type(0.0),     Type(10.0), true).sum();   //

  // Prior on sigma
  nlp -= dcauchy(sigma,   Type(0.0), Type(5.0), true);

 // Prior in Hyperparameters
  nlp -= dnorm(phi,     Type(0.0),     Type(2.0), true);   //
  nlp -= dnorm(kappa,   Type(0.0),     Type(2.0), true);   //
   

  

//=============================================================================================================
// Objective function is sum of negative log likelihood components
  //int n = y.size();
  Type nll = 0.0;	 
  

  matrix<Type> cov(dist_mat);
  for(int i=0; i<cov.rows(); i++)
    for(int j=0; j<cov.cols(); j++)
      cov(i,j) = matern(dist_mat(i,j), phi, kappa);
  
  nll += MVNORM_t<Type>(cov)(u);
  
  // Linear predictor "mu"
  vector<Type> mu = X*beta + u; 
  
// Likelihood
  nll -= sum(dnorm(y, mu, sigma, true)); // 
  
  
// Jacobian adjustment 
   nll -= logsigma + logphi + logkappa;
   
  
// Calculate joint negative log likelihood
   Type jnll = nll + nlp;
   
   Type rho = sqrt(8.0)/ kappa;
  
// Predictions  
   vector<Type> pred = mu;
  
//=====================================
// REPORT SECTION
//=====================================
  REPORT(mu);
  REPORT(sigma);
  REPORT(phi);
  REPORT(kappa);
  REPORT(rho);
  REPORT(cov);
  REPORT(u);
  REPORT(pred);
  return jnll;
}