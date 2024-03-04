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

template<class Type>
Type ldhalfnorm(Type x, Type var){
  return 0.5*log(2)-0.5*log(var*M_PI)+pow(x,2)/(2*var);
}



// helper function to make sparse SPDE precision matrix
// Inputs:
//    logkappa: log(kappa) parameter value
//    logtau: log(tau) parameter value
//    M0, M1, M2: these sparse matrices are output from R::INLA::inla.spde2.matern()$param.inla$M*
// template<class Type>
// SparseMatrix<Type> spde_Q(Type logkappa, Type logtau, SparseMatrix<Type> M0,
//                           SparseMatrix<Type> M1, SparseMatrix<Type> M2) {
//   SparseMatrix<Type> Q;
//   Type kappa2 = exp(2. * logkappa);
//   Type kappa4 = kappa2*kappa2;
//   Q = pow(exp(logtau), 2.)  * (kappa4*M0 + Type(2.0)*kappa2*M1 + M2);
//   return Q;
// }


//=====================================================
// Initialize the TMB model
//=====================================================
template<class Type>
Type objective_function<Type>::operator() ()
{
  
  
  //=========================
  //      DATA SECTION
  //=========================
  //DATA_INTEGER(model); 	          // type of spatial model
  //DATA_INTEGER(normalization);    // scalling or not the GMRF (considering the param of Q)
  DATA_VECTOR(y);		              // Response
  //DATA_MATRIX(X);                 // Fixed effects
  DATA_VECTOR(x1);                 // Fixed effects
  DATA_FACTOR(site);		          // Random effect index for observation i
  DATA_STRUCT(spde_mat, spde_t);  // Three matrices needed for representing the GMRF, see p. 8 in Lindgren et al. (2011)
  
  DATA_SPARSE_MATRIX(A);          // used to project from spatial mesh to data locations
  
  // DATA_SCALAR(tau0);
  // DATA_SCALAR(kappa0);
  
  
  //=========================
  //   PARAMETER SECTION
  //=========================
  
  // Fixed effects
  //PARAMETER_VECTOR(beta);
  PARAMETER(beta0);
  PARAMETER(beta1);
  
  PARAMETER(logsigma_e);
  PARAMETER(logtau);
  PARAMETER(logkappa);
  
  PARAMETER_VECTOR(u);	          // spatial effects
  
  
  // Transformed parameters
  Type sigma_e  = exp(logsigma_e);
  Type tau = exp(logtau);
  Type kappa = exp(logkappa);
  
  
  SparseMatrix<Type> Q = Q_spde(spde_mat, kappa);
  

  // ===================================
  //               Priors
  // ===================================
  Type nlp = 0.0;
  nlp -= dnorm(beta0,     Type(1.0),     Type(3.0), true);   //
  nlp -= dnorm(beta1,     Type(1.0),     Type(3.0), true);   //
  // 
  // // Prior on sigma
  nlp -= dcauchy(sigma_e,   Type(0.0), Type(5.0), true);
  
  
  // Prior in Hyperparameters
  nlp -= dnorm(tau,      Type(0.0),   Type(0.5), true);   //
  nlp -= dnorm(kappa,    Type(0.0),   Type(0.5), true);   //
  //nlp -= dlognorm(tau,      tau0,     Type(5.0), true);   //
  //nlp -= dlognorm(kappa,    kappa0,   Type(5.0), true);   //
  // nlp -= ldhalfnorm(tau,      Type(2.0));
  // nlp -= ldhalfnorm(kappa,    Type(2.0));
  
  //nlp -= dlognorm(tau,      Type(0.0),   Type(1.0), true);   //
  //nlp -= dlognorm(kappa,    Type(0.0),   Type(1.0), true);   //
  
  //=============================================================================================================
  // Objective function is sum of negative log likelihood components
  int n = y.size();	                 // number of observations 

  Type nll = 0.0;	                   // likelihood for the observations
  Type nll_u = 0.0;		               // likelihood for the spatial effect
  
  
  // Linear predictor
  vector<Type> mu(n);
  //projS = A * u;                     // Project S at mesh pts to data pts
  // if(model==1)
  //   mu  = beta0 + beta1*x1 + ((A*u) / tau);
  // if(model==2) 
  mu  = beta0 + beta1*x1 + A*u;
  // if(model==3) 
  //   mu  = beta0 + beta1*x1 + u(site);
  
  //if(normalization==1)
  nll_u += GMRF(Q)(u); // returns negative already
  //if(normalization==2) 
  //nll_u += SCALE(GMRF(Q), 1/ tau)(u); // returns negative already
  
  
  // Probability of data conditional on fixed effect values
  for(int i=0; i<n; i++){
    // And the contribution of the likelihood
    nll -= dnorm(y(i), mu(i), sigma_e, true);
  }
  
  
  // Simule data from mu
  vector<Type> y_sim(n);
  vector<Type> srf(n);
  srf = A*u;
  for( int i=0; i<n; i++){
    SIMULATE {
      y_sim(i) = rnorm(beta0 + beta1*x1(i) + srf(i), sigma_e);
      REPORT(y_sim)};
  }
  
  
  Type rho = sqrt(8.0) / kappa;
  Type sigma_s = 1.0 / sqrt (4.0 * M_PI* 2.0*tau * 2.0 * kappa);
  
  
  // Jacobian adjustment 
  nll -= logsigma_e + logtau + logkappa;
  
  
  // Calculate joint negative log likelihood
  Type jnll = nll + nll_u + nlp;
  
  
  //=====================================================================================================
  // Reporting
  REPORT(jnll);
  //REPORT(beta);
  REPORT(beta0);
  REPORT(beta1);
  REPORT(u);
  REPORT(mu);
  REPORT(logsigma_e);
  REPORT(logtau);
  REPORT(logkappa);
  REPORT(sigma_e);
  REPORT(tau);
  REPORT(kappa);
  REPORT(rho);		         
  REPORT(sigma_s);		         
  
  //=======================================================================================================
  // AD report (standard devations)
  //ADREPORT(beta)
  ADREPORT(beta0);
  ADREPORT(beta1);
  ADREPORT(logsigma_e);
  ADREPORT(logtau);
  ADREPORT(logkappa);
  ADREPORT(rho);
  ADREPORT(sigma_s);
  
  // Derived geospatial components
  ADREPORT(rho);		               // geostatistical range
  ADREPORT(sigma_s);		         
  return jnll;
}