
// include libraries
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
#include <string>
using namespace density;
using Eigen::SparseMatrix;
using namespace R_inla; //includes SPDE-spesific functions, e.g. Q_spde()


//// ---------- Custom likelihood functions, used be used in template
//// below. These are not built into TMB like dnorm and dgamma are.
//log-normal likelihood
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}
// Inverse Gaussian
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





//// ---------- The CPUE model
template<class Type>
Type objective_function<Type>::operator() ()
{
  
  
//==============================================================================
// DATA SECTION
//==============================================================================
  
  DATA_INTEGER(likelihood); 	// Likelihood flag
  DATA_VECTOR(y);		// Observed catches (can be zero)
  DATA_FACTOR(site);		// Random effect index for observation i
  DATA_STRUCT(spdeMatrices,spde_t); //Three matrices needed for representing the GMRF, see p. 8 in Lindgren et al. (2011)


//==============================================================================
// PARAMETER SECTION
//==============================================================================

  // Parameters
  PARAMETER(intercept);		// intercept
  PARAMETER(theta); 		// logit of zero_prob
  PARAMETER(logsigma);		// log of observation SD
  PARAMETER(logtau);	// spatial variance
  PARAMETER(logkappa);		// kind of the decorelation
  
  PARAMETER_VECTOR(u);		// spatial random effects
  
  
  //======================================
  // Transformed quantities and parameters
  //======================================
  Type zero_prob = 1 / (1 + exp(-theta));
  Type sigma = exp(logsigma);
  Type tau = exp(logtau);
  Type kappa = exp(logkappa);
  
  
  //Construct sparce precision matrix for latent field---
  SparseMatrix<Type> Q = Q_spde(spdeMatrices, kappa);
    


  //===================================
  //               Priors
  //===================================
  Type nlp = Type(0.0);                          // negative log prior  (priors)

  // Intercept
  nlp -= dnorm(intercept, Type(0.0),   Type(5.0), true);


  // Theta
  nlp -= dcauchy(theta,   Type(0.0), Type(2.0));
  
  // Variance component
  nlp -= dcauchy(sigma,   Type(0.0), Type(2.0));
  
  
  // Hyperpriors
  nlp -= dlognorm(tau,    Type(0.0),   Type(2.0), true);
  nlp -= dlognorm(kappa,  Type(0.0),   Type(2.0), true);
  

  

//======================================
// Initialize likelihood
//======================================
  Type nll_spatial=0;	 // 
  Type nll =0;	 //

  // Linear predictor
  int n = y.size();
  vector<Type> lin_pred(n);
  for(int i=0; i<n; i++){
    //pred(i) = intercept + beta_lat*lat(i) + beta_lon*lon(i) + u(site(i));
    lin_pred(i) = intercept + u(site(i));
  }
  // Probability of random effects (hyperdistribution). This replaces neg_log_density (MVN)
  // Spatial random effects likelihood
  //nll += SCALE(GMRF(Q), 1/tau)(u); // returns negative already
  nll_spatial += GMRF(Q)(u); // returns negative already


  // Probability of data conditional on fixed effect values
  vector<Type> log_lik(n);
  for(int i=0; i<n; i++){
    // If data is zero
    if(y(i)==0){
      nll -= log( zero_prob );
      //log_lik(i) = log(zero_prob);
    } else {
      // Positive data contribution
      nll -= log( 1-zero_prob );
      //log_lik(i) = log(1-zero_prob);
      // And the contribution of the likelihood
      if(likelihood==1)
	       nll -= dinvgauss(y(i), exp(lin_pred(i)), sigma, true);
         //log_lik(i) = dinvgauss(y(i), exp(lin_pred(i)), sigma, true);
      else if(likelihood==2)
	       nll -= dlognorm(y(i), lin_pred(i), sigma, true);
         //log_lik(i) = dlognorm(y(i), lin_pred(i), sigma, true);
      else {
      	std::cout << "Invalid likelihood specified" << std::endl;
        return 0;
      }
    }
  }
  
  //Type nll = -log_lik.sum(); // total NLL

  

  // Geospatial derived quantities, given parameters
  Type range = sqrt(8) / kappa;
  Type sigma_s = 1 / sqrt(4 * M_PI * 2*tau * 2*kappa);



// Predictions
   vector<Type> pred = lin_pred;
  
// Jacobian adjustment for transformed parameters
   nll -= logsigma + logtau + logkappa;   // add logalpha? how?
  
  // Calculate joint negative log likelihood
  Type jnll = nll + nll_spatial + nlp;

  // Reporting
  REPORT(zero_prob);
  REPORT(pred);
  REPORT(logsigma);
  REPORT(logtau);
  REPORT(logkappa);
  REPORT(u);
  REPORT(range);
  REPORT(sigma_s);
  return jnll;
}
