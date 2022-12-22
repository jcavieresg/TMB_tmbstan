// Estimating parameters in a Tweedie distribution.
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
  DATA_VECTOR(y);
  
  PARAMETER(mu);
  PARAMETER(phi);
  PARAMETER(p);
  
  
//==========================================
// Transformed parameters
  // Type phi  = exp(logphi);
  // Type p = exp(logp);
  
//==========================================  
// Priors
  Type nlp = Type(0.0);                                 // negative log prior  (priors)
  
  nlp-= dcauchy(mu,    Type(0.0), Type(5.0), true);
  
  // Parameter Tweddie distribution
  nlp-= dcauchy(phi,   Type(1.0),   Type(5.0));  
  
  
  Type nll = 0;
  for(int i=0; i<y.size(); i++)
    nll -= dtweedie(y(i), mu, phi, p, true);
  
  
  
// Jacobian adjustment for transformed parameters
   //nll -= logphi;   // add logalpha? how?
  
  // Calculate joint negative log likelihood
  Type jnll = nll + nlp;
  
  //return nll;
  return jnll;
}