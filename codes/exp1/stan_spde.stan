
//========================
// DATA
//========================
data {
  int<lower = 1> n;    // n obs
  int<lower = 1> p;    // n par
  int<lower = 1> row_spar;    // rows sparse matrix 
  int<lower = 1> col_spar;    // cols sparse matrix 
  
  vector[n] time;      //The response
  int notcens[n];      //indicator vector stating which persons that were censored
  matrix[n,p] X;       //Design matrix for fixed effects
  // matrix[row_spar, col_spar] M0;     // SPDE matrices from INLA
  // matrix[row_spar, col_spar] M1;
  // matrix[row_spar, col_spar] M2;
  matrix[row_spar, col_spar] Q;     // SPDE matrices from INLA
  matrix[n, col_spar] A;     //Matrix for interpolating points witin triangles 
}

//========================
// PARAMETERS
//========================
parameters {
  vector[p] beta;  
  real log_tau;
  real log_kappa;
  real log_omega;  
  vector[row_spar] x;   // spatial random effect  
 }


//========================
// T.PARAMETERS
//========================
transformed parameters {
  real tau = exp(log_tau);
  real kappa = exp(log_kappa);
  real omega = exp(log_omega);  // Parameter of Weibull distribution
//------------------------------------------
  vector[n] eta;
  //matrix[row_spar, col_spar] Q;
  vector[row_spar] zeroes = rep_vector(0, row_spar);
  vector[n] delta;
  delta = (A*x)/tau;
  eta = X*beta + delta;
  //Q = kappa^4*M0 + 2.0*kappa^2*M1 + M2;
   
}


//========================
// MODEL
//========================
model {
  for(i in 1:n){    
    real lambda = exp(eta[i]);
    real t_omega = pow(time[i],omega);
    real S = exp(-lambda*t_omega);            // Survival function
    real f = lambda*omega*t_omega/time[i]*S;  // Weibull density
    if(notcens[i]){
      target += log(f);
    }else{
      target += log(S); //The pasient survived until cencoring
    }
  }
  //---------------------------------------------
  x~multi_normal_prec(zeroes,Q);
  log_tau~normal(-2,.1);
  log_kappa~normal(2,.1);
}

//========================
// G.QUANTITIES
//========================
generated quantities{
  real range = sqrt(8)/kappa;   //Distance at which correlation has dropped to 0.1, see p. 4 in Lindgren et al. (2011)
}
