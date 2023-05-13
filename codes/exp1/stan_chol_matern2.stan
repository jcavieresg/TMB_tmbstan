
data {
  int<lower=0> n;           // number of observations
  vector[n] y;              // response variable
  int<lower=0> k;           // number of columns in the model matrix
  matrix[n, k] X;           // the model matrix
  int<lower=0> n_knots;     // number of interpolation points on the mesh
  //matrix[n_knots, n_knots] spde_mat[3]; // spde matrices
  
  matrix[n_knots, n_knots] M0;  // SPDE matrices from INLA
  matrix[n_knots, n_knots] M1;
  matrix[n_knots, n_knots] M2;
  
  
  matrix[n, n_knots] A;     // projection matrix
  
}


parameters {
  // Parameters to be sampled
  vector[k] beta;         // Vector of parameters for fixed covariates
  real<lower=0> sigma_e;  // sigma for the error in the observations "y"
  real<lower=0> kappa;    // Spatial parameters of the spde method
  real<lower=0> tau;      
  vector[n_knots] u;

}


model {
  
  vector[n] grf;
  
  {
  matrix[n_knots, n_knots] Q;  
  //Q = cholesky_decompose((tau ^ 2) * (spde_mat[1] * kappa ^ 4 + spde_mat[2] * 2 * kappa ^ 2 + spde_mat[3]))'; // inverse of Q
  Q = cholesky_decompose((tau^2) * (kappa^4*M0 + 2*kappa^2*M1 + M2))'; // inverse of Q
  
  grf = rep_vector(0,  n_knots) +  mdivide_left_tri_low(Q, u);
}

  beta ~ normal(0, 2);
  u ~ std_normal();
  sigma_e ~ cauchy(0, 5);
  kappa ~ lognormal(0, 0.5);
  tau ~ lognormal(0, 0.5);
  y ~ normal(X*beta + A*grf, sigma_e);
  
}





