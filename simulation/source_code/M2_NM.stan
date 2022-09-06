data{
  int n; // number of subjects
  int J; // number of longitudinal outcomes
  int nobs; // number of observations
  int nmiss[J]; // number of missingness
  int id_long[nobs]; // ID
  int L0; // number of eigenvalues for U_i(t)
  int L1; // number of eigenvalues for W_ij(t)
  int P; // number of cubic B-spline functions
  int P_surv; 
  
  real Y1[nobs]; // Y1
  real Y2[nobs]; // Y2
  real Y3[nobs];
  int miss_index1[nmiss[1]]; 
  int miss_index2[nmiss[2]]; 
  int miss_index3[nmiss[3]]; 
  real time[nobs]; // observed time
  
  matrix[n, P_surv] x; 
  
  real surv_time[n]; // survival time
  real status[n]; // censoring status
  matrix[nobs, P] b; // cubic B-spline matrix
}
parameters{
  vector[P] A1; // coefficients for mu1
  vector[P] A2; // coefficients for mu2
  vector[P] A3;

  vector<lower=0>[J] sigma; // sd of Y_{ij}(t)
  
  real logh0; 
  vector[P_surv] gamma_x;
   
  real Y1_imp[nmiss[1]]; 
  real Y2_imp[nmiss[2]];
  real Y3_imp[nmiss[3]];
}
transformed parameters{
  real Y1_full[nobs] = Y1; 
  real Y2_full[nobs] = Y2;
  real Y3_full[nobs] = Y3;
  
  real mu1[nobs];
  real mu2[nobs];
  real mu3[nobs];
  
  real h[n]; // hazard function
  real H[n]; // cumulative hazard function
  real LL[n]; // log survival likelihood
  
  Y1_full[miss_index1] = Y1_imp;
  Y2_full[miss_index2] = Y2_imp;
  Y3_full[miss_index3] = Y3_imp;
  
  for (i in 1:nobs){
    mu1[i] = b[i]*A1;
    mu2[i] = b[i]*A2;
    mu3[i] = b[i]*A3;
  }
  
  for (i in 1:n){
    h[i] = exp(logh0 + x[i]*gamma_x);
    H[i] = exp(logh0 + x[i]*gamma_x)*surv_time[i];
  
    LL[i] = status[i]*log(h[i]) + (-H[i]);
  }
}
model{
  A1 ~ normal(0, 10);
  A2 ~ normal(0, 10);
  A3 ~ normal(0, 10);
  
  sigma ~ inv_gamma(0.1, 0.1);
  
  logh0 ~ normal(0, 10);
  gamma_x ~ normal(0, 10);
  
  target+=normal_lpdf(Y1_full | mu1, sigma[1]);
  target+=normal_lpdf(Y2_full | mu2, sigma[2]);
  target+=normal_lpdf(Y3_full | mu3, sigma[3]);
  target+=LL;
}
