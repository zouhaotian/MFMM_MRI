data{
  int n; // number of subjects
  int J;
  int nobs; // number of observations
  int nmiss[J]; 
  int id_long[nobs]; // ID
  int L0; // number of eigenvalues for U_i(t)
  int L1; // number of eigenvalues for W_ij(t)
  int P; // number of cubic B-spline functions
  int P_surv; // number of covariates in survival model
  
  real Y1[nobs]; // Y1
  real Y2[nobs]; // Y2
  real Y3[nobs];
  real Y4[nobs];
  real Y5[nobs];
  int miss_index1[nmiss[1]];  // missing index for Y1
  int miss_index2[nmiss[2]];
  int miss_index3[nmiss[3]]; 
  int miss_index4[nmiss[4]];
  int miss_index5[nmiss[5]]; 
  real time[nobs]; // observed time
  
  matrix[n, P_surv] x;
  
  real surv_time[n]; // survival time
  real status[n]; // censoring status
  matrix[nobs, P] b; // cubic B-spline matrix
  matrix[nobs, L0] phi; // estimated eigenfunctions for U_i(t)
  matrix[nobs, L1] psi; // estimated eigenfunctions for W_ij(t)
  
  int<lower=0> Ltau; // number of knots for baseline hazard function
  real tau[Ltau]; // knots for time 
  matrix[n, Ltau+1] h_grid; // observed survival time covariate for survival outcome
  matrix[n, Ltau] h_index; // observed time spline coefficient index for survival outcome
}
parameters{
  vector[P] A1; // coefficients for mu1
  vector[P] A2; // coefficients for mu2
  vector[P] A3;
  vector[P] A4;
  vector[P] A5;
  
  vector[J-1] beta; // beta
  real<lower=0> omega[J]; // scale parameter for error
  
  vector[Ltau] logh0; 
  vector[P_surv] gamma_x; // coefficients for x[1]-x[4]
  
  real Y1_imp[nmiss[1]]; // imputed Y1 
  real Y2_imp[nmiss[2]];
  real Y3_imp[nmiss[3]];
  real Y4_imp[nmiss[4]];
  real Y5_imp[nmiss[5]];
}
transformed parameters{
  real Y1_full[nobs] = Y1; // Full Y1 (observed + imputed)
  real Y2_full[nobs] = Y2;
  real Y3_full[nobs] = Y3;
  real Y4_full[nobs] = Y4;
  real Y5_full[nobs] = Y5;
  
  real mu1[nobs]; // m_{i1}(t)
  real mu2[nobs];
  real mu3[nobs];
  real mu4[nobs];
  real mu5[nobs];
  
  real logh0_obs_t[n]; // observed baseline hazard for subject i
  real rs_w[n]; // risk score for subject i
  real h[n]; // hazard function
  matrix[n, Ltau] H; // cumulative hazard in interval (k, k+1)
  real LL[n]; // log survival likelihood
  
  Y1_full[miss_index1] = Y1_imp;
  Y2_full[miss_index2] = Y2_imp;
  Y3_full[miss_index3] = Y3_imp;
  Y4_full[miss_index4] = Y4_imp;
  Y5_full[miss_index5] = Y5_imp;
  
  for (i in 1:nobs){
    
    mu1[i] = b[i]*A1;
    mu2[i] = b[i]*A2;
    mu3[i] = b[i]*A3;
    mu4[i] = b[i]*A4;
    mu5[i] = b[i]*A5;
  }
  
  for (i in 1:n){
    
    logh0_obs_t[i] = h_index[i]*logh0;
    rs_w[i] = x[i]*gamma_x;
    
    h[i] = exp(logh0_obs_t[i] +  rs_w[i]);
    for (k in 1:Ltau){
      H[i, k] = exp(rs_w[i])*exp(logh0[k])*(h_grid[i, k+1] - h_grid[i, k]);
    }
    
    LL[i] = status[i]*log(h[i]) + sum(-H[i]);
  }
}
model{
  A1 ~ normal(0, 10); // prior 
  A2 ~ normal(0, 10);
  A3 ~ normal(0, 10);
  A4 ~ normal(0, 10);
  A5 ~ normal(0, 10);
  
  beta ~ normal(0, 10);
  omega ~ inv_gamma(0.1, 0.1);
  
  logh0 ~ normal(0, 10);
  gamma_x ~ normal(0, 10);
  
  target+=skew_normal_lpdf(Y1_full | mu1, omega[1], 0); // likelihood
  target+=skew_normal_lpdf(Y2_full | mu2, omega[2], 0);
  target+=skew_normal_lpdf(Y3_full | mu3, omega[3], 0);
  target+=skew_normal_lpdf(Y4_full | mu4, omega[4], 0);
  target+=skew_normal_lpdf(Y5_full | mu5, omega[5], 0);
  target+=LL; // survival log likelihood
}
