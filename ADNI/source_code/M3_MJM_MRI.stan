data{
  int n; // number of subjects
  int J; 
  int nobs;
  int nmiss[J];
  int id_long[nobs]; // ID
  int P_surv; // number of covariates for survival model
  int Lb;

  real time[nobs]; // observed time
  matrix[n, Lb] B1; // functional coefficients matrix
  matrix[n, P_surv] w;
  
  real Y1[nobs]; // response in longitudinal: NA is repleaced by -100
  real Y2[nobs];
  real Y3[nobs];
  real Y4[nobs];
  real Y5[nobs];
  int miss_index1[nmiss[1]]; 
  int miss_index2[nmiss[2]];
  int miss_index3[nmiss[3]]; 
  int miss_index4[nmiss[4]];
  int miss_index5[nmiss[5]]; 
  
  real surv_time[n]; // survival time 
  real status[n]; // survival status
  
  int Ltau; // number of knots for longitudinal time and baseline hazard function
  real tau[Ltau];
  matrix[n, Ltau+1] h_grid; // observed survival time covariate for survival outcome
  matrix[n, Ltau] h_index; // observed time spline coefficient index for survival outcome
  
  vector[2] zero;
  
  int s1[n];
  int e1[n];
}
parameters{
  vector[2] beta1;
  vector[2] beta2;
  vector[2] beta3;
  vector[2] beta4;
  vector[2] beta5;
  vector[Lb] BX1;
  vector[Lb] BX2;
  vector[Lb] BX3;
  vector[Lb] BX4;
  vector[Lb] BX5;
  vector[Lb] BW;
  
  vector[2] u[n];
  vector<lower=0>[2] sigma_u;
  real<lower=-1, upper=1> rho;
  
  real v2[2];
  real v3[2];
  real v4[2];
  real v5[2];
  
  real<lower=0> omega[J];
  
  vector[Ltau] logh0; 
  vector[P_surv] gamma;
  vector[J] alpha;
  
  real Y1_imp[nmiss[1]]; 
  real Y2_imp[nmiss[2]];
  real Y3_imp[nmiss[3]];
  real Y4_imp[nmiss[4]];
  real Y5_imp[nmiss[5]];
}
transformed parameters{
  cov_matrix[2] Sigma;
  
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
  
  real f1[n]; //sum of BX1*B_mat
  real f2[n];
  real f3[n];
  real f4[n];
  real f5[n];
  real fw[n];
  real g1[n];
  real g2[n];
  real g3[n];
  real g4[n];
  real g5[n];
  real rs_1[n]; // risk score for Y1
  real rs_2[n]; // risk score for Y2
  real rs_3[n]; // risk score for Y3
  real rs_4[n];
  real rs_5[n];
  real rs_w[n]; // risk score for survival
  
  real coef_t[n]; 
  real logh0_obs_t[n];
  real h[n];
  matrix[n, Ltau] Ht; // cumulative hazard in interval (k, k+1)
  real LL[n];
  
  Sigma[1, 1] = sigma_u[1]*sigma_u[1];
  Sigma[1, 2] = rho*sigma_u[1]*sigma_u[2];
  Sigma[2, 1] = Sigma[1, 2];
  Sigma[2, 2] = sigma_u[2]*sigma_u[2];
  
  Y1_full[miss_index1] = Y1_imp;
  Y2_full[miss_index2] = Y2_imp;
  Y3_full[miss_index3] = Y3_imp;
  Y4_full[miss_index4] = Y4_imp;
  Y5_full[miss_index5] = Y5_imp;
 
  for (i in 1:n){
    f1[i] = B1[i]*BX1;
    f2[i] = B1[i]*BX2;
    f3[i] = B1[i]*BX3;
    f4[i] = B1[i]*BX4;
    f5[i] = B1[i]*BX5;
    fw[i] = B1[i]*BW;
    g1[i] = beta1[1] + f1[i] + u[i, 1];
    g2[i] = beta2[1] + f2[i] + v2[1]*u[i, 1];
    g3[i] = beta3[1] + f3[i] + v3[1]*u[i, 1];
    g4[i] = beta4[1] + f4[i] + v4[1]*u[i, 1];
    g5[i] = beta5[1] + f5[i] + v5[1]*u[i, 1];
    rs_1[i] = alpha[1]*g1[i];
    rs_2[i] = alpha[2]*g2[i];
    rs_3[i] = alpha[3]*g3[i];
    rs_4[i] = alpha[4]*g4[i];
    rs_5[i] = alpha[5]*g5[i];
    rs_w[i] = w[i]*gamma;
    
    coef_t[i] = alpha[1]*(beta1[2]+u[i, 2]) + alpha[2]*(beta2[2]+v2[2]*u[i, 2]) + 
                alpha[3]*(beta3[2]+v3[2]*u[i, 2]) + alpha[4]*(beta4[2]+v4[2]*u[i, 2]) + 
                alpha[5]*(beta5[2]+v5[2]*u[i, 2]);
    
    logh0_obs_t[i] = h_index[i]*logh0;
  }
  
  for (i in 1:nobs){
    mu1[i] = g1[id_long[i]] + time[i]*beta1[2] + u[id_long[i], 2]*time[i];
    mu2[i] = g2[id_long[i]] + time[i]*beta2[2] + v2[2]*u[id_long[i], 2]*time[i];
    mu3[i] = g3[id_long[i]] + time[i]*beta3[2] + v3[2]*u[id_long[i], 2]*time[i];
    mu4[i] = g4[id_long[i]] + time[i]*beta4[2] + v4[2]*u[id_long[i], 2]*time[i];
    mu5[i] = g5[id_long[i]] + time[i]*beta5[2] + v5[2]*u[id_long[i], 2]*time[i];
  }

  for (i in 1:n){
    h[i] = exp(logh0_obs_t[i] + rs_w[i] + 
               rs_1[i] + rs_2[i] + rs_3[i] + rs_4[i] + rs_5[i] + 
               coef_t[i]*surv_time[i]);
    for (k in 1:Ltau){
      Ht[i, k] = exp(logh0[k] + rs_w[i] + rs_1[i] + rs_2[i] + rs_3[i] + rs_4[i] + rs_5[i])*
                (exp(coef_t[i]*h_grid[i, k+1]) - exp(coef_t[i]*h_grid[i, k]))/coef_t[i];
    }
    LL[i] = status[i]*log(h[i]) + sum(-Ht[i]);
  }
  
}
model{
  beta1 ~ normal(0, 10); // prior for beta
  beta2 ~ normal(0, 10);
  beta3 ~ normal(0, 10);
  beta4 ~ normal(0, 10);
  beta5 ~ normal(0, 10);
  
  BX1 ~ normal(0, 10);
  BX2 ~ normal(0, 10);
  BX3 ~ normal(0, 10);
  BX4 ~ normal(0, 10);
  BX5 ~ normal(0, 10);
  BW ~ normal(0, 10);
  
  u ~ multi_normal(zero, Sigma);
  sigma_u ~ cauchy(0, 2.5);
  rho ~ uniform(-1, 1);
  
  v2 ~ normal(0, 10);
  v3 ~ normal(0, 10);
  v4 ~ normal(0, 10);
  v5 ~ normal(0, 10);
  
  omega ~ inv_gamma(0.1, 0.1);
  
  logh0 ~ normal(0, 10);
  gamma ~ normal(0, 10);
  alpha ~ normal(0, 10);
  
  target+=skew_normal_lpdf(Y1_full | mu1, omega[1], 0);
  target+=skew_normal_lpdf(Y2_full | mu2, omega[2], 0);
  target+=skew_normal_lpdf(Y3_full | mu3, omega[3], 0);
  target+=skew_normal_lpdf(Y4_full | mu4, omega[4], 0);
  target+=skew_normal_lpdf(Y5_full | mu5, omega[5], 0);
  target+=LL;
}
generated quantities{
  real log_lik[n];
  for (i in 1:n){
    log_lik[i] = 0;
    for (k in s1[i]:e1[i]){
      log_lik[i] = log_lik[i] + skew_normal_lpdf(Y1_full[k] | mu1[k], omega[1], 0);
      log_lik[i] = log_lik[i] + skew_normal_lpdf(Y2_full[k] | mu2[k], omega[2], 0);
      log_lik[i] = log_lik[i] + skew_normal_lpdf(Y3_full[k] | mu3[k], omega[3], 0);
      log_lik[i] = log_lik[i] + skew_normal_lpdf(Y4_full[k] | mu4[k], omega[4], 0);
      log_lik[i] = log_lik[i] + skew_normal_lpdf(Y5_full[k] | mu5[k], omega[5], 0);
    } 
    log_lik[i] = log_lik[i] + LL[i]; 
  }
}
