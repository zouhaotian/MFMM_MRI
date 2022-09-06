data{
  int n; // number of subjects
  int J; 
  int nobs;
  int nmiss[J];
  int ID1[nobs]; 
  int ID2[nobs]; 
  int ID3[nobs];
  int Lb;
  
  int P_surv; // number of covariates for survival model
  int Ltau; // number of knots for longitudinal time and baseline hazard function
  real tau[Ltau];
  
  matrix[nobs, Ltau] time_spline1;
  matrix[nobs, Ltau] time_spline2;
  matrix[nobs, Ltau] time_spline3;
  matrix[n, Lb] B1; // functional coefficients matrix
  matrix[n, Lb] B2;
  matrix[n, Lb] B3;
  matrix[n, Lb] Bw;
  real time1[nobs]; 
  real time2[nobs]; 
  real time3[nobs]; 
  
  matrix[n, P_surv] w;
  
  real Y1[nobs]; // response in longitudinal: NA is repleaced by -100
  real Y2[nobs];
  real Y3[nobs];
  int miss_index1[nmiss[1]]; 
  int miss_index2[nmiss[2]];
  int miss_index3[nmiss[3]]; 
  
  real surv_time[n]; // survival time 
  real status[n]; // survival status
  matrix[n, Ltau+1] h_grid; // observed survival time covariate for survival outcome
  matrix[n, Ltau] h_index; // observed time spline coefficient index for survival outcome
  
  vector[2] zero;
}
parameters{
  vector[1+Ltau] beta1;
  vector[1+Ltau] beta2;
  vector[1+Ltau] beta3;
  vector[Lb] BX1;
  vector[Lb] BX2;
  vector[Lb] BX3;
  vector[Lb] BW;
  
  vector[2] u[n];
  vector<lower=0>[2] sigma_u;
  real<lower=-1, upper=1> rho;
  
  real v2[2];
  real v3[2];
  
  real<lower=0> omega[J];
  
  vector[Ltau] logh0; 
  vector[P_surv] gamma;
  vector[J] alpha;
  
  real Y1_imp[nmiss[1]]; 
  real Y2_imp[nmiss[2]];
  real Y3_imp[nmiss[3]];
}
transformed parameters{
  cov_matrix[2] Sigma;
  real Y1_full[nobs] = Y1; 
  real Y2_full[nobs] = Y2;
  real Y3_full[nobs] = Y3;
  real mu1[nobs]; // Mean of longitudinal response
  real mu2[nobs];
  real mu3[nobs];
  real f1[n]; //sum of BX1*B_mat
  real f2[n];
  real f3[n];
  real fw[n];
  real g1[n];
  real g2[n];
  real g3[n];
  real rs_1[n]; // risk score for Y1
  real rs_2[n]; // risk score for Y2
  real rs_3[n]; // risk score for Y3
  real rs_w[n]; // risk score for survival
  
  matrix[n, Ltau] coef_t; 
  vector[Ltau] const_t;
  real logh0_obs_t[n];
  real coef_obs_t[n];
  real const_obs_t[n];
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
  
  const_t[1] = 0;
  for (k in 2:Ltau){
    const_t[k] = const_t[k-1] - tau[k]*(alpha[1]*beta1[k] + alpha[2]*beta2[k] + alpha[3]*beta3[k]);
  }
 
  for (i in 1:n){
    f1[i] = B1[i]*BX1;
    f2[i] = B2[i]*BX2;
    f3[i] = B3[i]*BX3;
    fw[i] = Bw[i]*BW;
    g1[i] = beta1[1] + f1[i] + u[i, 1];
    g2[i] = beta2[1] + f2[i] + v2[1]*u[i, 1];
    g3[i] = beta3[1] + f3[i] + v3[1]*u[i, 1];
    rs_1[i] = alpha[1]*g1[i];
    rs_2[i] = alpha[2]*g2[i];
    rs_3[i] = alpha[3]*g3[i];
    rs_w[i] = w[i]*gamma + fw[i];
    
    coef_t[i, 1] = alpha[1]*(beta1[2]+u[i, 2]) + alpha[2]*(beta2[2]+v2[2]*u[i, 2]) + alpha[3]*(beta3[2]+v3[2]*u[i, 2]);
    for (k in 2:Ltau){
      coef_t[i, k] = coef_t[i, k-1] + alpha[1]*beta1[k+1] + alpha[2]*beta2[k+1] + alpha[3]*beta3[k+1];
    }
    
    logh0_obs_t[i] = h_index[i]*logh0;
    coef_obs_t[i] = h_index[i]*coef_t[i]';
    const_obs_t[i] = h_index[i]*const_t;
  }
  
  for (i in 1:nobs){
    mu1[i] = g1[ID1[i]] + time_spline1[i]*beta1[2:(Ltau+1)] + u[ID1[i], 2]*time1[i];
  }
  
  for (i in 1:nobs){
    mu2[i] = g2[ID2[i]] + time_spline2[i]*beta2[2:(Ltau+1)] + v2[2]*u[ID2[i], 2]*time2[i];
  }
  
  for (i in 1:nobs){
    mu3[i] = g3[ID3[i]] + time_spline3[i]*beta3[2:(Ltau+1)] + v3[2]*u[ID3[i], 2]*time3[i];
  }
  
  for (i in 1:n){
    h[i] = exp(logh0_obs_t[i] + rs_w[i] + const_obs_t[i] + 
               rs_1[i] + rs_2[i] + rs_3[i] + 
               coef_obs_t[i]*surv_time[i]);
    for (k in 1:Ltau){
      Ht[i, k] = exp(logh0[k] + rs_w[i] + const_t[k] + rs_1[i] + rs_2[i] + rs_3[i])*
                (exp(coef_t[i, k]*h_grid[i, k+1]) - exp(coef_t[i, k]*h_grid[i, k]))/coef_t[i, k];
    }
    LL[i] = status[i]*log(h[i]) + sum(-Ht[i]);
  }
  
}
model{
  beta1 ~ normal(0, 10); // prior for beta
  beta2 ~ normal(0, 10);
  beta3 ~ normal(0, 10);
  
  BX1 ~ normal(0, 10);
  BX2 ~ normal(0, 10);
  BX3 ~ normal(0, 10);
  BW ~ normal(0, 10);
  
  u ~ multi_normal(zero, Sigma);
  sigma_u ~ cauchy(0, 2.5);
  rho ~ uniform(-1, 1);
  
  v2 ~ normal(0, 10);
  v3 ~ normal(0, 10);
  
  omega ~ inv_gamma(0.1, 0.1);
  
  logh0 ~ normal(0, 10);
  gamma ~ normal(0, 10);
  alpha ~ normal(0, 10);
  
  target+=normal_lpdf(Y1_full | mu1, omega[1]);
  target+=normal_lpdf(Y2_full | mu2, omega[2]);
  target+=normal_lpdf(Y3_full | mu3, omega[3]);
  target+=LL;
}
