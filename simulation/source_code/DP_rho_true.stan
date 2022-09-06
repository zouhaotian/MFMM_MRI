data{
  int n; // number of subjects
  int J; 
  int nobs;
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
  real f1[n]; //sum of BX1*B_mat
  real f2[n];
  real f3[n];
  real fw[n];
  real time1[nobs]; 
  real time2[nobs]; 
  real time3[nobs]; 
  
  matrix[n, P_surv] w;
  
  real Y1[nobs]; // response in longitudinal: NA is repleaced by -100
  real Y2[nobs];
  real Y3[nobs];
  
  real surv_time[n]; // survival time 
  real surv_time2[n]; // t + delta t
  real status[n]; // survival status
  matrix[n, Ltau+1] h_grid; // observed survival time covariate for survival outcome
  matrix[n, Ltau] h_index; // observed time spline coefficient index for survival outcome
  matrix[n, Ltau+1] h_grid2;
  matrix[n, Ltau] h_index2;
  
  vector[2] zero;
  
  vector[1+Ltau] beta1;
  vector[1+Ltau] beta2;
  vector[1+Ltau] beta3;

  vector<lower=0>[2] sigma_u;
  real<lower=-1, upper=1> rho;
  
  real v2[2];
  real v3[2];
  
  real<lower=0> omega[J];
  
  vector[Ltau] logh0; 
  vector[P_surv] gamma;
  vector[J] alpha;
  
  vector[2] u[n];
}
parameters{
  real z;
}
transformed parameters{
  cov_matrix[2] Sigma;
  real mu1[nobs]; // Mean of longitudinal response
  real mu2[nobs];
  real mu3[nobs];
  real ll1[nobs];
  real ll2[nobs];
  real ll3[nobs];
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
  
  real logh0_obs_t2[n];
  real coef_obs_t2[n];
  real const_obs_t2[n];
  
  real h[n];
  real h2[n];
  matrix[n, Ltau] Ht; // cumulative hazard in interval (k, k+1)
  matrix[n, Ltau] Ht2; 
  real LL[n];
  
  Sigma[1, 1] = sigma_u[1]*sigma_u[1];
  Sigma[1, 2] = rho*sigma_u[1]*sigma_u[2];
  Sigma[2, 1] = Sigma[1, 2];
  Sigma[2, 2] = sigma_u[2]*sigma_u[2];
  
  const_t[1] = 0;
  for (k in 2:Ltau){
    const_t[k] = const_t[k-1] - tau[k]*(alpha[1]*beta1[k] + alpha[2]*beta2[k] + alpha[3]*beta3[k]);
  }
 
  for (i in 1:n){
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
    
    logh0_obs_t2[i] = h_index2[i]*logh0;
    coef_obs_t2[i] = h_index2[i]*coef_t[i]';
    const_obs_t2[i] = h_index2[i]*const_t;
  }
  
  for (i in 1:nobs){
    mu1[i] = g1[ID1[i]] + time_spline1[i]*beta1[2:(Ltau+1)] + u[ID1[i], 2]*time1[i];
    if (Y1[i]==-100) ll1[i] = 0; else ll1[i] = normal_lpdf(Y1[i] | mu1[i], omega[1]);
  }
  
  for (i in 1:nobs){
    mu2[i] = g2[ID2[i]] + time_spline2[i]*beta2[2:(Ltau+1)] + v2[2]*u[ID2[i], 2]*time2[i];
    if (Y2[i]==-100) ll2[i] = 0; else ll2[i] = normal_lpdf(Y2[i] | mu2[i], omega[2]);
  }
  
  for (i in 1:nobs){
    mu3[i] = g3[ID3[i]] + time_spline3[i]*beta3[2:(Ltau+1)] + v3[2]*u[ID3[i], 2]*time3[i];
    if (Y3[i]==-100) ll3[i] = 0; else ll3[i] = normal_lpdf(Y3[i] | mu3[i], omega[3]);  
  }
  
  for (i in 1:n){
    h[i] = exp(logh0_obs_t[i] + rs_w[i] + const_obs_t[i] + 
               rs_1[i] + rs_2[i] + rs_3[i] + 
               coef_obs_t[i]*surv_time[i]);
    h2[i] = exp(logh0_obs_t2[i] + rs_w[i] + const_obs_t2[i] + 
                rs_1[i] + rs_2[i] + rs_3[i] + 
                coef_obs_t2[i]*surv_time2[i]);
    for (k in 1:Ltau){
      Ht[i, k] = exp(logh0[k] + rs_w[i] + const_t[k] + rs_1[i] + rs_2[i] + rs_3[i])*
                (exp(coef_t[i, k]*h_grid[i, k+1]) - exp(coef_t[i, k]*h_grid[i, k]))/coef_t[i, k];
      Ht2[i, k] = exp(logh0[k] + rs_w[i] + const_t[k] + rs_1[i] + rs_2[i] + rs_3[i])*
                 (exp(coef_t[i, k]*h_grid2[i, k+1]) - exp(coef_t[i, k]*h_grid2[i, k]))/coef_t[i, k];

    }
    LL[i] = status[i]*log(h[i]) + sum(-Ht[i]);
  }
}
model{
  z ~ normal(0, 1);
}
generated quantities{
  real S[n];
  real S2[n];
  real cond_S[n];
  for (i in 1:n){
    S[i] = exp(sum(-Ht[i]));
    S2[i] = exp(sum(-Ht2[i]));
    cond_S[i] = S2[i]/S[i];
  }
}
