data{
  int n; // number of subjects in testing dataset
  int J;
  int nobs; // number of observations prior or equal to time T in testing dataset
  int id_long[nobs]; // ID
  int P_surv; // number of coefficients for x_surv matrix (survival)
  int Lb;
  
  real time[nobs]; // observed time
  matrix[n, Lb] B1;
  matrix[n, P_surv] w;
  
  real Y1[nobs]; // Y1
  real Y2[nobs]; // Y2
  real Y3[nobs];
  real Y4[nobs];
  real Y5[nobs];
  
  int<lower=0> Ltau; // number of knots for baseline hazard function
  real tau[Ltau]; // knots for time 
  vector[Ltau+1] h_grid; // observed survival time covariate for survival outcome
  
  vector[2] zero;

  vector[2] beta1; // estimated coefficients for mu1
  vector[2] beta2; // estimated coefficients for mu2
  vector[2] beta3;
  vector[2] beta4;
  vector[2] beta5;
  vector[Lb] BX1;
  vector[Lb] BX2;
  vector[Lb] BX3;
  vector[Lb] BX4;
  vector[Lb] BX5;
  vector[Lb] BW;
  real sigma_u[2];
  real rho;
  real v2[2];
  real v3[2];
  real v4[2];
  real v5[2];

  vector[J] omega;
  
  vector[Ltau] logh0; 
  vector[P_surv] gamma;
  vector[J] alpha;
}
parameters{
  vector[2] u[n];
}
transformed parameters{
  cov_matrix[2] Sigma;
  
  real mu1[nobs]; // mean for Y1
  real mu2[nobs]; // mean for Y2
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
  
  real ll1[nobs];
  real ll2[nobs];
  real ll3[nobs];
  real ll4[nobs];
  real ll5[nobs];
  
  real coef_t[n]; 
  matrix[n, Ltau] Ht; // cumulative hazard in interval (k, k+1)
  real LL[n];
  
  Sigma[1, 1] = sigma_u[1]*sigma_u[1];
  Sigma[1, 2] = rho*sigma_u[1]*sigma_u[2];
  Sigma[2, 1] = Sigma[1, 2];
  Sigma[2, 2] = sigma_u[2]*sigma_u[2];

  for (i in 1:n){
    f1[i] = B1[i]*BX1;
    f2[i] = B1[i]*BX2;
    f3[i] = B1[i]*BX3;
    f4[i] = B1[i]*BX3;
    f5[i] = B1[i]*BX3;
    fw[i] = B1[i]*BW;
    g1[i] = beta1[1] + f1[i] + u[i, 1];
    g2[i] = beta2[1] + v2[1]*u[i, 1];
    g3[i] = beta3[1] + v3[1]*u[i, 1];
    g4[i] = beta4[1] + v4[1]*u[i, 1];
    g5[i] = beta5[1] + v5[1]*u[i, 1];
    
    rs_1[i] = alpha[1]*g1[i];
    rs_2[i] = alpha[2]*g2[i];
    rs_3[i] = alpha[3]*g3[i];
    rs_4[i] = alpha[4]*g4[i];
    rs_5[i] = alpha[5]*g5[i];
    rs_w[i] = w[i]*gamma;
    
    coef_t[i] = alpha[1]*(beta1[2]+u[i, 2]) + alpha[2]*(beta2[2]+v2[2]*u[i, 2]) + 
                alpha[3]*(beta3[2]+v3[2]*u[i, 2]) + alpha[4]*(beta4[2]+v4[2]*u[i, 2]) + 
                alpha[5]*(beta5[2]+v5[2]*u[i, 2]);
  }
  
  for (i in 1:nobs){
    mu1[i] = g1[id_long[i]] + time[i]*beta1[2] + u[id_long[i], 2]*time[i];
    mu2[i] = g2[id_long[i]] + time[i]*beta2[2] + v2[2]*u[id_long[i], 2]*time[i];
    mu3[i] = g3[id_long[i]] + time[i]*beta3[2] + v3[2]*u[id_long[i], 2]*time[i];
    mu4[i] = g4[id_long[i]] + time[i]*beta4[2] + v4[2]*u[id_long[i], 2]*time[i];
    mu5[i] = g5[id_long[i]] + time[i]*beta5[2] + v5[2]*u[id_long[i], 2]*time[i];
    if (Y1[i]==-100) ll1[i] = 0; else ll1[i] = skew_normal_lpdf(Y1[i] | mu1[i], omega[1], 0);
    if (Y2[i]==-100) ll2[i] = 0; else ll2[i] = skew_normal_lpdf(Y2[i] | mu2[i], omega[2], 0);
    if (Y3[i]==-100) ll3[i] = 0; else ll3[i] = skew_normal_lpdf(Y3[i] | mu3[i], omega[3], 0);
    if (Y4[i]==-100) ll4[i] = 0; else ll4[i] = skew_normal_lpdf(Y4[i] | mu4[i], omega[4], 0);
    if (Y5[i]==-100) ll5[i] = 0; else ll5[i] = skew_normal_lpdf(Y5[i] | mu5[i], omega[5], 0);
  }
  
  for (i in 1:n){
    for (k in 1:Ltau){
      Ht[i, k] = exp(logh0[k] + rs_w[i] + rs_1[i] + rs_2[i] + rs_3[i] + rs_4[i] + rs_5[i])*
                (exp(coef_t[i]*h_grid[k+1]) - exp(coef_t[i]*h_grid[k]))/coef_t[i];
    }
    LL[i] = sum(-Ht[i]);
  }
}
model{
  u ~ multi_normal(zero, Sigma);
  
  target+=sum(ll1);
  target+=sum(ll2);
  target+=sum(ll3);
  target+=sum(ll4);
  target+=sum(ll5);
  target+=LL;
}
