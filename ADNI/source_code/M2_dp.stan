data{
  int n; // number of subjects in testing dataset
  int J;
  int nobs; // number of observations prior or equal to time T in testing dataset
  int id_long[nobs]; // ID
  int L0; // number of eigenvalues for U_i(t)
  int L1; // number of eigenvalues for W_ij(t)
  int Lm; // number of eigenfunctions for f_mi(v)
  int P; // number of cubic B-spline functions
  int P_surv;
  
  real Y1[nobs]; // Y1
  real Y2[nobs]; // Y2
  real Y3[nobs];
  real Y4[nobs];
  real Y5[nobs];
  real time[nobs]; // observed time
  matrix[n, P_surv] x;
  real surv_time[n]; // observed survival time
  matrix[nobs, P] b; // orthogonal cubic B-spline matrix
  matrix[nobs, L0] phi; // estimated eigenfunctions for U_i(t)
  matrix[nobs, L1] psi; // estimated eigenfunctions for W_ij(t)
  matrix[n, Lm] m_mat; // m_{ij} matrix
  matrix[L0, Lm] f_l; // \int_V phi_{ml}(v) psi_{mj}(v) dv
  real beta_m; 
  
  int<lower=0> Ltau; // number of knots for baseline hazard function
  real tau[Ltau]; // knots for time 
  vector[Ltau+1] h_grid; // observed survival time covariate for survival outcome
  
  int<lower=0> G; // number of quadrature points
  vector[G] w; // weights for quadrature
  vector[Ltau] constant; // 
  row_vector[G] phi1_interval[Ltau]; // phi_1 evaluated at G quadrature points for E intervals
  row_vector[G] phi2_interval[Ltau]; // phi_2 evaluated at G quadrature points for E intervals
  row_vector[G] psi1_interval[Ltau]; // psi_1 evaluated at G quadrature points for E intervals

  vector[P] A1; // estimated coefficients for mu1
  vector[P] A2; // estimated coefficients for mu2
  vector[P] A3;
  vector[P] A4;
  vector[P] A5;
  vector<lower=0>[L0] sqrt_d0; // estimated standard deviation for FPC scores xi_il
  vector<lower=0>[L1] sqrt_d1; // estimated standard deviation for FPC scores zeta_ijl

  vector[J-1] beta;  // estimated beta 
  vector[J] sigma;
  
  vector[Ltau] logh0; 
  vector[P_surv] gamma_x;
  real gamma0; // coefficients for U_i(t)
  real gamma1[J];  // coefficients for W_{ij}(t)
  vector[Lm] gamma_m; 
}
parameters{
  matrix[L0, n] xi; // FPC scores for U_i(t)
  matrix[L1, n] zeta_1; // FPC scores for W_i1(t)
  matrix[L1, n] zeta_2; // FPC scores for W_i2(t)
  matrix[L1, n] zeta_3; // FPC scores for W_i3(t)
  matrix[L1, n] zeta_4;
  matrix[L1, n] zeta_5;
}
transformed parameters{
  real mu1[nobs]; // mean for Y1
  real mu2[nobs]; // mean for Y2
  real mu3[nobs];
  real mu4[nobs];
  real mu5[nobs];
  matrix[n, Lm] mu_matrix = beta_m*(xi'*f_l); // beta_m * (i,j-th element is sum_l xi_{il}*f_{lj})
  matrix[n, Lm] xi_m = m_mat - mu_matrix; // FPC scores for f_mi(v)
  
  real ll1[nobs];
  real ll2[nobs];
  real ll3[nobs];
  real ll4[nobs];
  real ll5[nobs];
  
  vector[L0] tmp_xi; // xi[, i]: FPC score for U_i(t)
  vector[L1] tmp_zeta1; // zeta1[, i]: FPC scores for W_i1(t)
  vector[L1] tmp_zeta2; // zeta2[, i]: FPC scores for W_i2(t)
  vector[L1] tmp_zeta3; // zeta3[, i]: FPC scores for W_i3(t)
  vector[L1] tmp_zeta4;
  vector[L1] tmp_zeta5;
  
  matrix[n, Ltau] H; // cumulative hazard in interval (k, k+1)
  real LL[n]; // log survival likelihood
  row_vector[G] U; // for subject i, interval e, U = tmp_xi[1]*phi_surv[i, e] 
  row_vector[G] W1; // for subject i, interval e, W1 = tmp_zeta1*psi_surv[i, e]
  row_vector[G] W2;
  row_vector[G] W3;
  row_vector[G] W4;
  row_vector[G] W5;
  row_vector[G] f; // f = exp(\gamma_0*U_i(t) + \sum_{j=1}^J \gamma_{1j}*W_{ij}(t))

  for (i in 1:nobs){
    tmp_xi = col(xi, id_long[i]); // FPC scores xi_{i} for subject i
    
    tmp_zeta1 = col(zeta_1, id_long[i]); // FPC scores zeta_{i1} for subject i
    tmp_zeta2 = col(zeta_2, id_long[i]); // FPC scores zeta_{i2} for subject i
    tmp_zeta3 = col(zeta_3, id_long[i]);
    tmp_zeta4 = col(zeta_4, id_long[i]);
    tmp_zeta5 = col(zeta_5, id_long[i]);
    
    mu1[i] = b[i]*A1 + phi[i]*tmp_xi + psi[i]*tmp_zeta1;
    mu2[i] = b[i]*A2 + beta[1]*(phi[i]*tmp_xi + psi[i]*tmp_zeta2);
    mu3[i] = b[i]*A3 + beta[2]*(phi[i]*tmp_xi + psi[i]*tmp_zeta3);
    mu4[i] = b[i]*A4 + beta[3]*(phi[i]*tmp_xi + psi[i]*tmp_zeta4);
    mu5[i] = b[i]*A5 + beta[4]*(phi[i]*tmp_xi + psi[i]*tmp_zeta5);
    
    if (Y1[i]==-100) ll1[i] = 0; else ll1[i] = skew_normal_lpdf(Y1[i] | mu1[i], sigma[1], 0);
    if (Y2[i]==-100) ll2[i] = 0; else ll2[i] = skew_normal_lpdf(Y2[i] | mu2[i], sigma[2], 0);
    if (Y3[i]==-100) ll3[i] = 0; else ll3[i] = skew_normal_lpdf(Y3[i] | mu3[i], sigma[3], 0);
    if (Y4[i]==-100) ll4[i] = 0; else ll4[i] = skew_normal_lpdf(Y4[i] | mu4[i], sigma[4], 0);
    if (Y5[i]==-100) ll5[i] = 0; else ll5[i] = skew_normal_lpdf(Y5[i] | mu5[i], sigma[5], 0);
  }
  
  for (i in 1:n){
    tmp_xi = col(xi, i);
    tmp_zeta1 = col(zeta_1, i);
    tmp_zeta2 = col(zeta_2, i);
    tmp_zeta3 = col(zeta_3, i);
    tmp_zeta4 = col(zeta_4, i);
    tmp_zeta5 = col(zeta_5, i);
              
    for (k in 1:Ltau){
      U = tmp_xi[1]*phi1_interval[k] + tmp_xi[2]*phi2_interval[k];
      W1 = tmp_zeta1[1]*psi1_interval[k];
      W2 = tmp_zeta2[1]*psi1_interval[k];
      W3 = tmp_zeta3[1]*psi1_interval[k];
      W4 = tmp_zeta4[1]*psi1_interval[k];
      W5 = tmp_zeta5[1]*psi1_interval[k];
      f = exp(gamma0*U + gamma1[1]*W1 + gamma1[2]*W2 + gamma1[3]*W3 + gamma1[4]*W4 + gamma1[5]*W5);
      H[i, k] = exp(x[i]*gamma_x + logh0[k] + xi_m[i]*gamma_m)*constant[k]*(f*w);
    }
    
    LL[i] = sum(-H[i]);
  }
}
model{
  for (i in 1:L0){
    xi[i] ~ normal(0, sqrt_d0[i]);
  }
  
  for (i in 1:L1){
    zeta_1[i] ~ normal(0, sqrt_d1[i]);
    zeta_2[i] ~ normal(0, sqrt_d1[i]);
    zeta_3[i] ~ normal(0, sqrt_d1[i]);
    zeta_4[i] ~ normal(0, sqrt_d1[i]);
    zeta_5[i] ~ normal(0, sqrt_d1[i]);
  }
  
  target+=sum(ll1);
  target+=sum(ll2);
  target+=sum(ll3);
  target+=sum(ll4);
  target+=sum(ll5);
  target+=LL;
}
