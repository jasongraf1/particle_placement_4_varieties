// generated with brms 2.18.0
functions {
  /* zero-one-inflated beta log-PDF of a single response
   * Args:
   *   y: response value
   *   mu: mean parameter of the beta part
   *   phi: precision parameter of the beta part
   *   zoi: zero-one-inflation probability
   *   coi: conditional one-inflation probability
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real zero_one_inflated_beta_lpdf(real y, real mu, real phi,
                                    real zoi, real coi) {
     row_vector[2] shape = [mu * phi, (1 - mu) * phi];
     if (y == 0) {
       return bernoulli_lpmf(1 | zoi) + bernoulli_lpmf(0 | coi);
     } else if (y == 1) {
       return bernoulli_lpmf(1 | zoi) + bernoulli_lpmf(1 | coi);
     } else {
       return bernoulli_lpmf(0 | zoi) + beta_lpdf(y | shape[1], shape[2]);
     }
   }
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> K_phi;  // number of population-level effects
  matrix[N, K_phi] X_phi;  // population-level design matrix
  int<lower=1> K_zoi;  // number of population-level effects
  matrix[N, K_zoi] X_zoi;  // population-level design matrix
  int<lower=1> K_coi;  // number of population-level effects
  matrix[N, K_coi] X_coi;  // population-level design matrix
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  // data for group-level effects of ID 2
  int<lower=1> N_2;  // number of grouping levels
  int<lower=1> M_2;  // number of coefficients per level
  int<lower=1> J_2[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_2_1;
  // data for group-level effects of ID 3
  int<lower=1> N_3;  // number of grouping levels
  int<lower=1> M_3;  // number of coefficients per level
  int<lower=1> J_3[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_3_1;
  // data for group-level effects of ID 4
  int<lower=1> N_4;  // number of grouping levels
  int<lower=1> M_4;  // number of coefficients per level
  int<lower=1> J_4[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_4_phi_1;
  // data for group-level effects of ID 5
  int<lower=1> N_5;  // number of grouping levels
  int<lower=1> M_5;  // number of coefficients per level
  int<lower=1> J_5[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_5_phi_1;
  // data for group-level effects of ID 6
  int<lower=1> N_6;  // number of grouping levels
  int<lower=1> M_6;  // number of coefficients per level
  int<lower=1> J_6[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_6_phi_1;
  // data for group-level effects of ID 7
  int<lower=1> N_7;  // number of grouping levels
  int<lower=1> M_7;  // number of coefficients per level
  int<lower=1> J_7[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_7_zoi_1;
  // data for group-level effects of ID 8
  int<lower=1> N_8;  // number of grouping levels
  int<lower=1> M_8;  // number of coefficients per level
  int<lower=1> J_8[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_8_zoi_1;
  // data for group-level effects of ID 9
  int<lower=1> N_9;  // number of grouping levels
  int<lower=1> M_9;  // number of coefficients per level
  int<lower=1> J_9[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_9_zoi_1;
  // data for group-level effects of ID 10
  int<lower=1> N_10;  // number of grouping levels
  int<lower=1> M_10;  // number of coefficients per level
  int<lower=1> J_10[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_10_coi_1;
  // data for group-level effects of ID 11
  int<lower=1> N_11;  // number of grouping levels
  int<lower=1> M_11;  // number of coefficients per level
  int<lower=1> J_11[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_11_coi_1;
  // data for group-level effects of ID 12
  int<lower=1> N_12;  // number of grouping levels
  int<lower=1> M_12;  // number of coefficients per level
  int<lower=1> J_12[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_12_coi_1;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  int Kc_phi = K_phi - 1;
  matrix[N, Kc_phi] Xc_phi;  // centered version of X_phi without an intercept
  vector[Kc_phi] means_X_phi;  // column means of X_phi before centering
  int Kc_zoi = K_zoi - 1;
  matrix[N, Kc_zoi] Xc_zoi;  // centered version of X_zoi without an intercept
  vector[Kc_zoi] means_X_zoi;  // column means of X_zoi before centering
  int Kc_coi = K_coi - 1;
  matrix[N, Kc_coi] Xc_coi;  // centered version of X_coi without an intercept
  vector[Kc_coi] means_X_coi;  // column means of X_coi before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
  for (i in 2:K_phi) {
    means_X_phi[i - 1] = mean(X_phi[, i]);
    Xc_phi[, i - 1] = X_phi[, i] - means_X_phi[i - 1];
  }
  for (i in 2:K_zoi) {
    means_X_zoi[i - 1] = mean(X_zoi[, i]);
    Xc_zoi[, i - 1] = X_zoi[, i] - means_X_zoi[i - 1];
  }
  for (i in 2:K_coi) {
    means_X_coi[i - 1] = mean(X_coi[, i]);
    Xc_coi[, i - 1] = X_coi[, i] - means_X_coi[i - 1];
  }
}
parameters {
  vector[Kc] b;  // population-level effects
  real Intercept;  // temporary intercept for centered predictors
  vector[Kc_phi] b_phi;  // population-level effects
  real Intercept_phi;  // temporary intercept for centered predictors
  vector[Kc_zoi] b_zoi;  // population-level effects
  real Intercept_zoi;  // temporary intercept for centered predictors
  vector[Kc_coi] b_coi;  // population-level effects
  real Intercept_coi;  // temporary intercept for centered predictors
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // standardized group-level effects
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  vector[N_2] z_2[M_2];  // standardized group-level effects
  vector<lower=0>[M_3] sd_3;  // group-level standard deviations
  vector[N_3] z_3[M_3];  // standardized group-level effects
  vector<lower=0>[M_4] sd_4;  // group-level standard deviations
  vector[N_4] z_4[M_4];  // standardized group-level effects
  vector<lower=0>[M_5] sd_5;  // group-level standard deviations
  vector[N_5] z_5[M_5];  // standardized group-level effects
  vector<lower=0>[M_6] sd_6;  // group-level standard deviations
  vector[N_6] z_6[M_6];  // standardized group-level effects
  vector<lower=0>[M_7] sd_7;  // group-level standard deviations
  vector[N_7] z_7[M_7];  // standardized group-level effects
  vector<lower=0>[M_8] sd_8;  // group-level standard deviations
  vector[N_8] z_8[M_8];  // standardized group-level effects
  vector<lower=0>[M_9] sd_9;  // group-level standard deviations
  vector[N_9] z_9[M_9];  // standardized group-level effects
  vector<lower=0>[M_10] sd_10;  // group-level standard deviations
  vector[N_10] z_10[M_10];  // standardized group-level effects
  vector<lower=0>[M_11] sd_11;  // group-level standard deviations
  vector[N_11] z_11[M_11];  // standardized group-level effects
  vector<lower=0>[M_12] sd_12;  // group-level standard deviations
  vector[N_12] z_12[M_12];  // standardized group-level effects
}
transformed parameters {
  vector[N_1] r_1_1;  // actual group-level effects
  vector[N_2] r_2_1;  // actual group-level effects
  vector[N_3] r_3_1;  // actual group-level effects
  vector[N_4] r_4_phi_1;  // actual group-level effects
  vector[N_5] r_5_phi_1;  // actual group-level effects
  vector[N_6] r_6_phi_1;  // actual group-level effects
  vector[N_7] r_7_zoi_1;  // actual group-level effects
  vector[N_8] r_8_zoi_1;  // actual group-level effects
  vector[N_9] r_9_zoi_1;  // actual group-level effects
  vector[N_10] r_10_coi_1;  // actual group-level effects
  vector[N_11] r_11_coi_1;  // actual group-level effects
  vector[N_12] r_12_coi_1;  // actual group-level effects
  real lprior = 0;  // prior contributions to the log posterior
  r_1_1 = (sd_1[1] * (z_1[1]));
  r_2_1 = (sd_2[1] * (z_2[1]));
  r_3_1 = (sd_3[1] * (z_3[1]));
  r_4_phi_1 = (sd_4[1] * (z_4[1]));
  r_5_phi_1 = (sd_5[1] * (z_5[1]));
  r_6_phi_1 = (sd_6[1] * (z_6[1]));
  r_7_zoi_1 = (sd_7[1] * (z_7[1]));
  r_8_zoi_1 = (sd_8[1] * (z_8[1]));
  r_9_zoi_1 = (sd_9[1] * (z_9[1]));
  r_10_coi_1 = (sd_10[1] * (z_10[1]));
  r_11_coi_1 = (sd_11[1] * (z_11[1]));
  r_12_coi_1 = (sd_12[1] * (z_12[1]));
  lprior += normal_lpdf(b | 0, 1.5);
  lprior += student_t_lpdf(Intercept | 3, 0, 2.5);
  lprior += normal_lpdf(b_phi | 0, 2);
  lprior += student_t_lpdf(Intercept_phi | 3, 0, 2.5);
  lprior += normal_lpdf(b_zoi | 0, 1.5);
  lprior += logistic_lpdf(Intercept_zoi | 0, 1);
  lprior += normal_lpdf(b_coi | 0, 1.5);
  lprior += logistic_lpdf(Intercept_coi | 0, 1);
  lprior += student_t_lpdf(sd_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sd_2 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sd_3 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sd_4 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sd_5 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sd_6 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sd_7 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sd_8 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sd_9 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sd_10 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sd_11 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sd_12 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] phi = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] zoi = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] coi = rep_vector(0.0, N);
    mu += Intercept + Xc * b;
    phi += Intercept_phi + Xc_phi * b_phi;
    zoi += Intercept_zoi + Xc_zoi * b_zoi;
    coi += Intercept_coi + Xc_coi * b_coi;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_2_1[J_2[n]] * Z_2_1[n] + r_3_1[J_3[n]] * Z_3_1[n];
    }
    for (n in 1:N) {
      // add more terms to the linear predictor
      phi[n] += r_4_phi_1[J_4[n]] * Z_4_phi_1[n] + r_5_phi_1[J_5[n]] * Z_5_phi_1[n] + r_6_phi_1[J_6[n]] * Z_6_phi_1[n];
    }
    for (n in 1:N) {
      // add more terms to the linear predictor
      zoi[n] += r_7_zoi_1[J_7[n]] * Z_7_zoi_1[n] + r_8_zoi_1[J_8[n]] * Z_8_zoi_1[n] + r_9_zoi_1[J_9[n]] * Z_9_zoi_1[n];
    }
    for (n in 1:N) {
      // add more terms to the linear predictor
      coi[n] += r_10_coi_1[J_10[n]] * Z_10_coi_1[n] + r_11_coi_1[J_11[n]] * Z_11_coi_1[n] + r_12_coi_1[J_12[n]] * Z_12_coi_1[n];
    }
    mu = inv_logit(mu);
    phi = exp(phi);
    zoi = inv_logit(zoi);
    coi = inv_logit(coi);
    for (n in 1:N) {
      target += zero_one_inflated_beta_lpdf(Y[n] | mu[n], phi[n], zoi[n], coi[n]);
    }
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(z_1[1]);
  target += std_normal_lpdf(z_2[1]);
  target += std_normal_lpdf(z_3[1]);
  target += std_normal_lpdf(z_4[1]);
  target += std_normal_lpdf(z_5[1]);
  target += std_normal_lpdf(z_6[1]);
  target += std_normal_lpdf(z_7[1]);
  target += std_normal_lpdf(z_8[1]);
  target += std_normal_lpdf(z_9[1]);
  target += std_normal_lpdf(z_10[1]);
  target += std_normal_lpdf(z_11[1]);
  target += std_normal_lpdf(z_12[1]);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  // actual population-level intercept
  real b_phi_Intercept = Intercept_phi - dot_product(means_X_phi, b_phi);
  // actual population-level intercept
  real b_zoi_Intercept = Intercept_zoi - dot_product(means_X_zoi, b_zoi);
  // actual population-level intercept
  real b_coi_Intercept = Intercept_coi - dot_product(means_X_coi, b_coi);
}
