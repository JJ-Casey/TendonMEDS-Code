// generated with brms 2.20.4
functions {
 /* compute correlated group-level effects
  * Args:
  *   z: matrix of unscaled group-level effects
  *   SD: vector of standard deviation parameters
  *   L: cholesky factor correlation matrix
  * Returns:
  *   matrix of scaled group-level effects
  */
  matrix scale_r_cor(matrix z, vector SD, matrix L) {
    // r is stored in another dimension order than z
    return transpose(diag_pre_multiply(SD, L) * z);
  }
  
real tendonModel (real stretch, real stretchm1, real stretch2, real stretchmstretch2m1, real stretchLogStretch, real nu, real eta, real tau, real rho) {
  real phimu = exp(1.05309738 + 1.30927056 * nu);
  real phiE = exp(6.83672018 + 0.47191773 * eta);
  real a = exp(-3.80045123 + 0.64387023 * tau) + 1;
  real b = exp(-3.59771868 + 0.7310165 * rho) + a;
  real c = (a + b) / 2;

  real index1 = (stretch2 < a^2);
  real index2 = (stretch2 >= a^2) * (stretch2 < c^2);
  real index3 = (stretch2 >= c^2) * (stretch2 < b^2);
  real index4 = (stretch2 >= b^2);

  // A + Bl + Cl^2 + DlLogl
  real part1 = phimu * stretchmstretch2m1;
  real part2 = index2 * (-(a ^ 2) + 2 * a * log(a) * stretch + stretch2 + -2 * a * stretchLogStretch) / ((b - a) * (c - a));
  real part3 = index3 * ((c ^ 2 / ((c - a) * (b - c)) - (a ^ 2) / ((b - a) * (c - a))) +
                              (2 * a * log(a) / ((b - a) * (c - a)) - 2 * c * log(c) / ((b - c) * (c - a))) * stretch +
                              (-1 / ((b - a) * (b - c))) * stretch2 +
                              (2 * b / ((b - a) * (b - c))) * stretchLogStretch);
  real part4 = index4 * (-1 +
                              (2 * a * log(a) / ((b - a) * (c - a))
                                  + 2 * b * log(b) / ((b - c) * (b - a))
                                  - 2 * c * log(c) / ((b - c) * (c - a))) * stretch);

  return part1 + phiE * stretchm1 * (part2 + part3 + part4);
}

}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  vector<lower=0>[N] se;  // known sampling error
  int<lower=1> K_Nu;  // number of population-level effects
  matrix[N, K_Nu] X_Nu;  // population-level design matrix
  int<lower=1> K_eta;  // number of population-level effects
  matrix[N, K_eta] X_eta;  // population-level design matrix
  int<lower=1> K_tau;  // number of population-level effects
  matrix[N, K_tau] X_tau;  // population-level design matrix
  int<lower=1> K_rho;  // number of population-level effects
  matrix[N, K_rho] X_rho;  // population-level design matrix
  // covariates for non-linear functions
  vector[N] C_1;
  vector[N] C_2;
  vector[N] C_3;
  vector[N] C_4;
  vector[N] C_5;
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  array[N] int<lower=1> J_1;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_Nu_1;
  vector[N] Z_1_eta_2;
  vector[N] Z_1_tau_3;
  vector[N] Z_1_rho_4;
  int<lower=1> NC_1;  // number of group-level correlations
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  vector<lower=0>[N] se2 = square(se);
}
parameters {
  vector[K_Nu] b_Nu;  // regression coefficients
  vector[K_eta] b_eta;  // regression coefficients
  vector[K_tau] b_tau;  // regression coefficients
  vector[K_rho] b_rho;  // regression coefficients
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  matrix[M_1, N_1] z_1;  // standardized group-level effects
  cholesky_factor_corr[M_1] L_1;  // cholesky factor of correlation matrix
}
transformed parameters {
  real sigma = 0;  // dispersion parameter
  matrix[N_1, M_1] r_1;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_1] r_1_Nu_1;
  vector[N_1] r_1_eta_2;
  vector[N_1] r_1_tau_3;
  vector[N_1] r_1_rho_4;
  // compute actual group-level effects
  r_1 = scale_r_cor(z_1, sd_1, L_1);
  r_1_Nu_1 = r_1[, 1];
  r_1_eta_2 = r_1[, 2];
  r_1_tau_3 = r_1[, 3];
  r_1_rho_4 = r_1[, 4];
}
model {
  real lprior = 0;  // prior contributions to the log posterior
  // likelihood not including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] nlp_Nu = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] nlp_eta = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] nlp_tau = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] nlp_rho = rep_vector(0.0, N);
    // initialize non-linear predictor term
    vector[N] mu;
    nlp_Nu += X_Nu * b_Nu;
    nlp_eta += X_eta * b_eta;
    nlp_tau += X_tau * b_tau;
    nlp_rho += X_rho * b_rho;
    for (n in 1:N) {
      // add more terms to the linear predictor
      nlp_Nu[n] += r_1_Nu_1[J_1[n]] * Z_1_Nu_1[n];
    }
    for (n in 1:N) {
      // add more terms to the linear predictor
      nlp_eta[n] += r_1_eta_2[J_1[n]] * Z_1_eta_2[n];
    }
    for (n in 1:N) {
      // add more terms to the linear predictor
      nlp_tau[n] += r_1_tau_3[J_1[n]] * Z_1_tau_3[n];
    }
    for (n in 1:N) {
      // add more terms to the linear predictor
      nlp_rho[n] += r_1_rho_4[J_1[n]] * Z_1_rho_4[n];
    }
    for (n in 1:N) {
      // compute non-linear predictor values
      mu[n] = (tendonModel(C_1[n] , C_2[n] , C_3[n] , C_4[n] , C_5[n] , nlp_Nu[n] , nlp_eta[n] , nlp_tau[n] , nlp_rho[n]));
    }
    target += normal_lupdf(Y | mu, se);
  }
  // priors not including constants
  lprior += std_normal_lupdf(b_Nu);
  lprior += std_normal_lupdf(b_eta);
  lprior += std_normal_lupdf(b_tau);
  lprior += std_normal_lupdf(b_rho);
  lprior += student_t_lupdf(sd_1 | 3, 0, 1);
  lprior += lkj_corr_cholesky_lupdf(L_1 | 1);
  target += lprior;
  target += std_normal_lupdf(to_vector(z_1));
}
generated quantities {
  // compute group-level correlations
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
  vector<lower=-1,upper=1>[NC_1] cor_1;
  // extract upper diagonal of correlation matrix
  for (k in 1:M_1) {
    for (j in 1:(k - 1)) {
      cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }
}
