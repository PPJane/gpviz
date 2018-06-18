functions{
  matrix gp_pred_rng(real[] x1,
                     vector y1,
                     real[] x2,
                     real alpha,
                     real rho,
                     real sigma,
                     real delta) {
    int N1 = size(x1);
    int N2 = size(x2);
    vector[N2] f2_mu;
    vector[N2] f2;
    vector[N2] f2_prime_mu;
    vector[N2] f2_prime_var;
    matrix[3, N2] f2_stat;
    {      
      matrix[N1, N1] L_K;
      vector[N1] L_div_y1;
      matrix[N1, N2] k_x1_x2;
      matrix[N1, N2] x_diff;
      matrix[N1, N2] k_x1_x2_p;
      matrix[N1, N2] L_div_k_x1_x2;
      vector[N1] L_div_k_x1_x2_p;
      matrix[N2, N2] f2_cov;
      matrix[N2, N2] diag_delta;
      matrix[N1, N1] K;
      
      K = cov_exp_quad(x1, alpha, rho) + diag_matrix(rep_vector(1, N1))*square(sigma);
      L_K = cholesky_decompose(K);
      k_x1_x2 = cov_exp_quad(x1, x2, alpha, rho);
      L_div_k_x1_x2 = mdivide_left_tri_low(L_K, k_x1_x2);
      L_div_y1 = mdivide_left_tri_low(L_K, y1);
      
      f2_mu = L_div_k_x1_x2' * L_div_y1;
      f2_cov = cov_exp_quad(x2, alpha, rho) - L_div_k_x1_x2' * L_div_k_x1_x2;
      
      diag_delta = diag_matrix(rep_vector(delta, N2));
      f2 = multi_normal_rng(f2_mu, f2_cov + diag_delta);
      
      for (i in 1:N1){
        for (j in 1:N2){
          x_diff[i, j] = x1[i] - x2[j];
        }
      }

      k_x1_x2_p = ( 1 / square(rho)) * x_diff .* k_x1_x2;

      for (i in 1:N2){
        L_div_k_x1_x2_p = mdivide_left_tri_low(L_K, k_x1_x2_p[,i]);
        f2_prime_mu[i] = L_div_k_x1_x2_p' * L_div_y1;
        f2_prime_var[i] = sqrt( square(alpha/rho) - L_div_k_x1_x2_p' * L_div_k_x1_x2_p );
      }

      f2_stat[1] = to_row_vector(f2);
      f2_stat[2] = to_row_vector(f2_prime_mu);
      f2_stat[3] = to_row_vector(f2_prime_var);
    }
    return f2_stat;
  }
}

data{
  int<lower=1> N1;
  int<lower=1> N2;
  real x1[N1];
  vector[N1] y1;
  real x2[N2];
  real<lower=0> rho_alpha;
  real<lower=0> rho_beta;
  real<lower=0> alpha_mean;
  real<lower=0> alpha_sd;
  real<lower=0> sigma_mean;
  real<lower=0> sigma_sd;
}

transformed data{
  real delta = 1e-9;
  vector[N1] mu;
  for(n in 1:N1) mu[n] = 0;
}

parameters{
  real<lower=0> alpha;
  real<lower=0> rho;
  real<lower=0> sigma;
}

model{
  matrix[N1,N1] Sigma;
  matrix[N1,N1] L_S;
  
  Sigma = cov_exp_quad(x1, alpha, rho) + diag_matrix(rep_vector(1, N1))*square(sigma);
  
  L_S = cholesky_decompose(Sigma);
  y1 ~ multi_normal_cholesky(mu, L_S);
  
  rho ~ inv_gamma(rho_alpha, rho_beta);
  alpha ~ normal(alpha_mean, alpha_sd);
  sigma ~ normal(sigma_mean, sigma_sd);
}

generated quantities{
  matrix[3, N2] f2_stat;
  vector[N2] f2;
  vector[N2] fp2;
  vector[N2] f2_prime_mu;
  vector[N2] f2_prime_var;
  vector[N2] f2_prime;
  f2_stat = gp_pred_rng(x1, y1, x2, alpha, rho, sigma, delta);
  f2 = f2_stat[1]';
  f2_prime_mu = f2_stat[2]';
  f2_prime_var = f2_stat[3]';

  for (n2 in 1:N2){
    fp2[n2] = normal_rng(f2[n2], sigma);
    f2_prime[n2] = normal_rng(f2_prime_mu[n2], f2_prime_var[n2]);
  }
}
