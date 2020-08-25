rm(list=ls())

start_time <- Sys.time()

# Estimators code for cluster sims
# Our estimators and CI
# Save data to test Koopmeiners estimators later on

# Set sim seed number
# This gets changed by the .pbs file
args <- commandArgs(trailingOnly = TRUE)
n.sim = as.numeric(args[1])

#########################
# Functions
#########################
library(matrixStats)
library(MASS)

calculate_emp_ROC_t = function(Y_d, r_d, Y_d_bar, r_d_bar, t)
{
  r_f = get_likelihood_ratio(Y_d, Y_d_bar, t)
  
  ROC_t_stg1 = calculate_ROC_t(Y_d[1:(n_d*r_d)], Y_d_bar[1:(n_d_bar*r_d_bar)], t)
  var_ROC_t_stg1 = ROC_t_stg1 * (1 - ROC_t_stg1) / (n_d * r_d) + (r_f ^ 2) * t * (1 - t) / (n_d_bar * r_d_bar)
  
  ROC_t_stg2 = calculate_ROC_t(Y_d[(n_d * r_d + 1):n_d], Y_d_bar[(n_d_bar * r_d_bar + 1):n_d_bar], t)
  var_ROC_t_stg2 = ROC_t_stg2 * (1 - ROC_t_stg2) / (n_d * (1 - r_d)) + (r_f ^ 2) * t * (1 - t) / (n_d_bar * (1 - r_d_bar))
  
  ROC_t_all = calculate_ROC_t(Y_d, Y_d_bar, t)
  var_ROC_t_all = ROC_t_all * (1 - ROC_t_all) / n_d + (r_f ^ 2) * t * (1 - t) / n_d_bar
  
  list(ROC_t_stg1 = ROC_t_stg1, var_ROC_t_stg1 = var_ROC_t_stg1, ROC_t_stg2 = ROC_t_stg2, var_ROC_t_stg2 = var_ROC_t_stg2, ROC_t_all = ROC_t_all, var_ROC_t_all = var_ROC_t_all)
}

get_likelihood_ratio = function(Y_d, Y_d_bar, t)
{
  f_d = density(Y_d)
  f_d_bar = density(Y_d_bar)
  c_bar = quantile(Y_d_bar, 1 - t)
  
  x1 = f_d$x[f_d$x <= c_bar][length(f_d$y[f_d$x <= c_bar])]
  y1 = f_d$y[f_d$x <= c_bar][length(f_d$y[f_d$x <= c_bar])]
  x2 = f_d$x[f_d$x >= c_bar][1]
  y2 = f_d$y[f_d$x >= c_bar][1]
  
  f_d_c_bar = y2 - ((y2 - y1) / (x2 - x1)) * x2 + ((y2 - y1) / (x2 - x1)) * c_bar
  
  x1 = f_d_bar$x[f_d_bar$x <= c_bar][length(f_d_bar$y[f_d_bar$x <= c_bar])]
  y1 = f_d_bar$y[f_d_bar$x <= c_bar][length(f_d_bar$y[f_d_bar$x <= c_bar])]
  x2 = f_d_bar$x[f_d_bar$x >= c_bar][1]
  y2 = f_d_bar$y[f_d_bar$x >= c_bar][1]
  
  f_d_bar_c_bar = y2 - ((y2 - y1) / (x2 - x1)) * x2 + ((y2 - y1) / (x2 - x1)) * c_bar
  
  as.numeric(f_d_c_bar / f_d_bar_c_bar)
}

calculate_ROC_t = function(X_d, X_d_bar, t)
{
  c_bar = quantile(X_d_bar, 1 - t)
  mean(X_d > c_bar)
}

RB_NP_estimate = function(Y_d, r_d, Y_d_bar, r_d_bar, t, K1)
{
  # Nonparametric Estimates
  r_f = get_likelihood_ratio(Y_d[(n_d * r_d + 1):n_d], Y_d_bar[(n_d_bar * r_d_bar + 1):n_d_bar], t)
  
  proceed_stg2_bootstrap = 0
  while(sum(proceed_stg2_bootstrap) == 0)
  {
    temp = array(Y_d, c(n_d, K1))
    temp1 = array(Y_d_bar, c(n_d_bar, K1))
    
    Y_d_bootstrap = apply(temp, 2, sample, size = n_d)
    Y_d_bar_bootstrap = apply(temp1, 2, sample, size = n_d_bar)
    
    # Get ROC(t)_stg1 for each bootstrap sample
    c_all = apply(Y_d_bar_bootstrap[1:(n_d_bar * r_d_bar),], 2, quantile, probs = 1 - t)
    temp = t(array(c_all, c(K1, n_d * r_d)))
    temp1 = array(Y_d_bootstrap[1:(n_d * r_d),] > temp, c(n_d * r_d, K1))
    ROC_t_stg1_bootstrap = colMeans(temp1)  
    var_ROC_t_stg1_bootstrap = ROC_t_stg1_bootstrap * (1 - ROC_t_stg1_bootstrap) / (n_d * r_d) + (r_f ^ 2) * t * (1 - t) / (n_d_bar * r_d_bar)
    proceed_stg2_bootstrap = as.numeric((ROC_t_stg1_bootstrap + 1.96 * sqrt(var_ROC_t_stg1_bootstrap) >= h))
    
    # Get ROC(t)_stg2 for each bootstrap sample
    c_all = apply(Y_d_bar_bootstrap[(n_d_bar * r_d_bar + 1):n_d_bar,], 2, quantile, probs = 1 - t)
    temp = t(array(c_all, c(K1, n_d * (1 - r_d))))
    temp1 = array(Y_d_bootstrap[(n_d * r_d + 1) : n_d,] > temp, c(n_d * (1 - r_d), K1))
    ROC_t_stg2_bootstrap_NP = colMeans(temp1) 
  }
  ROC_t_RB_NP = mean(ROC_t_stg2_bootstrap_NP[proceed_stg2_bootstrap == 1])
  
  list(ROC_t_RB_NP = ROC_t_RB_NP)
}

calc_CC = function(Y_d, Y_d_bar, beta_hat)
{
  n_d = nrow(Y_d)
  n_d_bar = nrow(Y_d_bar)
  n = n_d + n_d_bar
  x = rbind(cbind(rep(1, n_d), Y_d), cbind(rep(1, n_d_bar), Y_d_bar))
  D = c(rep(1, n_d), rep(0, n_d_bar))
  u_i_hat = x %*% beta_hat
  u = as.numeric(quantile(u_i_hat[D == 0], 1 - t))
  p_u = exp(u) / (1 + exp(u))
  p_i_hat = exp(u_i_hat) / (1 + exp(u_i_hat))
  p_hat_bar = mean(p_i_hat)
  omega = t(x) %*% (x * array(p_i_hat * (1 - p_i_hat), c(n, length(beta_hat))))/n
  omega_inv = solve(omega)  
  d_i = sqrt(diag(x %*% omega_inv %*% t(x)))
  
  # Copas and Corbett original correction
  S_hat_u = (p_u / (p_hat_bar * (n ^ (3 / 2)))) * sum((d_i - 1/(p_hat_bar * (1 - p_hat_bar) * d_i)) * dnorm(sqrt(n) * (u_i_hat-u) / d_i))
  
  # My adjustment to Copas and Corbett correction
  temp1 = p_i_hat * (1 - p_i_hat) * d_i / (p_u * (1 - p_u))
  temp2 = p_i_hat / (p_u * p_hat_bar * d_i)
  temp3 = (1 - p_i_hat) / (d_i * (1 - p_u) * (1 - p_hat_bar))
  S_hat_u_2 = ((p_u / (p_hat_bar * (n ^ (3 / 2)))) * sum((temp1 - temp2 - temp3) * dnorm(sqrt(n) * (u_i_hat - u) / d_i))) / (3 / 2)
  list(S_hat_u = S_hat_u, S_hat_u_2 = S_hat_u_2)
}

RB_NP_estimate_v2 = function(Y_d, r_d, Y_d_bar, r_d_bar, r_f, t, K1)
{
  proceed_stg2_bootstrap_CC_2 = 0
  while(sum(proceed_stg2_bootstrap_CC_2) < 2)
  {
    temp = array(seq(1:n_d), c(n_d, K1))
    temp1 = array(seq(1:n_d_bar), c(n_d_bar, K1))
    
    Y_d_ID_bootstrap = apply(temp, 2, sample, size = n_d)
    Y_d_bar_ID_bootstrap = apply(temp1, 2, sample, size = n_d_bar)
    
    beta_hat_stg1 = array(NA, c(K1, 4))
    beta_hat_stg2 = array(NA, c(K1, 4))
    Y_d_bootstrap_stg1 = array(NA, c(n_d, K1))
    Y_d_bar_bootstrap_stg1 = array(NA, c(n_d_bar, K1))
    Y_d_bootstrap_stg2 = array(NA, c(n_d, K1))
    Y_d_bar_bootstrap_stg2 = array(NA, c(n_d_bar, K1))
    CC_correction_2_stg1 = rep(NA, K1)
    CC_correction_2_stg2 = rep(NA, K1)
    for(i in 1:K1)
    {
      sample_case_data_bootstrap = data.frame(D = rep(1, n_d), Y1 = Y_d[Y_d_ID_bootstrap[ ,i], 1], Y2 = Y_d[Y_d_ID_bootstrap[ ,i], 2], Y3 = Y_d[Y_d_ID_bootstrap[ ,i], 3])
      sample_control_data_bootstrap = data.frame(D = rep(0, n_d_bar), Y1 = Y_d_bar[Y_d_bar_ID_bootstrap[ ,i], 1], Y2 = Y_d_bar[Y_d_bar_ID_bootstrap[ ,i], 2], Y3 = Y_d_bar[Y_d_bar_ID_bootstrap[ ,i], 3])
      sample_data_bootstrap = rbind(sample_case_data_bootstrap, sample_control_data_bootstrap)
      
      sample_data_stg1 = rbind(sample_case_data_bootstrap[1:(n_d * r_d),], sample_control_data_bootstrap[1:(n_d_bar * r_d_bar),])
      # Fit the GLM
      glm.out = glm(D ~ Y1 + Y2 + Y3, family = binomial(logit), data = sample_data_stg1)
      # Extract coefficients
      beta_hat_stg1[i,] = coef(glm.out)     
      
      temp2 = cbind(sample_control_data_bootstrap[ ,2], sample_control_data_bootstrap[ ,3], sample_control_data_bootstrap[ ,4])
      temp3 = array(beta_hat_stg1[i,2:4], c(3, 1))
      Y_d_bar_bootstrap_stg1[ ,i] = temp2 %*% temp3
      
      temp4 = cbind(sample_case_data_bootstrap[ ,2], sample_case_data_bootstrap[ ,3], sample_case_data_bootstrap[ ,4])
      Y_d_bootstrap_stg1[ ,i] = temp4 %*% temp3
      
      sample_data_stg2 = rbind(sample_case_data_bootstrap[(n_d * r_d + 1):n_d,], sample_control_data_bootstrap[(n_d_bar * r_d_bar + 1):n_d_bar,])
      # Fit the GLM
      glm.out = glm(D ~ Y1 + Y2 + Y3, family = binomial(logit), data = sample_data_stg2)
      # Extract coefficients
      beta_hat_stg2[i,] = coef(glm.out)     
      
      temp3 = array(beta_hat_stg2[i, 2:4], c(3, 1))
      Y_d_bar_bootstrap_stg2[ ,i] = temp2 %*% temp3
      Y_d_bootstrap_stg2[ ,i] = temp4 %*% temp3
      
      Y_d_stg1_CC = cbind(sample_data_stg1$Y1, sample_data_stg1$Y2, sample_data_stg1$Y3)[sample_data_stg1$D == 1,]
      Y_d_bar_stg1_CC = cbind(sample_data_stg1$Y1, sample_data_stg1$Y2, sample_data_stg1$Y3)[sample_data_stg1$D == 0,]
      results_CC_correction_stg1 = calc_CC(Y_d_stg1_CC, Y_d_bar_stg1_CC, beta_hat_stg1[i,])
      CC_correction_2_stg1[i] = results_CC_correction_stg1$S_hat_u_2
      
      Y_d_stg2_CC = cbind(sample_data_stg2$Y1, sample_data_stg2$Y2, sample_data_stg2$Y3)[sample_data_stg2$D == 1,]
      Y_d_bar_stg2_CC = cbind(sample_data_stg2$Y1, sample_data_stg2$Y2, sample_data_stg2$Y3)[sample_data_stg2$D == 0,]
      results_CC_correction_stg2 = calc_CC(Y_d_stg2_CC, Y_d_bar_stg2_CC, beta_hat_stg2[i,])
      CC_correction_2_stg2[i] = results_CC_correction_stg2$S_hat_u_2
      
    }
    
    # Get ROC(t)_stg1 for each bootstrap sample
    c_all = apply(Y_d_bar_bootstrap_stg1[1:(n_d_bar * r_d_bar),], 2, quantile, probs = 1 - t)
    temp = t(array(c_all, c(K1, n_d * r_d)))
    temp1 = array(Y_d_bootstrap_stg1[1:(n_d * r_d),] > temp, c(n_d * r_d, K1))
    ROC_t_stg1_bootstrap = colMeans(temp1)  
    var_ROC_t_stg1_bootstrap = ROC_t_stg1_bootstrap * (1 - ROC_t_stg1_bootstrap) / (n_d * r_d) + (r_f ^ 2) * t * (1 - t) / (n_d_bar * r_d_bar)
    proceed_stg2_bootstrap_CC_2 = as.numeric(((ROC_t_stg1_bootstrap - CC_correction_2_stg1) + 1.96 * sqrt(var_ROC_t_stg1_bootstrap) >= h))
    
    # Get ROC(t)_stg2 for each bootstrap sample
    c_all = apply(Y_d_bar_bootstrap_stg2[(n_d_bar * r_d_bar + 1):n_d_bar,], 2, quantile, probs = 1 - t)
    temp = t(array(c_all, c(K1, n_d * (1 - r_d))))
    temp1 = array(Y_d_bootstrap_stg2[(n_d * r_d + 1):n_d,] > temp, c(n_d  *(1 - r_d), K1))
    ROC_t_stg2_bootstrap_NP = colMeans(temp1) 
    ROC_t_stg2_bootstrap_NP_CC_2 = ROC_t_stg2_bootstrap_NP - CC_correction_2_stg2
  }
  
  ROC_t_RB_NP_CC_2 = mean(ROC_t_stg2_bootstrap_NP_CC_2[proceed_stg2_bootstrap_CC_2 == 1])
  beta_hat_CC_2 = apply(beta_hat_stg2[proceed_stg2_bootstrap_CC_2 == 1,], 2, mean)
  
  list(ROC_t_RB_NP_CC_2 = ROC_t_RB_NP_CC_2, beta_hat_CC_2 = beta_hat_CC_2)
}

###########################
# SET PARAMETERS
###########################
ROC_t = 0.55
t = 0.2

# Specify correlation matrix
rho_A_B = 0.7
rho_A_C = 0.5
rho_B_C = 0.4
sigma_cases_controls = array(c(1, rho_A_B, rho_A_C, rho_A_B, 1, rho_B_C, rho_A_C, rho_B_C, 1),
                             c(3, 3))
inv_sigma_cases_controls = solve(sigma_cases_controls)

# Specify all but one element of the mean vector
mu_B = qnorm(0.50) - qnorm(t)
mu_C = qnorm(0.47) - qnorm(t)

# Use the quadratic formula to solve for mu_A
a_A = inv_sigma_cases_controls[1, 1]
b_A = 2 * (mu_B * inv_sigma_cases_controls[1, 2] + mu_C * inv_sigma_cases_controls[1, 3])
c_A = (
  mu_B * mu_B * inv_sigma_cases_controls[2, 2] + mu_C * mu_C * inv_sigma_cases_controls[3, 3] + 2 * mu_B * mu_C * inv_sigma_cases_controls[2, 3] - (qnorm(ROC_t) - qnorm(t)) ^ 2
)
mu_A = (-b_A + sqrt(b_A ^ 2 - 4 * a_A * c_A)) / (2 * a_A)

mu = c(mu_A, mu_B, mu_C)
beta_true = inv_sigma_cases_controls %*% array(mu, c(3, 1))
# ROC(t)
pnorm(sqrt(t(mu) %*% inv_sigma_cases_controls %*% mu) + qnorm(t))

# Set parameters
n_d = 200 # Sample size for cases
n_d_bar = 200 # Sample size for controls
n = n_d + n_d_bar # Total sample size
r_d = 0.5 # Proportion of cases in stage 1
r_d_bar = 0.5 # Proportion of controls in stage 1
h = 0.7 # Continuation criterion
alpha = 0.05

K1 = 500 # Number of bootstrap samples for estimation
K2 = 1000 # Number of bootstrap samples for confidence intervals in parametric RB

# Construct objects for use later
ROC_t_stg1_CC_2_boot <- numeric()
proceed_stg2_CC_2_boot <- numeric()
ROC_t_RB_NP_CC_2_boot <- numeric()
ROC_t_all_CC_2_boot <- numeric()
case_sample_data_boot <- matrix(nrow = n_d, ncol = 3)
control_sample_data_boot <- matrix(nrow = n_d_bar, ncol = 3)

##########################
# SIMULATE and OUTPUT DATA
##########################
set.seed(n.sim)

# Draw case and control data
case_sample_data <-
  mvrnorm(n_d, mu = c(mu_A, mu_B, mu_C), Sigma = sigma_cases_controls)
control_sample_data <-
  mvrnorm(n_d_bar, mu = c(0, 0, 0), Sigma = sigma_cases_controls)
D <- numeric() # Case vs. Control
D[1:n_d] <- 1
D[(n_d + 1):(n_d + n_d_bar)] <- 0
# Combine case and control samples into a single data frame
sample_data = data.frame(rbind(case_sample_data, control_sample_data), D)

##########################
# CALCULATE ESTIMATES
##########################
# Estimate parameters of logistic regression model using stage 1 data
D <- numeric() # Case vs. Control
D[1:n_d * r_d] <- 1
D[(n_d * r_d + 1):(n_d * r_d + n_d_bar * r_d_bar)] <- 0
# Combine case and control samples into a single data frame
sample_data_stg1 = data.frame(rbind(case_sample_data[1:(n_d * r_d), ], control_sample_data[1:(n_d_bar * r_d_bar), ]), D)
# Fit the GLM
glm.out = glm(D ~ X1 + X2 + X3, family = binomial(logit), data = sample_data_stg1)
# Extract coefficients
beta_hat_stg1 = coef(glm.out)

# Calculate combined biomarker for all data using linear combination from logistic regression in stage 1
combined_marker_hat_stg1 = cbind(sample_data$X1, sample_data$X2, sample_data$X3) %*% array(beta_hat_stg1[2:4], c(3, 1))
Y_d_stg1 = combined_marker_hat_stg1[c(1:n_d)]
Y_d_bar_stg1 = combined_marker_hat_stg1[c((n_d + 1):(n_d + n_d_bar))]

# Empirical Estimates based on stage 1
ROC_results_stg1 = calculate_emp_ROC_t(Y_d_stg1, r_d, Y_d_bar_stg1, r_d_bar, t)

# Stage 1 results
ROC_t_stg1_se = sqrt(ROC_results_stg1$var_ROC_t_stg1)

# Stage 1 results with Copas and Corbett correction
Y_d_stg1_CC = cbind(sample_data_stg1$X1,
                    sample_data_stg1$X2,
                    sample_data_stg1$X3)[c(1:n_d * r_d), ]
Y_d_bar_stg1_CC = cbind(sample_data_stg1$X1,
                        sample_data_stg1$X2,
                        sample_data_stg1$X3)[c((n_d * r_d + 1):(n_d * r_d + n_d_bar * r_d_bar)), ]

results_CC_correction_stg1 = calc_CC(Y_d_stg1_CC, Y_d_bar_stg1_CC, beta_hat_stg1)

r_f = get_likelihood_ratio(Y_d_stg1_CC, Y_d_bar_stg1_CC, t)

ROC_t_stg1_CC_2 = ROC_results_stg1$ROC_t_stg1 - results_CC_correction_stg1$S_hat_u_2
# Continuation criterion
proceed_stg2_CC_2 = as.numeric((ROC_t_stg1_CC_2 + 1.96 * ROC_t_stg1_se) >= h)

ROC_t_stg2_stg2_CC_2 = NA
Pepe_sd_stg2_CC_2 = NA
ROC_t_RB_NP_CC_2 = NA
Pepe_sd_RB_NP_CC_2 = NA
boot_se_RB_NP = NA

if (proceed_stg2_CC_2 == 1) {
  
  ##################################################
  # Estimates based on all data logisitic regression
  ##################################################
  # Estimate parameters of logistic regression model using all the data
  D <- numeric() # Case vs. Control
  D[1:n_d] <- 1
  D[(n_d + 1):(n_d + n_d_bar)] <- 0
  # Combine case and control samples into a single data frame
  sample_data = data.frame(rbind(case_sample_data, control_sample_data), D)
  # Fit the GLM
  glm.out = glm(D ~ X1 + X2 + X3, family = binomial(logit), data = sample_data)
  # Extract coefficients
  beta_hat_all = coef(glm.out)
  
  # Calculate combined biomarker for all data using linear combination from logistic regression
  combined_marker_hat_all = cbind(sample_data$X1, sample_data$X2, sample_data$X3) %*% array(beta_hat_all[2:4], c(3, 1))
  Y_d_all = combined_marker_hat_all[c(1:n_d)]
  Y_d_bar_all = combined_marker_hat_all[c((n_d + 1):(n_d + n_d_bar))]
  
  # Empirical Estimates based on stage 1
  ROC_results_all = calculate_emp_ROC_t(Y_d_all, r_d, Y_d_bar_all, r_d_bar, t)
  
  # Biased results that use all the data with Copas and Corbett correction
  Y_d_all_CC = cbind(sample_data$X1, sample_data$X2, sample_data$X3)[c(1:n_d), ]
  Y_d_bar_all_CC = cbind(sample_data$X1, sample_data$X2, sample_data$X3)[c((n_d + 1):(n_d + n_d_bar)), ]
  
  results_CC_correction_all = calc_CC(Y_d_all_CC, Y_d_bar_all_CC, beta_hat_all)
  
  r_f = get_likelihood_ratio(Y_d_all_CC, Y_d_bar_all_CC, t)
  
  ##################################################
  # Estimates based on stage 2 logisitic regression
  ##################################################
  # Estimate parameters of logistic regression model using stage 2 data
  D <- numeric() # Case vs. Control
  D[1:n_d * (1 - r_d)] <- 1
  D[(n_d * (1 - r_d) + 1):(n_d * (1 - r_d) + n_d_bar * (1 - r_d_bar))] <- 0
  sample_data_stg2 = data.frame(rbind(case_sample_data[(n_d * r_d + 1):n_d, ], control_sample_data[(n_d_bar * r_d_bar + 1):n_d_bar, ]), D)
  # Fit the GLM
  glm.out = glm(D ~ X1 + X2 + X3, family = binomial(logit), data = sample_data_stg2)
  # Extract coefficients
  beta_hat_stg2 = coef(glm.out)
  
  # Calculate combined biomarker for all data using linear combination from logistic regression in stage 2
  combined_marker_hat_stg2 = cbind(sample_data$X1, sample_data$X2, sample_data$X3) %*% array(beta_hat_stg2[2:4], c(3, 1))
  Y_d_stg2 = combined_marker_hat_stg2[sample_data$D == 1]
  Y_d_bar_stg2 = combined_marker_hat_stg2[sample_data$D == 0]
  
  # Empirical Estimates based on stage 2 logistic regression
  ROC_results_stg2 = calculate_emp_ROC_t(Y_d_stg2, r_d, Y_d_bar_stg2, r_d_bar, t)
  
  # Inefficient Stage 2 results results with Copas and Corbett correction
  Y_d_stg2_CC = cbind(sample_data_stg2$X1,
                      sample_data_stg2$X2,
                      sample_data_stg2$X3)[sample_data_stg2$D == 1, ]
  Y_d_bar_stg2_CC = cbind(sample_data_stg2$X1,
                          sample_data_stg2$X2,
                          sample_data_stg2$X3)[sample_data_stg2$D == 0, ]
  
  results_CC_correction_stg2 = calc_CC(Y_d_stg2_CC, Y_d_bar_stg2_CC, beta_hat_stg2)
  
  ROC_t_stg2_stg2_CC_2 = ROC_results_stg2$ROC_t_stg2 - results_CC_correction_stg2$S_hat_u_2
  
  r_f = get_likelihood_ratio(Y_d_stg2_CC, Y_d_bar_stg2_CC, t)
  Pepe_var_stg2_CC_2 = ROC_t_stg2_stg2_CC_2 * (1 - ROC_t_stg2_stg2_CC_2) / nrow(Y_d_stg2_CC) + (r_f ^ 2) * t * (1 - t) / nrow(Y_d_bar_stg2_CC)
  Pepe_sd_stg2_CC_2 = sqrt(Pepe_var_stg2_CC_2)
  
  ######################################
  # Modified Rao-Blackwell Algorithm
  ######################################
  Y_d = cbind(sample_data$X1, sample_data$X2, sample_data$X3)[sample_data$D == 1, ]
  Y_d_bar = cbind(sample_data$X1, sample_data$X2, sample_data$X3)[sample_data$D == 0, ]
  
  # Nonparametric Rao-Blackwell Estimator
  r_f = get_likelihood_ratio(Y_d_stg2[(n_d * r_d + 1):n_d], Y_d_bar_stg2[(n_d_bar * r_d_bar + 1):n_d_bar], t)
  ROC_RB_NP_results = RB_NP_estimate_v2(Y_d, r_d, Y_d_bar, r_d_bar, r_f, t, K1)
  
  ROC_t_RB_NP_CC_2 = ROC_RB_NP_results$ROC_t_RB_NP_CC_2
  
  Pepe_var_RB_NP_CC_2 = ROC_t_RB_NP_CC_2 * (1 - ROC_t_RB_NP_CC_2) / nrow(Y_d) + (r_f ^ 2) * t * (1 - t) / nrow(Y_d_bar)
  Pepe_sd_RB_NP_CC_2 = sqrt(Pepe_var_RB_NP_CC_2)
  
  for (j in c(1:n)) {
    print(j)
    
    case_boot_ix <- sample(c(1:n_d), size = n_d, replace = TRUE)
    control_boot_ix <-
      sample(c(1:n_d_bar), size = n_d_bar, replace = TRUE)
    
    for (k in c(1:n_d)) {
      case_sample_data_boot[k, ] <- case_sample_data[case_boot_ix[k], ]
    }
    
    for (k in c(1:n_d_bar)) {
      control_sample_data_boot[k, ] <-
        control_sample_data[control_boot_ix[k], ]
    }
    
    D <- numeric() # Case vs. Control
    D[1:n_d] <- 1
    D[(n_d + 1):n] <- 0
    sample_data_boot = data.frame(rbind(case_sample_data_boot, control_sample_data_boot), D)
    
    # Stage 1
    D <- numeric() # Case vs. Control
    D[1:n_d * r_d] <- 1
    D[(n_d * r_d + 1):(n_d * r_d + n_d_bar * r_d_bar)] <- 0
    sample_data_stg1_boot = data.frame(rbind(case_sample_data_boot[1:(n_d * r_d), ], control_sample_data_boot[1:(n_d_bar * r_d_bar), ]), D)
    # Fit the GLM
    glm.out = glm(D ~ X1 + X2 + X3, family = binomial(logit), data = sample_data_stg1_boot)
    # Extract coefficients
    beta_hat_stg1_boot = coef(glm.out)
    
    # Calculate combined biomarker for all data using linear combination from logistic regression in stage 1
    combined_marker_hat_stg1_boot = cbind(sample_data_boot$X1,
                                          sample_data_boot$X2,
                                          sample_data_boot$X3) %*% array(beta_hat_stg1_boot[2:4], c(3, 1))
    Y_d_stg1_boot = combined_marker_hat_stg1_boot[c(1:n_d)]
    Y_d_bar_stg1_boot = combined_marker_hat_stg1_boot[c((n_d + 1):n)]
    
    # Empirical Estimates based on stage 1
    ROC_results_stg1_boot = calculate_emp_ROC_t(Y_d_stg1_boot, r_d, Y_d_bar_stg1_boot, r_d_bar, t)
    
    # Stage 1 results
    Y_d_stg1_boot_CC = cbind(sample_data_stg1_boot$X1,
                             sample_data_stg1_boot$X2,
                             sample_data_stg1_boot$X3)[c(1:n_d * r_d), ]
    Y_d_bar_stg1_boot_CC = cbind(sample_data_stg1_boot$X1,
                                 sample_data_stg1_boot$X2,
                                 sample_data_stg1_boot$X3)[c((n_d * r_d + 1):(n_d * r_d + n_d_bar * r_d_bar)), ]
    
    results_CC_correction_stg1_boot = calc_CC(Y_d_stg1_boot_CC, Y_d_bar_stg1_boot_CC, beta_hat_stg1_boot)
    
    r_f = get_likelihood_ratio(Y_d_stg1_boot_CC, Y_d_bar_stg1_boot_CC, t)
    
    ROC_t_stg1_CC_2_boot[j] = ROC_results_stg1_boot$ROC_t_stg1 - results_CC_correction_stg1_boot$S_hat_u_2
    ROC_t_stg1_boot_se = sqrt(ROC_results_stg1_boot$var_ROC_t_stg1)
    proceed_stg2_CC_2_boot[j] = as.numeric((ROC_t_stg1_CC_2_boot[j] + 1.96 * ROC_t_stg1_boot_se) >= h)

    if (proceed_stg2_CC_2_boot[j] == 1) {
      ##################################################
      # Estimates based on all data logisitic regression
      ##################################################
      # Estimate parameters of logistic regression model using all the data
      D <- numeric() # Case vs. Control
      D[1:n_d] <- 1
      D[(n_d + 1):n] <- 0
      sample_data_boot = data.frame(rbind(case_sample_data_boot[, c(1:3)], control_sample_data_boot[, c(1:3)]),
                                    D)
      # Fit the GLM
      glm.out = glm(D ~ X1 + X2 + X3, family = binomial(logit), data = sample_data_boot)
      # Extract coefficients
      beta_hat_all_boot = coef(glm.out)
      
      # Calculate combined biomarker for all data using linear combination from logistic regression
      combined_marker_hat_all_boot = cbind(sample_data_boot$X1,
                                           sample_data_boot$X2,
                                           sample_data_boot$X3) %*% array(beta_hat_all_boot[2:4], c(3, 1))
      Y_d_all_boot = combined_marker_hat_all_boot[c(1:n_d)]
      Y_d_bar_all_boot = combined_marker_hat_all_boot[c((n_d + 1):(n_d + n_d_bar))]
      
      # Empirical Estimates based on all data
      ROC_results_all_boot = calculate_emp_ROC_t(Y_d_all_boot, r_d, Y_d_bar_all_boot, r_d_bar, t)
      
      # Biased results that use all the data with Copas and Corbett correction
      Y_d_all_boot_CC = cbind(sample_data_boot$X1,
                              sample_data_boot$X2,
                              sample_data_boot$X3)[c(1:n_d), ]
      Y_d_bar_all_boot_CC = cbind(sample_data_boot$X1,
                                  sample_data_boot$X2,
                                  sample_data_boot$X3)[c((n_d + 1):n), ]
      ROC_t_all_boot_se = sqrt(ROC_results_all_boot$var_ROC_t_all)
      
      results_CC_correction_all_boot = calc_CC(Y_d_all_boot_CC, Y_d_bar_all_boot_CC, beta_hat_all_boot)
      
      r_f = get_likelihood_ratio(Y_d_all_boot_CC, Y_d_bar_all_boot_CC, t)
      
      ROC_t_all_CC_2_boot[j] = ROC_results_all_boot$ROC_t_all - results_CC_correction_all_boot$S_hat_u_2
      
      ##################################################
      # Estimates based on stage 2 logisitic regression
      ##################################################
      # Estimate parameters of logistic regression model using stage 2 data
      D <- numeric() # Case vs. Control
      D[1:n_d * (1 - r_d)] <- 1
      D[(n_d * (1 - r_d) + 1):(n_d * (1 - r_d) + n_d_bar * (1 - r_d_bar))] <- 0
      sample_data_stg2_boot = data.frame(rbind(case_sample_data_boot[(n_d * r_d + 1):n_d, ], control_sample_data_boot[(n_d_bar * r_d_bar + 1):n_d_bar, ]),
                                         D)
      # Fit the GLM
      glm.out = glm(D ~ X1 + X2 + X3, family = binomial(logit), data = sample_data_stg2_boot)
      # Extract coefficients
      beta_hat_stg2_boot = coef(glm.out)
      
      # Calculate combined biomarker for all data using linear combination from logistic regression in stage 2
      combined_marker_hat_stg2_boot = cbind(sample_data_boot$X1,
                                            sample_data_boot$X2,
                                            sample_data_boot$X3) %*% array(beta_hat_stg2_boot[2:4], c(3, 1))
      Y_d_stg2_boot = combined_marker_hat_stg2_boot[sample_data_boot$D == 1]
      Y_d_bar_stg2_boot = combined_marker_hat_stg2_boot[sample_data_boot$D == 0]
      
      ######################################
      # Modified Rao-Blackwell Algorithm
      ######################################
      Y_d_RB_NP_boot = cbind(sample_data_boot$X1,
                             sample_data_boot$X2,
                             sample_data_boot$X3)[sample_data_boot$D == 1, ]
      Y_d_bar_RB_NP_boot = cbind(sample_data_boot$X1,
                                 sample_data_boot$X2,
                                 sample_data_boot$X3)[sample_data_boot$D == 0, ]
      
      # Nonparametric Rao-Blackwell Estimator
      r_f = get_likelihood_ratio(Y_d_stg2_boot[(n_d * r_d + 1):n_d], Y_d_bar_stg2_boot[(n_d_bar * r_d_bar + 1):n_d_bar], t)
      ROC_RB_NP_results = RB_NP_estimate_v2(Y_d_RB_NP_boot, r_d, Y_d_bar_RB_NP_boot, r_d_bar, r_f, t, K1)
      
      ROC_t_RB_NP_CC_2_boot[j] = ROC_RB_NP_results$ROC_t_RB_NP_CC_2
      
    }
  }
}

mean_ROC_t_RB_NP_boot_CC_2 <-
  mean(ROC_t_RB_NP_CC_2_boot, na.rm = TRUE)
ss_RB_NP_boot = sum((ROC_t_RB_NP_CC_2_boot - mean_ROC_t_RB_NP_boot_CC_2) ^ 2,
                    na.rm = TRUE)

# Calculate bootstrap standard error
n_boot = sum(proceed_stg2_CC_2_boot)
frac = 1 / (n_boot - 1)
boot_se_RB_NP = sqrt(frac * ss_RB_NP_boot)

######### CONFIDENCE INTERVALS #########

# Wald-Type
lower_Wald <-
  mean_ROC_t_RB_NP_boot_CC_2 - qnorm(.975) * boot_se_RB_NP
upper_Wald <-
  mean_ROC_t_RB_NP_boot_CC_2 + qnorm(.975) * boot_se_RB_NP
coverage_Wald <-
  as.numeric(lower_Wald <= ROC_t && upper_Wald >= ROC_t)

# Percentile Bootstrap

lower_percentile <-
  quantile(ROC_t_RB_NP_CC_2_boot, 0.025, na.rm = TRUE)
upper_percentile <-
  quantile(ROC_t_RB_NP_CC_2_boot, 0.975, na.rm = TRUE)
coverage_percentile <-
  as.numeric(lower_percentile <= ROC_t && upper_percentile >= ROC_t)

# BC

a_hat <- 0
z_0_hat <-
  qnorm(sum(as.numeric(ROC_t_RB_NP_CC_2_boot < ROC_t_RB_NP_CC_2), na.rm = TRUE) / sum(proceed_stg2_CC_2_boot))
alpha <- 0.025
alpha_1 <-
  pnorm(z_0_hat + (z_0_hat + qnorm(alpha)) / (1 - a_hat * (z_0_hat + qnorm(alpha))))
alpha_2 <-
  pnorm(z_0_hat + (z_0_hat + qnorm(1 - alpha)) / (1 - a_hat * (z_0_hat + qnorm(1 - alpha))))
lower_BC <- quantile(ROC_t_RB_NP_CC_2_boot, alpha_1, na.rm = TRUE)
upper_BC <- quantile(ROC_t_RB_NP_CC_2_boot, alpha_2, na.rm = TRUE)
coverage_BC <- as.numeric(lower_BC <= ROC_t && upper_BC >= ROC_t)

######################
## Hypothesis Tests ##
######################

# Our null hypothesis is H_0: ROC(0.2) <= rho_0
rho_0 <- 0.60
# We reject H_0 if our lower bound is greater than rho_0
# We create a dummy variable that will be 1 if the null rejected and 0 otherwise
reject_H_0_Wald <- as.numeric(lower_Wald > rho_0)
reject_H_0_percentile <- as.numeric(lower_percentile > rho_0)
reject_H_0_BC <- as.numeric(lower_BC > rho_0)

end_time <- Sys.time()

temp = diff(c(start_time, end_time))
units(temp) = "mins"
runtime <- as.numeric(temp)

################
# OUTPUT RESULTS
################
output_for_csv = cbind(n.sim, proceed_stg2_CC_2, ROC_t_RB_NP_CC_2, boot_se_RB_NP,
                       lower_Wald, upper_Wald, coverage_Wald, reject_Wald = reject_H_0_Wald,
                       lower_percentile, upper_percentile, coverage_percentile, reject_percentile = reject_H_0_percentile,
                       lower_BC, upper_BC, coverage_BC, reject_BC = reject_H_0_BC,
                       runtime)
write.csv(output_for_csv, file = paste("results_boot_CI_HT_55", "_", n.sim, ".csv", sep = ""))