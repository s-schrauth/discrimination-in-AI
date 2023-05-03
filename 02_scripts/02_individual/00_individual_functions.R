# Functions for 01_individual_model.Rmd

## Functions for Individual Analysis
### Functions for pair-wise distance calculations and comparison of instance & probability distances
metric_fun = function(x_high, x_low, y_high, y_low){
  quot1 = log(max(x_low/y_low, y_low/x_low))
  quot2 = log(max(x_high/y_high, y_high/x_high))
  
  output = max(quot1, quot2)
  
  return(output)
}

diff_dist = function(prob, data){
  
  # D(Mx,My) (probability distance)
  ## define vectors for probabilities of getting a high or low COMPAS score
  prob_high = as.vector(prob$high)
  prob_low = as.vector(prob$low)
  
  ## empty matrix for pair-wise probability distance
  Ddist = matrix(data = NA,
                 nrow = nrow(prob),
                 ncol = nrow(prob))
  
  ## loop over lower triangle matrix of Ddist and calculate pair-wise distances
  for (i in 1:nrow(prob)) {
    for (j in 1:i) {
      x_high = prob_high[i]
      x_low = prob_low[i]
      
      y_high = prob_high[j]
      y_low = prob_low[j]
      
      Ddist[i,j] = metric_fun(x_high, x_low, y_high, y_low)
    }
  }
  
  ## normalise distances to a range of 0 to 1 
  Ddist = rescale(Ddist, to = c(0,1))
  
  
  # d(x,y) (instance distance)
  ## calculate pair-wise gower distance
  ddist = data %>% 
    select(!compas) %>% 
    daisy(metric = "gower") %>% 
    as.matrix()
  
  ## create lower triangle matrix for compatibility
  ddist[upper.tri(ddist)] = NA
  
  colnames(ddist) = NULL
  rownames(ddist) = NULL
  
  
  # Difference Ddist & ddist
  ## count the cases where individual fairness is not given
  diff_dist = ifelse(Ddist > ddist, 1, 0)
  
  ## calculate relative proportion of those cases w.r.t all pair-wise distances
  n_unfair = sum(diff_dist, na.rm = TRUE)
  n_complete = (length(diff_dist)/2)-nrow(diff_dist)
  
  n_unfair_rel = n_unfair/n_complete
  n_unfair_rel
  
  # Output
  output = list(diff_dist = diff_dist,
                n_unfair_rel = n_unfair_rel)
  
  return(output)
}

### Function for basic logistic regression by hand 
### (with manually coded optimization algorithm analog to glm)
my_method <- function(x, y, 
                      family = binomial(link = "logit"), 
                      control, ...) {
  # Initialize the beta coefficients to zero
  beta = rep(0, ncol(x))
  
  # Set the maximum number of iterations to 100
  max_iter = 100
  
  # Set the convergence threshold to 1e-6
  tol = 1e-6
  
  # Iterate until convergence or maximum number of iterations reached
  for (iter in 1:max_iter) {
    # Calculate the predicted probabilities using the current beta coefficients
    eta = x %*% beta
    mu = family$linkinv(eta)
    
    # Calculate the working response
    w = mu * (1 - mu)
    z = eta + (y - mu) / w
    
    # Compute the gradient of the loss function
    grad = t(x) %*% (w * z)
    
    beta_new = beta + (solve(t(x) %*% diag(rep(w, nrow = nrow(x))) %*% x) %*% grad)
    
    # Check for convergence
    if (max(abs(beta_new - beta)) < tol) {
      break
    }
    
    if(iter%%10 == 0) print(paste("Iteration:", iter))
    
    # Update beta and continue iterating
    beta = beta_new
  }
  
  # calculate predicted values, covariance matrix, standard error of coefficients, z-values and p-values
  p_hat = plogis(x %*% beta)
  V = solve(t(x) %*% diag(c(p_hat * (1 - p_hat))) %*% x)
  se = sqrt(diag(V))
  z = beta / se
  p = 2 * (1 - pnorm(abs(z)))
  
  # Return the coefficients, standard errors, z-values, and p-values as a list
  return(list(coefficients = beta, se = se, z = z, p = p))
}

my_predict = function(summary, data){
  coeffs = as.vector(summary$estimate)
  data_new = data %>%
    select(!c(compas)) %>%
    mutate(intercept = 1) %>%
    relocate(intercept, .before = 0) %>%
    as.matrix()
  
  prob_high = plogis(data_new %*% coeffs)
  
  prob_low = 1 - prob_high
  
  df = data.frame(low = prob_low,
                  high = prob_high)
  
  return(df)
}

## Functions for In-Processing
### New loss function
my_penalty_function <- function(mu, x, metric) {
  
  prob = data.frame(low = mu,
                    high = 1-mu)
  
  data = x[,2:ncol(x)] %>% 
    as.data.frame() %>% 
    mutate(compas = NA) %>% 
    mutate(sex = factor(ifelse(sex_Female == 1, "Female", "Male"))) %>% 
    mutate(race = factor(case_when(
      race_African.American == 1 ~ "African American",
      race_Hispanic == 1 ~ "Hispanic",
      race_Other == 1 ~ "Other",
      .default = "Caucasian"
    )
    )
    ) %>% 
    mutate(charge_degree = factor(ifelse(charge_degree_M == 1, "M", "F"))) %>% 
    mutate(two_year_recid = factor(ifelse(two_year_recid_yes == 1, "yes", "no"))) %>% 
    select(sex,
           race,
           age,
           days_in_jail,
           juv_fel_count,
           juv_misd_count,
           juv_other_count,
           priors_count,
           charge_degree,
           two_year_recid,
           compas)
  
  dist = diff_dist(prob = prob, data = data)
  
  penalty = dist$n_unfair_rel
  
  return(penalty)
}

### Function for logistic regression with penalty by hand 
my_method_with_penalty <- function(x, y, 
                                   family = binomial(link = "logit"), 
                                   control, ...) {
  # Initialize the beta coefficients to zero
  beta = rep(0, ncol(x))
  
  # Set the maximum number of iterations to 100
  max_iter = 2 #####!!!!!!!!!!
  
  # Set the convergence threshold to 1e-6
  tol = 1e-6
  
  # Set the penalty weight to 0.5
  lambda = control$lambda
  
  metric = control$metric
  
  # Set the penalty function
  penalty_function = control$penalty_function
  
  # Iterate until convergence or maximum number of iterations reached
  for (iter in 1:max_iter) {
    # Calculate the predicted probabilities using the current beta coefficients
    eta = x %*% beta
    mu = family$linkinv(eta)
    
    # Calculate the working response
    w = mu * (1 - mu)
    z = eta + (y - mu) / w
    
    # Compute the gradient of the loss function
    grad = t(x) %*% (w * z)
    
    # Compute the penalty term gradient
    pen_grad = penalty_function(mu = mu, x = x)
    
    # Compute the total gradient and update the beta coefficients
    total_grad = grad + lambda * pen_grad
    
    beta_new = beta + (solve(t(x) %*% diag(rep(w, nrow = nrow(x))) %*% x) %*% total_grad)
    
    # Check for convergence
    if (max(abs(beta_new - beta)) < tol) {
      break
    }
    
    if(iter%%10 == 0) print(paste("Iteration:", iter)) 
    
    # Update beta and continue iterating
    beta = beta_new
  }
  
  # calculate predicted values, covariance matrix, standard error of coefficients, z-values and p-values
  p_hat = plogis(x %*% beta)
  V = solve(t(x) %*% diag(c(p_hat * (1 - p_hat))) %*% x)
  se = sqrt(diag(V))
  z = beta / se
  p = 2 * (1 - pnorm(abs(z)))
  
  # Return the coefficients, standard errors, z-values, and p-values as a list
  return(list(coefficients = beta, se = se, z = z, p = p))
}

