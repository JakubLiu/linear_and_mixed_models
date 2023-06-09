# this is a simplified mixed model, that is simpler because of the following assumptions:

#   1.) the variance of each element in the random effect vector is set to 1
#   2.) the covariance between any two elements of the random effect vector is set to 0
#   * because of assumptions 1. and 2. the covariance matrix is an identity matrix
#   3.) the 'n-th' element of the outcome (y) vector corresponds to the 'n-th' element in the random effects vector
#   * because of assumption 3. the random effects incidence matrix is an identity matrix

# QBA LIU, Jakub Prorok 2023 -----------------------------------------------------------------------------------------------------------------


mixed_model <- function(outcome, fixed_effect, additive_gen_var, resid_var){
  
  # initializing empty fixed effect incidence matrix
  X <- matrix(rep(0, length(fixed_effect)*length(unique(fixed_effect))), nrow = length(fixed_effect), ncol = length(unique(fixed_effect)))
  
  # populating the fixed effect incidence matrix
  for (i in 1:nrow(X)){
    for (j in 1:ncol(X)){
      if (fixed_effect[i] == unique(fixed_effect)[j]){
        X[i,j] <- 1
      }
    }
  }
  
  # creating the random effect identity matrix
  Z <- diag(length(outcome))
  
  # creating covariance matrix
  A <- diag(length(outcome))
  
  # solvign mixed model equation
  upper_left <- t(X)%*%X
  lower_left <- t(Z)%*%X
  upper_right <- t(X)%*%Z
  lower_right <- t(Z)%*%Z + solve(A)*(additive_gen_var/resid_var)    
  upper <- cbind(upper_left, upper_right)
  lower <- cbind(lower_left, lower_right)
  full_matrix <- rbind(upper, lower)
  upper2 <- t(X)%*%outcome
  lower2 <- Z%*%outcome
  full_vector <- rbind(upper2, lower2)
  results_all <- solve(full_matrix)%*%full_vector
  sep <- 1:length(unique(fixed_effect))
  results_fixed_effects <- data.frame(cbind(unique(fixed_effect),results_all[sep]))
  colnames(results_fixed_effects) <- c('Fixed effect level', 'Fixed effect value')
  results_random_effects <- data.frame(cbind(outcome, results_all[-sep]))
  colnames(results_random_effects) <- c('Random effect outcome', 'Random effect value')
  
  # output
  list(message = 'Single trait mixed model.', fixed_effects = results_fixed_effects, random_effects = results_random_effects)
  
}


# testing the model
y <- rnorm(1000) # outcome vector
b <- sample(c('Holstein', 'Brown-Swiss', 'Kangayam'), length(y), replace = TRUE)  # fixed effect vector (cattle breed)
mixed_model(outcome = y, fixed_effect = b, additive_gen_var = 0.02, resid_var = 0.04)
