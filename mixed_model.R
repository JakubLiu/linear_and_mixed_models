# this is a simplified mixed model that is simpler, because of the following assumptions:

#   1.) the variance of each element in the random effect vector is set to 1
#   2.) the pedigree file is simplified so that:
#       2.1.) parents do not have children with their children
#   * because of assumptions 1. and 2. the covariance matrix has 1 on it's main diagonal
#   3.) the 'n-th' element of the outcome (y) vector corresponds to the 'n-th' element in the random effects vector
#   * because of assumption 3. the random effects incidence matrix is an identity matrix
# QBA LIU, Jakub Prorok 2023 ---------------------------------------------------------------------------------------------------------------


# 1. creating input data ---------------------------------------------------------------------------------------------------------------
# creating simple pedigree file
animal_ID <- seq(1,100)                                             # 100 animal IDs
mother_ID <- sample(seq(1,10), length(animal_ID), replace = TRUE)   # 10 different mothers for the 100 animals (coded as 1-10)
father_ID <- sample(seq(11,20), length(animal_ID), replace = TRUE)  # 10 different fathers for the 100 animals (coded as 11-20)
ped_file <- data.frame(animal_ID, mother_ID, father_ID)

# creating response variable vector
y <- rnorm(nrow(ped_file))

# creating fixed effect vector
b <- sample(c('Holstein', 'Brown-Swiss', 'Bison'), nrow(ped_file), replace = TRUE)  # in this example the fixed effect is breed



# 2. the function -----------------------------------------------------------------------------------------------------------------------
mixed_model <- function(outcome, fixed_effect, pedigree, additive_gen_var, resid_var){
  
  library(MASS)
  
  # initializing empty fixed effect incidence matrixa
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
  A <- diag(nrow = nrow(pedigree), ncol = nrow(pedigree))
  for (i in 1:nrow(pedigree)){
    for (j in i:nrow(pedigree)){
      if (pedigree[i,2] == pedigree[j,2]){
        if (pedigree[i,3] == pedigree[j,3]){
          A[i,j] <- 1
          A[j,i] <- 1
        }
        else {
          A[i,j] <- 0.5
          A[j,i] <- 0.5
        }
      }
      if (pedigree[i,3] == pedigree[j,3]){
        if (pedigree[i,2] == pedigree[j,2]){
          A[i,j] <- 1
          A[j,i] <- 1
        }
        else{
          A[i,j] <- 0.5
          A[j,i] <- 0.5
        }
      }          
    }
  }
  
  
  # solvign mixed model equation
  upper_left <- t(X)%*%X
  lower_left <- t(Z)%*%X
  upper_right <- t(X)%*%Z
  lower_right <- t(Z)%*%Z + ginv(A)*(additive_gen_var/resid_var)    
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
  list(message = 'Single trait mixed model.', fixed_effects = results_fixed_effects, random_effects = results_random_effects, covariance_matrix = A)
  
}

# 3. testing the function on the input data --------------------------------------------------------------------------------------------
mixed_model(outcome = y, fixed_effect = b, pedigree = ped_file, additive_gen_var = 0.17, resid_var = 0.24)
