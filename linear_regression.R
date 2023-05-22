# this function performes a linear regression analysis by hand
# the user specifies the response variable and the predictor variable/variables aswell as the dataset
# QBA LIU, Jakub Prorok 2023 -----------------------------------------------------------------------------------------------------------------------

# creating example data _________________________________________________________________________________________________________________________
y <- c(140, 155, 159, 179, 192, 200, 212, 215)
X1 <- c(60, 62, 67, 70, 71, 72, 75, 78)
X2 <- c(22, 25, 24, 20, 15, 14, 14, 11)
X3 <- c(12,10,14,23,12,10,16,15)
my_data <- data.frame(y, X1, X2, X3)

# actual function _______________________________________________________________________________________________________________________________
linear_model <- function(response_var, predictor_var, data){
  
  # editing input data-------------------------------------------------------------------------------------------------
  predictors <- as.matrix(data[, predictor_var])
  response <- as.matrix(data[, response_var])
  first_col <- matrix(rep(1, nrow(data)), nrow = nrow(data), ncol = 1)
  design_matrix <- cbind(first_col, predictors)
  
  # calculating regression coefficients and intercept ----------------------------------------------------------
  beta_coefs <- solve(t(design_matrix)%*%design_matrix)%*%t(design_matrix)%*%response
  if (length(predictor_var) > 1){rownames(beta_coefs)[1] <- 'intercept'}
  else if (length(predictor_var) == 1){rownames(beta_coefs) <- c('intercept', predictor_var)}
  
  # calculating intercept standard error-----------------------------------------------------------------------------
  predicted_response <- design_matrix%*%beta_coefs
  residuals <- response - predicted_response    
  S2e = t(residuals)%*%residuals / (nrow(data) - ncol(design_matrix))
  Se = sqrt(S2e)
  std_errors = sqrt(diag(c(S2e) * solve(t(design_matrix)%*%design_matrix)))
  
  # calculating R_squared-----------------------------------------------------------------------------------------
  numerator <- sum(residuals^2)
  mean_response <- mean(response)
  denominator <- sum((response - mean_response)^2)
  R2 <- 1 - numerator / denominator
  
  # testing significance of predictor variables --------------------------------------------------------------------
  t <- 1:length(beta_coefs)
  for (i in 1:length(beta_coefs)){
    t[i] <- beta_coefs[i]/std_errors[i]
  }
  p_values <- 1:length(beta_coefs)
  for (i in 1:length(beta_coefs)){
    p_values[i] <- 2*pt(abs(t[i]), nrow(data) - ncol(predictors) - 1, lower.tail = FALSE)
  }
  decision <- 1:length(beta_coefs)
  for (i in 1:length(beta_coefs)){
    if (p_values[i] <= 0.05){decision[i] <- 'significant'}
    else{decision[i] <- 'non significant'}
  }
  out <- data.frame(cbind(beta_coefs, matrix(std_errors), matrix(t), matrix(p_values), matrix(decision)))
  colnames(out) <- c('coefficients', 'standard errors', 'T-statistic', 'p-value', 'significance')
  
  # testing for normality of residuals-------------------------------------------------------------------------------    
  if (shapiro.test(residuals)$p.value < 0.05){ warning('Warning, residuals are not normally distributed.')}
  
  # testing for multicolinearity ---------------------------------------------------------------------------------
  for (j in 1:ncol(predictors)){
    while (j+1 <= ncol(predictors)) {
      p_val <- var.test(predictors[,j], predictors[,j+1])$p.val
      if (p_val <= 0.05){
        warning('Warning, at least two independent variables are colinear.')
        break
      }
      j <- j + 1
    }
  }
  #---------------------------------------------------------------------------------------------------------------
  #return(design_matrix)
  if (length(predictor_var) == 1){list(message = 'Simple linear regression.', basic_statistics = out, R_squared = R2)}
  else if (length(predictor_var) > 1 ){list(message = 'Multiple linear regression.', basic_statistics = out, R_squared = R2)}
}

# testing the function __________________________________________________________________________________________________________________________
print(linear_model(response_var = 'y', predictor_var = c('X1', 'X2', 'X3'), data = my_data))
print(linear_model(response_var = 'y', predictor_var = c('X1'), data = my_data))