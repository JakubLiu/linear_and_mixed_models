# a function that calls upon custom written T tests and one way ANOVA
# the user can define which test needs to be performed and specify other metrics like the alternative hypothesis or the significance level
# there are 5 custom functions
#   - one sample T test
#   - two sample unpaired T test
#   - two sample paired T test
#   - one way ANOVA
#   - a wrapper function
# QBA LIU, Jakub Prorok 2023 ---------------------------------------------------------------------------------------------------------------------------------

# one sample t_test
one_sample_t_test <- function(sample, mu0, alternative){
  
  
  if (shapiro.test(sample)$p.val <= 0.05){
    warning('Caution, the provided sample doesnt meet the normality requirement.')
  }
  #Check if input sample are numeric vector with at least 2 values
  if (!is.numeric(sample) || length(sample) < 2) {
    stop("Sample should be a numeric vector with at least two values.")
  }
  
  #Check if given mu0 is a numeric value
  if (!is.numeric(mu0)) {
    stop("The hypothesized population mean (mu0) should be a numeric value.")
  }
  
  # Calculate sample statistics  
  mean_sample <- mean(sample)
  std_dev_sample <- sd(sample)
  n_sample <- length(sample)
  df <- n_sample - 1
  T <- sqrt(n_sample)*(as.numeric(mean_sample) - as.numeric(mu0))/std_dev_sample
  
  # Calculate p-value based on alternative hypothesis
  if (alternative == 1){     # two sided
    p_value <- pt(T, df)
  }
  else if (alternative == 2){  # greater
    p_value <- pt(T, df, lower.tail = FALSE)
  }
  else if (alternative == 3) {   # less
    p_value <- 2*pt(abs(T), df, lower.tail = FALSE)
  }
  
  # Return results as a list
  list(message = 'Results of one sample T-test:', test_statistic = T, p_value = round(p_value, 4))
}

# two sample unpaired t_test
two_sample_unpaired_T_test <- function(sample1 = 0, sample2 = 0, alternative){
  
  #Check if input samples are numeric vectors
  if (!is.numeric(sample1) || !is.numeric(sample2)){
    stop("Input samples must be numeric vectors.")
  }
  
  if (shapiro.test(sample1)$p.val <= 0.05 || shapiro.test(sample1)$p.val <= 0.05){
    warning('Caution, one of the provided samples does not meet the normality requirement.')
  }
  
  if (var.test(sample1, sample2)$p.val <= 0.05){
    warning('Caution, both samples do not have equal variances.')
  }
  
  if (length(sample1) < 2 || length(sample2) < 2){
    stop('Sample inputs must have at least two values')
  }
  
  # Calculate sample statistics 
  mean_sample1 <- mean(sample1)
  mean_sample2 <- mean(sample2)
  var_sample1 <- var(sample1)
  var_sample2 <- var(sample2)
  n_sample1 <- length(sample1)
  n_sample2 <- length(sample2)
  df <- n_sample1 + n_sample2 - 2
  T = (mean_sample1-mean_sample2)/sqrt((n_sample1-1)*var_sample1+(n_sample2-1)*var_sample2)*sqrt((n_sample1*n_sample2*df) / (n_sample1+n_sample2))
  
  # Calculate p-value based on alternative hypothesis
  if (alternative == 1){    # two sided
    p_value <- pt(T, df)
  }
  else if (alternative == 2){   # greater
    p_value <- pt(T, df, lower.tail = FALSE)
  }
  else if (alternative == 3) {    # less
    p_value <- 2*pt(abs(T), df, lower.tail = FALSE)
  }
  
  list(message = 'Results of two sample unpaired t test:', test_statistic = T, p_value = round(p_value, 4))
}

# two sample paired t_test
two_sample_paired_T_test <- function(sample1, sample2, alternative){
  
  if (!is.vector(sample1) || !is.vector(sample2)){
    stop("Sample inputs must be vectors")
  }
  if (length(sample1) != length(sample2)){
    stop("Sample inputs must have the same length")
  }
  if (!is.numeric(sample1) || !is.numeric(sample2)){
    stop("Sample inputs must be numeric")
  }
  if (length(sample1) < 2 || length(sample2) < 2){
    stop("Sample inputs must have at least two values")
  }
  
  if (shapiro.test(sample1)$p.val <= 0.05 || shapiro.test(sample2)$p.val <= 0.05){
    warning('Caution, one of the samples does not meet the normality requirement.')
  }
  
  if (var.test(sample1, sample2)$p.val <= 0.05){
    warning('Caution, the provided samples do not have equal variances.')
  }
  # Calculate sample statistics
  D <- sample1 - sample2
  mean_D <- mean(D)
  std_dev_D <- sd(D)
  n_D <- length(D)
  df <- n_D - 1
  T <- sqrt(n_D)*mean_D/std_dev_D
  
  # Calculate p-value based on alternative hypothesis
  if (alternative == 1){    # two sided
    p_value <- pt(T, df)
  }
  else if (alternative == 2){    # greater
    p_value <- pt(T, df, lower.tail = FALSE)
  }
  else if (alternative == 3) {   # less
    p_value <- 2*pt(abs(T), df, lower.tail = FALSE)
  }
  
  list(message = 'Results of two sample paired t test:', test_statistic = T, p_value = round(p_value, 4))
}

# one way ANOVA
one_way_anova <- function(data, sign_level){
  
  df <- data
  
  # testing normality of dependent variable levels
  
  for (j in 1:ncol(df)){
    if (shapiro.test(df[,j])$p.val <= 0.05){
      warning('Caution, at least one of the independent variable levels does not meet the normality requirement.')
    }
  }
  
  # testing for homogenity of variance
  variances <- 1:ncol(df)
  
  for (j in 1:ncol(df)){
    variances[j] <- var(df[,j])
  }
  
  variances_sorted <- sort(variances, decreasing = FALSE)
  F_stat <- variances_sorted[1] / tail(variances_sorted, n = 1)
  df_numerator <- nrow(df) - 1
  df_denominator <- nrow(df) - 1
  p_value <- pf(F_stat, df_numerator, df_denominator, lower.tail = FALSE)
  
  if (p_value <= 0.05){
    warning('Caution, at least one of the independent variable levels does not meet the normality requirement.')
  }
  
  total <- c()  # initialize vector for overall mean
  group_means <- c()  # initialize vector that holds the group means (this vector will be used when calculating SSR)
  
  for (j in 1:ncol(df)){
    
    holder_variable <- paste0("group_", j,"_mean")  # dynamically create variable names for each group
    assign(holder_variable, mean(df[,j])) # assign the group mean to the created variable
    group_means <- append(group_means, mean(df[,j]))  # populate the group means vector
    total <- append(total, df[,j])  # populate the overall mean vector
  }
  
  total_mean <- mean(total) # calculate the overall mean
  
  SSE <- 0
  
  for (j in 1:ncol(df)){            # iterating by df columns
    for (i in 1:length(df[,j])){    # iterating by rows of the given df column
      value <- (df[i,j] - mean(df[,j]))^2  # calculating every individual square or difference
      SSE <- SSE + value                   # summing the squares of differences
    }
  }
  
  SSR <- 0
  
  for (i in 1:length(group_means)){   # looping through the group means vector
    n <- length(df[,i])                         # size of group i
    value <- n*(group_means[i] - total_mean)^2  # single value that will add up to SSR
    SSR <- SSR + value                          # calculating SSR
  }
  
  SST = SSE + SSR
  n <- nrow(df)*ncol(df)
  k <- ncol(df)
  MSR <- SSR / (k-1)
  MSE <- SSE / (n-k)
  F <- MSR / MSE
  critical_value <- qf(sign_level, k-1, n-k, lower.tail = FALSE)
  p_val <- pf(F, k-1, n-k, lower.tail = FALSE)
  
  if (F <= critical_value){
    decision <- paste0('Failed to reject null hypothesis (significance level ', sign_level, ')')
  }
  else{
    decision <- paste0('Null hypothesis can be rejected (significance level ', sign_level, ')')
  }
  
  if (p_val <= 1){
    list(message = 'Results of one-way ANOVA:', MSR = MSR, MSE = MSE, F = F, p_value = p_val, decision = decision)
  }
  
  else if (p_val > 1){
    list(message = 'Results of one-way ANOVA:', MSR = MSR, MSE = MSE, F = F, p_value = 1, decision = decision)
  }
}
