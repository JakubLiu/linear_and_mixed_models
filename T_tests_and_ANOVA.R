# a function that calls upon custom written T tests and one way ANOVA
# the user can define which test needs to be performed and specify other metrics like the alternative hypothesis or the significance level
# there are 5 custom functions
#   - one sample T test
#   - two sample unpaired T test
#   - two sample paired T test
#   - one way ANOVA
#   - a wrapper function
# QBA LIU 2023 ---------------------------------------------------------------------------------------------------------------------------------


# one sample t test
one_sample_t_test <- function(sample, mu0, alternative){
  
  mean_sample <- mean(sample)
  std_dev_sample <- sd(sample)
  n_sample <- length(sample)
  df <- n_sample - 1
  T <- sqrt(n_sample)*(as.numeric(mean_sample) - as.numeric(mu0))/std_dev_sample
  
  if (alternative == 1){     # two sided
    p_value <- pt(T, df)
  }
  else if (alternative == 2){  # greater
    p_value <- pt(T, df, lower.tail = FALSE)
  }
  else if (alternative == 3) {   # less
    p_value <- 2*pt(abs(T), df, lower.tail = FALSE)
  }
  
  list(message = 'One sample T-test.', test_statistic = T, p_value = p_value)
}

# two sample unpaired t test
two_sample_unpaired_T_test <- function(sample1, sample2, alternative){
  
  mean_sample1 <- mean(sample1)
  mean_sample2 <- mean(sample2)
  var_sample1 <- var(sample1)
  var_sample2 <- var(sample2)
  n_sample1 <- length(sample1)
  n_sample2 <- length(sample2)
  df <- n_sample1 + n_sample2 - 2
  T = (mean_sample1-mean_sample2)/sqrt((n_sample1-1)*var_sample1+(n_sample2-1)*var_sample2)*sqrt((n_sample1*n_sample2*df) / (n_sample1+n_sample2))
  
  if (alternative == 1){    # two sided
    p_value <- pt(T, df)
  }
  else if (alternative == 2){   # greater
    p_value <- pt(T, df, lower.tail = FALSE)
  }
  else if (alternative == 3) {    # less
    p_value <- 2*pt(abs(T), df, lower.tail = FALSE)
  }
  
  list(message = 'Two sample unpaired T-test.', test_statistic = T, p_value = p_value)
}

# two sample paired t test
two_sample_paired_T_test <- function(sample1, sample2, alternative){
  
  D <- sample1 - sample2
  mean_D <- mean(D)
  std_dev_D <- sd(D)
  n_D <- length(D)
  df <- n_D - 1
  T <- sqrt(n_D)*mean_D/std_dev_D
  
  if (alternative == 1){    # two sided
    p_value <- pt(T, df)
  }
  else if (alternative == 2){    # greater
    p_value <- pt(T, df, lower.tail = FALSE)
  }
  else if (alternative == 3) {   # less
    p_value <- 2*pt(abs(T), df, lower.tail = FALSE)
  }
  
  list(message = 'Two sample paired T-test.', test_statistic = T, p_value = p_value)
}

# one way ANOVA
one_way_anova <- function(data, sign_level){
  
  df <- data
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
  
  if (F <= critical_value){
    decision <- paste0('Failed to reject null hypothesis (significance level ', sign_level, ')')
  }
  else{
    decision <- paste0('Null hypothesis can be rejected (significance level ', sign_level, ')')
  }
  
  list(message = 'One way ANOVA', MSR = MSR, MSE = MSE, F = F, decision = decision)
}

# wrapper
Statistics <- function(sample1 = 'empty', sample2 = 'empty', data = 'empty', sign_level = 0.05){
  
  test_version <- readline(prompt = "Please choose statistical test version.
    - one sample T-test ---> press 1
    - two sample unpaired T-test ---> press 2
    - two sample paired T-test ---> press 3
    - one way ANOVA ---> press 4
    : ")
  
  if (test_version != 4){
    input_alternative <- readline(prompt = "Please specify alternative hypothesis.
        - two sided ---> press 1
        - greater ---> press 2
        - less ---> press 3
        : ")
  }
  
  
  # ONE SAMPLE T TEST --------------------------------------------------------------------------------------------
  if (test_version == 1){
    input_mu0 <- readline(prompt = "Please specify mu0 value
    : ")
    input_mu0 <- as.numeric(input_mu0)
    
    one_sample_t_test(sample = sample1, mu0 = input_mu0, alternative = input_alternative)
  }
  # --------------------------------------------------------------------------------------------------------------
  # TWO SAMPLE UNPAIRED T TEST -----------------------------------------------------------------------------------
  else if (test_version == 2){
    two_sample_unpaired_T_test(sample1 = sample1, sample2 = sample2, alternative = input_alternative)
  }
  # --------------------------------------------------------------------------------------------------------------
  # TWO SAMPLE PAIRED T TEST -------------------------------------------------------------------------------------
  else if (test_version == 3){
    two_sample_paired_T_test(sample1 = sample1, sample2 = sample2, alternative = input_alternative)
  }
  # -------------------------------------------------------------------------------------------------------------
  # ONE WAY ANOVA -----------------------------------------------------------------------------------------------
  else if (test_version == 4){
    one_way_anova(data = data, sign_level = sign_level)
  }
  # ------------------------------------------------------------------------------------------------------------
}

# test for T tests
proba = rnorm(100, 23, 10)
proba2 = rnorm(100, 21, 8)
Statistics(sample1 = proba, sample2 = proba2)

# test for ANOVA
group1 <- c(85, 86, 88, 75, 78, 94, 98, 79, 71, 80)
group2 <- c(91, 92, 93, 85, 87, 84, 82, 88, 95, 96)
group3 <- c(79, 78, 88, 94, 92, 85, 83, 85, 82, 81)
my_df <- data.frame(group1, group2, group3)
Statistics(sample1 = proba, sample2 = proba2, data = my_df, sign_level = 0.01)