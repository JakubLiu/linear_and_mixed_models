# three functions that compute:
#   - Pearsons correlation coefficient
#   - Spearmans correlation coefficient
#   - Tau Kendalls correlation coefficient
# QBA LIU 2023, Jakub Prorok -------------------------------------------------------------------------------------------------------------------


# Pearson
pearson_correlation <- function(x,y){
  data <- data.frame(cbind(x,y))                                        # checking for missing values
  all_observations <- nrow(data)
  data <- data[!is.na(data$x) & !is.na(data$y), ]
  cleaned_observations <- nrow(data)
  diff <- all_observations - cleaned_observations
  message <- paste0('Omitted ', diff, ' missing value observation(s).')
  if (diff > 0){print(message)}
  if (shapiro.test(x)$p.value < 0.05 || shapiro.test(y)$p.value < 0.05)# testing normality assumption
  {warning('Warning, normality assumption not met.')}
  mean_x <- mean(data$x)
  mean_y <- mean(data$y)
  numerator <- sum((x-mean_x)*(y-mean(y)))
  denominator <- sqrt(sum((x-mean_x)^2) * sum((y-mean_y)^2))
  pearson_coef <- numerator / denominator                              # calculating coefficient
  T <- pearson_coef * sqrt((length(x)-2) / (1-(pearson_coef)^2))
  p_value <- 2*pt(abs(T), df = length(x)-2, lower.tail = FALSE)        # calculating significance of coefficient
  list(message = "Pearsons correlation coefficient.", coefficient = pearson_coef, p_value = p_value)
  
}

# Spearman
spearman_correlation <- function(x,y){
  #function to create ranks
  ranking <- function(vector_in){
    ranked <- 1:length(vector_in)
    for (i in sort(unique(vector_in))){
      x1 <- which(sort(vector_in) == i)
      x2 <- c(which(vector_in == i))    
      for (i in 1:length(x2)){
        suppressWarnings(ranked[x2[i]] <- rep(mean(x1), length(x1)))
      }
    }
    return(ranked)
  }
  # function to calculate Spearman's coefficient
  rank_x <- ranking(vector_in = x)                  # creating ranks
  rank_y <- ranking(vector_in = y)
  n <- length(x)
  numerator <- 6*sum((rank_x - rank_y)^2)
  denominator <- n*(n^2 - 1)
  spearman_coef <- 1 - numerator / denominator         # calculating coefficient
  T <- spearman_coef * sqrt((n-2)/(1-spearman_coef^2))
  p_value = 2*pt(abs(T), df = n-2, lower.tail = FALSE)  # calculating significance of coefficient
  list(message = 'Spearman correlation.', coefficient = spearman_coef, p_value = p_value)
  
}

# Kendall
kendall_correlation <- function(x, y){
  
  # function to create ranks
  ranking <- function(vector_in){
    ranked <- 1:length(vector_in)
    for (i in sort(unique(vector_in))){
      x1 <- which(sort(vector_in) == i)
      x2 <- c(which(vector_in == i))    
      for (i in 1:length(x2)){
        suppressWarnings(ranked[x2[i]] <- rep(mean(x1), length(x1)))
      }
    }
    return(ranked)
  }
  
  # calculating coefficient
  C <- rep(0, length(x))
  D <- rep(0, length(y))
  data <- data.frame(cbind(ranking(vector_in = x), ranking(vector_in = y), C, D)) # creating ranks
  colnames(data) <- c('vector_1', 'vector_2', 'C', 'D')
  data <- data[order(data$vector_1),] # ordering data in ascending order according to one of the vectors
  
  # calculating number of concordonate and disconcordonate pairs
  for (i in 1:nrow(data)){
    C_sum <- 0
    D_sum <- 0
    for (j in (i+1):nrow(data)){
      if (j <= nrow(data)){
        if (data$vector_2[i] < data$vector_2[j]){C_sum <- C_sum + 1}
        else {D_sum <- D_sum + 1}
      }
    }
    data$C[i] <- C_sum
    data$D[i] <- D_sum
  }
  
  data$C[nrow(data)] <- 0
  data$D[nrow(data)] <- 0
  Tau <- (sum(data$C) - sum(data$D)) / (sum(data$C) + sum(data$D)) # calculating Tau Kendalls coefficient
  n <- nrow(data)
  numerator <- 3*Tau*sqrt(n*(n-1))
  denominator <- sqrt(2*(2*n + 5))
  Z <- numerator / denominator
  p_value <- 2*pnorm(abs(Z), lower.tail=FALSE) # calculating significnce of coefficient
  list(message = 'Tau Kendalls correlation coefficient.', coefficient = Tau, p_value = p_value)
  
}