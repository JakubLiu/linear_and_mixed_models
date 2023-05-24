# this function allows the user to decide which correlation coefficient to calculate
# the wrapper function calls the subfunctions defined in 'correlations_sourcedore.R'
# QBA LIU 2023, Jakub Prorok ---------------------------------------------------------------------------------------------------------------

# path to source code
source("C:/Users/Lenovo/Desktop/STUDIA/BIOINFORMATYKA/SEMESTR_VI/MODELE_LINOWE_I_MIESZANE/projekt_zaliczeniowy/correlations_sourcecode.R")

# wrapper function
correlation <- function(vector1, vector2, method){
  if (method == 'Pearson'){pearson_correlation(x = vector1, y = vector2)}
  else if (method == 'Spearman'){spearman_correlation(x = vector1, y = vector2)}
  else if (method == 'Kendall'){kendall_correlation(x = vector1, y = vector2)}
  else {stop('Plese select one of the following methods: Pearson, Spearman, Kendall')}
}

# testing
vec1 <- rnorm(500)
vec2 <- rnorm(500)
correlation(vector1 = vec1, vector2 = vec2, method = 'Kendall')

