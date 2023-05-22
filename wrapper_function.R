# this function uses the function written in T_tests_and_ANOVA.R
# it enables the user to specify which test to use, aswell as specify the alternative hypothesis
# simple error handling is also included
# QBA LIU, Jakub Prorok 2023 -------------------------------------------------------------------------------------------------------------

# path to source code
source("C:/Users/Lenovo/Desktop/STUDIA/BIOINFORMATYKA/SEMESTR_VI/MODELE_LINOWE_I_MIESZANE/projekt_zaliczeniowy/program1/T_tests_ANOVA.r")

Statistics <- function(sample_in_1 = c(1), sample_in_2 = c(1), data_in, sign_level = 0.05){
  
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
    
    one_sample_t_test(sample = sample_in_1, mu0 = input_mu0, alternative = input_alternative)
  }
  # --------------------------------------------------------------------------------------------------------------
  # TWO SAMPLE UNPAIRED T TEST -----------------------------------------------------------------------------------
  else if (test_version == 2){
    
    if (length(sample_in_1) == 1 || length(sample_in_2) == 1){
      stop('Please provide both samples.')
    }
   
    two_sample_unpaired_T_test(sample1 = sample_in_1, sample2 = sample_in_2, alternative = input_alternative)
  }
  # --------------------------------------------------------------------------------------------------------------
  # TWO SAMPLE PAIRED T TEST -------------------------------------------------------------------------------------
  else if (test_version == 3){
    
      if (length(sample_in_1) == 1 || length(sample_in_2) == 1){
        stop('Please provide both samples.')
      }
      
      test <- length(sample_in_1) == length(sample_in_2)
      
      if (test == FALSE){
        stop('Both samples have to be of equal lengths.')
      }
    
    two_sample_paired_T_test(sample1 = sample_in_1, sample2 = sample_in_2, alternative = input_alternative)
  }
  # -------------------------------------------------------------------------------------------------------------
  # ONE WAY ANOVA -----------------------------------------------------------------------------------------------
  else if (test_version == 4){
    one_way_anova(data = data_in, sign_level = sign_level)
  }
  # ------------------------------------------------------------------------------------------------------------
}
