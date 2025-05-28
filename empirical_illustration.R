########################################
# Empirical example                    #
# Version: may-2025                    #
########################################

library(mirt)
library(psych)
library(readxl)
library(haven)
library(tidyverse)
source("simulation_functions.R")
load("empirical_illustration_data.RData")

# Data -------------------------------------------------------------------------

items <- read_excel("IPIP_NEO_ItemKey.xlsx")
items$item_num <- 1:nrow(items)
items$Pol <- as.numeric(items$Pol)

ipip300 <- read_por("IPIP300.POR")
ipip300 <- ipip300 %>%
  filter(if_all(11:ncol(.), ~ . != 0))
ipip300_usa_age <- ipip300 %>% filter(COUNTRY == "USA", AGE >= 19, AGE <= 25)

set.seed(1234)
random_sample <- sample(1:nrow(ipip300_usa_age), 1000)
ipip300_usa_age_s <- ipip300_usa_age[random_sample,]
ipip286_usa_age_s <- ipip300_usa_age_s[,items$NameDataset]

cutoff_mixed <- 1.125
cutoff_equal <- 0.75
use_SD <- TRUE
balance_dim_comb <- FALSE

# Graded Model Analysis --------------------------------------------------------

dimnames <- c("N", "E", "O", "A", "C")

Q <- matrix(0, nrow = nrow(items), ncol = 5)
colnames(Q) <- dimnames
for (i in 1:nrow(items)) {
  trait <- items$Trait[i]
  if (trait %in% dimnames) {
    Q[i, trait] <- 1
  }
}

model.mirt <- paste0(
  "N = ", paste(which(Q[, "N"] == 1), collapse = ", "), "\n",
  "E = ", paste(which(Q[, "E"] == 1), collapse = ", "), "\n",
  "O = ", paste(which(Q[, "O"] == 1), collapse = ", "), "\n",
  "A = ", paste(which(Q[, "A"] == 1), collapse = ", "), "\n",
  "C = ", paste(which(Q[, "C"] == 1), collapse = ", "), "\n",
  "COV = N*E*O*A*C"
)

sv <- mirt(ipip286_usa_age_s, 
           model = model.mirt, 
           itemtype = "graded", 
           pars = "values")

for (i in seq_len(nrow(items))) {
  item_name <- items$NameDataset[i]
  trait <- items$Trait[i]
  pol <- items$Pol[i]
  
  if (pol != -1) next
  
  a_col_index <- match(trait, dimnames)
  a_param_name <- paste0("a", a_col_index)
  
  sv_idx <- which(sv$item == item_name & sv$name == a_param_name)
  
  if (length(sv_idx) == 1) {
    sv$value[sv_idx] <- -abs(sv$value[sv_idx])
  }
}

mirt.est <- mirt(ipip286_usa_age_s, 
                 model = mirt.model(model.mirt), 
                 itemtype = "graded", 
                 pars = sv, 
                 method = "MHRM")
param <- coef(mirt.est, simplify = TRUE)$items

param_pol <- vector("numeric", length = nrow(param))
for (i in 1:nrow(param)) {
  item_values <- param[i, 1:5]
  non_zero_index <- which(item_values != 0)[1]
  
  if (item_values[non_zero_index] > 0) {
    param_pol[i] <- 1
  } else if (item_values[non_zero_index] < 0) {
    param_pol[i] <- -1
  }
}

items$Param_pol <- param_pol
items$Test <- param_pol == items$Pol
items <- cbind(items, param)
items$Trait_num <- recode(items$Trait, "N" = 1, "E" = 2, "O" = 3, "A" = 4, "C" = 5)
items$Scale <- rowSums(items[, c("a1", "a2", "a3", "a4", "a5")], na.rm = TRUE)
items$Threshold <- items$d3

bank <- data.frame(
  item_num = items$item_num,
  trait = items$Trait_num,
  pol = items$Pol,
  Scale = items$Scale,
  Threshold = items$Threshold,
  SD_rating = items$Mean_SD
)

Sigma <- as.matrix(coef(mirt.est, simplify = TRUE)$cov)

# FCQ 1 ------------------------------------------------------------------------
# With SDM, all homopolar blocks with a mixed-keyed bank, total of 70 blocks

n_blocks_total <- 70
degree_het <- 3
n_blocks_hetero <- 0

######################################## 
## GA                                  #
######################################## 

# Optimization
result_GA_1 <- optimization(data_item_chars = bank, 
                            algorithm = "GA", 
                            use_SD = use_SD, 
                            cutoff_mixed = cutoff_mixed, 
                            cutoff_equal = cutoff_equal, 
                            Sigma = Sigma, 
                            n_blocks_total = n_blocks_total, 
                            n_blocks_hetero = n_blocks_hetero, 
                            nCores = 7, 
                            balance_dim_comb = balance_dim_comb, 
                            which.combs = NULL, 
                            type_theta = "Quad_with_Sigma")

solution_GA_1 <- result_GA_1$solution
result_GA_1$num_polarity
result_GA_1$combinations

# Obtain empirical reliability
res_GA_1 <- res_empirical(solution = solution_GA_1, 
                          Sigma = Sigma)
res_thetas_GA_1 <- res_GA_1$res

######################################## 
# SA blueprint                         #
######################################## 

# Optimization
result_SA_bp_1 <- optimization(data_item_chars = bank, 
                               algorithm = "SA_blueprint", 
                               use_SD = use_SD, 
                               cutoff_mixed = cutoff_mixed, 
                               cutoff_equal = cutoff_equal, 
                               n_blocks_total = n_blocks_total, 
                               n_blocks_hetero = n_blocks_hetero)

solution_SA_bp_1 <- result_SA_bp_1$solution
result_SA_bp_1$num_polarity
result_SA_bp_1$combinations

# Obtain empirical reliability
res_SA_bp_1 <- res_empirical(solution = solution_SA_bp_1, 
                             Sigma = Sigma)
res_thetas_SA_bp_1 <- res_SA_bp_1$res

######################################## 
## SA with optimization of a parameter #
######################################## 

# Optimization
result_SA_param_1 <- optimization(data_item_chars = bank, 
                                  algorithm = "SA_with_param", 
                                  use_SD = use_SD, 
                                  cutoff_mixed = cutoff_mixed, 
                                  cutoff_equal = cutoff_equal, 
                                  n_blocks_total = n_blocks_total, 
                                  n_blocks_hetero = n_blocks_hetero)

solution_SA_param_1 <- result_SA_param_1$solution
result_SA_param_1$num_polarity
result_SA_param_1$combinations

# Obtain empirical reliability
res_SA_param_1 <- res_empirical(solution = solution_SA_param_1, 
                                Sigma = Sigma)
res_thetas_SA_param_1 <- res_SA_param_1$res

######################################## 
# BF                                   #
######################################## 

# Optimization
result_BF_1 <- optimization(data_item_chars = bank, 
                            algorithm = "BF", 
                            use_SD = use_SD, 
                            Sigma = Sigma, 
                            n_blocks_total = n_blocks_total, 
                            n_blocks_hetero = n_blocks_hetero, 
                            BF_iter = 100, 
                            BF_relia_method = "empirical")

solution_BF_1 <- result_BF_1$solution
result_BF_1$num_polarity
result_BF_1$combinations

# Obtain empirical reliability
res_BF_1 <- res_empirical(solution = solution_BF_1, 
                          Sigma = Sigma)
res_thetas_BF_1 <- res_BF_1$res

# FCQ 2 ------------------------------------------------------------------------
# With SDM, 20% heteropolar blocks with a mixed-keyed bank, total of 70 blocks

n_blocks_total <- 70
degree_het <- 3
n_blocks_hetero <- 14  # 20% heteropolar

######################################## 
## GA                                  #
######################################## 

# Optimization
result_GA_2 <- optimization(data_item_chars = bank, 
                            algorithm = "GA", 
                            use_SD = use_SD, 
                            cutoff_mixed = cutoff_mixed, 
                            cutoff_equal = cutoff_equal, 
                            Sigma = Sigma, 
                            n_blocks_total = n_blocks_total, 
                            n_blocks_hetero = n_blocks_hetero, 
                            nCores = 7, 
                            balance_dim_comb = balance_dim_comb, 
                            which.combs = NULL, 
                            type_theta = "Quad_with_Sigma")

solution_GA_2 <- result_GA_2$solution
result_GA_2$num_polarity
result_GA_2$combinations

# Obtain empirical reliability
res_GA_2 <- res_empirical(solution = solution_GA_2, 
                          Sigma = Sigma)
res_thetas_GA_2 <- res_GA_2$res

######################################## 
# SA blueprint                         #
########################################

# Optimization
result_SA_bp_2 <- optimization(data_item_chars = bank, 
                               algorithm = "SA_blueprint", 
                               use_SD = use_SD, 
                               cutoff_mixed = cutoff_mixed, 
                               cutoff_equal = cutoff_equal, 
                               n_blocks_total = n_blocks_total, 
                               n_blocks_hetero = n_blocks_hetero)

solution_SA_bp_2 <- result_SA_bp_2$solution
result_SA_bp_2$num_polarity
result_SA_bp_2$combinations

# Obtain empirical reliability
res_SA_bp_2 <- res_empirical(solution = solution_SA_bp_2, 
                             Sigma = Sigma)
res_thetas_SA_bp_2 <- res_SA_bp_2$res

######################################## 
## SA with optimization of a parameter #
######################################## 

# Optimization
result_SA_param_2 <- optimization(data_item_chars = bank, 
                                  algorithm = "SA_with_param", 
                                  use_SD = use_SD, 
                                  cutoff_mixed = cutoff_mixed, 
                                  cutoff_equal = cutoff_equal, 
                                  n_blocks_total = n_blocks_total, 
                                  n_blocks_hetero = n_blocks_hetero)

solution_SA_param_2 <- result_SA_param_2$solution
result_SA_param_2$num_polarity
result_SA_param_2$combinations

# Obtain empirical reliability
res_SA_param_2 <- res_empirical(solution = solution_SA_param_2, 
                                Sigma = Sigma)
res_thetas_SA_param_2 <- res_SA_param_2$res

######################################## 
# BF                                   #
######################################## 

# Optimization
result_BF_2 <- optimization(data_item_chars = bank, 
                            algorithm = "BF", 
                            use_SD = use_SD, 
                            Sigma = Sigma, 
                            n_blocks_total = n_blocks_total, 
                            n_blocks_hetero = n_blocks_hetero, 
                            BF_iter = 100, 
                            BF_relia_method = "empirical")

solution_BF_2 <- result_BF_2$solution
result_BF_2$num_polarity
result_BF_2$combinations

# Obtain empirical reliability
res_BF_2 <- res_empirical(solution = solution_BF_2, 
                          Sigma = Sigma)
res_thetas_BF_2 <- res_BF_2$res

# FCQ 3 ------------------------------------------------------------------------
# With SDM, all homopolar blocks with a mixed-keyed bank, total of 35 blocks

n_blocks_total <- 35
degree_het <- 3
n_blocks_hetero <- 0  

######################################## 
## GA                                  #
######################################## 

# Optimization
result_GA_3 <- optimization(data_item_chars = bank, 
                            algorithm = "GA", 
                            use_SD = use_SD, 
                            cutoff_mixed = cutoff_mixed, 
                            cutoff_equal = cutoff_equal, 
                            Sigma = Sigma, 
                            n_blocks_total = n_blocks_total, 
                            n_blocks_hetero = n_blocks_hetero, 
                            nCores = 7, 
                            balance_dim_comb = balance_dim_comb, 
                            which.combs = NULL, 
                            type_theta = "Quad_with_Sigma")

solution_GA_3 <- result_GA_3$solution
result_GA_3$num_polarity
result_GA_3$combinations

# Obtain empirical reliability
res_GA_3 <- res_empirical(solution = solution_GA_3, 
                          Sigma = Sigma)
res_thetas_GA_3 <- res_GA_3$res

######################################## 
# SA blueprint                         #
########################################

# Optimization
result_SA_bp_3 <- optimization(data_item_chars = bank, 
                               algorithm = "SA_blueprint", 
                               use_SD = use_SD, 
                               cutoff_mixed = cutoff_mixed, 
                               cutoff_equal = cutoff_equal, 
                               n_blocks_total = n_blocks_total, 
                               n_blocks_hetero = n_blocks_hetero)

solution_SA_bp_3 <- result_SA_bp_3$solution
result_SA_bp_3$num_polarity
result_SA_bp_3$combinations

# Obtain empirical reliability
res_SA_bp_3 <- res_empirical(solution = solution_SA_bp_3, 
                             Sigma = Sigma)
res_thetas_SA_bp_3 <- res_SA_bp_3$res

######################################## 
## SA with optimization of a parameter #
######################################## 

# Optimization
result_SA_param_3 <- optimization(data_item_chars = bank, 
                                  algorithm = "SA_with_param", 
                                  use_SD = use_SD, 
                                  cutoff_mixed = cutoff_mixed, 
                                  cutoff_equal = cutoff_equal, 
                                  n_blocks_total = n_blocks_total, 
                                  n_blocks_hetero = n_blocks_hetero)

solution_SA_param_3 <- result_SA_param_3$solution
result_SA_param_3$num_polarity
result_SA_param_3$combinations

# Obtain empirical reliability
res_SA_param_3 <- res_empirical(solution = solution_SA_param_3,  
                                Sigma = Sigma)
res_thetas_SA_param_3 <- res_SA_param_3$res

######################################## 
# BF                                   #
######################################## 

# Optimization
result_BF_3 <- optimization(data_item_chars = bank, 
                            algorithm = "BF", 
                            use_SD = use_SD, 
                            Sigma = Sigma, 
                            n_blocks_total = n_blocks_total, 
                            n_blocks_hetero = n_blocks_hetero, 
                            BF_iter = 100, 
                            BF_relia_method = "empirical")

solution_BF_3 <- result_BF_3$solution
result_BF_3$num_polarity
result_BF_3$combinations

# Obtain empirical reliability
res_BF_3 <- res_empirical(solution = solution_BF_3, 
                          Sigma = Sigma)
res_thetas_BF_3 <- res_BF_3$res

# FCQ 4 ------------------------------------------------------------------------
# With SDM, 20% heteropolar blocks with a mixed-keyed bank, total of 35 blocks

n_blocks_total <- 35
degree_het <- 3
n_blocks_hetero <- 7  # 20% heteropolar

######################################## 
## GA                                  #
######################################## 

# Optimization
result_GA_4 <- optimization(data_item_chars = bank, 
                            algorithm = "GA", 
                            use_SD = use_SD, 
                            cutoff_mixed = cutoff_mixed, 
                            cutoff_equal = cutoff_equal, 
                            Sigma = Sigma, 
                            n_blocks_total = n_blocks_total, 
                            n_blocks_hetero = n_blocks_hetero, 
                            nCores = 7, 
                            balance_dim_comb = balance_dim_comb, 
                            which.combs = NULL, 
                            type_theta = "Quad_with_Sigma")

solution_GA_4 <- result_GA_4$solution
result_GA_4$num_polarity
result_GA_4$combinations

# Obtain empirical reliability
res_GA_4 <- res_empirical(solution = solution_GA_4, 
                          Sigma = Sigma)
res_thetas_GA_4 <- res_GA_4$res

######################################## 
# SA blueprint                         #
########################################

# Optimization
result_SA_bp_4 <- optimization(data_item_chars = bank, 
                               algorithm = "SA_blueprint", 
                               use_SD = use_SD, 
                               cutoff_mixed = cutoff_mixed, 
                               cutoff_equal = cutoff_equal, 
                               n_blocks_total = n_blocks_total, 
                               n_blocks_hetero = n_blocks_hetero)

solution_SA_bp_4 <- result_SA_bp_4$solution
result_SA_bp_4$num_polarity
result_SA_bp_4$combinations

# Obtain empirical reliability
res_SA_bp_4 <- res_empirical(solution = solution_SA_bp_4, 
                             Sigma = Sigma)
res_thetas_SA_bp_4 <- res_SA_bp_4$res

######################################## 
## SA with optimization of a parameter #
######################################## 

# Optimization
result_SA_param_4 <- optimization(data_item_chars = bank, 
                                  algorithm = "SA_with_param", 
                                  use_SD = use_SD, 
                                  cutoff_mixed = cutoff_mixed, 
                                  cutoff_equal = cutoff_equal, 
                                  n_blocks_total = n_blocks_total, 
                                  n_blocks_hetero = n_blocks_hetero)

solution_SA_param_4 <- result_SA_param_4$solution
result_SA_param_4$num_polarity
result_SA_param_4$combinations

# Obtain empirical reliability
res_SA_param_4 <- res_empirical(solution = solution_SA_param_4, 
                                Sigma = Sigma)
res_thetas_SA_param_4 <- res_SA_param_4$res

######################################## 
# BF                                   #
######################################## 

# Optimization
result_BF_4 <- optimization(data_item_chars = bank, 
                            algorithm = "BF", 
                            use_SD = use_SD, 
                            Sigma = Sigma, 
                            n_blocks_total = n_blocks_total, 
                            n_blocks_hetero = n_blocks_hetero, 
                            BF_iter = 100, 
                            BF_relia_method = "empirical")

solution_BF_4 <- result_BF_4$solution
result_BF_4$num_polarity
result_BF_4$combinations

# Obtain empirical reliability
res_BF_4 <- res_empirical(solution = solution_BF_4, 
                          Sigma = Sigma)
res_thetas_BF_4 <- res_BF_4$res
