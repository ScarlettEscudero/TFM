########################################
# Simulation                           #
# Version: may-2025                    #
########################################

library(foreach)
source("simulation_functions_v3.R")
load("FCoptim_readyforsim.RData")

# Conditions ------------------------------------------------------------------
ratio <- c(4, 8) # ratio 4:1, and 8:1 # Please note that this factor was originally called item-to-block ratio, which is still called this way in the code, but is re-named length in the manuscript, as it determines the questionnaire length (ratio 4:1 = length 80 blocks, ratio 8:1 = 40 blocks)
SD <- c(T, F)
degree_het <- c(1, # 1 = positive items with moderate cor (r = 0.5), thus 0% heteropolar
                2, # 2 = positive items with high cor (r = 0.8), thus 0% heteropolar
                3, # 3 = positive and negative items with high cor (r = 0.9), and 0% heteropolar blocks
                4) # 4 = positive and negative items with high cor (r = 0.9), and 25% heteropolar blocks
replicas <- 1:10

conditions <- expand.grid(
  ratio = ratio, # this will determine the FC test length
  SD = SD, 
  degree_het = degree_het,
  replicas = replicas)

conditions$condition_id <- rep(1:16, 
                               nrow(conditions)/16)
rownames(conditions) <- NULL

set.seed(0635)
conditions$seed <- as.numeric(paste0(seq(1:nrow(conditions)), sample(1:10000, nrow(conditions))))

# Simulation -------------------------------------------------------------------
# Data needed
set.seed(1720)
J <- 320
K <- 5
Sigma <- matrix(c(1, 0.21, 0, 0.25, 0.53,
                  0.21, 1, 0.4, 0, 0.27,
                  0, 0.4, 1, 0, 0,
                  0.25, 0, 0, 1, 0.24,
                  0.53, 0.27, 0, 0.24, 1), nrow = 5, ncol = 5, byrow = TRUE)  
thetas_true <- mvtnorm::rmvnorm(5000, rep(0,K), sigma = Sigma) # true thetas

cutoff_mixed <- 0.75
cutoff_equal <- 0.5

# Prepare lists for storage
solutions_file <- "all_tests.rds"
if (file.exists(solutions_file)) {
  all_solutions <- readRDS(solutions_file)
} else {
  all_solutions <- list()
}

thetas_GA_list <- list()
thetas_SA_bp_list <- list()
thetas_SA_param_list <- list()
thetas_BF_list <- list()

all_marg_reli_GA_list <- list()
all_marg_reli_SA_bp_list <- list()
all_marg_reli_SA_param_list <- list()
all_marg_reli_BF_list <- list()

# Implement simulation
start <- Sys.time()
foreach(i = 1:nrow(conditions), 
        .inorder = F) %do% { 
            try({
              
              print(paste0(i, " out of ", nrow(conditions))) # Progress of the simulation
              
              start_cond <- Sys.time()
              
              conds <- conditions[i, ]
              set.seed(conds$seed)
              
              cond_name <- paste("cond", conds$ratio, conds$SD, conds$degree_het, conds$replicas, conds$condition_id, conds$seed, sep = "_")

              # Generate bank and get necessary conditions
              bank_i <- generate_data(conds, J = J)

              n_blocks_total <- J/(conds$ratio)

              if (conds$degree_het == 1){ 
                n_blocks_hetero <- 0
              } else if (conds$degree_het == 2){ 
                n_blocks_hetero <- 0
              } else if (conds$degree_het == 3){ 
                n_blocks_hetero <- 0
              }  else if (conds$degree_het == 4){
                n_blocks_hetero <- (n_blocks_total)*0.25
              }
              
              use_SD <- conds$SD
              balance_dim_comb <- FALSE
              
              # Optimization and obtain results of each algorithm
              
              ## GA
              ### Optimization
              start_GA <- Sys.time()
              result_GA <- optimization(data_item_chars = bank_i, algorithm = "GA", use_SD = use_SD, cutoff_mixed = cutoff_mixed, cutoff_equal = cutoff_equal, Sigma = Sigma, n_blocks_total = n_blocks_total, n_blocks_hetero = n_blocks_hetero, nCores = 7, balance_dim_comb = balance_dim_comb, which.combs = NULL, type_theta = "Quad_with_Sigma")
              end_GA <- Sys.time()
              duration_GA <-  paste(round(as.numeric(end_GA-start_GA),2), attr(end_GA-start_GA, "units"))

              solution_GA <- result_GA$solution
              comb_pol_GA_hom <- result_GA$num_polarity[1]
              names(comb_pol_GA_hom) <- paste0(names(comb_pol_GA_hom), "_GA")
              comb_pol_GA_het <- result_GA$num_polarity[2]
              names(comb_pol_GA_het) <- paste0(names(comb_pol_GA_het), "_GA")
              comb_trait_GA <- result_GA$combinations
              names(comb_trait_GA) <- paste0(names(comb_trait_GA), "_GA")
              count_traits_GA <- result_GA$count_traits
              names(count_traits_GA) <- paste0(names(count_traits_GA), "_GA")

              ### Obtain estimated thetas and reliability
              res_GA <- get_res(solution = solution_GA, thetas_true = thetas_true, Sigma = Sigma, method = "both")
              res_thetas_GA <- res_GA$res
              names(res_thetas_GA) <- paste0(names(res_thetas_GA), "_GA")
              
              #### Save estimated thetas and marginal reliability for all thetas
              thetas_est_GA <- res_GA$thetas_est_com
              thetas_GA_list[[cond_name]] <- thetas_est_GA
              if (file.exists("thetas_est/thetas_est_GA.rds")) {
                thetas_GA_list <- readRDS("thetas_est/thetas_est_GA.rds")
              } else {
                thetas_GA_list <- list()
              }
              thetas_GA_list[[cond_name]] <- thetas_est_GA
              saveRDS(thetas_GA_list, file = "thetas_est/thetas_est_GA.rds")
              
              all_marg_reli_GA <- res_GA$all_marg_reli
              all_marg_reli_GA_list[[cond_name]] <- all_marg_reli_GA
              if (file.exists("all_marg_reli/all_marg_reli_GA.rds")) {
                all_marg_reli_GA_list <- readRDS("all_marg_reli/all_marg_reli_GA.rds")
              } else {
                all_marg_reli_GA_list <- list()
              }
              all_marg_reli_GA_list[[cond_name]] <- all_marg_reli_GA
              saveRDS(all_marg_reli_GA_list, file = "all_marg_reli/all_marg_reli_GA.rds")
              
          
              
              # SA blueprint
              ## Optimization
              start_SA_bp <- Sys.time()
              result_SA_bp <- optimization(data_item_chars = bank_i, algorithm = "SA_blueprint", use_SD = use_SD, cutoff_mixed = cutoff_mixed, cutoff_equal = cutoff_equal, n_blocks_total = n_blocks_total, n_blocks_hetero = n_blocks_hetero)
              end_SA_bp <- Sys.time()
              duration_SA_bp <- paste(round(as.numeric(end_SA_bp-start_SA_bp),2), attr(end_SA_bp-start_SA_bp, "units"))
              
              solution_SA_bp <- result_SA_bp$solution
              comb_pol_SA_bp_hom <- result_SA_bp$num_polarity[1]
              names(comb_pol_SA_bp_hom) <- paste0(names(comb_pol_SA_bp_hom), "_SA_bp")
              comb_pol_SA_bp_het <- result_SA_bp$num_polarity[2]
              names(comb_pol_SA_bp_het) <- paste0(names(comb_pol_SA_bp_het), "_SA_bp")
              comb_trait_SA_bp <- result_SA_bp$combinations
              names(comb_trait_SA_bp) <- paste0(names(comb_trait_SA_bp), "_SA_bp")
              count_traits_SA_bp <- result_SA_bp$count_traits
              names(count_traits_SA_bp) <- paste0(names(count_traits_SA_bp), "_SA_bp")
              iter_num_SA_bp <- result_SA_bp$iter_count
              names(iter_num_SA_bp) <- paste0("iter_num", "_SA_bp")
              
              ### Obtain estimated thetas and reliability
              res_SA_bp <- get_res(solution = solution_SA_bp, thetas_true = thetas_true, Sigma = Sigma, method = "both")
              res_thetas_SA_bp <- res_SA_bp$res
              names(res_thetas_SA_bp) <- paste0(names(res_thetas_SA_bp), "_SA_bp")

              #### Save estimated thetas
              thetas_est_SA_bp <- res_SA_bp$thetas_est_com
              thetas_SA_bp_list[[cond_name]] <- thetas_est_SA_bp
              if (file.exists("thetas_est/thetas_est_SA_bp.rds")) {
                thetas_SA_bp_list <- readRDS("thetas_est/thetas_est_SA_bp.rds")
              } else {
                thetas_SA_bp_list <- list()
              }
              thetas_SA_bp_list[[cond_name]] <- thetas_est_SA_bp
              saveRDS(thetas_SA_bp_list, file = "thetas_est/thetas_est_SA_bp.rds")
              
              all_marg_reli_SA_bp <- res_SA_bp$all_marg_reli
              all_marg_reli_SA_bp_list[[cond_name]] <- all_marg_reli_SA_bp
              if (file.exists("all_marg_reli/all_marg_reli_SA_bp.rds")) {
                all_marg_reli_SA_bp_list <- readRDS("all_marg_reli/all_marg_reli_SA_bp.rds")
              } else {
                all_marg_reli_SA_bp_list <- list()
              }
              all_marg_reli_SA_bp_list[[cond_name]] <- all_marg_reli_SA_bp
              saveRDS(all_marg_reli_SA_bp_list, file = "all_marg_reli/all_marg_reli_SA_bp.rds")
              
              
              ## SA with optimization of a parameter
              ### Optimization
              start_SA_param <- Sys.time()
              result_SA_param <- optimization(data_item_chars = bank_i, algorithm = "SA_with_param", use_SD = use_SD, cutoff_mixed = cutoff_mixed, cutoff_equal = cutoff_equal, n_blocks_total = n_blocks_total, n_blocks_hetero = n_blocks_hetero)
              end_SA_param <- Sys.time()
              duration_SA_param <- paste(round(as.numeric(end_SA_param-start_SA_param),2), attr(end_SA_param-start_SA_param, "units"))
              
              solution_SA_param <- result_SA_param$solution
              comb_pol_SA_param_hom <- result_SA_param$num_polarity[1]
              names(comb_pol_SA_param_hom) <- paste0(names(comb_pol_SA_param_hom), "_SA_param")
              comb_pol_SA_param_het <- result_SA_param$num_polarity[2]
              names(comb_pol_SA_param_het) <- paste0(names(comb_pol_SA_param_het), "_SA_param")
              comb_trait_SA_param <- result_SA_param$combinations
              names(comb_trait_SA_param) <- paste0(names(comb_trait_SA_param), "_SA_param")
              count_traits_SA_param <- result_SA_param$count_traits
              names(count_traits_SA_param) <- paste0(names(count_traits_SA_param), "_SA_param")
              iter_num_SA_param <- result_SA_param$iter_count
              names(iter_num_SA_param) <- paste0("iter_num", "_SA_param")
              
              ### Obtain estimated thetas and reliability
              res_SA_param <- get_res(solution = solution_SA_param, thetas_true = thetas_true, Sigma = Sigma, method = "both")
              res_thetas_SA_param <- res_SA_param$res
              names(res_thetas_SA_param) <- paste0(names(res_thetas_SA_param), "_SA_param")
              
              #### Save estimated thetas
              thetas_est_SA_param <- res_SA_param$thetas_est_com
              thetas_SA_param_list[[cond_name]] <- thetas_est_SA_param
              if (file.exists("thetas_est/thetas_est_SA_param.rds")) {
                thetas_SA_param_list <- readRDS("thetas_est/thetas_est_SA_param.rds")
              } else {
                thetas_SA_param_list <- list()
              }
              thetas_SA_param_list[[cond_name]] <- thetas_est_SA_param
              saveRDS(thetas_SA_param_list, file = "thetas_est/thetas_est_SA_param.rds")
              
              all_marg_reli_SA_param <- res_SA_param$all_marg_reli
              all_marg_reli_SA_param_list[[cond_name]] <- all_marg_reli_SA_param
              if (file.exists("all_marg_reli/all_marg_reli_SA_param.rds")) {
                all_marg_reli_SA_param_list <- readRDS("all_marg_reli/all_marg_reli_SA_param.rds")
              } else {
                all_marg_reli_SA_param_list <- list()
              }
              all_marg_reli_SA_param_list[[cond_name]] <- all_marg_reli_SA_param
              saveRDS(all_marg_reli_SA_param_list, file = "all_marg_reli/all_marg_reli_SA_param.rds")
              
              
              ## BF
              ### Optimization
              start_BF <- Sys.time()
              result_BF <- optimization(data_item_chars = bank_i, algorithm = "BF", use_SD = use_SD, n_blocks_total = n_blocks_total, n_blocks_hetero = n_blocks_hetero, Sigma = Sigma, BF_iter = 100, BF_relia_method = "marginal", thetas_true = thetas_true, conds = conds)
              end_BF <- Sys.time()
              duration_BF <- paste(round(as.numeric(end_BF-start_BF),2), attr(end_BF-start_BF, "units"))

              solution_BF <- result_BF$solution
              comb_pol_BF_hom <- result_BF$num_polarity[1]
              names(comb_pol_BF_hom) <- paste0(names(comb_pol_BF_hom), "_BF")
              comb_pol_BF_het <- result_BF$num_polarity[2]
              names(comb_pol_BF_het) <- paste0(names(comb_pol_BF_het), "_BF")
              comb_trait_BF <- result_BF$combinations
              names(comb_trait_BF) <- paste0(names(comb_trait_BF), "_BF")
              count_traits_BF <- result_BF$count_traits
              names(count_traits_BF) <- paste0(names(count_traits_BF), "_BF")

              ### Obtain estimated thetas and reliability
              res_BF <- get_res(solution = solution_BF, thetas_true = thetas_true, Sigma = Sigma, method = "both") 
              res_thetas_BF <- res_BF$res
              names(res_thetas_BF) <- paste0(names(res_thetas_BF), "_BF")

              #### Save estimated thetas
              thetas_est_BF <- res_BF$thetas_est_com
              thetas_BF_list[[cond_name]] <- thetas_est_BF
              if (file.exists("thetas_est/thetas_est_BF.rds")) {
                thetas_BF_list <- readRDS("thetas_est/thetas_est_BF.rds")
              } else {
                thetas_BF_list <- list()
              }
              thetas_BF_list[[cond_name]] <- thetas_est_BF
              saveRDS(thetas_BF_list, file = "thetas_est/thetas_est_BF.rds")
              
              all_marg_reli_BF <- res_BF$all_marg_reli
              all_marg_reli_BF_list[[cond_name]] <- all_marg_reli_BF
              if (file.exists("all_marg_reli/all_marg_reli_BF.rds")) {
                all_marg_reli_BF_list <- readRDS("all_marg_reli/all_marg_reli_BF.rds")
              } else {
                all_marg_reli_BF_list <- list()
              }
              all_marg_reli_BF_list[[cond_name]] <- all_marg_reli_BF
              saveRDS(all_marg_reli_BF_list, file = "all_marg_reli/all_marg_reli_BF.rds")
              
              
              
              end_cond <- Sys.time()
              duration_cond <- end_cond - start_cond
              names(duration_cond) <- paste0("duration_cond")
              
              
              results <- c(
                comb_pol_GA_hom,
                comb_pol_GA_het,
                comb_trait_GA,
                count_traits_GA,
                comb_pol_SA_bp_hom,
                comb_pol_SA_bp_het,
                comb_trait_SA_bp,
                count_traits_SA_bp,
                comb_pol_SA_param_hom,
                comb_pol_SA_param_het,
                comb_trait_SA_param,
                count_traits_SA_param,
                comb_pol_BF_hom,
                comb_pol_BF_het,
                comb_trait_BF,
                count_traits_BF,
                iter_num_SA_bp, 
                iter_num_SA_param,
                duration_GA = duration_GA,
                duration_SA_bp = duration_SA_bp,
                duration_SA_param = duration_SA_param,
                duration_BF = duration_BF,
                res_thetas_GA,
                res_thetas_SA_bp,
                res_thetas_SA_param,
                res_thetas_BF,
                duration_cond
              )
              
              write.table(x = c(conds, results),
                          file = "sim_res.txt", append = TRUE,
                          quote = FALSE, sep = "\t", col.names = !file.exists("sim_res.txt"),
                          row.names = FALSE)
              
              all_solutions[[paste0("cond_", i)]] <- list(
                conds = conds,
                solution_GA = solution_GA,
                solution_SA_bp = solution_SA_bp,
                solution_SA_param = solution_SA_param,
                solution_BF = solution_BF)
              saveRDS(all_solutions, file = solutions_file)
              
            })
        }

end <- Sys.time()
end-start

#save(list = ls(), file = "sim_res_obj.RData") 
