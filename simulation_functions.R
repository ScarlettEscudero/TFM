source("lib/selectBlocks.R")
source("lib/NHBSAFunctions.R")
source("lib/NHBSA.R")
source("lib/blockAssemblyNHBSA.R")
library(autoFC)
library(mirt)
library(dplyr)
library(tidyr)
library(Rglpk)

inn_diff <- function(vec, cutoff = cutoff_equal) {
  vec <- unlist(vec)
  if (max(vec) - min(vec) > cutoff) {
    return(1)
  }
  else return(0)
}

range_m <- function(x) {
  return(max(x) - min(x))
}

create_Cmatrix <- function(data_item_chars, use_SD, homo, cutoff_mixed = cutoff_mixed, cutoff_equal = cutoff_equal) {
  J <- nrow(data_item_chars)
  Cmatrix <- matrix(1, nrow = J, ncol = J)
  
  if (use_SD == TRUE){
    if (homo == TRUE) {
      for (i in 1:J) {
        for (j in 1:J) {
          # Unidimensional restriction 
          if (data_item_chars[i, "trait"] == data_item_chars[j, "trait"]) {
            Cmatrix[i, j] <- 0  
          }
          sd_diff <- abs(data_item_chars[i, "SD_rating"] - data_item_chars[j, "SD_rating"])
          polarity_i <- data_item_chars[i, "pol"]
          polarity_j <- data_item_chars[j, "pol"]
          
          if (polarity_i == polarity_j) { # Same polarity (1 1 o -1 -1)
            if (sd_diff >= cutoff_equal) {
              Cmatrix[i, j] <- 0
            }
          } else if (polarity_i != polarity_j) {  # Different polarity (-1 1 o 1 -1)
              Cmatrix[i, j] <- 0
          }
        }
      }
    } else { # homo == FALSE
      for (i in 1:J) {
        for (j in 1:J) {
          # Unidimensional restriction 
          if (data_item_chars[i, "trait"] == data_item_chars[j, "trait"]) {
            Cmatrix[i, j] <- 0  
          }
          sd_diff <- abs(data_item_chars[i, "SD_rating"] - data_item_chars[j, "SD_rating"])
          polarity_i <- data_item_chars[i, "pol"]
          polarity_j <- data_item_chars[j, "pol"]
          
          if (polarity_i == polarity_j) { # Same polarity (1 1 o -1 -1)
            if (sd_diff >= cutoff_equal) {
              Cmatrix[i, j] <- 0
            }
          } else if (polarity_i != polarity_j) {  # Different polarity (-1 1 o 1 -1)
            if (sd_diff >= cutoff_mixed) {
              Cmatrix[i, j] <- 0
            }
          }
        }
      }
    }
  } else { # use_SD == FALSE
    if (homo == TRUE) {
      for (i in 1:J) {
        for (j in 1:J) {
          # Unidimensional restriction
          if (data_item_chars[i, "trait"] == data_item_chars[j, "trait"]) {
            Cmatrix[i, j] <- 0
          }
          
          polarity_i <- data_item_chars[i, "pol"]
          polarity_j <- data_item_chars[j, "pol"]
          
          if (polarity_i != polarity_j) {  # Different polarity (-1 1 o 1 -1)
            Cmatrix[i, j] <- 0
          }
        }
      }
    } else { # homo == FALSE
      for (i in 1:J) {
        for (j in 1:J) {
          # Unidimensional restriction
          if (data_item_chars[i, "trait"] == data_item_chars[j, "trait"]) {
            Cmatrix[i, j] <- 0
          }
        }
      }
    }
  }
  
  return(Cmatrix)
}

generate_data <- function(conds, J){

  J <- as.numeric(J)
  K <- 5
  J.k <- J / K
  item_pol_mixed <- rep(rep(c(-1, 1), each = J.k / 2), K)
  item_pol_pos <- rep(1, J)
  
  param_relations <- as.numeric(conds$degree_het)  
  
    # Build bank with degree_relation_sd_a_d 
    if (param_relations == 1){ # Condition all positive items with moderate correlation (r = 0.5) between SD-a (and moderate correlation between a-d)

      mean_Scale <- 1.5
      mean_Threshold <- -0.5
      mean_SD <- 4
      var_Scale <- 0.25  # SD = 0.5
      var_Threshold <- 0.64 # SD = 0.8 
      var_SD <- 0.25 # SD = 0.5
      
      cor_Scale_Threshold <- 0.3  # Correlation between Scale and Threshold
      cor_Scale_SD <- 0.5  # Correlation between Scale and SD
      
      covariance_Scale_Threshold <- cor_Scale_Threshold * sqrt(var_Scale) * sqrt(var_Threshold)
      covariance_Scale_SD <- cor_Scale_SD * sqrt(var_Scale) * sqrt(var_SD)
      
      cov_matrix <- matrix(
        c(var_Scale, covariance_Scale_Threshold, covariance_Scale_SD, covariance_Scale_Threshold, var_Threshold, 0, covariance_Scale_SD, 0, var_SD), 
        nrow = 3, byrow = TRUE)
      
      normal_multivariate <- tmvtnorm::rtmvnorm(J, mean = c(mean_Scale, mean_Threshold, mean_SD), 
                                                sigma = cov_matrix, 
                                                lower = c(0, -Inf, 1), 
                                                upper = c(Inf, Inf, 5))
      Scale <- normal_multivariate[,1]
      Threshold <- normal_multivariate[,2]
      SD <- normal_multivariate[,3]
      
      bank <- data.frame(
        item_num = seq(1:J),
        trait = rep(1:K, each = J.k),
        pol = item_pol_pos,
        Scale = Scale,
        Threshold = Threshold,
        SD_rating = SD
      )
      
    } else if (param_relations == 2){ # Condition all positive items with high correlation (r = 0.8) between SD-a (and moderate correlation between a-d)

      mean_Scale <- 1.5
      mean_Threshold <- -0.5
      mean_SD <- 4
      var_Scale <- 0.25  # SD = 0.5
      var_Threshold <- 0.64 # SD = 0.8 
      var_SD <- 0.25     # SD = 0.5
      
      cor_Scale_Threshold <- 0.3  # Correlation between Scale and Threshold
      cor_Scale_SD <- 0.8  # Correlation between Scale and SD
      
      covariance_Scale_Threshold <- cor_Scale_Threshold * sqrt(var_Scale) * sqrt(var_Threshold)
      covariance_Scale_SD <- cor_Scale_SD * sqrt(var_Scale) * sqrt(var_SD)
      
      cov_matrix <- matrix(
        c(var_Scale, covariance_Scale_Threshold, covariance_Scale_SD, covariance_Scale_Threshold, var_Threshold, 0, covariance_Scale_SD, 0, var_SD), 
        nrow = 3, byrow = TRUE)
      
      normal_multivariate <- tmvtnorm::rtmvnorm(J, mean = c(mean_Scale, mean_Threshold, mean_SD), 
                                                sigma = cov_matrix, 
                                                lower = c(0, -Inf, 1), 
                                                upper = c(Inf, Inf, 5))
      Scale <- normal_multivariate[,1]
      Threshold <- normal_multivariate[,2]
      SD <- normal_multivariate[,3]
      
      bank <- data.frame(
        item_num = seq(1:J),
        trait = rep(1:K, each = J.k),
        pol = item_pol_pos,
        Scale = Scale,
        Threshold = Threshold,
        SD_rating = SD
      )
      
    } else if (param_relations == 3 | 4){ # Condition positive and negative items with high correlation (r = 0.9) between SD-a (and moderate correlation between a-d)

      mean_Scale_pos <- 1.5
      mean_Scale_neg <- -1.5
      mean_Threshold_pos <- -0.5
      mean_Threshold_neg <- -1
      mean_SD_pos <- 4
      mean_SD_neg <- 2
      
      var_Scale <- 0.25  # SD = 0.5
      var_Threshold_pos <- 0.64 # SD = 0.8 
      var_Threshold_neg <- 0.64 # SD = 0.8
      var_SD <- 0.25     # SD = 0.5
      
      cor_Scale_Threshold_pos <- 0.3  # Internal correlation between Scale and Threshold for positive and negative items
      cor_Scale_Threshold_neg <- -0.3  # Internal correlation between negative Scale and Threshold 
      cor_Scale_SD_pos <- 0.2  # Internal correlation between positive Scale and SD 
      cor_Scale_SD_neg <- 0.2  # Internal correlation between negative Scale and SD 
      
      covariance_Scale_Threshold_pos <- cor_Scale_Threshold_pos * sqrt(var_Scale) * sqrt(var_Threshold_pos)
      covariance_Scale_Threshold_neg <- cor_Scale_Threshold_neg * sqrt(var_Scale) * sqrt(var_Threshold_neg)
      covariance_Scale_SD_pos <- cor_Scale_SD_pos * sqrt(var_Scale) * sqrt(var_SD)
      covariance_Scale_SD_neg <- cor_Scale_SD_neg * sqrt(var_Scale) * sqrt(var_SD)
      
      cov_matrix_pos <- matrix(
        c(var_Scale, covariance_Scale_Threshold_pos, covariance_Scale_SD_pos, covariance_Scale_Threshold_pos, var_Threshold_pos, 0, covariance_Scale_SD_pos, 0, var_SD),
        nrow = 3, byrow = TRUE)
      
      cov_matrix_neg <- matrix(
        c(var_Scale, covariance_Scale_Threshold_neg, covariance_Scale_SD_neg, covariance_Scale_Threshold_neg, var_Threshold_neg, 0, covariance_Scale_SD_neg, 0, var_SD),
        nrow = 3, byrow = TRUE)
      
      pos_items <- which(item_pol_mixed == 1)
      mixed_gaussian_pos <- tmvtnorm::rtmvnorm(length(pos_items), 
                                               mean = c(mean_Scale_pos, mean_Threshold_pos, mean_SD_pos), 
                                               sigma = cov_matrix_pos, 
                                               lower = c(-Inf, -3, 1), 
                                               upper = c(Inf, Inf, 5))
      
      neg_items <- which(item_pol_mixed == -1)
      mixed_gaussian_neg <- tmvtnorm::rtmvnorm(length(neg_items), 
                                                     mean = c(mean_Scale_neg, mean_Threshold_neg, mean_SD_neg), 
                                                     sigma = cov_matrix_neg, 
                                                     lower = c(-Inf, -3, 1), 
                                                     upper = c(Inf, Inf, 5))
      
      Scale_pos <- mixed_gaussian_pos[, 1]
      Scale_neg <- mixed_gaussian_neg[, 1]
      Threshold_pos <- mixed_gaussian_pos[, 2]
      Threshold_neg <- mixed_gaussian_neg[, 2]
      SD_pos <- mixed_gaussian_pos[, 3]
      SD_neg <- mixed_gaussian_neg[, 3]
      
      # Mix indices without replacement
      shuffled_pos <- sample(seq_along(pos_items))
      shuffled_neg <- sample(seq_along(neg_items))
      
      # Assign random values maintaining row correspondence
      Scale <- numeric(length(item_pol_mixed))
      Threshold <- numeric(length(item_pol_mixed))
      SD <- numeric(length(item_pol_mixed))
      
      Scale[pos_items] <- Scale_pos[shuffled_pos]
      Threshold[pos_items] <- Threshold_pos[shuffled_pos]
      SD[pos_items] <- SD_pos[shuffled_pos]
      
      Scale[neg_items] <- Scale_neg[shuffled_neg]
      Threshold[neg_items] <- Threshold_neg[shuffled_neg]
      SD[neg_items] <- SD_neg[shuffled_neg]
      

      bank <- data.frame(
        item_num = seq(1:J),
        trait = rep(1:K, each = J.k),
        pol = item_pol_mixed,
        Scale = Scale,
        Threshold = Threshold,
        SD_rating = SD
      )
 
    } 
  return(bank)
}

optimization <- function(data_item_chars, algorithm = "GA", use_SD = TRUE, cutoff_mixed = 1.25, cutoff_equal = 0.5, Sigma, n_blocks_total, n_blocks_hetero = 0, nCores = 5, balance_dim_comb = FALSE, which.combs, type_theta, BF_iter = 100, BF_relia_method = "marginal", thetas_true, conds = conds){
  
  data_item_chars <- as.data.frame(data_item_chars)
  n_blocks_total <- as.numeric(n_blocks_total)
  n_blocks_hetero <- as.numeric(n_blocks_hetero)
  n_blocks_homo <- as.numeric(n_blocks_total - n_blocks_hetero)
  cutoff_equal <- as.numeric(cutoff_equal)
  cutoff_mixed <- as.numeric(cutoff_mixed)
  unique_traits <- unique(data_item_chars$trait)
  nCores <- as.numeric(nCores)
  if (n_blocks_hetero == 0) {homo <- TRUE} else {homo <- FALSE}
  
  K <- length(unique(data_item_chars$trait)) # number of traits
  J <- nrow(data_item_chars) # item pool size
  a <- data_item_chars$Scale # discrimination parameters
  d <- data_item_chars$Threshold # threshold parameters
  
  
  A.spec <- matrix(0, J, K)
  traits <- data_item_chars$trait
  all_combinations <- as.data.frame(t(combn(unique_traits, 2))) %>%
    rename(trait.1 = V1, trait.2 = V2) %>% 
    mutate(trait_comb = paste(pmin(trait.1, trait.2), pmax(trait.1, trait.2), sep = "-")) %>%
    distinct(trait_comb)  
  
  for (i in 1:J) {
    A.spec[i, traits[i]] <- 1
  }
  
  A <- a*A.spec 
  
  colnames(A) <- paste0("a", 1:K)

  if (algorithm == "GA") {
    
    item_data <- cbind(id = 1:J, trait = traits, pol = data_item_chars$pol, A, d1 = d)
    
    Cmatrix <- create_Cmatrix(data_item_chars, use_SD = use_SD, homo = homo, cutoff_mixed = cutoff_mixed, cutoff_equal = cutoff_equal)
    
    n_combin <- nrow(t(combn(K,2)))
    
    if(balance_dim_comb == TRUE & is.null(which.combs)){
      which.combs <- list(dir_cons = rep((n_blocks_total-n_blocks_hetero)/n_combin,n_combin),
                          inv_cons = rep(n_blocks_hetero/n_combin,n_combin))
    } else if (balance_dim_comb == TRUE & !is.null(which.combs)){
      which.combs <- which.combs
    } else if (balance_dim_comb == FALSE & is.null(which.combs)) {
      which.combs <- list(dir_cons = n_blocks_total - n_blocks_hetero,
                          inv_cons = n_blocks_hetero)
    } else {
      which.combs <- list(dir_cons = n_blocks_total - n_blocks_hetero,
                          inv_cons = n_blocks_hetero)
    }
    
    result_GA <- blockAssemblyNHBSA(criterion = "AvgMargReli", simpMethod = "Boundary", nBlocks = n_blocks_total, nHetero = n_blocks_hetero, itemData = item_data, Sigma = Sigma, Cmatrix = Cmatrix, TIRT = FALSE, nCores = nCores, balance_dim_comb = FALSE, which.combs = which.combs, theta = type_theta)
    result_chars <- as.data.frame(result_GA$blocks)
    colnames(result_chars)[colnames(result_chars)=="id.1"] <- "item_num.1"
    colnames(result_chars)[colnames(result_chars)=="id.2"] <- "item_num.2"
    
    result_chars <- result_chars %>%
      dplyr::left_join(data_item_chars %>% dplyr::select(item_num, SD_rating), 
                       by = c("item_num.1" = "item_num")) %>%
      dplyr::rename(SD_rating.1 = SD_rating)
    result_chars <- result_chars %>%
      dplyr::left_join(data_item_chars %>% dplyr::select(item_num, SD_rating), 
                       by = c("item_num.2" = "item_num")) %>%
      dplyr::rename(SD_rating.2 = SD_rating)
    
    result_blocks <- result_chars[,c("item_num.1","item_num.2")]
    
    iter_count <- 0

  } else if (algorithm == "SA_blueprint") {
    if (use_SD == T) {
      if (homo == T) {
     
        iter_count <- 0
        max_iter <- 200
        repeat {
          result_SA <- autoFC::sa_pairing_generalized(block = make_random_block(nrow(data_item_chars), n_blocks_total*2, 2),
                                                      item_chars = data_item_chars[,c("trait","SD_rating","pol")],
                                                      r = 0.999,
                                                      FUN = c("facfun","inn_diff","var"),
                                                      weights = c(10000,-1000000,-100000))
          
          proposed_items <- data_item_chars[c(t(result_SA$block_final)),]
          
          proposed_items_blocks <- data.frame(
            item_num.1 = proposed_items$item_num[seq(1, nrow(proposed_items)-1, by=2)],
            item_num.2 = proposed_items$item_num[seq(2, nrow(proposed_items), by=2)],
            trait.1 = proposed_items$trait[seq(1, nrow(proposed_items)-1, by=2)],
            trait.2 = proposed_items$trait[seq(2, nrow(proposed_items), by=2)],
            pol.1 = proposed_items$pol[seq(1, nrow(proposed_items)-1, by=2)],
            pol.2 = proposed_items$pol[seq(2, nrow(proposed_items), by=2)],
            SD_rating.1 = proposed_items$SD_rating[seq(1, nrow(proposed_items)-1, by=2)],
            SD_rating.2 = proposed_items$SD_rating[seq(2, nrow(proposed_items), by=2)],
            Scale.1 = proposed_items$Scale[seq(1, nrow(proposed_items)-1, by=2)],
            Scale.2 = proposed_items$Scale[seq(2, nrow(proposed_items), by=2)],
            d1.1 = proposed_items$Threshold[seq(1, nrow(proposed_items)-1, by=2)],
            d1.2 = proposed_items$Threshold[seq(2, nrow(proposed_items), by=2)])
          
          num_heteropolares_int <- proposed_items_blocks %>%
            filter((pol.1 == 1 & pol.2 == -1) |
                     (pol.1 == -1 & pol.2 == 1)) %>% nrow(); num_heteropolares_int
          
          iter_count <- iter_count + 1
          if (num_heteropolares_int == 0 || iter_count >= max_iter) break
        }

        result_blocks <- as.data.frame(result_SA$block_final)
        colnames(result_blocks) <- c("item_num.1","item_num.2")
        
        result_chars <- result_blocks %>%
          left_join(data_item_chars %>% rename_with(~ paste0(.x, ".1"), -item_num), by = c("item_num.1" = "item_num")) %>%
          left_join(data_item_chars %>% rename_with(~ paste0(.x, ".2"), -item_num), by = c("item_num.2" = "item_num"))
        
      } else { #homo == F

        pairs <- expand.grid(item.1 = data_item_chars$item_num, item.2 = data_item_chars$item_num)
        pairs <- pairs[pairs$item.1 < pairs$item.2, ]

        pairs <- pairs %>%
          left_join(data_item_chars, by = c("item.1" = "item_num")) %>%
          rename(trait.1 = trait, pol.1 = pol, SD_rating.1 = SD_rating, Scale.1 = Scale, d1.1 = Threshold) %>%
          left_join(data_item_chars, by = c("item.2" = "item_num")) %>%
          rename(trait.2 = trait, pol.2 = pol, SD_rating.2 = SD_rating, Scale.2 = Scale, d1.2 = Threshold)

        valid_pairs <- pairs %>%
          mutate(
            homopolar = (pol.1 == pol.2),
            SD_restriction = abs(SD_rating.1 - SD_rating.2) <= cutoff_mixed) %>%
          filter(SD_restriction, trait.1 != trait.2, !homopolar)  # Only heteropolar y multidimensional

        num_pairs <- nrow(valid_pairs)

        unique_items <- unique(c(valid_pairs$item.1, valid_pairs$item.2)) # Number of unique items (because of the SD restriction, not all 320 items are in valid_pairs)
        length_unique_items <- length(unique_items)

        item_index <- match(unique_items, sort(unique_items))

        mat <- matrix(0, nrow = length_unique_items + 11, ncol = num_pairs) # Constraints for length_unique_items + 11 additional constraints

        for (j in 1:length_unique_items) {
          item_j <- unique_items[j]
          mat[j, ] <- ifelse(valid_pairs$item.1 == item_j | valid_pairs$item.2 == item_j, 1, 0)
        }

        mat[length_unique_items + 1,] <- ifelse(valid_pairs$trait.1 == 1 | valid_pairs$trait.2 == 1, 1, 0) # Restriction row for items with trait 1 (for lower constraint)
        mat[length_unique_items + 2,] <- ifelse(valid_pairs$trait.1 == 2 | valid_pairs$trait.2 == 2, 1, 0) # Restriction row for items with trait 2 (for lower constraint)
        mat[length_unique_items + 3,] <- ifelse(valid_pairs$trait.1 == 3 | valid_pairs$trait.2 == 3, 1, 0) # Restriction row for items with trait 3 (for lower constraint)
        mat[length_unique_items + 4,] <- ifelse(valid_pairs$trait.1 == 4 | valid_pairs$trait.2 == 4, 1, 0) # Restriction row for items with trait 4 (for lower constraint)
        mat[length_unique_items + 5,] <- ifelse(valid_pairs$trait.1 == 5 | valid_pairs$trait.2 == 5, 1, 0) # Restriction row for items with trait 5 (for lower constraint)
        mat[length_unique_items + 6,] <- ifelse(valid_pairs$trait.1 == 1 | valid_pairs$trait.2 == 1, 1, 0) # Restriction row for items with trait 1 (for higher constraint)
        mat[length_unique_items + 7,] <- ifelse(valid_pairs$trait.1 == 2 | valid_pairs$trait.2 == 2, 1, 0) # Restriction row for items with trait 2 (for higher constraint)
        mat[length_unique_items + 8,] <- ifelse(valid_pairs$trait.1 == 3 | valid_pairs$trait.2 == 3, 1, 0) # Restriction row for items with trait 3 (for higher constraint)
        mat[length_unique_items + 9,] <- ifelse(valid_pairs$trait.1 == 4 | valid_pairs$trait.2 == 4, 1, 0) # Restriction row for items with trait 4 (for higher constraint)
        mat[length_unique_items + 10,] <- ifelse(valid_pairs$trait.1 == 5 | valid_pairs$trait.2 == 5, 1, 0) # Restriction row for items with trait 5 (for higher constraint)

        mat[length_unique_items + 11,] <- ifelse((valid_pairs$pol.1 == 1 & valid_pairs$pol.2 == -1) | (valid_pairs$pol.1 == -1 & valid_pairs$pol.2 == 1), 1, 0)

        dir <- c(rep("<=", length_unique_items), rep(">=",5), rep("<=",5), "==")
        rhs <- c(rep(1, length_unique_items), rep(((n_blocks_hetero/5)*2)*0.8,5), rep(((n_blocks_hetero/5)*2)*1.2,5), n_blocks_hetero)

        obj <- runif(num_pairs)

        sol <- Rglpk_solve_LP(obj, mat, dir, rhs, max = TRUE, types = rep("B", num_pairs))

        selected_blocks <- valid_pairs[sol$solution == 1, ]
        colnames(selected_blocks)[1:2] <- c("item_num.1", "item_num.2")

        selected_blocks <- selected_blocks %>%
          select(item_num.1, item_num.2, trait.1, trait.2,
                 pol.1, pol.2, SD_rating.1, SD_rating.2,
                 Scale.1, Scale.2,
                 d1.1, d1.2)

        picked_items <- selected_blocks %>%
          pivot_longer(
            cols = starts_with("item_num"),
            names_to = "group",
            values_to = "item_num"
          ) %>%
          mutate(group = as.numeric(gsub("item_num.", "", group))) %>%
          arrange(row_number(), group) %>%
          mutate(
            trait = ifelse(group == 1, trait.1, trait.2),
            pol = ifelse(group == 1, pol.1, pol.2),
            SD_rating = ifelse(group == 1, SD_rating.1, SD_rating.2),
            Scale = ifelse(group == 1, Scale.1, Scale.2),
            Threshold = ifelse(group == 1, d1.1, d1.2)
          ) %>%
          select(item_num, trait, pol, Scale, Threshold, SD_rating)

        item_info_rest <- data_item_chars %>%
          filter(!item_num %in% picked_items$item_num)
        
        iter_count <- 0
        max_iter <- 200
        repeat {
          result_SA <- autoFC::sa_pairing_generalized(block = make_random_block(nrow(item_info_rest), (n_blocks_total-n_blocks_hetero)*2, 2),
                                                      item_chars = item_info_rest[,c("trait","SD_rating","pol")],
                                                      r = 0.999,
                                                      FUN = c("facfun","inn_diff","var"),
                                                      weights = c(10000,-1000000,-100000)) 
          
          rest_items <- item_info_rest[c(t(result_SA$block_final)),]
          
          rest_items_blocks <- data.frame(
            item_num.1 = rest_items$item_num[seq(1, nrow(rest_items)-1, by=2)],
            item_num.2 = rest_items$item_num[seq(2, nrow(rest_items), by=2)],
            trait.1 = rest_items$trait[seq(1, nrow(rest_items)-1, by=2)],
            trait.2 = rest_items$trait[seq(2, nrow(rest_items), by=2)],
            pol.1 = rest_items$pol[seq(1, nrow(rest_items)-1, by=2)],
            pol.2 = rest_items$pol[seq(2, nrow(rest_items), by=2)],
            SD_rating.1 = rest_items$SD_rating[seq(1, nrow(rest_items)-1, by=2)],
            SD_rating.2 = rest_items$SD_rating[seq(2, nrow(rest_items), by=2)],
            Scale.1 = rest_items$Scale[seq(1, nrow(rest_items)-1, by=2)],
            Scale.2 = rest_items$Scale[seq(2, nrow(rest_items), by=2)],
            d1.1 = rest_items$Threshold[seq(1, nrow(rest_items)-1, by=2)],
            d1.2 = rest_items$Threshold[seq(2, nrow(rest_items), by=2)])
          
          num_heteropolares_int <- rest_items_blocks %>%
            filter((pol.1 == 1 & pol.2 == -1) |
                     (pol.1 == -1 & pol.2 == 1)) %>% nrow()
          
          iter_count <- iter_count + 1
          if (num_heteropolares_int == 0 || iter_count >= max_iter) break
        }

        result_blocks <- as.data.frame(matrix(c(picked_items$item_num, rest_items$item_num), ncol = 2, byrow = TRUE))
        colnames(result_blocks) <- c("item_num.1", "item_num.2")
        
        result_chars <- result_blocks %>%
          left_join(data_item_chars %>% rename_with(~ paste0(.x, ".1"), -item_num), by = c("item_num.1" = "item_num")) %>%
          left_join(data_item_chars %>% rename_with(~ paste0(.x, ".2"), -item_num), by = c("item_num.2" = "item_num"))
        
      }
    } else { # use_SD == F
      if (homo == T) {
        iter_count <- 0
        max_iter <- 200
        repeat {
          result_SA <- autoFC::sa_pairing_generalized(block = make_random_block(nrow(data_item_chars), n_blocks_total*2, 2),
                                                      item_chars = data_item_chars[,c("trait","pol")],
                                                      r = 0.999,
                                                      FUN = c("facfun","var"),
                                                      weights = c(10000,-100000))
          
          proposed_items <- data_item_chars[c(t(result_SA$block_final)),]
          
          proposed_items_blocks <- data.frame(
            item_num.1 = proposed_items$item_num[seq(1, nrow(proposed_items)-1, by=2)],
            item_num.2 = proposed_items$item_num[seq(2, nrow(proposed_items), by=2)],
            trait.1 = proposed_items$trait[seq(1, nrow(proposed_items)-1, by=2)],
            trait.2 = proposed_items$trait[seq(2, nrow(proposed_items), by=2)],
            pol.1 = proposed_items$pol[seq(1, nrow(proposed_items)-1, by=2)],
            pol.2 = proposed_items$pol[seq(2, nrow(proposed_items), by=2)],
            SD_rating.1 = proposed_items$SD_rating[seq(1, nrow(proposed_items)-1, by=2)],
            SD_rating.2 = proposed_items$SD_rating[seq(2, nrow(proposed_items), by=2)],
            Scale.1 = proposed_items$Scale[seq(1, nrow(proposed_items)-1, by=2)],
            Scale.2 = proposed_items$Scale[seq(2, nrow(proposed_items), by=2)],
            d1.1 = proposed_items$Threshold[seq(1, nrow(proposed_items)-1, by=2)],
            d1.2 = proposed_items$Threshold[seq(2, nrow(proposed_items), by=2)])
          
          num_heteropolares_int <- proposed_items_blocks %>%
            filter((pol.1 == 1 & pol.2 == -1) |
                     (pol.1 == -1 & pol.2 == 1)) %>% nrow()
          
          iter_count <- iter_count + 1
          if (num_heteropolares_int == 0 || iter_count >= max_iter) break
        }
        result_blocks <- as.data.frame(result_SA$block_final)
        colnames(result_blocks) <- c("item_num.1","item_num.2")
        
        result_chars <- result_blocks %>%
          left_join(data_item_chars %>% rename_with(~ paste0(.x, ".1"), -item_num), by = c("item_num.1" = "item_num")) %>%
          left_join(data_item_chars %>% rename_with(~ paste0(.x, ".2"), -item_num), by = c("item_num.2" = "item_num"))

      } else { #homo == F

        pairs <- expand.grid(item.1 = data_item_chars$item_num, item.2 = data_item_chars$item_num)
        pairs <- pairs[pairs$item.1 < pairs$item.2, ]
        
        pairs <- pairs %>%
          left_join(data_item_chars, by = c("item.1" = "item_num")) %>%
          rename(trait.1 = trait, pol.1 = pol, SD_rating.1 = SD_rating, Scale.1 = Scale, d1.1 = Threshold) %>%
          left_join(data_item_chars, by = c("item.2" = "item_num")) %>%
          rename(trait.2 = trait, pol.2 = pol, SD_rating.2 = SD_rating, Scale.2 = Scale, d1.2 = Threshold)
        
        valid_pairs <- pairs %>%
          mutate(
            homopolar = (pol.1 == pol.2)) %>%
          filter(trait.1 != trait.2, !homopolar)  
        
        num_pairs <- nrow(valid_pairs)
        
        unique_items <- unique(c(valid_pairs$item.1, valid_pairs$item.2)) # Number of unique items (because of the SD restriction, not all 320 items are in valid_pairs)
        length_unique_items <- length(unique_items)
        
        item_index <- match(unique_items, sort(unique_items))
        
        mat <- matrix(0, nrow = length_unique_items + 11, ncol = num_pairs) # Constraints for length_unique_items + 11 additional constraints
        
        for (j in 1:length_unique_items) {
          item_j <- unique_items[j]
          mat[j, ] <- ifelse(valid_pairs$item.1 == item_j | valid_pairs$item.2 == item_j, 1, 0)
        }
        
        mat[length_unique_items + 1,] <- ifelse(valid_pairs$trait.1 == 1 | valid_pairs$trait.2 == 1, 1, 0) # Restriction row for items with trait 1 (for lower constraint)
        mat[length_unique_items + 2,] <- ifelse(valid_pairs$trait.1 == 2 | valid_pairs$trait.2 == 2, 1, 0) # Restriction row for items with trait 2 (for lower constraint)
        mat[length_unique_items + 3,] <- ifelse(valid_pairs$trait.1 == 3 | valid_pairs$trait.2 == 3, 1, 0) # Restriction row for items with trait 3 (for lower constraint)
        mat[length_unique_items + 4,] <- ifelse(valid_pairs$trait.1 == 4 | valid_pairs$trait.2 == 4, 1, 0) # Restriction row for items with trait 4 (for lower constraint)
        mat[length_unique_items + 5,] <- ifelse(valid_pairs$trait.1 == 5 | valid_pairs$trait.2 == 5, 1, 0) # Restriction row for items with trait 5 (for lower constraint)
        mat[length_unique_items + 6,] <- ifelse(valid_pairs$trait.1 == 1 | valid_pairs$trait.2 == 1, 1, 0) # Restriction row for items with trait 1 (for higher constraint)
        mat[length_unique_items + 7,] <- ifelse(valid_pairs$trait.1 == 2 | valid_pairs$trait.2 == 2, 1, 0) # Restriction row for items with trait 2 (for higher constraint)
        mat[length_unique_items + 8,] <- ifelse(valid_pairs$trait.1 == 3 | valid_pairs$trait.2 == 3, 1, 0) # Restriction row for items with trait 3 (for higher constraint)
        mat[length_unique_items + 9,] <- ifelse(valid_pairs$trait.1 == 4 | valid_pairs$trait.2 == 4, 1, 0) # Restriction row for items with trait 4 (for higher constraint)
        mat[length_unique_items + 10,] <- ifelse(valid_pairs$trait.1 == 5 | valid_pairs$trait.2 == 5, 1, 0) # Restriction row for items with trait 5 (for higher constraint)
        
        mat[length_unique_items + 11,] <- ifelse((valid_pairs$pol.1 == 1 & valid_pairs$pol.2 == -1) | (valid_pairs$pol.1 == -1 & valid_pairs$pol.2 == 1), 1, 0)
        
        dir <- c(rep("<=", length_unique_items), rep(">=",5), rep("<=",5), "==")
        rhs <- c(rep(1, length_unique_items), rep(((n_blocks_hetero/5)*2)*0.8,5), rep(((n_blocks_hetero/5)*2)*1.2,5), n_blocks_hetero)
        
        obj <- runif(num_pairs)
        
        sol <- Rglpk_solve_LP(obj, mat, dir, rhs, max = TRUE, types = rep("B", num_pairs))
        
        selected_blocks <- valid_pairs[sol$solution == 1, ]
        colnames(selected_blocks)[1:2] <- c("item_num.1", "item_num.2")
        
        selected_blocks <- selected_blocks %>%
          select(item_num.1, item_num.2, trait.1, trait.2,
                 pol.1, pol.2, SD_rating.1, SD_rating.2,
                 Scale.1, Scale.2,
                 d1.1, d1.2)
        
        picked_items <- selected_blocks %>%
          pivot_longer(
            cols = starts_with("item_num"),
            names_to = "group",
            values_to = "item_num"
          ) %>%
          mutate(group = as.numeric(gsub("item_num.", "", group))) %>%
          arrange(row_number(), group) %>%
          mutate(
            trait = ifelse(group == 1, trait.1, trait.2),
            pol = ifelse(group == 1, pol.1, pol.2),
            SD_rating = ifelse(group == 1, SD_rating.1, SD_rating.2),
            Scale = ifelse(group == 1, Scale.1, Scale.2),
            Threshold = ifelse(group == 1, d1.1, d1.2)
          ) %>%
          select(item_num, trait, pol, Scale, Threshold, SD_rating)
        
        item_info_rest <- data_item_chars %>%
          filter(!item_num %in% picked_items$item_num)
        
        iter_count <- 0
        max_iter <- 200
        repeat {
        result_SA <- autoFC::sa_pairing_generalized(block = make_random_block(nrow(item_info_rest), (n_blocks_total-n_blocks_hetero)*2, 2),
                                                    item_chars = item_info_rest[,c("trait","pol")],
                                                    r = 0.999,
                                                    FUN = c("facfun","var"),
                                                    weights = c(10000,-100000)) 
        
        rest_items <- item_info_rest[c(t(result_SA$block_final)),]
        rest_items_blocks <- data.frame(
          item_num.1 = rest_items$item_num[seq(1, nrow(rest_items)-1, by=2)],
          item_num.2 = rest_items$item_num[seq(2, nrow(rest_items), by=2)],
          trait.1 = rest_items$trait[seq(1, nrow(rest_items)-1, by=2)],
          trait.2 = rest_items$trait[seq(2, nrow(rest_items), by=2)],
          pol.1 = rest_items$pol[seq(1, nrow(rest_items)-1, by=2)],
          pol.2 = rest_items$pol[seq(2, nrow(rest_items), by=2)],
          SD_rating.1 = rest_items$SD_rating[seq(1, nrow(rest_items)-1, by=2)],
          SD_rating.2 = rest_items$SD_rating[seq(2, nrow(rest_items), by=2)],
          Scale.1 = rest_items$Scale[seq(1, nrow(rest_items)-1, by=2)],
          Scale.2 = rest_items$Scale[seq(2, nrow(rest_items), by=2)],
          d1.1 = rest_items$Threshold[seq(1, nrow(rest_items)-1, by=2)],
          d1.2 = rest_items$Threshold[seq(2, nrow(rest_items), by=2)])
        
        num_heteropolares_int <- rest_items_blocks %>%
          filter((pol.1 == 1 & pol.2 == -1) |
                   (pol.1 == -1 & pol.2 == 1)) %>% nrow()
        
        iter_count <- iter_count + 1
        if (num_heteropolares_int == 0 || iter_count >= max_iter) break
        }
        
        result_blocks <- as.data.frame(matrix(c(picked_items$item_num, rest_items$item_num), ncol = 2, byrow = TRUE))
        colnames(result_blocks) <- c("item_num.1", "item_num.2")
        
        result_chars <- result_blocks %>%
          left_join(data_item_chars %>% rename_with(~ paste0(.x, ".1"), -item_num), by = c("item_num.1" = "item_num")) %>%
          left_join(data_item_chars %>% rename_with(~ paste0(.x, ".2"), -item_num), by = c("item_num.2" = "item_num"))
      }
    }
    
    item_chars_1 <- data_item_chars %>%
      rename_with(~ paste0(.x, ".1"), -item_num)
    item_chars_2 <- data_item_chars %>%
      rename_with(~ paste0(.x, ".2"), -item_num)
    
    result_chars <- result_blocks %>%
      left_join(item_chars_1, by = c("item_num.1" = "item_num")) %>%
      left_join(item_chars_2, by = c("item_num.2" = "item_num"))
    
    result_chars <- result_chars %>%
      mutate(
        d1.1 = Threshold.1,
        d1.2 = Threshold.2) %>%
      mutate(
        a1.1 = ifelse(trait.1 == 1, Scale.1, 0),
        a2.1 = ifelse(trait.1 == 2, Scale.1, 0),
        a3.1 = ifelse(trait.1 == 3, Scale.1, 0),
        a4.1 = ifelse(trait.1 == 4, Scale.1, 0),
        a5.1 = ifelse(trait.1 == 5, Scale.1, 0),
        a1.2 = ifelse(trait.2 == 1, Scale.2, 0),
        a2.2 = ifelse(trait.2 == 2, Scale.2, 0),
        a3.2 = ifelse(trait.2 == 3, Scale.2, 0),
        a4.2 = ifelse(trait.2 == 4, Scale.2, 0),
        a5.2 = ifelse(trait.2 == 5, Scale.2, 0))
    
    result_chars <- result_chars %>%
      select(item_num.1, item_num.2, 
             trait.1, trait.2,
             pol.1, pol.2, SD_rating.1, SD_rating.2,
             d1.1, d1.2,
             a1.1, a2.1, a3.1, a4.1, a5.1,
             a1.2, a2.2, a3.2, a4.2, a5.2)
    
  } else if (algorithm == "SA_with_param") {
    if (use_SD == T) {
      if (homo == T) {
        iter_count <- 0
        max_iter <- 200
        repeat {
          result_SA <- autoFC::sa_pairing_generalized(block = make_random_block(nrow(data_item_chars), n_blocks_total*2, 2),
                                                      item_chars = data_item_chars[,c("trait","SD_rating","pol","Scale")],
                                                      r = 0.999,
                                                      FUN = c("facfun","inn_diff","var", "var"),
                                                      weights = c(10000,-1000000,-100000,10000))
          
          proposed_items <- data_item_chars[c(t(result_SA$block_final)),]
          
          proposed_items_blocks <- data.frame(
            item_num.1 = proposed_items$item_num[seq(1, nrow(proposed_items)-1, by=2)],
            item_num.2 = proposed_items$item_num[seq(2, nrow(proposed_items), by=2)],
            trait.1 = proposed_items$trait[seq(1, nrow(proposed_items)-1, by=2)],
            trait.2 = proposed_items$trait[seq(2, nrow(proposed_items), by=2)],
            pol.1 = proposed_items$pol[seq(1, nrow(proposed_items)-1, by=2)],
            pol.2 = proposed_items$pol[seq(2, nrow(proposed_items), by=2)],
            SD_rating.1 = proposed_items$SD_rating[seq(1, nrow(proposed_items)-1, by=2)],
            SD_rating.2 = proposed_items$SD_rating[seq(2, nrow(proposed_items), by=2)],
            Scale.1 = proposed_items$Scale[seq(1, nrow(proposed_items)-1, by=2)],
            Scale.2 = proposed_items$Scale[seq(2, nrow(proposed_items), by=2)],
            d1.1 = proposed_items$Threshold[seq(1, nrow(proposed_items)-1, by=2)],
            d1.2 = proposed_items$Threshold[seq(2, nrow(proposed_items), by=2)])
          
          num_heteropolares_int <- proposed_items_blocks %>%
            filter((pol.1 == 1 & pol.2 == -1) |
                     (pol.1 == -1 & pol.2 == 1)) %>% nrow()
          
          iter_count <- iter_count + 1
          if (num_heteropolares_int == 0 || iter_count >= max_iter) break
        }
        
        result_blocks <- as.data.frame(result_SA$block_final)
        colnames(result_blocks) <- c("item_num.1","item_num.2")
        
        result_chars <- result_blocks %>%
          left_join(data_item_chars %>% rename_with(~ paste0(.x, ".1"), -item_num), by = c("item_num.1" = "item_num")) %>%
          left_join(data_item_chars %>% rename_with(~ paste0(.x, ".2"), -item_num), by = c("item_num.2" = "item_num"))
        
      } else { #homo == F
        
        pairs <- expand.grid(item.1 = data_item_chars$item_num, item.2 = data_item_chars$item_num)
        pairs <- pairs[pairs$item.1 < pairs$item.2, ]

        pairs <- pairs %>%
          left_join(data_item_chars, by = c("item.1" = "item_num")) %>%
          rename(trait.1 = trait, pol.1 = pol, SD_rating.1 = SD_rating, Scale.1 = Scale, d1.1 = Threshold) %>%
          left_join(data_item_chars, by = c("item.2" = "item_num")) %>%
          rename(trait.2 = trait, pol.2 = pol, SD_rating.2 = SD_rating, Scale.2 = Scale, d1.2 = Threshold)

        valid_pairs <- pairs %>%
          mutate(
            homopolar = (pol.1 == pol.2),
            SD_restriction = abs(SD_rating.1 - SD_rating.2) <= cutoff_mixed) %>%
          filter(SD_restriction, trait.1 != trait.2, !homopolar)  # Only heteropolar and multidimensional

        num_pairs <- nrow(valid_pairs)

        unique_items <- unique(c(valid_pairs$item.1, valid_pairs$item.2)) # Number of unique items (because of the SD restriction, not all 320 items are in valid_pairs)
        length_unique_items <- length(unique_items)

        item_index <- match(unique_items, sort(unique_items))

        mat <- matrix(0, nrow = length_unique_items + 11, ncol = num_pairs) # Constraints for length_unique_items + 11 additional constraints

        for (j in 1:length_unique_items) {
          item_j <- unique_items[j]
          mat[j, ] <- ifelse(valid_pairs$item.1 == item_j | valid_pairs$item.2 == item_j, 1, 0)
        }

        mat[length_unique_items + 1,] <- ifelse(valid_pairs$trait.1 == 1 | valid_pairs$trait.2 == 1, 1, 0) # Restriction row for items with trait 1 (for lower constraint)
        mat[length_unique_items + 2,] <- ifelse(valid_pairs$trait.1 == 2 | valid_pairs$trait.2 == 2, 1, 0) # Restriction row for items with trait 2 (for lower constraint)
        mat[length_unique_items + 3,] <- ifelse(valid_pairs$trait.1 == 3 | valid_pairs$trait.2 == 3, 1, 0) # Restriction row for items with trait 3 (for lower constraint)
        mat[length_unique_items + 4,] <- ifelse(valid_pairs$trait.1 == 4 | valid_pairs$trait.2 == 4, 1, 0) # Restriction row for items with trait 4 (for lower constraint)
        mat[length_unique_items + 5,] <- ifelse(valid_pairs$trait.1 == 5 | valid_pairs$trait.2 == 5, 1, 0) # Restriction row for items with trait 5 (for lower constraint)
        mat[length_unique_items + 6,] <- ifelse(valid_pairs$trait.1 == 1 | valid_pairs$trait.2 == 1, 1, 0) # Restriction row for items with trait 1 (for higher constraint)
        mat[length_unique_items + 7,] <- ifelse(valid_pairs$trait.1 == 2 | valid_pairs$trait.2 == 2, 1, 0) # Restriction row for items with trait 2 (for higher constraint)
        mat[length_unique_items + 8,] <- ifelse(valid_pairs$trait.1 == 3 | valid_pairs$trait.2 == 3, 1, 0) # Restriction row for items with trait 3 (for higher constraint)
        mat[length_unique_items + 9,] <- ifelse(valid_pairs$trait.1 == 4 | valid_pairs$trait.2 == 4, 1, 0) # Restriction row for items with trait 4 (for higher constraint)
        mat[length_unique_items + 10,] <- ifelse(valid_pairs$trait.1 == 5 | valid_pairs$trait.2 == 5, 1, 0) # Restriction row for items with trait 5 (for higher constraint)

        mat[length_unique_items + 11,] <- ifelse((valid_pairs$pol.1 == 1 & valid_pairs$pol.2 == -1) | (valid_pairs$pol.1 == -1 & valid_pairs$pol.2 == 1), 1, 0)

        dir <- c(rep("<=", length_unique_items), rep(">=",5), rep("<=",5), "==")
        rhs <- c(rep(1, length_unique_items), rep(((n_blocks_hetero/5)*2)*0.8,5), rep(((n_blocks_hetero/5)*2)*1.2,5), n_blocks_hetero)

        obj <- runif(num_pairs)

        sol <- Rglpk_solve_LP(obj, mat, dir, rhs, max = TRUE, types = rep("B", num_pairs))

        selected_blocks <- valid_pairs[sol$solution == 1, ]
        colnames(selected_blocks)[1:2] <- c("item_num.1", "item_num.2")

        selected_blocks <- selected_blocks %>%
          select(item_num.1, item_num.2, trait.1, trait.2,
                 pol.1, pol.2, SD_rating.1, SD_rating.2,
                 Scale.1, Scale.2,
                 d1.1, d1.2)

        picked_items <- selected_blocks %>%
          pivot_longer(
            cols = starts_with("item_num"),
            names_to = "group",
            values_to = "item_num"
          ) %>%
          mutate(group = as.numeric(gsub("item_num.", "", group))) %>%
          arrange(row_number(), group) %>%
          mutate(
            trait = ifelse(group == 1, trait.1, trait.2),
            pol = ifelse(group == 1, pol.1, pol.2),
            SD_rating = ifelse(group == 1, SD_rating.1, SD_rating.2),
            Scale = ifelse(group == 1, Scale.1, Scale.2),
            Threshold = ifelse(group == 1, d1.1, d1.2)
          ) %>%
          select(item_num, trait, pol, Scale, Threshold, SD_rating)

        item_info_rest <- data_item_chars %>%
          filter(!item_num %in% picked_items$item_num)
        
        
        iter_count <- 0
        max_iter <- 200
        repeat {
        result_SA <- autoFC::sa_pairing_generalized(block = make_random_block(nrow(item_info_rest), (n_blocks_total-n_blocks_hetero)*2, 2),
                                                    item_chars = item_info_rest[,c("trait","SD_rating","pol","Scale")],
                                                    r = 0.999,
                                                    FUN = c("facfun","inn_diff","var","var"),
                                                    weights = c(10000,-1000000,-100000,10000)) 
        
        rest_items <- item_info_rest[c(t(result_SA$block_final)),]
        rest_items_blocks <- data.frame(
          item_num.1 = rest_items$item_num[seq(1, nrow(rest_items)-1, by=2)],
          item_num.2 = rest_items$item_num[seq(2, nrow(rest_items), by=2)],
          trait.1 = rest_items$trait[seq(1, nrow(rest_items)-1, by=2)],
          trait.2 = rest_items$trait[seq(2, nrow(rest_items), by=2)],
          pol.1 = rest_items$pol[seq(1, nrow(rest_items)-1, by=2)],
          pol.2 = rest_items$pol[seq(2, nrow(rest_items), by=2)],
          SD_rating.1 = rest_items$SD_rating[seq(1, nrow(rest_items)-1, by=2)],
          SD_rating.2 = rest_items$SD_rating[seq(2, nrow(rest_items), by=2)],
          Scale.1 = rest_items$Scale[seq(1, nrow(rest_items)-1, by=2)],
          Scale.2 = rest_items$Scale[seq(2, nrow(rest_items), by=2)],
          d1.1 = rest_items$Threshold[seq(1, nrow(rest_items)-1, by=2)],
          d1.2 = rest_items$Threshold[seq(2, nrow(rest_items), by=2)])
        
        num_heteropolares_int <- rest_items_blocks %>%
          filter((pol.1 == 1 & pol.2 == -1) |
                   (pol.1 == -1 & pol.2 == 1)) %>% nrow()
        
        iter_count <- iter_count + 1
        if (num_heteropolares_int == 0 || iter_count >= max_iter) break
        }
        
        result_blocks <- as.data.frame(matrix(c(picked_items$item_num, rest_items$item_num), ncol = 2, byrow = TRUE))
        colnames(result_blocks) <- c("item_num.1", "item_num.2")
        
        result_chars <- result_blocks %>%
          left_join(data_item_chars %>% rename_with(~ paste0(.x, ".1"), -item_num), by = c("item_num.1" = "item_num")) %>%
          left_join(data_item_chars %>% rename_with(~ paste0(.x, ".2"), -item_num), by = c("item_num.2" = "item_num"))
      }
    } else { # use_SD == F
      if (homo == T) {
        iter_count <- 0
        max_iter <- 200
        repeat {
          result_SA <- autoFC::sa_pairing_generalized(block = make_random_block(nrow(data_item_chars), n_blocks_total*2, 2),
                                                      item_chars = data_item_chars[,c("trait","pol","Scale")],
                                                      r = 0.999,
                                                      FUN = c("facfun","var","var"),
                                                      weights = c(10000,-100000,10000))
          
          proposed_items <- data_item_chars[c(t(result_SA$block_final)),]
          
          proposed_items_blocks <- data.frame(
            item_num.1 = proposed_items$item_num[seq(1, nrow(proposed_items)-1, by=2)],
            item_num.2 = proposed_items$item_num[seq(2, nrow(proposed_items), by=2)],
            trait.1 = proposed_items$trait[seq(1, nrow(proposed_items)-1, by=2)],
            trait.2 = proposed_items$trait[seq(2, nrow(proposed_items), by=2)],
            pol.1 = proposed_items$pol[seq(1, nrow(proposed_items)-1, by=2)],
            pol.2 = proposed_items$pol[seq(2, nrow(proposed_items), by=2)],
            SD_rating.1 = proposed_items$SD_rating[seq(1, nrow(proposed_items)-1, by=2)],
            SD_rating.2 = proposed_items$SD_rating[seq(2, nrow(proposed_items), by=2)],
            Scale.1 = proposed_items$Scale[seq(1, nrow(proposed_items)-1, by=2)],
            Scale.2 = proposed_items$Scale[seq(2, nrow(proposed_items), by=2)],
            d1.1 = proposed_items$Threshold[seq(1, nrow(proposed_items)-1, by=2)],
            d1.2 = proposed_items$Threshold[seq(2, nrow(proposed_items), by=2)])
          
          num_heteropolares_int <- proposed_items_blocks %>%
            filter((pol.1 == 1 & pol.2 == -1) |
                     (pol.1 == -1 & pol.2 == 1)) %>% nrow()
          
          iter_count <- iter_count + 1
          if (num_heteropolares_int == 0 || iter_count >= max_iter) break
        }

        result_blocks <- as.data.frame(result_SA$block_final)
        colnames(result_blocks) <- c("item_num.1","item_num.2")
        
        result_chars <- result_blocks %>%
          left_join(data_item_chars %>% rename_with(~ paste0(.x, ".1"), -item_num), by = c("item_num.1" = "item_num")) %>%
          left_join(data_item_chars %>% rename_with(~ paste0(.x, ".2"), -item_num), by = c("item_num.2" = "item_num"))
        
      } else { #homo == F
        
        pairs <- expand.grid(item.1 = data_item_chars$item_num, item.2 = data_item_chars$item_num)
        pairs <- pairs[pairs$item.1 < pairs$item.2, ]
        
        pairs <- pairs %>%
          left_join(data_item_chars, by = c("item.1" = "item_num")) %>%
          rename(trait.1 = trait, pol.1 = pol, SD_rating.1 = SD_rating, Scale.1 = Scale, d1.1 = Threshold) %>%
          left_join(data_item_chars, by = c("item.2" = "item_num")) %>%
          rename(trait.2 = trait, pol.2 = pol, SD_rating.2 = SD_rating, Scale.2 = Scale, d1.2 = Threshold)
        
        valid_pairs <- pairs %>%
          mutate(
            homopolar = (pol.1 == pol.2)) %>%
          filter(trait.1 != trait.2, !homopolar)  
        
        num_pairs <- nrow(valid_pairs)
        
        unique_items <- unique(c(valid_pairs$item.1, valid_pairs$item.2)) # Number of unique items (because of the SD restriction, not all 320 items are in valid_pairs)
        length_unique_items <- length(unique_items)
        
        item_index <- match(unique_items, sort(unique_items))
        
        mat <- matrix(0, nrow = length_unique_items + 11, ncol = num_pairs) # Constraints for length_unique_items + 11 additional constraints
        
        for (j in 1:length_unique_items) {
          item_j <- unique_items[j]
          mat[j, ] <- ifelse(valid_pairs$item.1 == item_j | valid_pairs$item.2 == item_j, 1, 0)
        }
        
        mat[length_unique_items + 1,] <- ifelse(valid_pairs$trait.1 == 1 | valid_pairs$trait.2 == 1, 1, 0) # Restriction row for items with trait 1 (for lower constraint)
        mat[length_unique_items + 2,] <- ifelse(valid_pairs$trait.1 == 2 | valid_pairs$trait.2 == 2, 1, 0) # Restriction row for items with trait 2 (for lower constraint)
        mat[length_unique_items + 3,] <- ifelse(valid_pairs$trait.1 == 3 | valid_pairs$trait.2 == 3, 1, 0) # Restriction row for items with trait 3 (for lower constraint)
        mat[length_unique_items + 4,] <- ifelse(valid_pairs$trait.1 == 4 | valid_pairs$trait.2 == 4, 1, 0) # Restriction row for items with trait 4 (for lower constraint)
        mat[length_unique_items + 5,] <- ifelse(valid_pairs$trait.1 == 5 | valid_pairs$trait.2 == 5, 1, 0) # Restriction row for items with trait 5 (for lower constraint)
        mat[length_unique_items + 6,] <- ifelse(valid_pairs$trait.1 == 1 | valid_pairs$trait.2 == 1, 1, 0) # Restriction row for items with trait 1 (for higher constraint)
        mat[length_unique_items + 7,] <- ifelse(valid_pairs$trait.1 == 2 | valid_pairs$trait.2 == 2, 1, 0) # Restriction row for items with trait 2 (for higher constraint)
        mat[length_unique_items + 8,] <- ifelse(valid_pairs$trait.1 == 3 | valid_pairs$trait.2 == 3, 1, 0) # Restriction row for items with trait 3 (for higher constraint)
        mat[length_unique_items + 9,] <- ifelse(valid_pairs$trait.1 == 4 | valid_pairs$trait.2 == 4, 1, 0) # Restriction row for items with trait 4 (for higher constraint)
        mat[length_unique_items + 10,] <- ifelse(valid_pairs$trait.1 == 5 | valid_pairs$trait.2 == 5, 1, 0) # Restriction row for items with trait 5 (for higher constraint)
        
        mat[length_unique_items + 11,] <- ifelse((valid_pairs$pol.1 == 1 & valid_pairs$pol.2 == -1) | (valid_pairs$pol.1 == -1 & valid_pairs$pol.2 == 1), 1, 0)
        
        dir <- c(rep("<=", length_unique_items), rep(">=",5), rep("<=",5), "==")
        rhs <- c(rep(1, length_unique_items), rep(((n_blocks_hetero/5)*2)*0.8,5), rep(((n_blocks_hetero/5)*2)*1.2,5), n_blocks_hetero)
        
        obj <- runif(num_pairs)
        
        sol <- Rglpk_solve_LP(obj, mat, dir, rhs, max = TRUE, types = rep("B", num_pairs))
        
        selected_blocks <- valid_pairs[sol$solution == 1, ]
        colnames(selected_blocks)[1:2] <- c("item_num.1", "item_num.2")
        
        selected_blocks <- selected_blocks %>%
          select(item_num.1, item_num.2, trait.1, trait.2,
                 pol.1, pol.2, SD_rating.1, SD_rating.2,
                 Scale.1, Scale.2,
                 d1.1, d1.2)
        
        picked_items <- selected_blocks %>%
          pivot_longer(
            cols = starts_with("item_num"),
            names_to = "group",
            values_to = "item_num"
          ) %>%
          mutate(group = as.numeric(gsub("item_num.", "", group))) %>%
          arrange(row_number(), group) %>%
          mutate(
            trait = ifelse(group == 1, trait.1, trait.2),
            pol = ifelse(group == 1, pol.1, pol.2),
            SD_rating = ifelse(group == 1, SD_rating.1, SD_rating.2),
            Scale = ifelse(group == 1, Scale.1, Scale.2),
            Threshold = ifelse(group == 1, d1.1, d1.2)
          ) %>%
          select(item_num, trait, pol, Scale, Threshold, SD_rating)
        
        item_info_rest <- data_item_chars %>%
          filter(!item_num %in% picked_items$item_num)
        
        
        iter_count <- 0
        max_iter <- 200
        repeat {
        result_SA <- autoFC::sa_pairing_generalized(block = make_random_block(nrow(item_info_rest), (n_blocks_total-n_blocks_hetero)*2, 2),
                                                    item_chars = item_info_rest[,c("trait","pol","Scale")],
                                                    r = 0.999,
                                                    FUN = c("facfun","var","var"),
                                                    weights = c(10000,-100000,10000)) 
        
        rest_items <- item_info_rest[c(t(result_SA$block_final)),]
        rest_items_blocks <- data.frame(
          item_num.1 = rest_items$item_num[seq(1, nrow(rest_items)-1, by=2)],
          item_num.2 = rest_items$item_num[seq(2, nrow(rest_items), by=2)],
          trait.1 = rest_items$trait[seq(1, nrow(rest_items)-1, by=2)],
          trait.2 = rest_items$trait[seq(2, nrow(rest_items), by=2)],
          pol.1 = rest_items$pol[seq(1, nrow(rest_items)-1, by=2)],
          pol.2 = rest_items$pol[seq(2, nrow(rest_items), by=2)],
          SD_rating.1 = rest_items$SD_rating[seq(1, nrow(rest_items)-1, by=2)],
          SD_rating.2 = rest_items$SD_rating[seq(2, nrow(rest_items), by=2)],
          Scale.1 = rest_items$Scale[seq(1, nrow(rest_items)-1, by=2)],
          Scale.2 = rest_items$Scale[seq(2, nrow(rest_items), by=2)],
          d1.1 = rest_items$Threshold[seq(1, nrow(rest_items)-1, by=2)],
          d1.2 = rest_items$Threshold[seq(2, nrow(rest_items), by=2)])
        
        num_heteropolares_int <- rest_items_blocks %>%
          filter((pol.1 == 1 & pol.2 == -1) |
                   (pol.1 == -1 & pol.2 == 1)) %>% nrow()
        
        iter_count <- iter_count + 1
        if (num_heteropolares_int == 0 || iter_count >= max_iter) break
        }
        
        result_blocks <- as.data.frame(matrix(c(picked_items$item_num, rest_items$item_num), ncol = 2, byrow = TRUE))
        colnames(result_blocks) <- c("item_num.1", "item_num.2")
        
        result_chars <- result_blocks %>%
          left_join(data_item_chars %>% rename_with(~ paste0(.x, ".1"), -item_num), by = c("item_num.1" = "item_num")) %>%
          left_join(data_item_chars %>% rename_with(~ paste0(.x, ".2"), -item_num), by = c("item_num.2" = "item_num"))
      }
    }
    
    item_chars_1 <- data_item_chars %>%
      rename_with(~ paste0(.x, ".1"), -item_num)
    item_chars_2 <- data_item_chars %>%
      rename_with(~ paste0(.x, ".2"), -item_num)
    
    result_chars <- result_blocks %>%
      left_join(item_chars_1, by = c("item_num.1" = "item_num")) %>%
      left_join(item_chars_2, by = c("item_num.2" = "item_num"))
    
    result_chars <- result_chars %>%
      mutate(
        d1.1 = Threshold.1,
        d1.2 = Threshold.2) %>%
      mutate(
        a1.1 = ifelse(trait.1 == 1, Scale.1, 0),
        a2.1 = ifelse(trait.1 == 2, Scale.1, 0),
        a3.1 = ifelse(trait.1 == 3, Scale.1, 0),
        a4.1 = ifelse(trait.1 == 4, Scale.1, 0),
        a5.1 = ifelse(trait.1 == 5, Scale.1, 0),
        a1.2 = ifelse(trait.2 == 1, Scale.2, 0),
        a2.2 = ifelse(trait.2 == 2, Scale.2, 0),
        a3.2 = ifelse(trait.2 == 3, Scale.2, 0),
        a4.2 = ifelse(trait.2 == 4, Scale.2, 0),
        a5.2 = ifelse(trait.2 == 5, Scale.2, 0))
    
    result_chars <- result_chars %>%
      select(item_num.1, item_num.2, 
             trait.1, trait.2,
             pol.1, pol.2, SD_rating.1, SD_rating.2,
             d1.1, d1.2,
             a1.1, a2.1, a3.1, a4.1, a5.1,
             a1.2, a2.2, a3.2, a4.2, a5.2)
   } else if (algorithm == "BF"){
    reliability_results <-  vector("list", BF_iter)
    result_chars_list <- vector("list", BF_iter)
    
    theta_uncorr <- mirt:::QMC_quad(500, K)
    Sigma_chol_matrix <- chol(Sigma)
    theta <- theta_uncorr %*% Sigma_chol_matrix
    
    prior_Info <- array(solve(Sigma), dim = c(K, K, nrow(theta))) 
    euc_Weights <- rep(1, nrow(theta))
    euc_Weights <- euc_Weights/sum(euc_Weights)
    
    for (i in 1:BF_iter){
      if (use_SD == T){
        pairs <- expand.grid(item.1 = data_item_chars$item_num, item.2 = data_item_chars$item_num)
        pairs <- pairs[pairs$item.1 < pairs$item.2, ]
        
        pairs <- pairs %>%
          left_join(data_item_chars, by = c("item.1" = "item_num")) %>%
          rename(trait.1 = trait, pol.1 = pol, SD_rating.1 = SD_rating, Scale.1 = Scale, d1.1 = Threshold) %>%
          left_join(data_item_chars, by = c("item.2" = "item_num")) %>%
          rename(trait.2 = trait, pol.2 = pol, SD_rating.2 = SD_rating, Scale.2 = Scale, d1.2 = Threshold)
        
        valid_pairs <- pairs %>%
          mutate(
            homopolar = (pol.1 == pol.2),  
            SD_restriction = ifelse(
              homopolar,  
              abs(SD_rating.1 - SD_rating.2) <= cutoff_equal,  # Restriction for homopolar blocks
              abs(SD_rating.1 - SD_rating.2) <= cutoff_mixed)) %>%   # Restriction for heteropolar blocks
          filter(SD_restriction, trait.1 != trait.2) # Restriction for only multidimensional blocks
      } else { #use_SD == F
        pairs <- expand.grid(item.1 = data_item_chars$item_num, item.2 = data_item_chars$item_num)
        pairs <- pairs[pairs$item.1 < pairs$item.2, ]
        
        pairs <- pairs %>%
          left_join(data_item_chars, by = c("item.1" = "item_num")) %>%
          rename(trait.1 = trait, pol.1 = pol, SD_rating.1 = SD_rating, Scale.1 = Scale, d1.1 = Threshold) %>%
          left_join(data_item_chars, by = c("item.2" = "item_num")) %>%
          rename(trait.2 = trait, pol.2 = pol, SD_rating.2 = SD_rating, Scale.2 = Scale, d1.2 = Threshold)
        
        valid_pairs <- pairs %>%
          filter(trait.1 != trait.2)
      }
      
      num_pairs <- nrow(valid_pairs)
      
      unique_items <- unique(c(valid_pairs$item.1, valid_pairs$item.2)) #number of unique items (because of the SD restriction, not all items are in valid_pairs)
      length_unique_items <- length(unique_items)
      
      mat <- matrix(0, nrow = length_unique_items + 13, ncol = num_pairs) 
      
      for (j in 1:length_unique_items) {
        mat[j,] <- ifelse(valid_pairs$item.1 == j | valid_pairs$item.2 == j, 1, 0)
      } # Restriction rows for each item
      
      mat[length_unique_items + 1,] <- ifelse(valid_pairs$trait.1 == 1 | valid_pairs$trait.2 == 1, 1, 0) # Restriction row for items with trait 1 (for lower distribution)
      mat[length_unique_items + 2,] <- ifelse(valid_pairs$trait.1 == 2 | valid_pairs$trait.2 == 2, 1, 0) # Restriction row for items with trait 2 (for lower distribution)
      mat[length_unique_items + 3,] <- ifelse(valid_pairs$trait.1 == 3 | valid_pairs$trait.2 == 3, 1, 0) # Restriction row for items with trait 3 (for lower distribution)
      mat[length_unique_items + 4,] <- ifelse(valid_pairs$trait.1 == 4 | valid_pairs$trait.2 == 4, 1, 0) # Restriction row for items with trait 4 (for lower distribution)
      mat[length_unique_items + 5,] <- ifelse(valid_pairs$trait.1 == 5 | valid_pairs$trait.2 == 5, 1, 0) # Restriction row for items with trait 5 (for lower distribution)
      mat[length_unique_items + 6,] <- ifelse(valid_pairs$trait.1 == 1 | valid_pairs$trait.2 == 1, 1, 0) # Restriction row for items with trait 1 (for higher distribution)
      mat[length_unique_items + 7,] <- ifelse(valid_pairs$trait.1 == 2 | valid_pairs$trait.2 == 2, 1, 0) # Restriction row for items with trait 2 (for higher distribution)
      mat[length_unique_items + 8,] <- ifelse(valid_pairs$trait.1 == 3 | valid_pairs$trait.2 == 3, 1, 0) # Restriction row for items with trait 3 (for higher distribution)
      mat[length_unique_items + 9,] <- ifelse(valid_pairs$trait.1 == 4 | valid_pairs$trait.2 == 4, 1, 0) # Restriction row for items with trait 4 (for higher distribution)
      mat[length_unique_items + 10,] <- ifelse(valid_pairs$trait.1 == 5 | valid_pairs$trait.2 == 5, 1, 0) # Restriction row for items with trait 5 (for higher distribution)
      mat[length_unique_items + 11,] <- ifelse((valid_pairs$pol.1 == 1 & valid_pairs$pol.2 == -1) | (valid_pairs$pol.1 == -1 & valid_pairs$pol.2 == 1), 1, 0)
      mat[length_unique_items + 12,] <- ifelse((valid_pairs$pol.1 == 1 & valid_pairs$pol.2 == 1) | (valid_pairs$pol.1 == -1 & valid_pairs$pol.2 == -1), 1, 0)
      mat[length_unique_items + 13,] <- 1  # Restriction row for not permiting the repetition of items, since valid_pairs is already a filter for it. Thus, all values are 1 because we consider all have valid blocks
      
      dir <- c(rep("<=", length_unique_items), rep(">=",5), rep("<=",5), "==", "==", "==")
      rhs <- c(rep(1, length_unique_items), rep(((n_blocks_total/5)*2)*0.8,5), rep(((n_blocks_total/5)*2)*1.2,5), n_blocks_hetero, n_blocks_homo, n_blocks_total)
      
      obj <- runif(num_pairs)
      
      sol <- Rglpk_solve_LP(obj, mat, dir, rhs, max = TRUE, types = rep("B", num_pairs))
      
      result_chars <- valid_pairs[sol$solution == 1, ]
      colnames(result_chars)[1:2] <- c("item_num.1", "item_num.2")
      
      result_chars <- result_chars %>%
        mutate(
          a1.1 = ifelse(trait.1 == 1, Scale.1, 0),
          a2.1 = ifelse(trait.1 == 2, Scale.1, 0),
          a3.1 = ifelse(trait.1 == 3, Scale.1, 0),
          a4.1 = ifelse(trait.1 == 4, Scale.1, 0),
          a5.1 = ifelse(trait.1 == 5, Scale.1, 0),
          a1.2 = ifelse(trait.2 == 1, Scale.2, 0),
          a2.2 = ifelse(trait.2 == 2, Scale.2, 0),
          a3.2 = ifelse(trait.2 == 3, Scale.2, 0),
          a4.2 = ifelse(trait.2 == 4, Scale.2, 0),
          a5.2 = ifelse(trait.2 == 5, Scale.2, 0)) %>%
        select(item_num.1, item_num.2, trait.1, trait.2,
               pol.1, pol.2, SD_rating.1, SD_rating.2,
               d1.1, d1.2,
               a1.1, a2.1, a3.1, a4.1, a5.1,
               a1.2, a2.2, a3.2, a4.2, a5.2)
      
      result_chars_list[[i]] <- result_chars
      
    # Obtain estimated thetas and reliability for all candidates
        
        if (BF_relia_method == "empirical"){
          res_thetas_BF <- try(res_empirical(solution = result_chars, Sigma = Sigma))
          reliability_results[[i]] <- if (!inherits(res_thetas_BF, "try-error")) as.list(res_thetas_BF) else list(cor_est_true_general = NA)
          
        } else if (BF_relia_method == "marginal"){
          result_chars_matrix <- as.matrix(result_chars)
          storage.mode(result_chars_matrix) <- "numeric"
          
          selected_BlockInfo <- blocksInfo(result_chars_matrix, Sigma, theta, TIRT = FALSE)
          res_relia <- AvgMargReli(selected_BlockInfo, prior_Info, euc_Weights)
          reliability_results[[i]] <- res_relia
        }
    }
    
     # Select the best candidate
      if (BF_relia_method == "empirical"){
        relia_values <- sapply(reliability_results, function(x) if (!is.null(x) && "cor_est_true_general" %in% names(x)) x["cor_est_true_general"] else NA)
        best_BF_solution <- which.max(relia_values)
        result_chars <- result_chars_list[[best_BF_solution]]
        
        flat_all_reliability <- lapply(reliability_results, function(x) {
          if (is.null(x)) {
            return(data.frame(matrix(ncol = 0, nrow = 1)))  
          } else {
            return(as.data.frame(t(x))) 
          }
        })
        
        res_relia_BF <- flat_all_reliability[[best_BF_solution]]
        
      } else if (BF_relia_method == "marginal"){
        relia_values <- sapply(reliability_results, function(x) if (is.null(x)) NA else x)
        best_BF_solution <- which.max(relia_values)
        
        result_chars <- result_chars_list[[best_BF_solution]]
        
        flat_all_reliability <- lapply(reliability_results, function(x) {
          if (is.null(x)) {
            return(data.frame(matrix(ncol = 0, nrow = 1)))  
          } else {
            return(as.data.frame(t(x))) 
          }
        })
        
        res_relia_BF <- try(get_res(solution = result_chars, thetas_true = thetas_true, Sigma = Sigma, method = "both")) # method = "marginal"
      }
      
    
    result_blocks <- as.data.frame(result_chars[,c("item_num.1", "item_num.2")])
    
    df_reliability_results <- do.call(cbind, flat_all_reliability)
    
     write.table(x = data.frame(conds, c(df_reliability_results)), 
                 file = "sim_res_BF_allRes.txt", 
                 append = TRUE, row.names = FALSE, col.names = !file.exists("sim_res_BF_allRes.txt")) # This file has all the BF solutions
    
    iter_count <- 0
  }
  
  # Obtain all results characteristics 
  comb <- result_chars %>%
    group_by(trait_comb = paste(pmin(trait.1, trait.2), pmax(trait.1, trait.2), sep = "-")) %>%
    summarise(count = n(), .groups = "drop") %>%
    right_join(all_combinations, by = "trait_comb") %>%  
    replace_na(list(count = 0)) %>%  
    pivot_wider(names_from = trait_comb, values_from = count, values_fill = list(count = 0)) 
  names(comb) <- c("dim1_dim2", "dim1_dim3", "dim1_dim4", "dim1_dim5", "dim2_dim3", "dim2_dim4", "dim2_dim5", "dim3_dim4", "dim3_dim5", "dim4_dim5")
  
  count_traits <- table(c(result_chars$trait.1, result_chars$trait.2))
  names(count_traits) <- c("count_dim1","count_dim2","count_dim3","count_dim4","count_dim5")
  
  num_homopolares <- result_chars %>% 
    filter((pol.1 == 1 & pol.2 == 1) |
             (pol.1 == -1 & pol.2 == -1)) %>% nrow()
  
  num_heteropolares <- result_chars %>%
    filter((pol.1 == 1 & pol.2 == -1) |
             (pol.1 == -1 & pol.2 == 1)) %>% nrow()
  comb_polarity <- data.frame(homopolar = num_homopolares, heteropolar = num_heteropolares)
  
  # Outputs
    results <- list(solution = result_chars, blocks = result_blocks, combinations = comb, count_traits = count_traits, num_polarity = comb_polarity, iter_count = iter_count)

  
  return(results)
}

get_res <- function(solution, thetas_true = thetas_true, Sigma = Sigma, method = "marginal") { 
  
  blocks <- as.data.frame(cbind(solution$item_num.1, solution$item_num.2))
  colnames(blocks) <- c("id.1", "id.2")
  ndim <- K
  n_blocks_total <- nrow(blocks)
  solution <- solution %>%
    mutate(a.1 = rowSums(select(., starts_with("a") & ends_with(".1"))),
           a.2 = rowSums(select(., starts_with("a") & ends_with(".2"))))
  
  A.matrix.1 <- solution[,c("a1.1", "a2.1", "a3.1", "a4.1", "a5.1")]
  A.matrix.2 <- solution[,c("a1.2", "a2.2", "a3.2", "a4.2", "a5.2")]
  A.matrix <- as.matrix(A.matrix.2 - A.matrix.1) 
  colnames(A.matrix) <- NULL
  d1 <- -solution$a.1*solution$d1.1
  d2 <- -solution$a.2*solution$d1.2
  d <- d2-d1
  
  if (method == "empirical"){ 
    data.FC <- mirt::simdata(a = A.matrix, 
                             d = matrix(d, ncol = 1), 
                             Theta = thetas_true, 
                             itemtype = "2PL")
    
    mod.specs <- rep(NA, K)
    for(k in 1:K){
      mod.specs[k] <- paste0("F", k," = ", paste(c(which(solution$trait.1 == k), which(solution$trait.2 == k)), collapse = ", "))
    }
    mod.specs[ncol(Sigma)+1] <- paste0("COV = ", paste(paste0("F",1:ncol(Sigma)), collapse = "*"))
    mirt.model_FC <- mirt.model(paste(mod.specs, collapse = "\n"))
    
    # Estimated model
    mod_est <- mirt(data.FC, model = mirt.model_FC, itemtype = "2PL", method = "MHRM", SE = FALSE, verbose = FALSE, technical = list(NCYCLES = 500))
    #param_est <- as.data.frame(coef(mod_est, simplify = TRUE)$items)
    
    thetas_est_com <- fscores(mod_est, response.pattern = data.FC, QMC = T, full.scores = TRUE) # estimated thetas/factor scores 
    thetas_est <- thetas_est_com[,1:5]    
    thetas_est_se <- thetas_est_com[,6:10]
    #empirical_relia <- empirical_rxx(thetas_est)
    
    # Obtain comparative measures
    
    cor_est_true_F1 <- round(cor(thetas_est[,1], thetas_true[,1])^2, 4) 
    cor_est_true_F2 <- round(cor(thetas_est[,2], thetas_true[,2])^2, 4) 
    cor_est_true_F3 <- round(cor(thetas_est[,3], thetas_true[,3])^2, 4) 
    cor_est_true_F4 <- round(cor(thetas_est[,4], thetas_true[,4])^2, 4) 
    cor_est_true_F5 <- round(cor(thetas_est[,5], thetas_true[,5])^2, 4)
    cor_est_true_general <- rowMeans(cbind(cor_est_true_F1,
                                           cor_est_true_F2, 
                                           cor_est_true_F3, 
                                           cor_est_true_F4, 
                                           cor_est_true_F5))
    
    cor_est_true_notsq_F1 <- round(cor(thetas_est[,1], thetas_true[,1]), 4) 
    cor_est_true_notsq_F2 <- round(cor(thetas_est[,2], thetas_true[,2]), 4) 
    cor_est_true_notsq_F3 <- round(cor(thetas_est[,3], thetas_true[,3]), 4) 
    cor_est_true_notsq_F4 <- round(cor(thetas_est[,4], thetas_true[,4]), 4) 
    cor_est_true_notsq_F5 <- round(cor(thetas_est[,5], thetas_true[,5]), 4)
    
    rmse_general <- round(sqrt((sum((thetas_true - thetas_est)^2))/(nrow(thetas_true)*ncol(thetas_true))), 4) 
    rmse_F1 <- round(sqrt((sum((thetas_true[,1] - thetas_est[,1])^2))/(nrow(thetas_true))), 4)
    rmse_F2 <- round(sqrt((sum((thetas_true[,2] - thetas_est[,2])^2))/(nrow(thetas_true))), 4)
    rmse_F3 <- round(sqrt((sum((thetas_true[,3] - thetas_est[,3])^2))/(nrow(thetas_true))), 4)
    rmse_F4 <- round(sqrt((sum((thetas_true[,4] - thetas_est[,4])^2))/(nrow(thetas_true))), 4)
    rmse_F5 <- round(sqrt((sum((thetas_true[,5] - thetas_est[,5])^2))/(nrow(thetas_true))), 4)
    rmse_means <- rowMeans(cbind(rmse_F1, rmse_F2, rmse_F3, rmse_F4, rmse_F5))
    
    bias_general <- round((sum(thetas_true - thetas_est))/(nrow(thetas_true)*ncol(thetas_true)), 4)
    bias_F1 <- round((sum(thetas_true[,1] - thetas_est[,1]))/(nrow(thetas_true)), 4)
    bias_F2 <- round((sum(thetas_true[,2] - thetas_est[,2]))/(nrow(thetas_true)), 4)
    bias_F3 <- round((sum(thetas_true[,3] - thetas_est[,3]))/(nrow(thetas_true)), 4)
    bias_F4 <- round((sum(thetas_true[,4] - thetas_est[,4]))/(nrow(thetas_true)), 4)
    bias_F5 <- round((sum(thetas_true[,5] - thetas_est[,5]))/(nrow(thetas_true)), 4)
    bias_means <- rowMeans(cbind(bias_F1, bias_F2, bias_F3, bias_F4, bias_F5))
    
    
    res <- c(
      cor_est_true_F1 = cor_est_true_F1,
      cor_est_true_F2 = cor_est_true_F2,
      cor_est_true_F3 = cor_est_true_F3,
      cor_est_true_F4 = cor_est_true_F4,
      cor_est_true_F5 = cor_est_true_F5,
      cor_est_true_general = cor_est_true_general,
      cor_est_true_notsq_F1 = cor_est_true_notsq_F1,
      cor_est_true_notsq_F2 = cor_est_true_notsq_F2,
      cor_est_true_notsq_F3 = cor_est_true_notsq_F3,
      cor_est_true_notsq_F4 = cor_est_true_notsq_F4,
      cor_est_true_notsq_F5 = cor_est_true_notsq_F5,
      rmse_general = rmse_general,
      rmse_F1 = rmse_F1,
      rmse_F2 = rmse_F2,
      rmse_F3 = rmse_F3,
      rmse_F4 = rmse_F4,
      rmse_F5 = rmse_F5,
      rmse_means = rmse_means,
      bias_general = bias_general,
      bias_F1 = bias_F1,
      bias_F2 = bias_F2,
      bias_F3 = bias_F3,
      bias_F4 = bias_F4,
      bias_F5 = bias_F5,
      bias_means = bias_means
    ) 
  } else if (method == "marginal") {
    prior_Info <- array(solve(Sigma), dim = c(ncol(thetas_true), ncol(thetas_true), nrow(thetas_true)))
    euc_Weights <- rep(1, nrow(thetas_true))
    euc_Weights <- euc_Weights/sum(euc_Weights)
    
    result_chars_matrix <- as.matrix(solution)
    storage.mode(result_chars_matrix) <- "numeric"
    
    selected_BlockInfo <- blocksInfo(result_chars_matrix, Sigma, thetas_true, TIRT = FALSE)
    res <- AvgMargReli(selected_BlockInfo, prior_Info, euc_Weights)
    names(res) <- "AvgMargReli"
    
  } else if (method == "both"){
    data.FC <- mirt:::simdata(a = A.matrix, 
                              d = matrix(d, ncol = 1), 
                              Theta = thetas_true, 
                              itemtype = "2PL")
    
    mod.specs <- rep(NA, K)
    for(k in 1:K){
      mod.specs[k] <- paste0("F", k," = ", paste(c(which(solution$trait.1 == k), which(solution$trait.2 == k)), collapse = ", "))
    }
    mod.specs[ncol(Sigma)+1] <- paste0("COV = ", paste(paste0("F",1:ncol(Sigma)), collapse = "*"))
    mirt.model_FC <- mirt.model(paste(mod.specs, collapse = "\n"))
    
    sv <- mirt(data.FC, model = mirt.model_FC, itemtype = "2PL", pars = "values")
    item_names <- unique(sv$item)
    
    A <- matrix(0, nrow(solution), 5)
    for(i in 1:nrow(solution)){
      A[i,solution$trait.1[i]] <-  as.numeric(solution$pol.1[i])
      A[i,solution$trait.2[i]] <- as.numeric(solution$pol.2[i])
    }
    
    for(j in 1:nrow(solution)) {
      item1 <- item_names[solution$item_num.1[j]]
      item2 <- item_names[solution$item_num.2[j]]
      
      for(k in 1:5) {
        idx1 <- which(sv$item == item1 & sv$name == paste0("a", k))
        if (length(idx1) == 1) {
          sv[idx1, "value"] <- sv[idx1, "value"] * A[j, k]
        }
        idx2 <- which(sv$item == item2 & sv$name == paste0("a", k))
        if (length(idx2) == 1) {
          sv[idx2, "value"] <- sv[idx2, "value"] * A[j, k]
        }
      }
    }
    
    # Estimated model
    mod_est <- mirt(data.FC, model = mirt.model_FC, itemtype = "2PL", method = "MHRM", SE = FALSE, verbose = FALSE, technical = list(NCYCLES = 500))
    #param_est <- as.data.frame(coef(mod_est, simplify = TRUE)$items)
    
    thetas_est_com <- fscores(mod_est, response.pattern = data.FC, QMC = T, full.scores = TRUE) # estimated thetas/factor scores 
    thetas_est <- thetas_est_com[,1:5]
    thetas_est_se <- thetas_est_com[,6:10]
    #empirical_relia <- empirical_rxx(thetas_est)
    
    # Obtain comparative measures
    
    cor_est_true_F1 <- round(cor(thetas_est[,1], thetas_true[,1])^2, 4) 
    cor_est_true_F2 <- round(cor(thetas_est[,2], thetas_true[,2])^2, 4) 
    cor_est_true_F3 <- round(cor(thetas_est[,3], thetas_true[,3])^2, 4) 
    cor_est_true_F4 <- round(cor(thetas_est[,4], thetas_true[,4])^2, 4) 
    cor_est_true_F5 <- round(cor(thetas_est[,5], thetas_true[,5])^2, 4)
    cor_est_true_general <- rowMeans(cbind(cor_est_true_F1,
                                           cor_est_true_F2, 
                                           cor_est_true_F3, 
                                           cor_est_true_F4, 
                                           cor_est_true_F5))
    
    cor_est_true_notsq_F1 <- round(cor(thetas_est[,1], thetas_true[,1]), 4) 
    cor_est_true_notsq_F2 <- round(cor(thetas_est[,2], thetas_true[,2]), 4) 
    cor_est_true_notsq_F3 <- round(cor(thetas_est[,3], thetas_true[,3]), 4) 
    cor_est_true_notsq_F4 <- round(cor(thetas_est[,4], thetas_true[,4]), 4) 
    cor_est_true_notsq_F5 <- round(cor(thetas_est[,5], thetas_true[,5]), 4)
    
    rmse_general <- round(sqrt((sum((thetas_true - thetas_est)^2))/(nrow(thetas_true)*ncol(thetas_true))), 4) 
    rmse_F1 <- round(sqrt((sum((thetas_true[,1] - thetas_est[,1])^2))/(nrow(thetas_true))), 4)
    rmse_F2 <- round(sqrt((sum((thetas_true[,2] - thetas_est[,2])^2))/(nrow(thetas_true))), 4)
    rmse_F3 <- round(sqrt((sum((thetas_true[,3] - thetas_est[,3])^2))/(nrow(thetas_true))), 4)
    rmse_F4 <- round(sqrt((sum((thetas_true[,4] - thetas_est[,4])^2))/(nrow(thetas_true))), 4)
    rmse_F5 <- round(sqrt((sum((thetas_true[,5] - thetas_est[,5])^2))/(nrow(thetas_true))), 4)
    rmse_means <- rowMeans(cbind(rmse_F1, rmse_F2, rmse_F3, rmse_F4, rmse_F5))
    
    bias_general <- round((sum(thetas_true - thetas_est))/(nrow(thetas_true)*ncol(thetas_true)), 4)
    bias_F1 <- round((sum(thetas_true[,1] - thetas_est[,1]))/(nrow(thetas_true)), 4)
    bias_F2 <- round((sum(thetas_true[,2] - thetas_est[,2]))/(nrow(thetas_true)), 4)
    bias_F3 <- round((sum(thetas_true[,3] - thetas_est[,3]))/(nrow(thetas_true)), 4)
    bias_F4 <- round((sum(thetas_true[,4] - thetas_est[,4]))/(nrow(thetas_true)), 4)
    bias_F5 <- round((sum(thetas_true[,5] - thetas_est[,5]))/(nrow(thetas_true)), 4)
    bias_means <- rowMeans(cbind(bias_F1, bias_F2, bias_F3, bias_F4, bias_F5))
    
    prior_Info <- array(solve(Sigma), dim = c(ncol(thetas_true), ncol(thetas_true), nrow(thetas_true)))
    euc_Weights <- rep(1, nrow(thetas_true))
    euc_Weights <- euc_Weights/sum(euc_Weights)
    
    result_chars_matrix <- as.matrix(solution)
    storage.mode(result_chars_matrix) <- "numeric"
    
    selected_BlockInfo <- blocksInfo(result_chars_matrix, Sigma, thetas_true, TIRT = FALSE)
    res_marg <- AvgMargReli_s(selected_BlockInfo, prior_Info, euc_Weights)
    res_AvgMargReli <- res_marg$avg_reli
    names(res_AvgMargReli) <- "AvgMargReli"
    all_marg_reli <- res_marg$all_reli
    
    res <- c(
      cor_est_true_F1 = cor_est_true_F1,
      cor_est_true_F2 = cor_est_true_F2,
      cor_est_true_F3 = cor_est_true_F3,
      cor_est_true_F4 = cor_est_true_F4,
      cor_est_true_F5 = cor_est_true_F5,
      cor_est_true_general = cor_est_true_general,
      cor_est_true_notsq_F1 = cor_est_true_notsq_F1,
      cor_est_true_notsq_F2 = cor_est_true_notsq_F2,
      cor_est_true_notsq_F3 = cor_est_true_notsq_F3,
      cor_est_true_notsq_F4 = cor_est_true_notsq_F4,
      cor_est_true_notsq_F5 = cor_est_true_notsq_F5,
      rmse_general = rmse_general,
      rmse_F1 = rmse_F1,
      rmse_F2 = rmse_F2,
      rmse_F3 = rmse_F3,
      rmse_F4 = rmse_F4,
      rmse_F5 = rmse_F5,
      rmse_means = rmse_means,
      bias_general = bias_general,
      bias_F1 = bias_F1,
      bias_F2 = bias_F2,
      bias_F3 = bias_F3,
      bias_F4 = bias_F4,
      bias_F5 = bias_F5,
      bias_means = bias_means,
      res_AvgMargReli
    ) 

  }
  
  return(list(res = res, thetas_est_com = as.data.frame(thetas_est_com), all_marg_reli = all_marg_reli))
}


AvgMargReli_s <- function(selectedBlockInfo, priorInfo, eucWeights){
  
  nTraits <- ncol(priorInfo[,,1])
  testInfo = apply(selectedBlockInfo,c(1:3),sum)
  
  BayesTestInfo = array(sapply(1:dim(testInfo)[3],function(i){
    testInfo[,,i]+priorInfo[,,i]
  }),dim=dim(testInfo))
  
  ThetaVars <- apply(BayesTestInfo,3,function(x) (diag(solve(x))))
  
  reli <- 1-t(eucWeights)%*%t(ThetaVars)
  avg_reli <- mean(1-t(eucWeights)%*%t(ThetaVars))
  all_reli <- colMeans(1 - ThetaVars)
  
  return(list(reli = reli, avg_reli = avg_reli, all_reli = all_reli))
}

res_empirical <- function(solution, Sigma = Sigma) { 
  
  blocks <- as.data.frame(cbind(solution$item_num.1, solution$item_num.2))
  colnames(blocks) <- c("id.1", "id.2")
  K <- nrow(Sigma)
  n_blocks_total <- nrow(blocks)
  solution <- solution %>%
    mutate(a.1 = rowSums(select(., starts_with("a") & ends_with(".1"))),
           a.2 = rowSums(select(., starts_with("a") & ends_with(".2"))))
  
  A.matrix.1 <- solution[,c("a1.1", "a2.1", "a3.1", "a4.1", "a5.1")]
  A.matrix.2 <- solution[,c("a1.2", "a2.2", "a3.2", "a4.2", "a5.2")]
  A.matrix <- as.matrix(A.matrix.2 - A.matrix.1) 
  colnames(A.matrix) <- NULL
  d1 <- -solution$a.1*solution$d1.1
  d2 <- -solution$a.2*solution$d1.2
  d <- d2-d1
  
  #set.seed(789)
  data.FC <- mirt::simdata(a = A.matrix, 
                           d = matrix(d, ncol = 1), 
                           N = 5000,
                           itemtype = "2PL",
                           sigma = Sigma)
  
  mod.specs <- rep(NA, K)
  for(k in 1:K){
    mod.specs[k] <- paste0("F", k," = ", paste(c(which(solution$trait.1 == k), which(solution$trait.2 == k)), collapse = ", "))
  }
  mod.specs[ncol(Sigma)+1] <- paste0("COV = ", paste(paste0("F",1:ncol(Sigma)), collapse = "*"))
  mirt.model_FC <- mirt.model(paste(mod.specs, collapse = "\n"))
  
  # Estimated model
  mod_est <- mirt(data.FC, model = mirt.model_FC, itemtype = "2PL", method = "MHRM", SE = FALSE, verbose = FALSE, technical = list(NCYCLES = 500))
  #param_est <- as.data.frame(coef(mod_est, simplify = TRUE)$items)
  
  thetas_est_com <- fscores(mod_est, response.pattern = data.FC, QMC = T, full.scores = TRUE) # estimated thetas/factor scores 
  empirical_relia <- empirical_rxx(thetas_est_com)
  mean_empirical_relia <- mean(empirical_relia)
  names(mean_empirical_relia) <- "mean_empirical_relia"
  
  res <- c(
    empirical_relia, 
    mean_empirical_relia
  ) 
  
  return(list(res = res))
}