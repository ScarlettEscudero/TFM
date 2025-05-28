#NHBSAFunctions.R

######################################################################################################
###################################         initpop         ##########################################
######################################################################################################

initpop <- function (opt,pop,x){
  
  # Function: pop = initpop(opt, pop, varargin)
  # Description: Initialize population.
  # Syntax:
  #   pop = initpop(opt, pop, varargin)
  #     (default) Create a random initial population by permuting vector varargin.
  #
  # Parameters: 
  #   pop : an empty population
  #
  #*************************************************************************
  initpopPermute <- function(opt,pop,x){
    
    nVar = opt$numVar;
    
    popsize = opt$popsize;
    
    
    for( i in 1:popsize){
      var <- sample(x,nVar)
      pop[[i]]$var = var
    }
    return(pop)
  }
  
  #*************************************************************************
  # 1. Identify parameters
  #*************************************************************************
  method <- NULL
  if(is.data.frame(x)){
    method <- "intro"
  } else if(is.numeric(x)){
    method = 'permute'
  } else if(is.null(method)){
    stop("Invalid pop initialization argument.")
  }
  
  
  #*************************************************************************
  # 2. Initialize population
  #*************************************************************************
  
  if(method == 'permute'){
    pop <- initpopPermute(opt, pop, x)
  } else if(method == "intro"){
    for(i in 1:opt$popsize){
      pop[[i]]$var = as.numeric(x[i,])
    }
  }
  
  return(pop)
}

######################################################################################################
###################################        blocksInfo        ##########################################
######################################################################################################

lik2TIRT <- function(blockData, Sigma){
  
  nTraits <- ncol(Sigma)
  
  dpar <- blockData[,"d1.2"]-blockData[,"d1.1"]
  
  A.1 <- blockData[,paste0("a",1:nTraits,".1")]
  A.2 <- blockData[,paste0("a",1:nTraits,".2")]
  
  Lambda.1<-(A.1/1.702)/sqrt(1 + diag((A.1/1.702) %*% Sigma %*% t((A.1/1.702))))
  Lambda.2<-(A.2/1.702)/sqrt(1 + diag((A.2/1.702) %*% Sigma %*% t((A.2/1.702))))
  
  psi.1 <- diag(diag(1 - Lambda.1 %*% Sigma %*% t(Lambda.1)))
  psi.2 <- diag(diag(1 - Lambda.2 %*% Sigma %*% t(Lambda.2)))
  
  sigma.1 <- Lambda.1 %*% Sigma %*% t(Lambda.1) + psi.1
  sigma.2 <- Lambda.2 %*% Sigma %*% t(Lambda.2) + psi.2
  
  psi <- psi.1 + psi.2
  
  Lambda.fc <- Lambda.2-Lambda.1
  
  sigma <- Lambda.fc %*% Sigma %*% t(Lambda.fc) + psi
  
  # Precompute the scaling vector
  scaling_vector <- 1 / sqrt(diag(sigma))
  
  # Apply scaling to each row of Lambda.fc
  Lambda.fc.std <- Lambda.fc * scaling_vector
  
  uniqueness.sqrt.std <- sqrt(1-diag(Lambda.fc.std %*% Sigma %*% t(Lambda.fc.std)))
  
  A.fc <- 1.702*Lambda.fc.std/uniqueness.sqrt.std 
  d.fc <- dpar/uniqueness.sqrt.std
  
  pars.TIRT <- cbind(A.fc, d.fc)
  colnames(pars.TIRT) <- c(paste0("a",1:nTraits), "d")
  
  return(pars.TIRT)
}

######################################################################################################
###################################        blocksInfo        ##########################################
######################################################################################################


blocksInfo <- function(blockData, Sigma, theta, TIRT = TRUE){
  
  nTraits <- ncol(Sigma)
  
  nBlocks <- nrow(blockData)
  
  if(TIRT){
    pars.TIRT <- lik2TIRT(blockData, Sigma)
    
    dpar <- pars.TIRT[,"d"]
    
    A <- pars.TIRT[,paste0("a",1:nTraits)]
    
  } else {
    
    dpar <- blockData[,"d1.2"]-blockData[,"d1.1"]
    
    A <- matrix(blockData[,paste0("a",1:nTraits,".2")]-blockData[,paste0("a",1:nTraits,".1")], nrow=nBlocks)
    
  }
  
  I <- array(NA,dim = c(nTraits,nTraits,nBlocks))
  for(i in 1:nBlocks){
    I[,,i] <- A[i,]%*%t(A[i,])
  }
  
  P  <- t(1/(1+exp(-(A%*%t(theta)+dpar))))
  PQ <- P*(1-P);
  
  blockInfo <- array(NA,dim=c(nTraits,nTraits,nrow(theta),nBlocks))
  for(i in 1:nBlocks){
    for(j in 1:nrow(theta)){
      blockInfo[,,j,i] = I[,,i]*PQ[j,i]
    }
  }
  return(blockInfo)
}

######################################################################################################
###################################       Aoptimality      ##########################################
######################################################################################################

Aoptimality <- function(selectedBlockInfo, priorInfo){
  
  nTraits <- ncol(priorInfo[,,1])
  testInfo = apply(selectedBlockInfo,c(1:3),sum)
  
  BayesTestInfo = array(sapply(1:dim(testInfo)[3],function(i){
    testInfo[,,i]+priorInfo[,,i]
  }),dim=dim(testInfo))
  
  A <- apply(BayesTestInfo,3,function(x) sum(diag(solve(x))))
  
  return(A)
}

######################################################################################################
###############################  Average Marginal Reliability  #######################################
######################################################################################################

AvgMargReli <- function(selectedBlockInfo, priorInfo, eucWeights){
  
  nTraits <- ncol(priorInfo[,,1])
  testInfo = apply(selectedBlockInfo,c(1:3),sum)
  
  BayesTestInfo = array(sapply(1:dim(testInfo)[3],function(i){
    testInfo[,,i]+priorInfo[,,i]
  }),dim=dim(testInfo))
  
  ThetaVars <- apply(BayesTestInfo,3,function(x) (diag(solve(x))))
  
  A <- mean(1-t(eucWeights)%*%t(ThetaVars))
  
  return(A)
}
######################################################################################################
###################################       Doptimality      ##########################################
######################################################################################################

Doptimality <- function(selectedBlockInfo, priorInfo){
  
  nTraits <- ncol(priorInfo[,,1])
  testInfo = apply(selectedBlockInfo,c(1:3),sum)
  
  BayesTestInfo = array(sapply(1:dim(testInfo)[3],function(i){
    testInfo[,,i]+priorInfo[,,i]
  }),dim=dim(testInfo))
  
  D <- apply(BayesTestInfo,3,det)
  
  return(D)
}

######################################################################################################
###################################       Eoptimality      ##########################################
######################################################################################################

Eoptimality <- function(selectedBlockInfo, priorInfo){
  
  nTraits <- ncol(priorInfo[,,1])
  testInfo = apply(selectedBlockInfo,c(1:3),sum)
  
  BayesTestInfo = array(sapply(1:dim(testInfo)[3],function(i){
    testInfo[,,i]+priorInfo[,,i]
  }),dim=dim(testInfo))
  
  E <- apply(BayesTestInfo,3,function(x) min(eigen(x)$values))
  
  return(E)
}

######################################################################################################
###################################        minVar       ##########################################
######################################################################################################

minVar <- function(selectedBlockInfo, priorInfo){
  
  nTraits <- ncol(priorInfo[,,1])
  testInfo = apply(selectedBlockInfo,c(1:3),sum)
  
  BayesTestInfo = array(sapply(1:dim(testInfo)[3],function(i){
    testInfo[,,i]+priorInfo[,,i]
  }),dim=dim(testInfo))
  
  mV <- apply(BayesTestInfo,3,function(x) (diag(solve(x))))
  
  return(mV)
}

######################################################################################################
###################################          objfun         ##########################################
######################################################################################################

objfun <- function(x,varargin){
  # Objective function : Test problem 'BLOCK ASSEMBLY'.
  #*************************************************************************
  itemData   <- varargin$itemData
  Sigma      <- varargin$Sigma
  criterion  <- varargin$criterion
  itemIdx    <- varargin$itemIdx
  priorInfo   <- varargin$priorInfo
  theta      <- varargin$theta
  eucWeights <- varargin$eucWeights
  simpMethod <- varargin$simpMethod
  TIRT <- varargin$TIRT
  
  selected_blocks<-selectBlocks(x,itemData,itemIdx)
  
  selectedBlockInfo <- blocksInfo(selected_blocks,Sigma,theta, TIRT)
  
  if(criterion == 'Aoptimality'){
    
    #Aoptimality
    #*************************************************************************
    obj = Aoptimality(selectedBlockInfo, priorInfo)
    
  } else if(criterion == 'AvgMargReli'){
    
    #Average Marginal Reliability 
    #*************************************************************************
    obj = AvgMargReli(selectedBlockInfo, priorInfo, eucWeights)
    
  } else if(criterion == 'Eoptimality'){
    
    #Eoptimality
    #*************************************************************************
    obj = -Eoptimality(selectedBlockInfo, priorInfo)
    
  } else if(criterion == 'Doptimality'){
    
    #Doptimality
    #*************************************************************************
    obj = -Doptimality(selectedBlockInfo, priorInfo)
    
  } else if(criterion == 'minVar'){
    
    #Doptimality
    #*************************************************************************
    obj = minVar(selectedBlockInfo, priorInfo)
    
  }
  
  if(simpMethod == "Boundary"){
    
    y <- max(obj)
    
  } else if(simpMethod == "Average"){
    
    y <- sum(obj*eucWeights)
    
  }
  
  return(list(y=y))
}

######################################################################################################
###################################         evaluate        ##########################################
######################################################################################################

evaluate <- function(opt, pop, state, varargin){
  # Function: [pop, state] = evaluate(opt, pop, state, varargin)
  # Description: Evaluate the objective functions of each individual in the
  #   population.
  #
  #*************************************************************************
  
  evalIndividual <- function(indi, varargin){
    # Function: [indi, evalTime] = evalIndividual(indi, varargin)
    # Description: Evaluate one objective function.
    
    idx <- NULL
    
    start.time <- Sys.time()
    
    objres <- objfun(indi$var, varargin)
    
    y <- objres$y
    # cons <- objres$cons
    
    end.time <- Sys.time()
    
    evalTime <- difftime(end.time,start.time)
    
    # Save the objective values and constraint violations
    indi$obj <- y;
    
    return(list(indi = indi, evalTime = evalTime))
  }
  
  N = length(pop);
  allTime = rep(NA,N);  # allTime : use to calculate average evaluation times
  
  
  if(is.null(opt$nCores)){
    #*************************************************************************
    # Evaluate objective function in serial
    #*************************************************************************
    for(i in 1:N){
      evalIndi <- evalIndividual(pop[[i]], varargin)
      pop[[i]] <- evalIndi$indi
      allTime[i] <- evalIndi$evalTime
      
    }
  } else {
    #*************************************************************************
    # Evaluate objective function in parallel
    #*************************************************************************
    comb <- function(x, ...) {
      lapply(seq_along(x),
             function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
    }
    
    evalIndi = foreach(i = 1:N,
                       .combine='comb', .multicombine=TRUE, .export = envir,
                       .init=list(list(), list())) %dopar% {
                         try({
                           evalIndividual(pop[[i]], varargin)
                         })
                       }
    
    pop <- evalIndi[[1]]
    allTime <- unlist(evalIndi[[2]])
  }
  
  #*************************************************************************
  # Statistics
  #*************************************************************************
  state$avgEvalTime   = sum(allTime) / length(allTime);
  state$evaluateCount = state$evaluateCount + length(pop);
  
  
  return(list(pop = pop, state = state))
}

######################################################################################################
###################################        Node Histogram         ##########################################
######################################################################################################

NH <- function(nVar,parent,fraction,pop,condProbs,Bratio,tabvars,nBlocks,varargin){
  
  nTraits <- varargin$nTraits
  itemData <- varargin$itemData
  itemIdx <- varargin$itemIdx
  nHetero <- varargin$nHetero
  nPool <- nrow(itemData)
  
  condProbs2 <- condProbs
  trait_combinations <- t(combn(nTraits,2))
  n_combin <- nrow(trait_combinations)
  
  dir_cons <- varargin$which.combs$dir_cons
  inv_cons <- varargin$which.combs$inv_cons
  
  tour  <- parent$var
  failed <- 0
  check_list <- 0
  offspring = rep(NA,nVar)
  
  
  if (varargin$balance_dim_comb == TRUE){
    n_direct_by_trait <- rep(0,n_combin)
    n_inverse_by_trait <- rep(0,n_combin)
    
    if(fraction > 0){
      cuts <- sort(sample(nVar,fraction))
      firstCut <- sample(1:(fraction-1),1)
      
      template <- cuts[firstCut]:cuts[firstCut+1]
      
      # copy tour into offspring and add mirror values/indices
      for (k in template){
        probs <- condProbs2[,k]
        
        paste(sort(c(itemData[k,"trait"], itemData[tour[k],"trait"])), collapse = ".")
        
        which.comb <- as.numeric(factor(paste(sort(c(itemData[k,"trait"], itemData[tour[k],"trait"])), collapse = "."),
                                        levels = apply(cbind(trait_combinations[,1], trait_combinations[,2]), 1, function(x) paste(sort(x), collapse = "."))))
        if((itemData[k,3]+itemData[tour[k],3]) != 0){
          if(probs[tour[k]]==1 & (n_direct_by_trait[which.comb] < dir_cons[which.comb])){
            offspring[k] = tour[k]
            offspring[offspring[k]] = k
            
            n_direct_by_trait[which.comb] <- n_direct_by_trait[which.comb] + 1
            condProbs2[k,] <- condProbs2[,k] <- 0
            condProbs2[offspring[k],] <- condProbs2[,offspring[k]] <- 0
            
            for(combs in which((n_direct_by_trait - dir_cons)>=0)){
              condProbs2[itemData[,2] == trait_combinations[combs,2] & itemData[,3] == 1,itemData[,2] == trait_combinations[combs,1] & itemData[,3] == 1] <- 0
              condProbs2[itemData[,2] == trait_combinations[combs,1] & itemData[,3] == 1,itemData[,2] == trait_combinations[combs,2] & itemData[,3] == 1] <- 0
              condProbs2[itemData[,2] == trait_combinations[combs,2] & itemData[,3] == -1,itemData[,2] == trait_combinations[combs,1] & itemData[,3] == -1] <- 0
              condProbs2[itemData[,2] == trait_combinations[combs,1] & itemData[,3] == -1,itemData[,2] == trait_combinations[combs,2] & itemData[,3] == -1] <- 0
            }
          }
        } else if((itemData[k,3]+itemData[tour[k],3])==0){
          if(probs[tour[k]]==1 & (n_inverse_by_trait[which.comb] < inv_cons[which.comb])){
            offspring[k] = tour[k]
            offspring[offspring[k]] = k
            
            n_inverse_by_trait[which.comb] <- n_inverse_by_trait[which.comb] + 1
            condProbs2[k,] <- condProbs2[,k] <- 0
            condProbs2[offspring[k],] <- condProbs2[,offspring[k]] <- 0
            
            for(combs in which((n_inverse_by_trait - inv_cons)>=0)){
              condProbs2[itemData[,2] == trait_combinations[combs,2] & itemData[,3] == 1,itemData[,2] == trait_combinations[combs,1] & itemData[,3] == -1] <- 0
              condProbs2[itemData[,2] == trait_combinations[combs,1] & itemData[,3] == 1,itemData[,2] == trait_combinations[combs,2] & itemData[,3] == -1] <- 0
              condProbs2[itemData[,2] == trait_combinations[combs,2] & itemData[,3] == -1,itemData[,2] == trait_combinations[combs,1] & itemData[,3] == 1] <- 0
              condProbs2[itemData[,2] == trait_combinations[combs,1] & itemData[,3] == -1,itemData[,2] == trait_combinations[combs,2] & itemData[,3] == 1] <- 0
              
            }
          }
        }
      }
    }
    
    while(any(is.na(offspring))){
      
      if(any(rowSums(condProbs2)>0)){
        k <- ifelse(sum(rowSums(condProbs2)>0)>1,sample(which(rowSums(condProbs2)>0),1),which(rowSums(condProbs2)>0))
        tab <- tabvars[,k]
        epsi <- length(pop)/(nPool/(2*nBlocks))/sum(condProbs[k,])*Bratio
        probs <- (tab+epsi)*condProbs2[k,]
        
        probs[offspring[!is.na(offspring)]] <- probs[which(!is.na(offspring))] <- 0
        
        if(sum(probs)>0){
          
          probs <- probs/sum(probs)      
          
          offspring[k] <- which(rmultinom(n=1,size=1,prob=probs)==1)
          offspring[offspring[k]] = k
        } else {
          offspring[k] = k
        }
        
        if((itemData[k,"pol"]+itemData[offspring[k],"pol"]) != 0){
          n_direct_by_trait[as.numeric(factor(paste(sort(c(itemData[k,"trait"], itemData[offspring[k],"trait"])), collapse = "."),
                                              levels = apply(cbind(trait_combinations[,1], trait_combinations[,2]), 1, function(x) paste(sort(x), collapse = "."))))] <- 
            n_direct_by_trait[as.numeric(factor(paste(sort(c(itemData[k,"trait"], itemData[offspring[k],"trait"])), collapse = "."),
                                                levels = apply(cbind(trait_combinations[,1], trait_combinations[,2]), 1, function(x) paste(sort(x), collapse = "."))))] + 1
        } else if((itemData[k,"pol"]+itemData[offspring[k],"pol"]) == 0){
          n_inverse_by_trait[as.numeric(factor(paste(sort(c(itemData[k,"trait"], itemData[offspring[k],"trait"])), collapse = "."),
                                               levels = apply(cbind(trait_combinations[,1], trait_combinations[,2]), 1, function(x) paste(sort(x), collapse = "."))))] <- 
            n_inverse_by_trait[as.numeric(factor(paste(sort(c(itemData[k,"trait"], itemData[offspring[k],"trait"])), collapse = "."),
                                                 levels = apply(cbind(trait_combinations[,1], trait_combinations[,2]), 1, function(x) paste(sort(x), collapse = "."))))] + 1
        }
        
        print(n_direct_by_trait)
        print(n_inverse_by_trait)
        condProbs2[k,] <- condProbs2[,k] <- 0
        condProbs2[offspring[k],] <- condProbs2[,offspring[k]] <- 0
        
        for(combs in which((n_direct_by_trait - dir_cons)>=0)){
          condProbs2[itemData[,"trait"] == trait_combinations[combs,2] & itemData[,"pol"] == 1,itemData[,"trait"] == trait_combinations[combs,1] & itemData[,"pol"] == 1] <- 0
          condProbs2[itemData[,"trait"] == trait_combinations[combs,1] & itemData[,"pol"] == 1,itemData[,"trait"] == trait_combinations[combs,2] & itemData[,"pol"] == 1] <- 0
          condProbs2[itemData[,"trait"] == trait_combinations[combs,2] & itemData[,"pol"] == -1,itemData[,"trait"] == trait_combinations[combs,1] & itemData[,"pol"] == -1] <- 0
          condProbs2[itemData[,"trait"] == trait_combinations[combs,1] & itemData[,"pol"] == -1,itemData[,"trait"] == trait_combinations[combs,2] & itemData[,"pol"] == -1] <- 0
          
        }
        
        for(combs in which((n_inverse_by_trait - inv_cons)>=0)){
          condProbs2[itemData[,"trait"] == trait_combinations[combs,2] & itemData[,"pol"] == 1,itemData[,"trait"] == trait_combinations[combs,1] & itemData[,"pol"] == -1] <- 0
          condProbs2[itemData[,"trait"] == trait_combinations[combs,1] & itemData[,"pol"] == 1,itemData[,"trait"] == trait_combinations[combs,2] & itemData[,"pol"] == -1] <- 0
          condProbs2[itemData[,"trait"] == trait_combinations[combs,2] & itemData[,"pol"] == -1,itemData[,"trait"] == trait_combinations[combs,1] & itemData[,"pol"] == 1] <- 0
          condProbs2[itemData[,"trait"] == trait_combinations[combs,1] & itemData[,"pol"] == -1,itemData[,"trait"] == trait_combinations[combs,2] & itemData[,"pol"] == 1] <- 0
        }
      } else if(all(rowSums(condProbs2)==0)){
        offspring[is.na(offspring)] <- which(is.na(offspring))
      }
    }
  } else { #varargin$balance_dim_comb == FALSE 
    n_direct_by_trait <- 0
    n_inverse_by_trait <- 0
    
    if(fraction > 0){
      cuts <- sort(sample(nVar,fraction))
      firstCut <- sample(1:(fraction-1),1)
      
      template <- cuts[firstCut]:cuts[firstCut+1]
      
      # copy tour into offspring and add mirror values/indices
      for (k in template){
        probs <- condProbs2[,k]
        
        if((itemData[k,3]+itemData[tour[k],3]) != 0){
          if(probs[tour[k]]==1 & (n_direct_by_trait< dir_cons)){
            offspring[k] = tour[k]
            offspring[offspring[k]] = k
            
            n_direct_by_trait <- n_direct_by_trait + 1
            condProbs2[k,] <- condProbs2[,k] <- 0
            condProbs2[offspring[k],] <- condProbs2[,offspring[k]] <- 0
            
            for(combs in which((n_direct_by_trait - dir_cons)>=0)){
              condProbs2[itemData[,3] == -1,itemData[,3] == -1] <- 0
              condProbs2[itemData[,3] ==  1,itemData[,3] ==  1] <- 0
            }
          }
        } else if((itemData[k,3]+itemData[tour[k],3])==0){
          if(probs[tour[k]]==1 & (n_inverse_by_trait < inv_cons)){
            offspring[k] = tour[k]
            offspring[offspring[k]] = k
            
            n_inverse_by_trait <- n_inverse_by_trait + 1
            condProbs2[k,] <- condProbs2[,k] <- 0
            condProbs2[offspring[k],] <- condProbs2[,offspring[k]] <- 0
            
            for(combs in which((n_inverse_by_trait - inv_cons)>=0)){
              condProbs2[itemData[,3] == 1,itemData[,3] == -1] <- 0
              condProbs2[itemData[,3] == -1,itemData[,3] == 1] <- 0
              
            }
          }
        }
      }
    }
    
    while(any(is.na(offspring))){
      
      if(any(rowSums(condProbs2)>0)){
        k <- ifelse(sum(rowSums(condProbs2)>0)>1,sample(which(rowSums(condProbs2)>0),1),which(rowSums(condProbs2)>0))
        tab <- tabvars[,k]
        epsi <- length(pop)/(nPool/(2*nBlocks))/sum(condProbs[k,])*Bratio
        probs <- (tab+epsi)*condProbs2[k,]
        
        probs[offspring[!is.na(offspring)]] <- probs[which(!is.na(offspring))] <- 0
        
        if(sum(probs)>0){
          
          probs <- probs/sum(probs)      
          
          offspring[k] <- which(rmultinom(n=1,size=1,prob=probs)==1)
          offspring[offspring[k]] = k
        } else {
          offspring[k] = k
        }
        
        if((itemData[k,"pol"]+itemData[offspring[k],"pol"]) != 0){
          n_direct_by_trait <- 
            n_direct_by_trait + 1
        } else if((itemData[k,"pol"]+itemData[offspring[k],"pol"]) == 0){
          n_inverse_by_trait <- 
            n_inverse_by_trait + 1
        }
        
        print(n_direct_by_trait)
        print(n_inverse_by_trait)
        condProbs2[k,] <- condProbs2[,k] <- 0
        condProbs2[offspring[k],] <- condProbs2[,offspring[k]] <- 0
        
        for(combs in which((n_direct_by_trait - dir_cons)>=0)){
          condProbs2[itemData[,3] == -1,itemData[,3] == -1] <- 0
          condProbs2[itemData[,3] ==  1,itemData[,3] == 1] <- 0
        }
        
        for(combs in which((n_inverse_by_trait - inv_cons)>=0)){
          condProbs2[itemData[,3] == 1,itemData[,3] == -1] <- 0
          condProbs2[itemData[,3] == -1,itemData[,3] == 1] <- 0
        }
      } else if(all(rowSums(condProbs2)==0)){
        offspring[is.na(offspring)] <- which(is.na(offspring))
      }
    }
  }
  
  return(offspring)
}

######################################################################################################
###################################      crossoverOp        ##########################################
######################################################################################################

crossoverOp <- function(opt, pop, state, fraction, varargin){
  
  # Function: pop = crossoverOp(opt, pop, state)
  # Description: Crossover operator. All of the individuals would be do crossover, but
  #   only "crossoverFraction" of design variables of an individual would changed.
  #*************************************************************************
  
  crsNH <- function(parent, fraction, pop, condProbs, varargin,tabvars,opt){
    # Function: [child3, child2] = crsPMX(parent1, parent2, fraction)
    # Description: (For integrer coding) PMX mirrored crossover
    #
    #*************************************************************************
    
    child = parent
    
    nVar = length(parent$var);
    
    mirror.vars <- NH(nVar,parent,fraction,pop,condProbs,opt$Bratio,tabvars,opt$nBlocks, varargin)
    
    child$var <- mirror.vars
    
    return(child)
  }
  
  #*************************************************************************
  # 1. Check for the parameters
  #*************************************************************************
  # determine the crossover method
  
  crossOp <- opt$crossOp
  nBlocks <- opt$nBlocks
  
  nVar = opt$numVar;
  Bratio <- opt$Bratio
  
  if(crossOp=="NH"){
    condProbs <- opt$condProbs
    vars <- t(sapply(1:length(pop),function(i) pop[[i]]$var))
    tabvars <- matrix(0,nVar,nVar)
    
    for(k in 1:nVar){
      bla <- table(vars[,k])
      tabvars[as.numeric(names(bla)),k] <- bla
    }
    
    #Serial implementation
    if(is.null(opt$nCores)){
      newind <- sample(1:length(pop))
      for (ind in 1:length(pop)){ # Popsize should be even number
        # Create children
        child  <- crsNH(pop[[ind]], fraction, pop, condProbs, varargin, tabvars, opt)
        
        pop[[newind[ind]]] <- child
      }
    } else {
      #Parallel implementation
      newind <- sample(1:length(pop))
      pop[newind] <- foreach(j = 1:length(pop),
                             .export = c(envir)) %dopar% {
                               try({
                                 crsNH(pop[[j]], fraction, pop, condProbs, varargin, tabvars,opt)
                               })
                             }
      # foreach(j = 1:(opt$nCores*2)) %dopar% {
      #   rm(list = ls())
      #   gc()
      # }
    }
    
  }
  
  return(pop)
}

envir <- ls(environment())