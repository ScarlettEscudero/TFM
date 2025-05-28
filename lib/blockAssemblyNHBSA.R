library(mvtnorm)

source("lib/selectBlocks.R")
source("lib/NHBSAFunctions.R")
source("lib/NHBSA.R")

blockAssemblyNHBSA <- function(criterion, simpMethod, nBlocks, nHetero, itemData, Sigma, Cmatrix, priorInfo = NULL, TIRT = TRUE, nCores = 2, balance_dim_comb, which.combs = NULL, theta){
  
  nTraits <- ncol(Sigma)
  nQuads <-  500 # number of quadrature points for computing the objective function.
  nPool <- nrow(itemData)
  nGen <- 30 # maximum number of generations to converge.
  popsize <- nPool # size of the populations in each generation.
  Bratio <- .0002 # random error added to crossover
  nCutpoints <- 2 # number of splits to paste equal genes from parent to child
  
  possibleBlocks <- combn(1:nPool,2)
  
  nTraitComb <- nrow(t(combn(nTraits,2)))
  
  #if(nHetero/nTraitComb != round(nHetero/nTraitComb) & is.null(which.combs)) {stop("Argument which.combs must be defined.")} else {
  #  if(sum(which.combs$inv_cons) != nHetero) {stop("Argument which.combs does not properly corrrespond to nHetero.")}
  #}
  
  
  trait_combinations <- t(combn(nTraits,2))
  n_combin <- nrow(trait_combinations)

 if(balance_dim_comb == TRUE & !is.null(which.combs)){
   which.combs <- list(dir_cons = rep((nBlocks-nHetero)/n_combin,n_combin),
                       inv_cons = rep(nHetero/n_combin,n_combin))
 } else if (balance_dim_comb == FALSE & is.null(which.combs)) {
   which.combs <- list(dir_cons = nBlocks - nHetero,
                       inv_cons = nHetero)
 } else if (balance_dim_comb == FALSE & !is.null(which.combs)) {
   which.combs <- which.combs
 } else {
   which.combs <- list(dir_cons = nBlocks - nHetero,
                       inv_cons = nHetero)
 }
  
  if(theta == "Quad_without_Sigma"){
    theta <- mirt:::QMC_quad(nQuads,nTraits)
  } else if (theta == "Quad_with_Sigma") {
    theta_uncorr <- mirt:::QMC_quad(nQuads,nTraits)
    Sigma_chol_matrix <- chol(Sigma)
    theta <- theta_uncorr %*% Sigma_chol_matrix
  }
  else {
    theta <- theta
  }
  
  condProbs <- Cmatrix
  
  opt <- NULL
  opt$popsize <- popsize                   # population size
  opt$maxGen  <- nGen                # max generation
  opt$nBlocks <- nBlocks
  opt$numVar <- nrow(itemData)        # number of design variables
  opt$initfun <-  c(1:nrow(itemData))
  opt$crossoverFraction <- nCutpoints
  opt$condProbs <- condProbs
  opt$crossOp <- "NH"
  opt$Bratio <- Bratio
  opt$nCores <- nCores
  
  itemIdx <- 1:nPool
  #eucWeights <- dmvnorm(theta, sigma = Sigma)
  #eucWeights <- eucWeights/sum(eucWeights)
  eucWeights <- rep(1, nrow(theta))
  eucWeights <- eucWeights/sum(eucWeights)
  
  if(is.null(priorInfo)){
    priorInfo <- array(solve(Sigma), dim = c(nTraits, nTraits, nQuads))
  }
  
  varargin <- list(itemData = itemData, Sigma = Sigma, criterion = criterion,
                   itemIdx = itemIdx, nBlocks = nBlocks, nHetero = nHetero, nTraits = nTraits,
                   priorInfo = priorInfo, theta = theta, eucWeights = eucWeights, simpMethod = simpMethod, which.combs = which.combs, TIRT = TIRT, balance_dim_comb = balance_dim_comb)
  
  result <- nhbsa(opt, varargin)                # begin the optimization!
  
  lastPop <- result$pop
  opt <- result$opt
  
  lastGen <- result$state$currentGen
  lastGenChange <- result$state$repSolution
  
  lastObj <- t(sapply(lastPop,function(i){
    c(i$obj)
  }))
  
  which.best <- which.max(lastObj)
  solution <- lastPop[[which.best]]$var
  
  totalTime <- result$state$totalTime
  
  selectedBlocks <- as.matrix(selectBlocks(solution,itemData,itemIdx))
  
  selectedBlockInfoList <- blocksInfo(selectedBlocks,Sigma,theta, TIRT)
  
  selectedBlockInfo <- selectedBlockInfoList
  
  objVal <- objfun(solution, varargin)$y
  
  reliab <- 1-colSums(t(minVar(selectedBlockInfoList, priorInfo = priorInfo))*eucWeights)
  
  write.table(selectedBlocks,
              file = paste0('out/SelectedBlocks_',
                            criterion,'_',simpMethod,'_',nBlocks,'_',nHetero,'_',nPool,'.dat'),
              quote=F,row.names = F,col.names = T)
  
  write.table(solution,
              file = paste0('out/Solution_',
                            criterion,'_',simpMethod,'_',nBlocks,'_',nHetero,'_',nPool,'.dat'),
              quote=F,row.names = F,col.names = F)
  
  return(list(blocks = selectedBlocks, reliability = reliab))
}

