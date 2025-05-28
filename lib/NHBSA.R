library("foreach")
library("doParallel")
library("snow")
library("abind")

nhbsa <- function(opt, varargin){
  
  start.Totaltime <- Sys.time()
  
  #*************************************************************************
  nVar    = opt$numVar;
  nObj    = opt$numObj;
  popsize = opt$popsize;
  fraction = opt$crossoverFraction
  opt2 <- opt; opt2$Bratio <- 10000
  nCores <- opt$nCores
  
  if(!is.null(nCores)){
    cl <- makeCluster(nCores,type="SOCK", outfile="")
    registerDoParallel(cl,cores=nCores)
  }
  
  pop <- rep(list(list(var = rep(0,nVar), obj = 0)),popsize)
  
  varreg <- NULL
  # state: optimization state of one generation
  
  state <- list(currentGen = 1, evaluateCount = 0, totalTime = 0, avgEvalTime = 0, repSolution = 0)
  
  #result <- list(pops = rep(list(pop),opt$maxGen), states = rep(list(state),opt$maxGen), opt = opt)
  # each row is the population of one generation
  # each row is the optimizaiton state of one generation
  # use for output
  
  # global variables
  #if STOP_NSGA!=0, then stop the optimizaiton
  assign("STOP_NSGA",0,.GlobalEnv)
  
  
  #*************************************************************************
  # initialize the P0 population
  #*************************************************************************
  
  ngen = 1;
  pop1 <- initpop(opt, pop, opt$initfun);
  if(!is.data.frame(opt$initfun)) pop1 <- crossoverOp(opt2, pop1, state, fraction = 0, varargin = varargin)
  
  pop1.evaluate <- evaluate(opt, pop1, state, varargin)
  
  pop1 <- pop1.evaluate$pop
  state <- pop1.evaluate$state
  
  # state
  state$currentGen = ngen;
  state$totalTime = as.numeric(difftime(Sys.time(),start.Totaltime,units = "sec"))
  
  vars <- sapply(pop1,function(i) i$var)
  vars.coll <- apply(vars, 2, paste, collapse = ".")
  objs <- sapply(pop1,function(i) i$obj)
  old.objs <- sapply(pop1,function(i) i$obj)
  old.solution <- pop1[[which.max(old.objs)]]$var
  
  varreg <- abind(varreg,vars, along = 3)
  
  write.table(x = t(c(ngen,mean(objs),min(objs),as.numeric(difftime(Sys.time(),start.Totaltime,units = "sec")),old.solution)),
              file = paste0('out/AllSolutions_',
                            varargin$criterion,'_',varargin$simpMethod,'_',opt$nBlocks,'_',varargin$nHetero,'_',opt$numVar,'.dat'),
              quote=F,row.names = F,col.names = F)
  
  #*************************************************************************
  # NHBSA iteration
  #*************************************************************************
  
  while( ngen < opt$maxGen & STOP_NSGA==0){
    
    ngen = ngen+1;
    state$currentGen <- ngen;
    
    
    
    if((ngen/50)==round(ngen/50)){
      
      if(!is.null(nCores)) stopCluster(cl)
      
      if(!is.null(nCores)){
        cl <- makeCluster(nCores,type="SOCK")
        registerDoParallel(cl,cores=nCores)
      }
    }
    
    pop2 <- crossoverOp(opt, pop1, state, fraction = opt$crossoverFraction, varargin = varargin)
    
    vars2.coll <- sapply(pop2,function(i) paste(i$var, collapse = "."))
    
    match.pop2.pop1 <- match(vars2.coll, vars.coll)
    # print(sum(is.na(match.pop2.pop1)))
    
    pop2.evaluate <- tryCatch({
      evaluate(opt, pop2[is.na(match.pop2.pop1)], state, varargin)
    }, error = function(e){})
    
    pop2 <- c(pop2.evaluate$pop, pop1[match.pop2.pop1[!is.na(match.pop2.pop1)]])

    N <- length(pop1)
    
    obj1     = matrix(unlist(sapply(pop1,function(n){ n$obj})),N,byrow=T)
    obj2     = matrix(unlist(sapply(pop2,function(n){ n$obj})),N,byrow=T)
    
    
    new.pop <- c(pop1,pop2)[order(c(obj1,obj2), decreasing = TRUE)[1:N]]
    new.obj <- c(obj1,obj2)[order(c(obj1,obj2), decreasing = TRUE)[1:N]]
    
    which.best <- which.max(new.obj)
    new.solution <- new.pop[[which.best]]$var
    
    pop1 <- new.pop
    
    state$totalTime = as.numeric(difftime(Sys.time(),start.Totaltime,units = "sec"))
    
    if(any(old.solution != new.solution)){
      state$repSolution <- ngen
    } 
    
    write.table(t(c((ngen-1)*popsize+which.best,mean(new.obj),max(new.obj),state$totalTime,new.solution)),
                file = paste0('out/AllSolutions_',
                              varargin$criterion,'_',varargin$simpMethod,'_',opt$nBlocks,'_',varargin$nHetero,'_',opt$numVar,'.dat'),
                quote=F,row.names = F,col.names = F,append=T)
    
    old.solution <- new.solution

    vars <- sapply(pop1,function(i) i$var)
    vars.coll <- apply(vars, 2, paste, collapse = ".")
    varsdiv <- apply(vars,1,function(i) length(table(i)))
    
    varreg <- abind(varreg, vars, along = 3)

    if(all(varsdiv==1)){
      STOP_NSGA = 1;
    }
    old.objs <- new.obj
  }
  
  if(!is.null(nCores)) stopCluster(cl)
  
  saveRDS(object = varreg, file = paste0('out/AllVectors_',
                                         varargin$criterion,'_',varargin$simpMethod,'_',opt$nBlocks,'_',varargin$nHetero,'_',opt$numVar,'.rds'))
  
  return(list(pop = pop1, state = state, opt = opt))
  
}
