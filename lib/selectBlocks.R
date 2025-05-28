#selectBlocks.R

selectBlocks <- function(x,itemData,itemIdx){
  
  nTraits = max(itemData[,"trait"])

  pairedBlocks <- NULL
  for(i in 1:length(x)){
    pairedBlocks <- rbind(pairedBlocks,c(t(itemData[c(itemIdx[i],x[i]),]),
      abs(itemData[itemIdx[i],"trait"] - itemData[x[i],"trait"])>0))
  }

  bidim <- pairedBlocks[,ncol(pairedBlocks)]
  
  selectedBlocks<-pairedBlocks[bidim == 1, -ncol(pairedBlocks)]
  colnames(selectedBlocks) <- c(paste0(colnames(itemData),".1"),paste0(colnames(itemData),".2"))
  
  selected_blocks<-selectedBlocks[!duplicated(data.frame(list(do.call(pmin,data.frame(selectedBlocks[,c("id.1","id.2")])),do.call(pmax,data.frame(selectedBlocks[,c("id.1","id.2")]))))),]
  
  return(selected_blocks)
  
}
