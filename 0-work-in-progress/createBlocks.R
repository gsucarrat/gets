##dgp:
set.seed(123)
y <- rnorm(20)
x <- coredata(iim(y))

##function that create blocks (same algorithm as in isat() and
##blocksFun()):
createBlocks <- function(y, x, no.of.blocks=NULL, max.block.size=30,
  ratio.threshold=0.8, keep=NULL)
{
  ##initiate:
  NROWx <- NROW(x)
  NCOLx <- NCOL(x)

  ##determine no. of blocks:
  if( is.null(no.of.blocks) ){
    blockratio.value <- NCOLx/(ratio.threshold*NCOLx)
    blocksize.value <-
      NCOLx/min(NROWx*ratio.threshold, max.block.size)
    no.of.blocks <- max(2,blockratio.value,blocksize.value)
    no.of.blocks <- ceiling(no.of.blocks)
    no.of.blocks <- min(NCOLx, no.of.blocks) #ensure blocks < NCOL
  }

  ##make partitions:
  blocksize <- ceiling(NCOLx/no.of.blocks)
  partitions.t2 <- blocksize
  for(j in 1:no.of.blocks){
    if( blocksize*j <= NCOLx ){
      partitions.t2[j] <- blocksize*j
    }
  }
  ##check if last block contains last regressor:
  if( partitions.t2[ length(partitions.t2) ] < NCOLx ){
    partitions.t2 <- c(partitions.t2, NCOLx)
  }
  blocksadj <- length(partitions.t2)
  partitions.t1 <- partitions.t2 + 1
  partitions.t1 <- c(1, partitions.t1[ -blocksadj ])

  ##make the list w/blocks:
  blocks <- list()
  for(i in 1:blocksadj){
    blocks[[i]] <- partitions.t1[i]:partitions.t2[i]
  }

  ##add keep entries to each block:
  if( !is.null(keep) ){
    for(i in 1:length(blocks)){
      blocks[[i]] <- sort(union(keep, blocks[[i]]))
    }
  }

  #return:
  return(blocks)

} #end createBlocks()