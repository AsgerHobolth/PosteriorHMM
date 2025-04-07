##-------------------------------------------------------------
## General functions for Finite Markov Chain Imbedding with 2 states:
## - FMCI for the number of state 1 to state 2 transitions
## - FMCI for the number of state 2 positions
## - FMCI for the run length in state 2
## - FMCI for the longest run in state 2
## For each summary statistics the code consists of two functions:
## One function for the transition probability matrix
## One function for calculating the distribution of the summary 
## statistics (using an initial probability vector and the transition 
## probability matrix)
##---------------------------------------------------------------
## Finite Markov chain imbedding:
## Number of transitions from state 1 to state 2
##---------------------------------------------------------------
## Setting up the Markov chain:
## Transition probability matrix for FMCI
## Name: n12TransMatFct
## Input: 
## a (prb of staying in state 1)
## b (prb of staying in state 2)
## nMax (max number of transitions)
## Output: 
## Transition probability matrix for number of transitions
n12TransMatFct <- function( a,b,nMax ){
  TransMat <- matrix( 0,nrow=(2+2*nMax),ncol=(2+2*nMax) )
  ## Transitions from odd numbered states
  for (i in 0:(nMax-1)){
    TransMat[(1+2*i),(1+2*i)] <- b
    TransMat[(1+2*i),(2+2*i)] <- 1-b
  }
  ## Transitions from even numbered states
  for (i in 1:nMax){
    TransMat[(2*i),(2*i)] <- a
    TransMat[(2*i),(1+2*i)] <- 1-a
  }
  ## Final transitions
  TransMat[(1+2*nMax),(2+2*nMax)] <- 1
  TransMat[(2+2*nMax),(2+2*nMax)] <- 1
  return( TransMat )
}
##---------------------------------------------------------------
## Distribution for the number of state 1 to state 2 transitions
## Name: n12PrbFct 
## Input: pi1,pi2,aVec,bVec
## nMax is the maximum number of transitions considered
## Output: 
## Vector of probabilities for the number of transitions
n12PrbFct <- function( pi1,pi2,aVec,bVec,nMax ){
  bgn <- c(pi2,pi1,rep(0,2*nMax))
  prb <- bgn
  for (i in 1:length(aVec)){
    prb <- prb%*%n12TransMatFct( a=aVec[i],b=bVec[i],nMax=nMax )
  }
  prbOut <- colSums(matrix(prb,nrow=2))
  return(prbOut)
}
##-----------------------------------------------------------------
## Finite Markov chain imbedding: 
## Number of positions in state 2
##-----------------------------------------------------------------
## Setting up the Markov chain
## Transition matrix probability for FMCI
## Name: n2TransMatFct
## Input:
## a is the probability of staying in state 1
## b is the probability of staying in state 2
## nMax is the max number of positions in state 2 
## Output:
## Transition probability for the number of positions in state 2
n2TransMatFct <- function( a,b,nMax ){
  TransMat <- matrix(0,nrow=(2+2*nMax),ncol=(2+2*nMax))
  ## Transitions from odd numbered states
  for (i in 0:(nMax-1)){
    TransMat[(1+2*i),(2+2*i)] <- 1-b
    TransMat[(1+2*i),(3+2*i)] <- b
  }
  ## Transitions from even numbered states
  for (i in 1:nMax){
    TransMat[(2*i),(2*i)] <- a
    TransMat[(2*i),(1+2*i)] <- 1-a
  }
  ## Final transitions
  TransMat[(1+2*nMax),(2+2*nMax)] <- 1
  TransMat[(2+2*nMax),(2+2*nMax)] <- 1
  return( TransMat )
}
##-----------------------------------------------
## Distribution for the number of positions in state 2
## Name: n2PrbFct
## Input:
## pi1,pi2,aVec,bVec
## nMax is the max number of positions in state 2 
## Output:
## Vector of probabilities for the number of positions in state 2
n2PrbFct <- function( pi1,pi2,aVec,bVec,nMax ){
  bgn <- c( 0,pi1,pi2,rep(0,2*nMax-1) )
  prb <- bgn
  for (i in 1:length(aVec)){
    prb <- prb%*%n2TransMatFct( a=aVec[i],b=bVec[i],nMax=nMax)
  }
  prbOut <- colSums(matrix(prb,nrow=2))
  return(prbOut)
}
##-------------------------------------------------------------
## Finite Markov chain imbedding:
## Exact run length
##-------------------------------------------------------------
## Setting up the Markov chain
## Transition probability for FMCI
## Name:
## r2TransMatFct
## Input: 
## a is the probability of staying in state 1
## b is the probability of staying in state 2
## mxRuns is the max number of runs
## k is the run length
## Output: 
## Transition probability matrix for exactly k runs
r2TransMatFct <- function( a,b,mxRuns,k){
  blockMat <- matrix(0,nrow=(k+2),ncol=(k+2))
  blockMat[1,1:2] <- c(b,1-b)
  blockMat[2,2:3] <- c(a,1-a)
  if (k>1){
    blockMat[3:(k+1),2] <- rep(1-b,(k-1))
    for (i in (3:(k+1))){
      blockMat[i,(i+1)] <- b
    }
  }
  blockMat[(k+2),1] <- b
  ## Fill out the (mxRuns+1) blocks
  sz <- (mxRuns+1)*(k+2)+1
  TransMat <- matrix(0,nrow=sz,ncol=sz)
  for (i in 0:mxRuns){
    bgnIx <- 1+i*(k+2)
    endIx <- (i+1)*(k+2)
    TransMat[bgnIx:endIx,bgnIx:endIx] <- blockMat
  }
  for (i in 1:mxRuns){
    TransMat[i*(k+2),i*(k+2)+2] <- 1-b
  }
  TransMat[(mxRuns+1)*(k+2),sz] <- 1-b
  TransMat[sz,sz] <- 1
  return( TransMat )
}
##-------------------------------------------
## Mean number of runs of exactly length k  
## Name: r2MeanFct
## Input:
## pi1,pi2,aVec,bVec,mxRuns,k
## Output: 
## Mean number of runs of exactly length k
r2MeanFct <- function( pi1,pi2,aVec,bVec,mxRuns,k){
  sz <- (mxRuns+1)*(k+2)+1
  bgn <- c(0,pi1,pi2,rep(0,sz-3))
  prb <- bgn
  for (i in 1:length(aVec)){
    prb <- prb%*%
    r2TransMatFct( a=aVec[i],b=bVec[i],mxRuns=mxRuns,k=k)
  }
  ## Mean number of success runs of size exactly k
  wgtVec <- rep(0,sz)
  for (i in 0:mxRuns){
    bgnIx <- 1+i*(k+2)
    endIx <- (i+1)*(k+2)
    wgtVec[bgnIx:endIx] <- i
  }
  wgtVec[sz] <- mxRuns+1
  meank <- sum( prb*wgtVec )
  return( meank )
}
##-------------------------------------------------------------
## Finite Markov chain imbedding:
## Longest run in state 2 
##-------------------------------------------------------------
## Setting up the Markov chain
## Transition probability for FMCI
## Name:
## L2TransMatFct
## Input: a,b,mxRun
## Output: 
## Transition probability matrix for longest run 
BlockMatFct <- function( a,b,k ){
  BlockMat <- matrix(0,nrow=(k+2),ncol=(k+2))
  BlockMat[1,c(1,k+2)] <- c(1-b,b)
  BlockMat[2,1:2] <- c(a,1-a)
  for (i in 3:(k+2)){
    BlockMat[i,c(1,i)] <- c(1-b,b)
  }
  return(BlockMat)
}
L2TransMatFct <- function( a,b,mxRun ){
  ## Sizes of blocks
  BlockSz <- c(1,3:(mxRun+2),1)
  CumSumBlockSz <- cumsum( BlockSz )
  ## Sizes of transition matrix
  sz <- sum( BlockSz )
  TransMat <- matrix(0,nrow=sz,ncol=sz)
  TransMat[1,1:2] <- c(a,1-a)
  TransMat[sz,sz] <- 1
  ## Fill out the mxRun blocks of the transition matrix
  for (j in 1:mxRun){
    bgnIx <- CumSumBlockSz[j]+1
    endIx <- CumSumBlockSz[j+1]
    TransMat[bgnIx:endIx,(bgnIx+1):(endIx+1)] <- BlockMatFct( a,b,j )
  }
  return( TransMat )
}
##---------------------------------------------------------------
## Distribution for the longest run in state 2
## Name: L2PrbFct 
## Input: pi1,pi2,aVec,bVec,mxRun
## mxRun is the longest run considered
## Output: 
## Vector of probabilities for the longest run
L2PrbFct <- function( pi1,pi2,aVec,bVec,mxRun){
  BlockSz <- c(1,3:(mxRun+2),1)
  CumSumBlockSz <- cumsum( BlockSz )
  ## Size of initial probability vector 
  sz <- sum( BlockSz )
  bgn <- c( pi1,pi2,rep(0,(sz-2)) )
  prb <- bgn
  for (i in 1:length(aVec)){
    prb <- prb%*%L2TransMatFct( a=aVec[i],b=bVec[i],mxRun=mxRun )
  }
  CumSumPrb <- cumsum(prb)
  prbOut <- diff(c(0,CumSumPrb[CumSumBlockSz]))
  return(prbOut)
}
