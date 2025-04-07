# Name: HMMPoisExpectationsFct.R
# Author: Asger Hobolth
# Purpose:
# Calculates expected transition count matrix
# Calculates likelihood
# Calculates posterior probability
#--------------------------------------------------------
# Input:
# InitProb: Initial probabilities
# TransProb: Transition probabilities;
#            Probabilities between hidden states
#            nHS times nHS matrix
# Lambda: Emission probabilities for Poisson 
#         Vector of length nHS
# ObsSeq: Observed sequence
#
# Output:
# $TransCnt: Expected transition count matrix
# $PostProb: Posterior probability matrix
# $Lk: Log likelihood
# $BackwardLik: Backward matrix
#---------------------------------------------------------
HMMPoisExpectationsFct <- function(InitProb,TransProb,Lambda,ObsSeq){
  len <- length(ObsSeq)
  nHS <- nrow(TransProb)
  #-------------------
  # Forward algorithm 
  #-------------------
  # Define ForwardLik matrix
  ForwardLik <- matrix(0,nrow=len,ncol=nHS)
  # Start condition
  ForwardLik[1,] <- InitProb*dpois(ObsSeq[1],Lambda)
  # Determine ForwardLik by recursion
  for (k in 2:len){
    for (j in 1:nHS){
      ForwardLik[k,j] <- dpois(ObsSeq[k],Lambda[j])*
                         sum(TransProb[,j]*
                             ForwardLik[k-1,])
    }
  }
  ForwardLikVal <- sum(ForwardLik[len,])
  #cat("Likelihood from Forward algorithm:",ForwardLikVal,"\n")
  #--------------------
  # Backward algorithm
  #--------------------
  # Define BackwardLik
  BackwardLik <- matrix(0,nrow=len,ncol=nHS)
  # Start condition
  BackwardLik[len,] <- rep(1,nHS)
  # Determine BackwardLik by recursion
  for (k in (len-1):1){
    for (j in 1:nHS){
      BackwardLik[k,j] <- sum(TransProb[j,1:nHS]*
                              dpois(ObsSeq[k+1],Lambda)*  
                              BackwardLik[k+1,])
    }
  }
  BackwardLikVal <- sum(InitProb*
                        dpois(ObsSeq[1],Lambda)*
                        BackwardLik[1,])
  #cat("Likelihood from Backward algorithm:",BackwardLikVal,"\n")
  ##-----------------------
  ## Posterior probability
  ##-----------------------
  PostProb <- exp(log(BackwardLik)+log(ForwardLik)-log(BackwardLikVal))
  ##----------------------------
  ## Expected transition counts
  ##----------------------------
  TransCnt <- matrix(0,nrow=nHS,ncol=nHS)
  for (k in 1:nHS){
    for (l in 1:nHS){
      Probkl <- ForwardLik[1:(len-1),k]*BackwardLik[2:len,l]/
        BackwardLikVal*TransProb[k,l]*dpois(ObsSeq[2:len],Lambda[l])
      TransCnt[k,l] <- sum(Probkl)
    }
  }
  output <- list()
  output$TransCnt <- TransCnt
  output$PostProb <- PostProb
  output$Lk <- BackwardLikVal
  output$Backward <- BackwardLik
  return(output)
}

