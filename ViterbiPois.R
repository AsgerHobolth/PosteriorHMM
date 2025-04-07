## Name: ViterbiPoisFct
## Author: Asger Hobolth
## Purpose:
## Calculates the path through the hidden states with the 
## highest probability using the Viterbi algorithm.
## This path is also called the Best Hidden Sequence
##--------------------------------------------------------
## Input:
## InitProb: Initial probabilities
##           (vector of length nHS)
##           [nHS=number of hidden states]
## TransProb: Transition probabilities:
##            Probabilities between hidden states
##            (nHS times nHS matrix)
## Lambda: Emission probabilities
##           Vector of rates of length nHS 
## ObsSeq: Observed sequence
##
## Output:
## list consisting of 
## BackTrack: Best Hidden Sequence
## MaxValue: Corresponds to EG's delta 
## MaxArrow: Corresponds to EG's psi
##---------------------------------------------------------
ViterbiPoisFct <- function(InitProb,TransProb,Lambda,ObsSeq){
  Len <- length(ObsSeq)
  nHS <- nrow(TransProb)
  ## Define MaxValue and MaxArrow
  MaxValue <- matrix(0,nrow=Len,ncol=nHS)
  MaxArrow <- matrix(0,nrow=Len,ncol=nHS)
  ## Initialization
  MaxValue[1,] <- InitProb*dpois(ObsSeq[1],Lambda)
  ## Determine MaxValue and MaxArrow by recursion
  for (k in 2:Len){
    for (j in 1:nHS){
      Val <- TransProb[,j]*
        rep(dpois(ObsSeq[k],Lambda[j]),nHS)*
          MaxValue[k-1,]
      MaxValue[k,j] <- max(Val)
      MaxArrow[k,j] <- which.max(Val)
    }
  }  
  ## Backtrack 
  BackTrack <- rep(0,Len)
  BackTrack[Len] <- which.max(MaxValue[Len,])
  for (i in Len:2) {
    BackTrack[i-1] <- MaxArrow[i,BackTrack[i]]
  }
  out <- list()
  out$BackTrack <- BackTrack
  out$MaxValue <- MaxValue
  out$MaxArrow <- MaxArrow
  return(out)
}



