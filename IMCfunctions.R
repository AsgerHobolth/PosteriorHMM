##-------------------------------------------------------------
## General functions for 
## (i) Finding the inhomogeneous Markov chain: pi1,pi2,aVec,bVec  
## (ii) Simulating from an inhomogeneous Markov chain
## (iii) Summarizing the simulations in empirical distributions
##-------------------------------------------------------------
## (i) Finding the inhomogeneous Markov chain
##-------------------------------------------------------------
## Input:
## Parameters of HMM: InitProb, TransProb, Lambda
## Observed sequence: ObsSeq
## Output:
## Parameters of inhomogeneous Markov chain for the hidden 
## state sequence: 
## Initial distribution (pi1, pi2)
## Vectors of probabilities for staying (aVec, bVec)
prmIMC <- function(IntPrb,TransPrb,Lam,ObsSeq){
  HMMRes <- HMMPoisExpectationsFct(InitProb=IntPrb,TransProb=TransPrb,
                                   Lambda=Lam,ObsSeq=ObsSeq)
  ## Posterior initial distribution
  eta1 <- HMMRes$Backward[1,1]*dpois(ObsSeq[1],Lam[1])*IntPrb[1]/
    HMMRes$Lk
  eta2 <- HMMRes$Backward[1,2]*dpois(ObsSeq[1],Lam[2])*IntPrb[2]/
    HMMRes$Lk
  ## Posterior probabilities of staying
  Len <- length(ObsSeq)
  aSeq <- HMMRes$Backward[2:Len,1]*TransPrb[1,1]*dpois(ObsSeq[2:Len],Lam[1])/
    HMMRes$Backward[1:(Len-1),1]
  bSeq <- HMMRes$Backward[2:Len,2]*TransPrb[2,2]*dpois(ObsSeq[2:Len],Lam[2])/
    HMMRes$Backward[1:(Len-1),2]
  out <- list()
  out$eta1 <- eta1 ; out$eta2 <- eta2
  out$aSeq <- aSeq ; out$bSeq <- bSeq
  return( out )
}
##-------------------------------------------------------------
## (ii) Function for simulating from an inhomogeneous Markov chain
##-------------------------------------------------------------
## Input:
## pi1,pi2,aVec,bVec
## Output:
## Random sample from the Markov chain
rIMC <- function( pi1,pi2,aVec,bVec ){
  y <- sample(x=1:2,size=1,prob=c(pi1,pi2))
  x <- y
  for (i in 1:length(aVec)){
    xNew <- ifelse(x==1,
                   sample(1:2,1,prob=c(aVec[i],1-aVec[i])),
                   sample(2:1,1,prob=c(bVec[i],1-bVec[i])))
    x <- xNew
    y <- c(y,xNew)
  }
  return(y)
}