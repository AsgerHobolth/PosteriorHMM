##------------------------------------------------------
## Analysis of Fetal Lamb Movements
##-----------------------------------------------------
## Data is from Guttorp (1995) Exercise d9 page 124
##-----------------------------------------------------
ObsSeq <- c(0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,
            0,1,0,0,0,0,2,2,0,0,0,0,1,0,0,
            1,1,0,0,1,1,1,0,0,1,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,
            1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,1,0,0,0,0,0,7,3,2,3,2,4,
            0,0,0,0,1,0,0,0,0,0,0,0,1,0,2,
            0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,
            2,1,0,0,1,0,0,0,1,0,1,1,0,0,0,
            1,0,0,1,0,0,0,1,2,0,0,0,1,0,1,
            1,0,1,0,0,2,0,1,2,1,1,2,1,0,1,
            1,0,0,1,1,0,0,0,1,1,1,0,4,0,0,
            2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
## Plot data
plot(ObsSeq,pch=19,cex=.7,xlab="",
     ylab="Number of movements",main="Fetal Lamb Movements")
##-------------------------------------------
## We analyse the data using a HMM with two
## hidden states and Poisson emissions
## State 1 is a state with low activity and 
## state 2 is a state with high activity.
## The parameters are thus the initial probabilities for
## the Markov chain to begin in state 1 or state 2, the
## probabilities for staying in state 1 or 2 (a and b),
## and the low and high Poisson rates (lam1 and lam2)
##----------------------------------------------------
## MLE for the six parameters using the 
## Baum-Welch (EM) algorithm:
source("HMMPoisExpectations.R")
IntPrb <- c(0.5,0.5)
a <- 0.8 ; b <- 0.1 
TransProb <- matrix(c(a,1-a,1-b,b),byrow=TRUE,nrow=2)
lam1 <- 0.3 ; lam2 <- 3
Lambda <- c(lam1,lam2)
HMMRes <- HMMPoisExpectationsFct(InitProb=IntPrb,TransProb=TransProb,
                                 Lambda=Lambda,ObsSeq=ObsSeq)
## Run EM algorithm
nIter <- 40
for (i in 1:nIter){
  ## Updating initial probability
  IntPrb <- HMMRes$PostProb[1,]
  ## Updating rates
  lam1Est <- sum(ObsSeq*HMMRes$PostProb[,1])/sum(HMMRes$PostProb[,1])
  lam2Est <- sum(ObsSeq*HMMRes$PostProb[,2])/sum(HMMRes$PostProb[,2])
  Lambda <- c(lam1Est,lam2Est)
  ## Updating transitions
  aEst <- HMMRes$TransCnt[1,1]/sum(HMMRes$TransCnt[1,])
  bEst <- HMMRes$TransCnt[2,2]/sum(HMMRes$TransCnt[2,])
  TransProb <- matrix(c(aEst,1-aEst,1-bEst,bEst),byrow=TRUE,nrow=2)
  HMMRes <- HMMPoisExpectationsFct(InitProb=IntPrb,TransProb=TransProb,
                                   Lambda=Lambda,ObsSeq=ObsSeq)
  cat("Iteration",i,HMMRes$Lk,"\n")
}
cat("Maximum likelihood estimates:",
    "IntPrb:",IntPrb,
    "lam1Est:",lam1Est,"lam2Est:",lam2Est,"aEst:",aEst,"bEst:",bEst,"\n")
## Plot posterior probabilities
plot(1:length(ObsSeq),HMMRes$PostProb[,2],pch=19,col="blue",cex=.5)
##----------------------------------------------------------------------
## Final parameters and tables
##----------------------------------------------------------------------
## Initial Probability
IntPrb <- c(1,0)
## Transitions
a <- 0.989 ; b <- 0.703 
TransPrb <- matrix(c(a,1-a,1-b,b),byrow=TRUE,nrow=2)
## Emissions
lam1 <- 0.278 ; lam2 <- 3.217
Lam <- c(lam1,lam2)
HMMRes <- HMMPoisExpectationsFct(InitProb=IntPrb,TransProb=TransPrb,
                                 Lambda=Lam,ObsSeq=ObsSeq)
##--------------------------------------------
## Inhomogeneous Markov chain: Transitions 
##--------------------------------------------
source("IMCfunctions.R")
resIMC <- prmIMC(IntPrb,TransProb,Lambda,ObsSeq)
Len <- length(ObsSeq)
## Plot inhomogeneous transition probabilities
plot(2:Len,resIMC$aSeq,pch=19,cex=.5,col="tomato",ylim=c(0,1),type="s",xlab="",
     main="Inhomogeneous posterior probability of staying in state",
     lwd=3,cex.lab=1.2,ylab="Probability of staying")
points(2:Len,resIMC$bSeq,pch=19,cex=.5,col="darkgreen",type="s",lwd=2)
##-------------------------------------------------------------------
## Simulate nSim hidden state sequences from the posterior
nSim <- 1000 
rSimMat <- matrix(0,nrow=nSim,ncol=Len)
for (sim in 1:nSim){
  rSimMat[sim,] <- rIMC( resIMC$eta1,resIMC$eta2,resIMC$aSeq,resIMC$bSeq )
}
##---------------------------------
## Visualization of simulations
##---------------------------------
plot(0,0,xlim=c(0,Len),ylim=c(0,nSim+100),
     xlab="",ylab="Simulation number",col="white",cex.lab=1.2,
     main="Decoding and posterior hidden state simulations")
for (sim in 1:nSim){
  SimSeq <- which(rSimMat[sim,]>1.5)
  rect(SimSeq,sim,SimSeq+1,sim+1,col="orange",border=NA)
}
## Viterbi decoding
source("ViterbiPois.R")
HMMVit <- ViterbiPoisFct(InitProb=c(0.5,0.5),TransProb=TransProb,
                         Lambda=Lambda,ObsSeq=ObsSeq)
VitDec <- which(HMMVit$BackTrack>1.5)
rect(VitDec,1085,VitDec+1,1115,col="green",border=NA)
mtext("Viterbi",side=2,at=1110,cex=1.1,las=1)
## Posterior decoding
PostDec <- which(HMMRes$PostProb[,2]>0.5)
rect(PostDec,1035,PostDec+1,1065,col="blue",border=NA)
mtext("Posterior",side=2,at=1060,cex=1.1,las=1)
##----------------------------------------------------------------
## Posterior probability and Empirical frequencies
##----------------------------------------------------------------
plot(colMeans(rSimMat-1),pch=19,cex=.5,col="orange",ylim=c(0,1),type="s",
     main="Posterior probability of state",
     xlab="Time (intervals of 5 seconds)",
     lwd=4,cex.lab=1.2,ylab="Probability of state")
points(1:Len,HMMRes$PostProb[,2],pch=19,cex=.5,col="forestgreen",type="s",lwd=2)
legend("topleft",c("Exact","Empirical"),bty="n",cex=1.1,
       lty=1,col=c("forestgreen","orange"),lwd=c(3,4))
##----------------------------------------------------------
##-------------------- FMCI --------------------------------
##----------------------------------------------------------
## Number of jumps from state 1 to state 2
##----------------------------------------------------------
source("FMCIfunctions.R")
## Simulation-based number of jumps
n12EmpVec <- rep(0,nSim)
for (sim in 1:nSim){
  n12EmpVec[sim] <- 
    sum( substring(paste(rSimMat[sim,],collapse=""),1:(Len-1),2:Len)=="12" )
}
n12max <- max(n12EmpVec)+2
## Analytical number of jumps
n12PrbVec <- n12PrbFct( pi1=resIMC$eta1,pi2=resIMC$eta2,
                        aVec=resIMC$aSeq,bVec=resIMC$bSeq,
                        nMax=n12max )
## Plot empirical and analytical number of jumps
plot(0:(max(n12EmpVec)+1),tabulate(n12EmpVec+1,nbins=(n12max))/nSim,type="h",
     cex.lab=1.2,cex.main=1.2,lwd=6,col="orange",
     xlab="Number of jumps",ylab="Probability",
     main="Distribution of number of jumps",xlim=c(0,n12max))
points((0:n12max),n12PrbVec,col="forestgreen",pch=19,cex=1,type="b",lwd=2)
legend("topright",c("Analytical (FMCI)","Empirical (simulations)"),bty="n",
       cex=1.2,lty=1,col=c("forestgreen","orange"),lwd=c(3,4))
##----------------------------------------------------------------
## Number of positions in state 2
##----------------------------------------------------------------
## Simulation-based number of positions in state 2
n2EmpVec <- rowSums(rSimMat-1)
## Analytical number of state 2 positions
nMax <- max( n2EmpVec )+2
n2PrbVec <- n2PrbFct( pi1=resIMC$eta1,pi2=resIMC$eta2,
                      aVec=resIMC$aSeq,bVec=resIMC$bSeq,
                      nMax=nMax )
## Plot the two distributions
plot(0:(nMax-1),tabulate(n2EmpVec+1,nbins=nMax)/nSim,type="h",
     xlab="Number of state 2 positions in hidden Markov chain",
     ylab="Probability",lwd=6,xlim=c(0,nMax),
     main="Distribution of number of state 2 positions",cex.lab=1.2,col="orange")
points(0:nMax,n2PrbVec,col="forestgreen",pch=19,cex=1,type="b",lwd=2)
legend("topright",c("Analytical (FMCI)","Empirical (simulations)"),
       col=c("forestgreen","orange"),lwd=c(3,4),lty=1,bty="n",cex=1.2)
##----------------------------------------------------------------
## Run lengths in state 2
##----------------------------------------------------------------
## Empirical number of run lengths in state 2
x <- rle(rSimMat[1,])
rL <- x$lengths[which(x$values==2)]
for (sim in 2:nSim){
  x <- rle(rSimMat[sim,])
  rL <- c(rL,x$lengths[which(x$values==2)])
}
mxrL <- max(rL)+2
NEmpRuns <- tabulate(rL,nbins=mxrL)/nSim
## Expected analytical run length number 
NMeanRuns <- rep(0,mxrL)
for (k in 1:mxrL){
  mxRuns <- 1+2*ceiling(NEmpRuns[k])
  NMeanRuns[k] <- 
    r2MeanFct( pi1=resIMC$eta1,pi2=resIMC$eta2,
               aVec=resIMC$aSeq,bVec=resIMC$bSeq,
               mxRuns=mxRuns,k=k)
}
## Plot the two distributions
plot(1:mxrL,NEmpRuns,type="h",col="orange",
     xlim=c(1,mxrL+1),lwd=4,
     main="Expected number of hidden state 2 exact run length",
     xlab="Run length",ylab="Expected number of exact runs")
points(1:mxrL,NMeanRuns,col="forestgreen",pch=19,cex=1.2,type="b",lwd=1.5)
legend("topright",c("Analytical (FMCI)","Empirical (simulations)"),
       col=c("forestgreen","orange"),lwd=c(3,4),lty=1,bty="n",cex=0.9)
##----------------------------------------------------------------
## Longest run in state 2
##----------------------------------------------------------------
## Empirical number of longest run in state 2
EmpL2 <- rep(0,nSim)
for (sim in 1:nSim){
  x <- rle( rSimMat[sim,] )
  EmpL2[sim] <- max( x$lengths[which(x$values==2)] )
}
mxL2 <- max(EmpL2)+3
NEmpL2 <- tabulate(EmpL2,nbins=mxL2)/nSim
## Analytical longest run
L2Prb <- L2PrbFct( pi1=resIMC$eta1,pi2=resIMC$eta2,
                   aVec=resIMC$aSeq,bVec=resIMC$bSeq,
                   mxRun=mxL2)
## Plot the two distributions
plot(1:mxL2,NEmpL2,type="h",col="orange",
     xlim=c(1,mxL2+1),lwd=4,ylim=c(0,1),
     main="Distribution of longest run in state 2",
     xlab="Longest run",ylab="Probability")
points(1:(mxL2+1),L2Prb[2:(mxL2+2)],col="forestgreen",pch=19,cex=0.9,type="b",lwd=1.5)
legend("topright",c("Analytical (FMCI)","Empirical (simulations)"),
       col=c("forestgreen","orange"),lwd=c(3,4),lty=1,bty="n",cex=0.9)


