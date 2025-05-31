# Randomly simulates a determinantally-thinned Poisson point process. 
#
# A determinantally-thinned Poisson point process is essentially a discrete
# determinantal point process whose underlying state space is a single 
# realization of a Poisson point process defined on some bounded continuous space. 
#
# For more details, see the paper by Blaszczyszyn and Keeler[1].
#
# Author: H.P. Keeler, Inria/ENS, Paris, and University of Melbourne, Melbourne, 2018.
#
#References:
#[1] Blaszczyszyn and Keeler, Determinantal thinning of point processes with network learning applications, 2018.

rm(list=ls(all=TRUE)); #clear all variables
graphics.off(); #close all figures
cat("\014"); #clear screen

#START -- Parameters -- START
#Poisson point process parameters
lambda=50; #intensity (ie mean density) of the Poisson process

#choose kernel
choiceKernel=2;#1 for Gaussian (ie squared exponetial ); 2 for Cauchy
sigma=1;# parameter for Gaussian and Cauchy kernel
alpha=1;# parameter for Cauchy kernel

#Simulation window parameters
xMin=0;xMax=1;yMin=0;yMax=1;
xDelta=xMax-xMin;yDelta=yMax-yMin; #rectangle dimensions
areaTotal=xDelta*yDelta; #area of rectangle
#END -- Parameters -- END

#Simulate a Poisson point process on a rectangle
numbPoints=rpois(1,areaTotal*lambda);#Poisson number of points
xx=xDelta*runif(numbPoints)+xMin;#x coordinates of Poisson points
yy=xDelta*runif(numbPoints)+yMin;#y coordinates of Poisson points

# START -- CREATE L matrix -- START 
sizeL=numbPoints;
#Calculate Gaussian or kernel kernel based on grid x/y values
#all squared distances of x/y difference pairs
xxDiff=(outer(xx,rep(1,sizeL))-outer(rep(1,sizeL),xx));
yyDiff=(outer(yy,rep(1,sizeL))-outer(rep(1,sizeL),yy));
rrDiffSquared=(xxDiff^2+yyDiff^2);
if (choiceKernel==1){
  #Gaussian/squared exponential kernel
  L=lambda*exp(-(rrDiffSquared)/sigma^2);
} else{ 
  if (choiceKernel==2){
    #Cauchy kernel
    L=lambda/(1+rrDiffSquared/sigma^2)^(alpha+1/2);
  } else{
    stop('choiceKernel has to be equal to 1 or 2.');
  }
}
# END-- CREATE L matrix -- # END

# START Simulating/sampling DPP
#Eigen decomposition
tmp=eigen(L); eigenValuesL=tmp$values; eigenVectorsL=tmp$vectors
eigenVectorsL[,2]=tmp$vectors[,3]; eigenVectorsL[,3]=tmp$vectors[,2]# REMOVE later

eigenValuesK <- eigenValuesL / (1+eigenValuesL); #eigenvalues of K
indexEigen <- which(runif(sizeL) <= eigenValuesK); #index of successful Bernoulli trials

#number of points in the DPP realization
numbPointsDPP<-length(indexEigen);
#retrieve eigenvectors corresponding to successful Bernoulli trials
spaceV <- eigenVectorsL[,indexEigen]; #subspace V
indexDPP <- rep(0,numbPointsDPP); #vector for final DPP configuration
numbPointsRemain=numbPointsDPP; #number of remaining points

if (numbPointsDPP>0){
  #Loop through for all points
  for (ii in 1:numbPointsDPP){
    #Compute probabilities for each point i
    Prob_i <- rowSums(spaceV^2); #sum across rows
    Prob_i <- Prob_i / sum(Prob_i); #normalize
    
    #Choose a new point using PMF Prob_i
    indexCurrent <- min(which(cumsum(Prob_i)>runif(1)));
    indexDPP[ii]<-indexCurrent;
    numbPointsRemain= numbPointsRemain-1; #update
    
    if (ii <numbPointsDPP){ 
      #Choose a vector to eliminate
      jj<-min(which(spaceV[indexCurrent,]!=0));    
      columnVj <- spaceV[,jj];
      spaceV <- spaceV[,-jj];
      
      spaceV=matrix(spaceV,sizeL,numbPointsRemain); #reshape matrix
      
      #Update matrix V
      spaceV <- spaceV - outer(columnVj,(spaceV[indexCurrent,]/columnVj[indexCurrent])); #remove Vj component from the space
      #subspaceV<-orth(subspaceV); # Orthogonalize V using SVD
      tempQR=qr(spaceV,k=ii);
      spaceV<-qr.Q(tempQR); #Orthonormalize using Householder method
    }
      
  }
}
indexDPP <- sort(indexDPP); #sort points
#END - Simulating/sampling DPP - END

#Plotting
#Plot Poisson point process
plot(xx,yy,xlab='x', ylab='y', col="black",pch=1,cex=3);
colorAll=colors(); #list all colors
colorRand=colorAll[sample(length(colorAll),1)] #randomly choose one color
#Plot determinantally-thinned Poisson point process
points(xx[indexDPP],yy[indexDPP],col=colorRand,pch=16,cex=2);
