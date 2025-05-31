# This file fits a determinatally-thinned point process to a 
# (dependently-)thinned-point process based on the method outlined in the 
# paper by Blaszczyszyn and Keeler[1], which is essentially the method 
# developed by Kulesza and Taskar[2] in Section 4.1.1.
#
# This is the second file (of three files) to run to reproduce results similar 
# to those presented in the paper by Blaszczyszyn and Keeler[1].
#
# The data used for fitting (or training) is stored in the file Subset.mat, 
# which is generated with the MATLAB file SubsetGenerate.m; see 
#
# https://github.com/hpaulkeeler/DetPoisson_MATLAB
# 
# The fitting paramters are stored locally in the file SubsetFitParam.npz 
#
# This code was originally written by H.P. Keeler in MATLAB; see 
# https://github.com/hpaulkeeler/DetPoisson_MATLAB
# 
# Author: H.P. Keeler, Inria/ENS, Paris, and University of Melbourne,
# Melbourne, 2019
#
# References:
# [1] Blaszczyszyn and Keeler, Determinantal thinning of point processes
# with network learning applications, 2018.
# [2] Kulesza and Taskar, "Determinantal point processes for machine 
# learning",Now Publisers, 2012

import numpy as np; #NumPy package for arrays, random number generation, etc
import matplotlib.pyplot as plt #for plotting
from matplotlib import collections  as mc #for plotting line segments
from scipy.io import loadmat #for reading mat files
from scipy.optimize import minimize #For optimizing
from scipy.stats import poisson #for the Poisson probability mass function

from funNeighbourL import funNeighbourL

plt.close("all"); # close all figures

T=100;#Number of training/learning samples

###START Load up values from MATLAB .mat file START###
dataMATLAB=loadmat('Subset.mat');
lambda0=np.double(dataMATLAB['lambda']);#intensity of underlying Poisson PP
xx0=np.double(dataMATLAB['xx0']); 
yy0=np.double(dataMATLAB['yy0']); 
areaSample=np.double(dataMATLAB['areaSample']); #area of sample window
rSub=np.double(dataMATLAB['rSub']); #radius of matern or triangular process
lambdaSub=np.double(dataMATLAB['lambdaSub']); #intensity of subset PP
windowSample=dataMATLAB['windowSample'][0]; #vector describing window dims
choiceModel=np.int(dataMATLAB['choiceModel']); #model number (ie 1,2 or 3)
labelModel=str(dataMATLAB['labelModel'][0]); #name/label of model
booleDisk=np.int(dataMATLAB['booleDisk'])!=0; #if simulation window is disk
#x/y values of all underlying Poisson PPs
ppStructTemp=dataMATLAB['ppStructPoisson']; 
numbTrain=ppStructTemp.size; #total number of simulations
#extract data for underlying Poisson point processes
xxList=[np.concatenate(ppStructTemp[ss][0][0]) for ss in range(numbTrain)];
yyList=[np.concatenate(ppStructTemp[ss][0][1]) for ss in range(numbTrain)];
ppXYPoisson=[(xxList[ss],yyList[ss])for ss in range(numbTrain)];

nList=[np.int(ppStructTemp[ss][0][2]) for ss in range(numbTrain)];
nArray=np.array(nList);
#extract data for subset point processes
indexSubTemp=dataMATLAB['indexCellSub'];  
indexListSub=[np.array(np.concatenate(indexSubTemp[ss][0])-1,dtype=int) for ss in range(numbTrain)];
#NOTE: need to subtract one from MATLAB indices as Python indexing starts at zero.

###END Load up values from MATLAB .mat file END###

###START -- Model fitting parameters -- START
#The parameters can be changed, but these were observed to work well.
if any(choiceModel==np.array([1,2])):
    #Fitting parameters for Matern hard-core (I/II) point process
    N=1; #number of neighbours for distances -- use N=1 or N=2
    M=0;#Must be N=M-1 or M=0.
    booleOptSigma=True; #Set to true to also optimize sigma value
    choiceKernel=1; #1 for Gaussian kernel, 2 for Cauchy kernel
    sigma=1;#sigma for Gaussian or Cauchy kernel
    #Sigma value is ignored if booleOptSigma=true
else:
    # Fitting parameters for triangle point process
    N=2; #number of neighbours for distances -- use N=1 or N=2
    M=1;#Must be N=M-1 or M=0.
    booleOptSigma=False; #Set to true to also optimize sigma value
    choiceKernel=1; #1 for Gaussian kernel, 2 for Cauchy kernel
    sigma=0;#sigma for Gaussian or Cauchy kernel 
    #Sigma value is ignored if booleOptSigma=true

###END -- Model fitting parameters -- END

##Probability of a Poisson realization with two few points
probLessM=np.sum(poisson.pmf(np.arange(N+1),lambda0));
if any(nArray<=N):
    raise SystemExit('Underlying Poisson realization needs at least N points');

#total number of possible training/learning samples
numbTrain=nArray.size; 
if T>numbTrain:
    raise SystemExit('Not enough training samples ie T>numbSim.');

#Deterministic (ie gradient) optimization method
if booleOptSigma:
    thetaGuess=np.ones(N+M+2); #Take initial guess for sigma values
else:
    thetaGuess=np.ones(N+M+1); #Take initial guess for theta values

#Function definitions for log-likelihood.
def funLikelihood_data(T,ppXYPoisson,indexListSub,choiceKernel,lambda0,sigma,theta,booleOptSigma,N,M):
    if booleOptSigma: 
        #sets sigma to one of the theta parameters to optimize
        sigma=theta[-1:]; 
        theta=theta[:-1];     
    #initialize  vector    
    logLikelihoodVector=np.zeros(T);

    #Loop through all training/learning samples
    for tt in range(T):    
        xx=ppXYPoisson[tt][0];yy=ppXYPoisson[tt][1];
        indexSub=indexListSub[tt]; #index for sub point process
        
        #Create L matrix (ie for Phi) based on nearest neighbours
        L=funNeighbourL(xx,yy,lambda0,choiceKernel,sigma,theta,N,M);
               
        #Create sub L matrix (ie for Psi)
        subL=L[np.ix_(indexSub,indexSub)]; 
        
        logLikelihoodVector[tt]=(np.log(np.linalg.det(subL))-np.log(np.linalg.det(L+np.eye(L.shape[0]))));    
    #END for-loop                 
    logLikelihood=np.sum(logLikelihoodVector);
    return logLikelihood

#function to maximize. See above for funLikelihood_Data function
def funMax_theta(theta):
    return funLikelihood_data(T,ppXYPoisson,indexListSub,choiceKernel,lambda0,sigma,theta,booleOptSigma,N,M);

#define function to be minimized
def funMin(theta):
    return (-1*funMax_theta(theta)); 

#Minimize function -- may take a while. 
resultsOpt=minimize(funMin,thetaGuess, method='BFGS',options={'disp': True});
thetaMax=resultsOpt.x; 

if booleOptSigma:
    sigma=thetaMax[-1:]; #retrieve sigma values from theta vector
    thetaMax=thetaMax[:-1]; 

print('sigma = ', sigma);    
print('thetaMax', thetaMax);
choiceModelFitted=choiceModel; #record which model was used for fitting

#save model fitting parameters in a .npz file
np.savez('SubsetFitParam',thetaMax=thetaMax,T=T,sigma=sigma,N=N,M=M,\
         choiceModelFitted=choiceModelFitted,booleOptSigma=booleOptSigma,\
         choiceKernel=choiceKernel);
