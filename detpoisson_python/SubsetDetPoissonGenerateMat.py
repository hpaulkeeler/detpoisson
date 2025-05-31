# This simulates determinatally-thinned point processes that have been 
# fitted to thinned-point process based on the method outlined in the paper 
# by Blaszczyszyn and Keeler[1], which is essentially the method developed 
# by Kulesza and Taskar[2]. 
#
# This is the third file (of three files) to run to reproduce results similar 
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
# REQUIREMENTS:
# Uses Statistics (and Machine learning) Toolbox for random variable.
#
# Author: H.P. Keeler, Inria/ENS, Paris, and University of Melbourne,
# Melbourne, 2019
#
# References:
# [1] Blaszczyszyn and Keeler, Determinantal thinning of point processes
# with network learning applications, 2018.
# [2] Kulesza and Taskar, "Determinantal point processes for machine
# learning",Now Publisers, 2012
# [3] Chiu, Stoyan, Kendall, and Mecke, "Stochastic geometry and its
# applications", Wiley.
# [4] Baddeley, Rubak and Turner, "Spatial point patterns: Methodology and
# applications with R, 2016.
# [5] Shirai and Takahashi, "Random point fields associated with certain
# Fredholm determinants I -- fermion, poisson and boson point", 2003.12


import numpy as np; #NumPy package for arrays, random number generation, etc
import matplotlib.pyplot as plt #for plotting
from matplotlib import collections  as mc #for plotting line segments
from scipy.io import loadmat #for reading mat files
from scipy.optimize import minimize #For optimizing
from scipy.stats import poisson #for the Poisson probability mass function

from funNeighbourL import funNeighbourL
from funSimSimpleLDPP import funSimSimpleLDPP


plt.close("all"); # close all figures

numbSim=10**3;

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
numbSimTotal=ppStructTemp.size; #total number of simulations 
#extract data for underlying Poisson point processes
xxList=[np.concatenate(ppStructTemp[ss][0][0]) for ss in range(numbSimTotal)];
yyList=[np.concatenate(ppStructTemp[ss][0][1]) for ss in range(numbSimTotal)];
ppXYPoisson=[(xxList[ss],yyList[ss])for ss in range(numbSimTotal)];

nList=[np.int(ppStructTemp[ss][0][2]) for ss in range(numbSimTotal)];
nArray=np.array(nList);
#extract data for subset point processes
indexSubTemp=dataMATLAB['indexCellSub'];  
indexListSub=[np.array(np.concatenate(indexSubTemp[ss][0])-1,dtype=int) for ss in range(numbSimTotal)];
#NOTE: need to subtract one from MATLAB indices as Python indexing starts at zero.

###END Load up values from MATLAB .mat file END###

###START Load up values from Python SubsetFitParam.npz file START###
fileVarsFitted=np.load('SubsetFitParam.npz'); #load file
thetaMax=fileVarsFitted['thetaMax'];
T=fileVarsFitted['T'];
sigma=fileVarsFitted['sigma'];
N=fileVarsFitted['N'];
M=fileVarsFitted['M'];
choiceModelFitted=fileVarsFitted['choiceModelFitted'];
booleOptSigma=fileVarsFitted['booleOptSigma'];
choiceKernel=fileVarsFitted['choiceKernel'];
###END Load up values from Python SubsetFitParam.npz file END###

if (numbSim+T>numbSimTotal):
    raise SystemExit('Need to create more realizations with SubsetGenerate.m');
else:
    #Look at unused realizations (ie the ones not used for fitting)
    #select a random subset of unused realizations
    ttValuesPerm=np.arange(T,numbSim+T);
    np.random.shuffle(ttValuesPerm);
    ttValues=ttValuesPerm[np.arange(numbSim)];


#initialize  variables for collecting statistics
numbSub=np.zeros(numbSim); #number of points in subset point process
numbDPP=np.zeros(numbSim); #number of points in detrminantal point process

for ss in range(numbSim):
    tt=ttValues[ss]; 
    xxPoisson=ppXYPoisson[tt][0];
    yyPoisson=ppXYPoisson[tt][1];
    indexSub=indexListSub[tt]; #index for sub point process    
    numbSub[ss]=indexSub.size; 
    
    #generate L matrix based on Poisson point realization
    L=funNeighbourL(xxPoisson,yyPoisson,lambda0,choiceKernel,sigma,thetaMax,N,M);
    #Eigen decomposition
    eigenValuesL, eigenVectorsL=np.linalg.eig(L);
        
    #Simulate next DPP generation
    indexDPP=funSimSimpleLDPP(eigenVectorsL,eigenValuesL);
    numbDPP[ss]=indexDPP.size;
    
lambdaEmpDPP=np.mean(numbDPP)/areaSample #empirical intensity of DPP
print('lambdaEmpDPP = ',lambdaEmpDPP);
lambdaEmpSub=np.mean(numbSub)/areaSample #empirical intensity of subset PP
print('lambdaEmpSub = ',lambdaEmpSub);

#Plotting 
#x/y values of subset point process
xxSub=xxPoisson[indexSub]; yySub=yyPoisson[indexSub];
#x/y values of determinantal point process
xxDPP=xxPoisson[indexDPP]; yyDPP=yyPoisson[indexDPP];

markerSize=12; #marker size for the Poisson points
#Plot Poisson point process
plt.plot(xxPoisson,yyPoisson,'ko',markerfacecolor="None",markersize=markerSize);
#Plot subset point process
plt.plot(xxSub,yySub,'rx',markersize=markerSize/2);
#Plot determinantally-thinned Poisson point process
plt.plot(xxDPP,yyDPP,'b+',markersize=markerSize);


plt.xlabel('x'); plt.ylabel('y');
plt.legend(('Poisson Process',labelModel,'Determinantal Poisson'));
plt.axis('equal');
#plt.axis('off');



