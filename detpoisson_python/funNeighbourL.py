# function funNeighbourL(xx,yy,lambda0,choiceKernel,sigma,theta,N,M), returns L
# This function file creates an L(-ensemble-)matrix, as detailed in the
# paper by Blaszczyszyn and Keeler[1](Section IV).
#
# The quality features or covariates (q_x in the above paper) are based on
# the nearest neighbour distances. The similiarirty matrix (S in the paper)
# which creates replusion among the points, can be formed from either
# Gaussian or Cauchy kernel function.
#
# This Python code was originally written by H.P. Keeler in MATLAB; see 
# https://github.com/hpaulkeeler/DetPoisson_MATLAB
#
# INPUTS:
# xx and yy are the x and y values of the underlying discrete state space,
# which is usually a realization of a point process on a bounded continuous
# 2-D space such as a square or a disk.
#
# lambda0 is the point intensity/density (ie average number of points per
# unit area) of the point process, which is used to rescale the distances.
#
# choiceKernel is a variable that takes value 1 (for Gaussian) or 2 (for
# Cauchy) to select the kernel function.
#
# sigma is a parameter of kernel function.
#
# theta is a fitting parameters for the quality features/covariates.
#
# N is the number of neighbouring points.
#
# M is the number of distances between neighbour points. M is optional, but
# when used, M must be equal to zero or N-1.
#
# OUTPUTS:
# An L-ensemble kernel matrix for a determinantal point process on a
# discrete space; see [1] for details.
#
# Author: H.P. Keeler, Inria/ENS, Paris, and University of Melbourne,
# Melbourne, 2019.
#
# References:
# [1] Blaszczyszyn and Keeler, Determinantal thinning of point processes
# with network learning applications, 2018.


import numpy as np; #NumPy package for arrays, random number generation, etc
import matplotlib.pyplot as plt #for plotting
from scipy.io import loadmat #for reading mat files

from sklearn.neighbors import NearestNeighbors

plt.close("all"); # close all figures

### TEMP 
n=6;
#xx,yy=np.meshgrid(np.linspace(0,1,n),np.linspace(0,2,n));
#xx=np.linspace(0,1,n); 
xx=np.array([0,0.1,0.3,0.6,1,2]);
yy=np.linspace(0,3,n)
#yy=np.zeros(n);
lambda0=3;choiceKernel=1;sigma=1;
theta=np.array([1,2,3,4,5,6]);
N=3;
M=2;


def funNeighbourL(xx,yy,lambda0,choiceKernel,sigma,theta,N,M):
    #Check that M is the right value if it exists
    if ('M' in locals()):
        if (M!=N-1) and (M>0):
            raise SystemExit('M must be equal to N-1 or zero.');
        
            
    theta=theta.reshape(-1); #theta needs to be a vector
    xx=xx.reshape(-1); yy=yy.reshape(-1); #xx/yy need to be vectors
    meanD=1/np.sqrt(lambda0); #rescaling constant
    xx=xx/meanD;yy=yy/meanD; #rescale distances
    
    alpha=1; #an additional parameter for the Cauchy kernel
    
    #START - Creation of L matrix
    sizeL=xx.size; #width/height of L (ie cardinality of state space)
    
    #START -- Create q (ie quality feature/covariage) vector
    #zeroth term
    thetaFeature=theta[0]*np.ones(sizeL);
    if N>0:
        #1st to N th terms
        booleMatrix=(np.ones(sizeL)-np.eye(sizeL))==1; #removes distances to self
        indexMatrix=np.tile(range(sizeL),(sizeL,1));
        
        #find every pair combination
        xxMatrix=np.tile(xx,(sizeL-1,1)); #repeat vector xx sizeL-1 times
        yyMatrix=np.tile(yy,(sizeL-1,1)); #repeat vector yy sizeL-1 times
        
        #find rows/cols for accessing x/y values        
        rowTemp=np.tile(range(sizeL-1),(1,sizeL));    
        colTemp=indexMatrix[booleMatrix];     
        #access x/y values and reshape vectors
        xxTemp=xxMatrix[rowTemp,colTemp].reshape((sizeL-1,sizeL),order='F');    
        yyTemp=yyMatrix[rowTemp,colTemp].reshape((sizeL-1,sizeL),order='F');    
        
        #differences in distances for all pairs   
        xxNearDiff=xxMatrix-xxTemp;yyNearDiff=yyMatrix-yyTemp;       
        #calculate nearest distance for all pairs    
        distNearTemp=np.hypot(xxNearDiff,yyNearDiff);
        
        #sort distances (across cols)    
        indexNear=np.argsort(distNearTemp,axis=0); #sorting index    
        distNearTemp=distNearTemp[indexNear.reshape(-1),np.tile(range(sizeL),sizeL-1)];    
        distNear=distNearTemp.reshape((sizeL-1,sizeL));     
        
        #replicate parameter vector
        thetaMatrix=np.tile(theta[1:],(sizeL,1));
        #distance between nearest neighbours neighbours
        theta_distNear=thetaMatrix[:,range(N)].T*distNear[range(N),:];
        #add contribution from parameters and features
        thetaFeature=thetaFeature+np.sum(theta_distNear,axis=0);
        
        #Run for distances between nearest neighbours
        #N+1 to M th terms
        if ('M' in locals()) and (M>0):
            #find rows/cols of nearest neighbours
            colBetweenTemp=np.tile(range(sizeL),sizeL-1);
            rowBetween=np.column_stack((indexNear.reshape(-1),colBetweenTemp));    
            indexBetweenTemp=[np.ravel_multi_index(rowBetween[ii],(sizeL-1,sizeL)) for ii in range(sizeL*(sizeL-1))];           
            colBetween=colTemp[indexBetweenTemp];
            
            #access x/y values and reshape vectors                
            xxTemp=xxTemp[indexNear.reshape(-1),np.tile(range(sizeL),sizeL-1)];
            xxTemp=xxTemp.reshape((sizeL-1,sizeL));
            yyTemp=yyTemp[indexNear.reshape(-1),np.tile(range(sizeL),sizeL-1)];
            yyTemp=yyTemp.reshape((sizeL-1,sizeL));
            
            #differences in distances for all pairs   
            xxNearDiff=(xxTemp[0:M,:]-xxTemp[1:M+1,:]);        
            yyNearDiff=(yyTemp[0:M,:]-yyTemp[1:M+1,:]);        
            #distance between nearest neighbours neighbours
            distBetween=np.hypot(xxNearDiff,yyNearDiff);
            
            #dot product of parameters and features
            theta_distBetween=thetaMatrix[:,range(N,N+M)].T*distBetween;
            #add contribution from parameters and features
            thetaFeature=thetaFeature +np.sum(theta_distBetween,axis=0);
        
    
    qVector=np.exp(thetaFeature); #find q vector (ie feature/covariate values)
    
    #END -- Create q vector
    
    #START - Create similarity matrix S
    if sigma!=0:
        #all squared distances of x/y difference pairs    
        xxDiff=np.outer(xx, np.ones((sizeL,)))-np.outer( np.ones((sizeL,)),xx);
        yyDiff=np.outer(yy, np.ones((sizeL,)))-np.outer( np.ones((sizeL,)),yy);    
        rrDiffSquared=(xxDiff**2+yyDiff**2);
        if choiceKernel==1:
            ##Gaussian kernel
            SMatrix=np.exp(-(rrDiffSquared)/sigma**2);
        else:
            ##Cauchy kernel
            SMatrix=1/(1+rrDiffSquared/sigma**2)**(alpha+1/2);      
    else:
        SMatrix=np.eye(sizeL);
    
    
    #END - Create similarity matrix S with Gaussian kernel
    
    #START Create L matrix
    qMatrix=np.tile(qVector,(np.shape(SMatrix)[0],1)); #q diagonal matrix
    L=(qMatrix.T)*SMatrix*qMatrix;
    #END Create L matrix    
    return L;
