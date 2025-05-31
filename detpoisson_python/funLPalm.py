# This function calculates the L matrix for a Palm distribution
# conditioned on points existing at locations of the statespace indexed by 
# indexPalm. The method is discussed in the paper by  Blaszczyszyn and 
# Keeler (in Section B of the appendix), which is based on Proposition 1.2
# in the paper by Borodin and Rains[1], but an equivalent result appears 
# in the paper by Shirai and Takahashi[3]; see Theorem 6.5 and 
# Corolloary 6.6 in [3].
#
# This Python code was originally written by H.P. Keeler in MATLAB; see 
# https://github.com/hpaulkeeler/DetPoisson_MATLAB
#
# INPUTS:
# L = A square L(-matrix-)kernel, which must be (semi-)positive-definite. 
# indexPalm = an index set for the conditioned point, where all the points 
# of the underlying statespace correspond to the rows (or columns) of L.
#
# LPalm= The (reduced) Palm version of the L matrix, which is a square 
# matrix with dimension of size(L,1)-length(indexPalm).
#
# Author: H.P. Keeler, Inria/ENS, Paris, and University of Melbourne,
# Melbourne, 2019.
#
# References:
# [1] Blaszczyszyn and Keeler, "Determinantal thinning of point processes 
# with network learning applications", 2018.
# [2] Borodin and Rains, "Eynard-Mehta theorem, Schur process, and their 
# Pfaffian analogs", 2005
# [3] Shirai and Takahashi, "Random point fields associated with certain 
# Fredholm determinants I -- fermion, poisson and boson point", 2003.

import numpy as np; #NumPy package for arrays, random number generation, etc
import matplotlib.pyplot as plt #for plotting
from matplotlib import collections  as mc #for plotting line segments

#plt.close("all"); # close all figures
#L=np.random.uniform(0,1,(4,4));
#L=np.array([[3, 2, 1], [4, 5,6], [9, 8,7]]);
#indexPalm=np.array([0,1]); #indexing starts at zero

def funLPalm(L,indexPalm):
    sizeL=L.shape[1];
    #if max(indexPalm)>sizeL:
        #error('The index for the Palm points larger than matrix L.');
    indexAll=np.arange(sizeL);    #index
    indexRemain=np.setdiff1d(indexAll,indexPalm);#index of remaining points/locations
    identTemp=np.eye(sizeL); #identity matrix
    identTemp[indexPalm,indexPalm]=0;#some ones set to zero
    invLBoro=np.matmul(np.eye(sizeL),np.linalg.inv(identTemp+L));
    identTemp=np.eye(indexRemain.size);

    #reduce matrix in two steps
    #invLBoroTemp=invLBoro[indexRemain,:]; #access rows
    #invLBoroTemp=invLBoroTemp[:,indexRemain]; #access columns
    
    #reduce matrix in one step    
    invLBoroTemp=invLBoro[np.ix_(indexRemain,indexRemain)];    
    
    #Borodin and Rains construction -- see equation (24) in [1]
    LPalm=np.matmul(identTemp,np.linalg.inv(invLBoroTemp))-identTemp;
    return LPalm

