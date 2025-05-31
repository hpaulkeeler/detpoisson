# Randomly simulates a determinantally-thinned Poisson point process. It
# then tests the empirical results against analytic ones.
#
# A determinantally-thinned Poisson point process is essentially a discrete
# determinantal point process whose underlying state space is a single 
# realization of a Poisson point process defined on some bounded continuous 
# space. 
#
# For more details, see the paper by Blaszczyszyn and Keeler[1].
#
# Author: H.P. Keeler, Inria/ENS, Paris, and University of Melbourne, 
# Melbourne, 2018.
#
#References:
#[1] Blaszczyszyn and Keeler, Determinantal thinning of point processes with 
#network learning applications, 2018.

#import relevant libraries
import numpy as np

from funSimSimpleLDPP import funSimSimpleLDPP
from funLtoK import funLtoK

##set random seed for reproducibility
#np.random.seed(1);

numbSim=10**4; #number of simulations

#START -- Parameters -- START
#Poisson point process parameters
lambda0=10; #intensity (ie mean density) of the Poisson process

#choose kernel
choiceKernel=1;#1 for Gaussian (ie squared exponetial ); 2 for Cauchy
sigma=1;# parameter for Gaussian and Cauchy kernel
alpha=1;# parameter for Cauchy kernel

#Simulation window parameters
xMin=0;xMax=1;yMin=0;yMax=1;
xDelta=xMax-xMin;yDelta=yMax-yMin; #rectangle dimensions
areaTotal=xDelta*yDelta; #area of rectangle
#END -- Parameters -- END

#Simulate a Poisson point process
numbPoints = np.random.poisson(lambda0*areaTotal);#Poisson number of points
xx = np.random.uniform(xMin,xMax,numbPoints);#x coordinates of Poisson points
yy = np.random.uniform(yMin,yMax,numbPoints);#y coordinates of Poisson points

# START -- CREATE L matrix -- START 
sizeL=numbPoints;
#Calculate Gaussian or kernel kernel based on grid x/y values
#all squared distances of x/y difference pairs
xxDiff=np.outer(xx, np.ones((numbPoints,)))-np.outer( np.ones((numbPoints,)),xx);
yyDiff=np.outer(yy, np.ones((numbPoints,)))-np.outer( np.ones((numbPoints,)),yy)
rrDiffSquared=(xxDiff**2+yyDiff**2);

if choiceKernel==1:
    #Gaussian/squared exponential kernel
    L=lambda0*np.exp(-(rrDiffSquared)/sigma**2);
    
elif choiceKernel==2:
        #Cauchy kernel
    L=lambda0/(1+rrDiffSquared/sigma**2)**(alpha+1/2); 
        
else:        
    raise Exception('choiceKernel has to be equal to 1 or 2.');
     
# END-- CREATE L matrix -- # END

#START Testing DPP simulation START#
#Eigen decomposition
eigenValL, eigenVectL=np.linalg.eig(L);

#run simulations with tests
probX_i_Emp=np.zeros(numbPoints); #initialize variables
indexTest=np.arange(2); #choose a subset of [0 numbPoints-1]
probTestEmp=0; #initialize variables
#loop through for each simulation
for ss in range(numbSim):
    #run determinantal simuation
    indexDPP=funSimSimpleLDPP(eigenVectL,eigenValL); #returns index 
    probX_i_Emp[indexDPP]=probX_i_Emp[indexDPP]+1;        
    
    countTemp=0; #initialize count
    for ii in range(len(indexTest)):
        #check that each point of test subset appears
        countTemp=countTemp+any(indexDPP==indexTest[ii]);
    
    probTestEmp=probTestEmp+(countTemp==len(indexTest));
    
    
#empirically estimate the probabilities of each point appearing
probX_i_Emp=probX_i_Emp/numbSim
print('probX_iEmp = ',probX_i_Emp);

#calculate exactly the probabilities of each point appearing
K=funLtoK(L);
probX_i_Exact=np.diag(K);
print('probX_i_Exact = ',probX_i_Exact);

#empirically estimate the probabilities of test subset appearing
probTestEmp=probTestEmp/numbSim
print('probTestEmp = ',probTestEmp)

#calculate exactly the probabilities of test subset appearing
probTestExact=np.linalg.det(K[:,indexTest][indexTest,:])
print('probTestExact = ',probTestExact);

#END Testing DPP simulation END#