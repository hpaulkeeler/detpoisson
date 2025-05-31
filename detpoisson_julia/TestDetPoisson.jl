#Simulate a determinantally-thinned Poisson point process on a rectangle
#Author: H. Paul Keeler, 2019.

#Note: Need the .+ for adding a scalar to an array
#Also need . for sqrt, exp, cos, sin etc and assigning scalars to arrays
#Index and Boolean arrrays need to be vectors eg v=zeros(n),  NOT v=zeros(n,1)

#clearconsole(); #for clearing Julia REPL console

using Distributions #for random simulations
using PyPlot #uses Python for plotting
using LinearAlgebra;

#include some helper functions
include("funSimSimpleLDPP.jl"); #for simulations
include("funLtoK.jl"); #for converting a non-singular L matrix to a singular matrix

#set random seed for reproducibility
#rng(1);

numbSim=10^4; #number of simulations

###START -- Parameters -- ###START
#Poisson point process parameters
lambda=10; #intensity (ie mean density) of the Poisson process

#choose kernel
choiceKernel=1; #1 for Gaussian (ie squared exponetial );2 for Cauchy
sigma=1;# parameter for Gaussian and Cauchy kernel
alpha=1;# parameter for Cauchy kernel

#Simulation window parameters
xMin=0;xMax=1;yMin=0;yMax=1;
xDelta=xMax-xMin;yDelta=yMax-yMin; #rectangle dimensions
areaTotal=xDelta*yDelta; #area of rectangle
###END -- Parameters -- ###END

### START -- Simulate a Poisson point process on a rectangle ###START
#Simulate Poisson point process
numbPoints=rand(Poisson(areaTotal*lambda)); #Poisson number of points
xx=xDelta*(rand(numbPoints)).+xMin;#x coordinates of Poisson points
yy=xDelta*(rand(numbPoints)).+yMin;#y coordinates of Poisson points
### END -- Simulate a Poisson point process on a rectangle --END ###

#numbPoints=5; xx=-2:1:2; yy=3:7 #TEMP

### START -- CREATE L matrix -- START ###
#all squared distances of x/y difference pairs
xxDiff=kron(xx,ones(1,numbPoints))-kron(ones(numbPoints,1),transpose(xx));
yyDiff=kron(yy,ones(1,numbPoints))-kron(ones(numbPoints,1),transpose(yy));
rrDiffSquared=(xxDiff.^2+yyDiff.^2);
if choiceKernel==1
    ##Gaussian/squared exponential kernel
    L=lambda*exp.(-(rrDiffSquared)/sigma^2);
elseif choiceKernel==2
    ##Cauchy kernel
    L=lambda./(1+rrDiffSquared/sigma^2).^(alpha+1/2);
else
    println("choiceKernel has to be equal to 1 or 2.");
end

### END-- CREATE L matrix -- ### END
L=Symmetric(L); #convert to Symmetric marix

### START Testing DPP simulation START ###
#Retrieve eigenvalues and eigenvectors
eigenVectL=eigvecs(L); #eigen decomposition -- vectors
eigenValL=eigvals(L); #eigen decomposition -- values

#run simulations with tests
global probX_i_Emp=zeros(numbPoints); #initialize variables
indexTest=2:3; #choose a subset of [1 numbPoints]
global probTestEmp=0; #initialize variables
#loop through for each simulation
for ss=1:numbSim
    #run determinantal simuation
    indexDPP=funSimSimpleLDPP(eigenVectL,eigenValL); #returns index
    global probX_i_Emp[indexDPP]=probX_i_Emp[indexDPP].+1;

    countTemp=0; #initialize count
    for ii=1:length(indexTest)
        #check that each point of test subset appears
        countTemp=countTemp+any(indexDPP.==indexTest[ii]);
    end
    global probTestEmp=probTestEmp.+(countTemp.==length(indexTest));
end

#empirically estimate the probabilities of each point appearing
probX_i_Emp=probX_i_Emp./numbSim
println("probX_i_Emp = ", probX_i_Emp);


#calculate exactly the probabilities of each point appearing
K=funLtoK(L);
probX_i_Exact=diag(K)
println("probX_i_Exact = ", probX_i_Exact);

#empirically estimate the probabilities of test subset appearing
probTestEmp=probTestEmp./numbSim
println("probTestEmp = ", probTestEmp);

#calculate exactly the probabilities of test subset appearing
probTestExact=det(K[indexTest,indexTest])
println("probTestExact = ", probTestExact);

###END Testing DPP simulation END###
