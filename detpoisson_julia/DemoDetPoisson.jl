#Simulate a determinantally-thinned Poisson point process on a rectangle
#Author: H. Paul Keeler, 2019.

#Note: Need the .+ for adding a scalar to an array
#Also need . for sqrt, exp, cos, sin etc and assigning scalars to arrays
#Index and Boolean arrrays need to be vectors eg v=zeros(n),  NOT v=zeros(n,1)

#clearconsole(); #for clearing Julia REPL console

using Distributions #for random simulations
using LinearAlgebra;
using PyPlot #uses Python for plotting
PyPlot.close("all");  # close all PyPlot figures

#set random seed for reproducibility
#Random.seed!(1234)

###START -- Parameters -- ###START
#Poisson point process parameters
lambda=50; #intensity (ie mean density) of the Poisson process

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

#START - Sampling/simulating DPP - START
# #Retrieve eigenvalues and eigenvectors
eigenVectL=eigvecs(L); #eigen decomposition -- vectors
eigenValL=(eigvals(L)); #eigen decomposition -- values
eigenVectK = (eigenVectL);
eigenValK = eigenValL./(1 .+eigenValL); #eigenvalues of K

#Bernoulli trials (ie coin flips) to determine number of points
booleEigen = (rand(length(eigenValK)) .<= eigenValK);

#number of points in the DPP realization
numbPointsDPP=sum(booleEigen);
#retrieve eigenvectors corresponding to successful Bernoulli trials
global spaceV = eigenVectK[:,booleEigen]; #subspace V
indexDPP = zeros(Int8,numbPointsDPP); #index for final DPP configuration

global numbPointsRemain=numbPointsDPP;
#Loop through for all points
for ii=1:numbPointsDPP
    #Compute probabilities for each point i
    global Prob_i=vec(sum(spaceV.^2,dims=2)); #sum across rows
    Prob_i=Prob_i./sum(Prob_i); #normalize

    #Choose a point (from 1 to numbPoints) using (prob mass function) Prob_i
    uRand=rand(1);
    indexCurrent= findall(x -> x> uRand[1], cumsum(Prob_i))[1];
    indexDPP[ii]=indexCurrent; #update index

    if ii<numbPointsDPP
        #Choose a vector to remove
        jj=findall(x-> x>0, abs.(spaceV[indexCurrent,:]))[1];
        columnVj=spaceV[:,jj]; #j-th column of V
        booleKeep_j=(jj).!=(1:numbPointsRemain);
        global spaceV=spaceV[:,booleKeep_j]; #remove column
        #probably a better way for removing columns

        #Update matrix V by removing Vj component from the space
        spaceV=spaceV-kron(columnVj, transpose(spaceV[indexCurrent,:]./columnVj[indexCurrent]));

        #Orthonormalize
        decompQR=qr(spaceV);
        spaceV=Matrix(decompQR.Q);
        #need Matrix function for properly sized matrix
    end

    global numbPointsRemain= numbPointsRemain-1; #update
end
indexDPP=sort(indexDPP);

#Plotting
#Plot Poisson point process
PyPlot.scatter(xx,yy, edgecolor="k", facecolor="none");
PyPlot.xlabel("x"); plt.ylabel("y");
#random color vector
vectorColor=rand(3);
#Plot determinantally-thinned Poisson point process
PyPlot.scatter(xx[indexDPP],yy[indexDPP],edgecolor="none",facecolor=vectorColor);

println("Program has ended.");
