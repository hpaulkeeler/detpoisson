
function funSimSimpleLDPP(eigenVectL,eigenValL)

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
return indexDPP=sort(indexDPP);
end
