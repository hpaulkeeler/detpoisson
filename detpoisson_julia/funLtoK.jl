# K=funLtoK(L)
# The function funLtoK(L) converts a kernel L matrix into a (normalized)
# kernel K matrix. The K matrix has to be semi-positive definite.
#
# Code available here:
# Keeler, 2018, https://github.com/hpaulkeeler/DetPoisson_Julia
#
# #TEST
# B=[3, 2, 1; 4, 5,6; 9, 8,7];
# L=B'*B
# L =
#
#   106    98    90
#    98    93    88
#    90    88    86
# K=funLtoK(L)
# K =
#
#     0.7602    0.3348   -0.0906
#     0.3348    0.3320    0.3293
#    -0.0906    0.3293    0.7492

function funLtoK(L)
    eigenVectLK=eigvecs(L); #eigen decomposition -- vectors
    eigenValL=(eigvals(L)); #eigen decomposition -- values
    eigenValK = eigenValL./(1 .+eigenValL); #eigenvalues of K
    eigenValK=Diagonal(eigenValK); #eigenvalues of L as diagonal matrix
    K=eigenVectLK*eigenValK*(eigenVectLK'); #recombine from eigen components
    K=real(K); #make sure all values are real
    return K;
end
