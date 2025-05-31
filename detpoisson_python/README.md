# DetPoissonPython

Randomly simulates/samples a determinantally-thinned Poisson point process on a rectangle. I believe this is a new type of point process, originally proposed by Blaszczyszyn and Keeler in the paper[1]: 

https://arxiv.org/abs/1810.08672

A determinantally-thinned (Poisson) point process is essentially a discrete determinantal point process whose underlying state space is a single realization of a (Poisson) point process defined on some bounded continuous space. This is a repulsive point process, where the repulsion depends on the kernel and average density of points. For more details, see the paper by Blaszczyszyn and Keeler[1].

An obvious question is whether a determinantally-thinned Poisson point process is *also* a determinantal point process? The answer, we believe, is no, but it's not obvious. 

If you use this code in a publication, please cite the aforementioned paper by Blaszczyszyn and Keeler[1]. Unless stated otherwise, H.P. Keeler wrote this Python code, which is based on MATLAB also written by H.P. Keeler. For further details, see https://github.com/hpaulkeeler/DetPoisson_MATLAB

To simulate/sample the (discrete) determinantal point process, I modified the code in sample_dpp.py from this repository:

https://github.com/mbp28/determinantal-point-processes

## Other code repositories

I originally wrote all code in R and in MATLAB, which both have a very similar structure; see:  

https://github.com/hpaulkeeler/DetPoisson_R 

https://github.com/hpaulkeeler/DetPoisson_MATLAB

It should be noted that there are a number of repositories with Python code for simulating/sampling (discrete) determinantal point processes. A very comprehensive one with various simulation/sampling methods (for both discrete and continuous point processes) is the following: 

https://github.com/guilgautier/DPPy

Other repositories include:

https://github.com/javiergonzalezh/dpp

https://github.com/mbp28/determinantal-point-processes

https://github.com/ChengtaoLi/dpp

https://github.com/mehdidc/dpp

## Author
H.P. Keeler, Inria/ENS, Paris, and University of Melbourne, Melbourne, 2018.

## References
[1] Blaszczyszyn and Keeler, Determinantal thinning of point processes with network learning applications, 2018.
