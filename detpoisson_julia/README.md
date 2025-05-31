# DetPoisson_Julia

Randomly simulates/samples a determinantally-thinned Poisson point process on a rectangle. I believe this is a new type of point process, originally proposed by Blaszczyszyn and Keeler in the paper[1]: 

https://arxiv.org/abs/1810.08672

Run the file: DemoDetPoisson.jl

A determinantally-thinned (Poisson) point process is essentially a discrete determinantal point process whose underlying state space is a single realization of a (Poisson) point process defined on some bounded continuous space. This is a repulsive point process, where the repulsion depends on the kernel and average density of points. For more details, see the paper by Blaszczyszyn and Keeler[1].

An obvious question is whether a determinantally-thinned Poisson point process is *also* a determinantal point process? The answer, we believe, is no, but it's not obvious. 

If you use this code in a publication, please cite the aforementioned paper by Blaszczyszyn and Keeler[1]. Unless stated otherwise, H.P. Keeler wrote this Python code, which is based on MATLAB also written by H.P. Keeler. For further details, see https://github.com/hpaulkeeler/DetPoisson_MATLAB

## Other repositories

I originally wrote all the code in R and in MATLAB, which both have a very similar structure; see:  

https://github.com/hpaulkeeler/DetPoisson_R 

https://github.com/hpaulkeeler/DetPoisson_MATLAB

I have also written it in Python, with the simulation part copied from other code. For details see:

https://github.com/hpaulkeeler/DetPoisson_Python

After (re)writing my simple code in Julia, I noticed this repository with Julia code that does more advanced sampling/simulating of determinantal point processes:

https://github.com/alshedivat/DeterminantalPointProcesses.jl

## References

[1] Blaszczyszyn and Keeler, Determinantal thinning of point processes with network learning applications, 2018.
