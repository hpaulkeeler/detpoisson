# DetPoisson_R
Randomly simulates a determinantally-thinned Poisson point process on a rectangle. I believe this is a new type of point process, originally proposed by Blaszczyszyn and Keeler in the paper[1]: 

https://arxiv.org/abs/1810.08672

A determinantally-thinned (Poisson) point process is essentially a discrete determinantal point process whose underlying state space is a single realization of a (Poisson) point process defined on some (bounded) continuous space. This is a repulsive point process, where the repulsion depends on the kernel and average density of points. For more details, see the paper by Blaszczyszyn and Keeler[1].

An obvious question is whether a determinantally-thinned Poisson point process is *also* a determinantal point process? The answer, we believe, is no, but it's from obvious. 

References:
[1] Blaszczyszyn and Keeler, Determinantal thinning of point processes with network learning applications, 2018.
