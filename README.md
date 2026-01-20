# detpoisson

Simulating and fitting determinantally-thinned Poisson point processes.

## Overview

A determinantally-thinned Poisson point process is a discrete determinantal point process whose underlying state space is a single realization of a Poisson point process defined on some bounded continuous space. This is a repulsive point process, where the repulsion depends on the kernel and average density of points.

This new type of point process was originally proposed in the paper by Blaszczyszyn and Keeler:

> B. Blaszczyszyn and H.P. Keeler, "Determinantal thinning of point processes with network learning applications," 2018.
> https://arxiv.org/abs/1810.08672

An obvious question is whether a determinantally-thinned Poisson point process is *also* a determinantal point process? The answer, we believe, is no, but it's far from obvious.

## Repository Structure

```
detpoisson/
├── detpoisson_matlab/     # MATLAB implementation (full workflow)
├── detpoisson_python/     # Python implementation
├── detpoisson_julia/      # Julia implementation
├── detpoisson_r/          # R implementation
├── detpoisson-article.pdf # Paper with full details
└── LICENSE                # Apache 2.0
```

## Quick Start

### Demonstration

Run `DemoDetPoisson` in your preferred language to simulate a single determinantally-thinned Poisson point process:

**MATLAB:**
```matlab
cd detpoisson_matlab
DemoDetPoisson
```

**Python:**
```bash
cd detpoisson_python
python DemoDetPoisson.py
```

**Julia:**
```julia
include("detpoisson_julia/DemoDetPoisson.jl")
```

**R:**
```r
source("detpoisson_r/DemoDetPoisson.R")
```

## Reproducing Results from the Paper

The MATLAB files can reproduce the numerical results from the paper. The workflow fits determinantally-thinned Poisson point processes to dependently-thinned Poisson point processes such as Matern hard-core point processes.

### Three-Step Workflow

1. **Generate training data:** Run `SubsetGenerate.m` to create Matern I/II or Triangle hard-core point processes. Saves results to `Subset.mat`.

2. **Fit the model:** Run `SubsetDetPoissonFit.m` to fit a determinantal thinning model using maximum likelihood. Saves parameters to `SubsetFitParam.mat`.

3. **Validate the fit:** Run `SubsetDetPoissonGenerate.m` to compare the fitted model against the original process using nearest-neighbour and contact distributions.

To reproduce the exact results from the paper, set the random seed to one by uncommenting:
```matlab
seedRand=1; rng(seedRand)
```
in `SubsetGenerate.m` and `SubsetDetPoissonGenerate.m`.

## Key Concepts

### L-ensembles and K-kernels

The determinantal point process is defined via an L-ensemble matrix where:
- **L matrix**: Positive semi-definite matrix defining the point process
- **K kernel**: Related to L by K = L(I + L)^(-1), with eigenvalues in [0,1]

### Quality and Similarity

The L matrix decomposes as L[x,y] = q_x * S[x,y] * q_y where:
- **Quality (q_x)**: Measures the "goodness" of point x (higher = more likely to be retained)
- **Similarity (S[x,y])**: Measures repulsion between points (higher = less likely both appear)

### Sampling Algorithm

The code uses the Kulesza-Taskar algorithm:
1. Eigendecompose L to get eigenvalues and eigenvectors
2. Transform eigenvalues: λ_K = λ_L / (1 + λ_L)
3. Perform Bernoulli trials with eigenvalues as success probabilities
4. Iteratively select points and orthonormalize the remaining subspace

## Applications

- **Wireless network modeling**: Base station layouts exhibiting repulsion
- **Sleep/power schemes**: Coordinated transmitter on-off patterns
- **Transmission scheduling**: Medium access control with geometric dependencies

## References

1. B. Blaszczyszyn and H.P. Keeler, "Determinantal thinning of point processes with network learning applications," 2018. https://arxiv.org/abs/1810.08672

2. A. Kulesza and B. Taskar, "Determinantal point processes for machine learning," Foundations and Trends in Machine Learning, 2012.

## License

This project is licensed under the Apache License 2.0. See [LICENSE](LICENSE) for details.
