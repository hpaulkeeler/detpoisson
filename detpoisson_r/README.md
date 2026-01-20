# detpoisson_r

R implementation for simulating determinantally-thinned Poisson point processes.

## Overview

A determinantally-thinned Poisson point process is a discrete determinantal point process whose underlying state space is a single realization of a Poisson point process defined on some bounded continuous space. This is a repulsive point process, where the repulsion depends on the kernel and average density of points.

For more details, see the paper by Blaszczyszyn and Keeler:
https://arxiv.org/abs/1810.08672

An obvious question is whether a determinantally-thinned Poisson point process is *also* a determinantal point process? The answer, we believe, is no, but it's far from obvious.

## Files

| File | Description |
|------|-------------|
| `DemoDetPoisson.R` | Demonstration of simulating a determinantally-thinned Poisson process |

## Quick Start

```r
source("DemoDetPoisson.R")
```

## Features

The demonstration script:
1. Simulates a Poisson point process on a unit square
2. Constructs an L-matrix using Gaussian or Cauchy kernel
3. Performs eigendecomposition and Bernoulli sampling
4. Applies the Kulesza-Taskar algorithm to select points
5. Visualizes the original Poisson process and the determinantal subset

## Requirements

- R 3.0+
- Base R (no additional packages required)

## Reference

B. Blaszczyszyn and H.P. Keeler, "Determinantal thinning of point processes with network learning applications," 2018. https://arxiv.org/abs/1810.08672

See also: [detpoisson_matlab](../detpoisson_matlab), [detpoisson_python](../detpoisson_python), [detpoisson_julia](../detpoisson_julia)
