# detpoisson_julia

Julia implementation for simulating determinantally-thinned Poisson point processes.

## Overview

A determinantally-thinned Poisson point process is a discrete determinantal point process whose underlying state space is a single realization of a Poisson point process defined on some bounded continuous space. This is a repulsive point process, where the repulsion depends on the kernel and average density of points.

For more details, see the paper by Blaszczyszyn and Keeler:
https://arxiv.org/abs/1810.08672

## Files

| File | Description |
|------|-------------|
| `DemoDetPoisson.jl` | Basic demonstration of simulating a determinantally-thinned Poisson process |
| `TestDetPoisson.jl` | Tests DPP simulation against analytical probabilities |
| `funLtoK.jl` | Convert L-ensemble matrix to normalized K-kernel |
| `funSimSimpleLDPP.jl` | Simulate DPP using Kulesza-Taskar algorithm |

## Quick Start

```julia
include("DemoDetPoisson.jl")
```

Or for testing:

```julia
include("TestDetPoisson.jl")
```

## Requirements

- Julia 1.0+
- LinearAlgebra (standard library)
- Distributions
- Plots (for visualization)

## Other Julia DPP Libraries

For more advanced DPP implementations, see:
- [DeterminantalPointProcesses.jl](https://github.com/alshedivat/DeterminantalPointProcesses.jl)

## Reference

B. Blaszczyszyn and H.P. Keeler, "Determinantal thinning of point processes with network learning applications," 2018. https://arxiv.org/abs/1810.08672

See also: [detpoisson_matlab](../detpoisson_matlab), [detpoisson_python](../detpoisson_python), [detpoisson_r](../detpoisson_r)
