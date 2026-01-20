# detpoisson_python

Python implementation for simulating and fitting determinantally-thinned Poisson point processes.

## Overview

A determinantally-thinned Poisson point process is a discrete determinantal point process whose underlying state space is a single realization of a Poisson point process defined on some bounded continuous space. This is a repulsive point process, where the repulsion depends on the kernel and average density of points.

For more details, see the paper by Blaszczyszyn and Keeler:
https://arxiv.org/abs/1810.08672

## Files

### Demonstration and Testing

| File | Description |
|------|-------------|
| `DemoDetPoisson.py` | Basic demonstration of simulating a determinantally-thinned Poisson process |
| `TestDetPoisson.py` | Tests DPP simulation against analytical probabilities |

### Fitting Workflow

| File | Description |
|------|-------------|
| `SubsetDetPoissonFitMat.py` | Fit determinantal model to MATLAB-generated training data |
| `SubsetDetPoissonGenerateMat.py` | Generate statistics using fitted parameters |

### Core Functions

| File | Description |
|------|-------------|
| `funLtoK.py` | Convert L-ensemble matrix to normalized K-kernel |
| `funSimSimpleLDPP.py` | Simulate DPP using Kulesza-Taskar algorithm |
| `funNeighbourL.py` | Create L-matrix using nearest-neighbour features |
| `funLPalm.py` | Compute Palm distribution of L-matrix |

## Quick Start

```bash
# Basic demonstration
python DemoDetPoisson.py

# Test against analytical results
python TestDetPoisson.py
```

## Requirements

- NumPy
- SciPy
- Matplotlib

## Other Python DPP Libraries

For more advanced DPP implementations, see:
- [DPPy](https://github.com/guilgautier/DPPy) - Comprehensive sampling for discrete and continuous DPPs
- [dpp](https://github.com/javiergonzalezh/dpp)
- [determinantal-point-processes](https://github.com/mbp28/determinantal-point-processes)

## Reference

B. Blaszczyszyn and H.P. Keeler, "Determinantal thinning of point processes with network learning applications," 2018. https://arxiv.org/abs/1810.08672

See also: [detpoisson_matlab](../detpoisson_matlab), [detpoisson_julia](../detpoisson_julia), [detpoisson_r](../detpoisson_r)
