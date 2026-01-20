# detpoisson_matlab

MATLAB implementation for simulating and fitting determinantally-thinned Poisson point processes.

## Overview

This is the most complete implementation, containing the full workflow for reproducing the numerical results from the paper by Blaszczyszyn and Keeler. It includes code for:
- Simulating determinantally-thinned Poisson point processes
- Generating Matern hard-core and Triangle point processes
- Fitting determinantal models using maximum likelihood
- Validating fits using nearest-neighbour and contact distributions

## Files

### Demonstration and Testing

| File | Description |
|------|-------------|
| `DemoDetPoisson.m` | Basic demonstration of simulating a determinantally-thinned Poisson process |
| `TestDetPoisson.m` | Tests DPP simulation against analytical probabilities |

### Full Workflow (Paper Reproduction)

| File | Description |
|------|-------------|
| `SubsetGenerate.m` | Step 1: Generate Matern I/II or Triangle hard-core processes |
| `SubsetDetPoissonFit.m` | Step 2: Fit determinantal model via maximum likelihood |
| `SubsetDetPoissonGenerate.m` | Step 3: Validate fit using distribution functions |

### Core Functions

| File | Description |
|------|-------------|
| `funLtoK.m` | Convert L-ensemble matrix to normalized K-kernel |
| `funSimSimpleLDPP.m` | Simulate DPP using eigendecomposition |
| `funNeighbourL.m` | Create L-matrix using nearest-neighbour features |
| `funLPalm.m` | Compute Palm distribution of L-matrix |

## Quick Start

```matlab
% Basic demonstration
DemoDetPoisson

% Test against analytical results
TestDetPoisson
```

## Reproducing Paper Results

```matlab
% Step 1: Generate training data (Matern II or Triangle process)
SubsetGenerate  % Creates Subset.mat

% Step 2: Fit determinantal model
SubsetDetPoissonFit  % Creates SubsetFitParam.mat

% Step 3: Validate the fit
SubsetDetPoissonGenerate  % Generates comparison plots
```

For exact paper reproduction, uncomment `seedRand=1; rng(seedRand)` in Steps 1 and 3.

## Requirements

- MATLAB (tested on R2018a+)
- Optimization Toolbox (for `fminsearch` in fitting)

## Reference

B. Blaszczyszyn and H.P. Keeler, "Determinantal thinning of point processes with network learning applications," 2018. https://arxiv.org/abs/1810.08672

See also: [detpoisson_python](../detpoisson_python), [detpoisson_julia](../detpoisson_julia), [detpoisson_r](../detpoisson_r)
