# Fast projection onto the top-k-sum constraint

## Introduction
We consider the Euclidean projection onto the top-k-sum: 
```math
    \min_{\boldsymbol{x} \in \mathbb{R}^n}\ \frac{1}{2}\| \boldsymbol{x} - \boldsymbol{a} \|^2 \quad
    \text{subject to}\ x \in \mathcal{T}_{(k)}^r:=\left\{\boldsymbol{x} \in \mathbb{R}^n: \mathrm{T}_{(k)}(\boldsymbol{x}):=\sum_{i=1}^k \vec{x}_i \leq r\right\}，
```
where $r \in \mathbb{R}$, $k \in \{1,2,\dots,n\}$ and sequence $\boldsymbol{a} \in \mathbb{R}^n$ are given, $\vec{x}_1 \geq \vec{x}_2 \geq \dots \geq \vec{x}_n$ are the components of $\boldsymbol{x}$ in nonincreasing order. The related paper for this repository can be seen at <https://arxiv.org/abs/2512.10255>.
> 
## Repository Structure
The repository contains the following directories. 
- ```src/```: implementations of the algorithms, including 
    1. EIPS: an efficient intersection point searching algorithm callable by ```project_topksum_EIPS!```;
    2. ESGS, PLCP, GRID, Gurobi Solver. These methods are implemented in [1]. 
- ```experiment/```: scripts for running the experiments and processing the results. 
- ```plot/```: generated plots for the experiments.
- ```initial_point_selecting/```: results of the Experiment *Initial points selecting*. 
- ```complexity/``` and ```complexity_tau_r=1/```: results of the Experiment *Complexity*. 
- ```time_compare/```: results of the Experiment *Time comparing*. 

## Quick Start
A minimal example is provided in ```test.jl``` and ```demo.jl```.
- ```test.jl```: an example calling all methods and showing the difference among the methods. 
- ```demo.jl```: scripts for running EIPS for different $r$ and $k$. 


## Reference
[1] Roth,J.,Cui,Y.:Top-k-sum.jl. <https://doi.org/10.5281/zenodo.14189753>
