## General information

This repository contains the implementation of paper [Submodular Maximization Under the Intersection of Matroid and Knapsack Constraints](https://arxiv.org/abs/2307.09487). It includes experiments on the movie recommendation and weighted max-cut applications and corresponding parametric sensitivity analysis. Part of codes use the open-source toolkit from Chris Harshaw's work in https://github.com/crharshaw/SubmodularGreedy.jl.


## Algorithms

```
greedy.jl: The code of Greedy is function: greedy(pq::PriorityQueue, f_diff, ind_add_oracle; knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing}=nothing, density_ratio::AbstractFloat=0.0, epsilon::AbstractFloat=0.0, opt_size_ub::Integer=length(pq), verbose::Bool=false) in this file.

repeated-greedy.jl: The code of Repeated Greedy is function: repeated_greedy(gnd::Array{<:Integer}, f_diff, ind_add_oracle; num_sol::Integer=0, k::Integer=0, knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing}=nothing, monotone::Bool=false, epsilon::Float64=0.0, opt_size_ub::Integer=maximum(gnd), verbose_lvl::Integer=1) in this file.

simultaneous-greedys.jl: The code of DENSITYSEARCH is function: simultaneous_greedys(gnd::Array{<:Integer}, f_diff, ind_add_oracle; num_sol::Integer=0, k::Integer=0, knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing}=nothing, extendible::Bool=false, monotone::Bool=false, epsilon::AbstractFloat=0.0, opt_size_ub::Integer=length(gnd), verbose_lvl::Integer=1) in this file. Note that the param 'knapsack_constraints' must not be 'nothing' or this algorithm will deteriorate to the normal Simultaneous Greedy Algorithm.

fantom.jl: The code of FANTOM is function: fantom(gnd::Array{<:Integer},f_diff,ind_add_oracle;knapsack_constraints::Union{Array{<:AbstractFloat,2},Nothing}=nothing,epsilon::AbstractFloat=0.5,k::Integer) in this file.

sprout.jl: The code of the binary search part of SPROUT++ is function: main_part_sprout(param_c::Integer,fA_value::AbstractFloat, gnd::Array{<:Integer}, f_diff, ind_add_oracle; num_sol::Integer=0, k::Integer=0, knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing}=nothing, extendible::Bool=True, monotone::Bool=false, epsilon::AbstractFloat=0.0, opt_size_ub::Integer=length(gnd), verbose_lvl::Integer=1,mu::AbstractFloat=1.0) in this file. For different problems, we construct the SPROUT++ algorithm by call this function and perform the partial enumeration and other techniques in the corresponding setting files.
```

## Settings

```
SubmodularGreedy.jl: This is a Julia file to export the algorithms.

helper-funs.jl: The fundamental algorithmic tools for test algorithms, i.e., Greedy, Repeated Greedy, FANTOM, DENSITYSEARCHSGS, SPROUT++.
```

## Citation

If you find our work helpful, please cite:

```bibtex
@article{gu2023submodular,
  author       = {Yu{-}Ran Gu and
                  Chao Bian and
                  Chao Qian},
  title        = {Submodular Maximization under the Intersection of Matroid and Knapsack
                  Constraints},
  journal      = {CoRR},
  volume       = {abs/2307.09487},
  year         = {2023},
}
```
