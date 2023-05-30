------------------------------------------------------------------------------------------
	                   Readme for the SPROUT++ algorithm Package	          
------------------------------------------------------------------------------------------

A. General information:

This package includes the Julia code of the SPROUT++ algorithm for Submodular Maximization Under the Intersection of Matroid and Knapsack Constraints. Part of these codes use the open source toolkit from chris harshaw's work in https://github.com/crharshaw/SubmodularGreedy.jl.

*************************************************************


1. First part includes the following algorithmic files:

helper-funs.jl: The fundamental algorithmic tools for test algorithms, i.e., Greedy, Repeated Greedy, FANTOM, DENSITYSEARCHSGS, SPROUT++.

greedy.jl: The code of Greedy is function: greedy(pq::PriorityQueue, f_diff, ind_add_oracle; knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing}=nothing, density_ratio::AbstractFloat=0.0, epsilon::AbstractFloat=0.0, opt_size_ub::Integer=length(pq), verbose::Bool=false) in this file.

repeated-greedy.jl: The code of Repeated Greedy is function: repeated_greedy(gnd::Array{<:Integer}, f_diff, ind_add_oracle; num_sol::Integer=0, k::Integer=0, knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing}=nothing, monotone::Bool=false, epsilon::Float64=0.0, opt_size_ub::Integer=maximum(gnd), verbose_lvl::Integer=1) in this file.

simultaneous-greedys.jl: The code of DENSITYSEARCH is function: simultaneous_greedys(gnd::Array{<:Integer}, f_diff, ind_add_oracle; num_sol::Integer=0, k::Integer=0, knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing}=nothing, extendible::Bool=false, monotone::Bool=false, epsilon::AbstractFloat=0.0, opt_size_ub::Integer=length(gnd), verbose_lvl::Integer=1) in this file. Note that the param 'knapsack_constraints' must not be 'nothing' or this algorithm will deteriorate to the normal Simultaneous Greedy Algorithm.

fantom.jl: The code of FANTOM is function: fantom(gnd::Array{<:Integer},f_diff,ind_add_oracle;knapsack_constraints::Union{Array{<:AbstractFloat,2},Nothing}=nothing,epsilon::AbstractFloat=0.5,k::Integer) in this file.

sprout.jl: The code of the binary search part of SPROUT++ is function: main_part_sprout(param_c::Integer,fA_value::AbstractFloat, gnd::Array{<:Integer}, f_diff, ind_add_oracle; num_sol::Integer=0, k::Integer=0, knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing}=nothing, extendible::Bool=True, monotone::Bool=false, epsilon::AbstractFloat=0.0, opt_size_ub::Integer=length(gnd), verbose_lvl::Integer=1,mu::AbstractFloat=1.0) in this file. For different problems, we construct the SPROUT++ algorithm by call this function and perform the partial enumeration and other techniques in the corresponding setting files.

SubmodularGreedy.jl: This is a Julia file to export the algorithms.


2. Second part includes the following setting files:

movie_recommendation_2knapsacks.jl: The setting of the experiment of "varying knapsack budget" and "varying knapsack budget" on the movie recommendation problem. In this file we consider two knapsacks c_1 and c_2. We vary the knapsack budget through the param 'budget_param' from 0.64 to 1.0 when the param 'card_limit' is fixed to 10. Likewise, we vary the maximum number of allowed movies through the param 'card_limit' from 4 to 10 when the param 'budget_param' is fixed to 1.0.

movie_recommendation_3knapsacks.jl: The setting of the experiment of "varying knapsack budget" and "varying knapsack budget" on the movie recommendation problem. In this file we consider three knapsacks c_1, c_2 and c_3.  We vary the knapsack budget through the param 'budget_param' from 0.64 to 1.0 when the param 'card_limit' is fixed to 10. Likewise, we vary the maximum number of allowed movies through the param 'card_limit' from 4 to 10 when the param 'budget_param' is fixed to 1.0.

weighted_max_cut.jl: The setting of the experiment of "varying knapsack budget" and "varying knapsack budget" on the max-cut problem. We vary the knapsack budget through the param 'budget_param' from 0.64 to 1.0 when the param 'card_limit' is fixed to 10. Likewise, we vary the maximum number of allowed nodes through the param 'card_limit' from 4 to 10 when the param 'budget_param' is fixed to 1.0.

parametric_sensitivity_analysis_tc.jl: The setting of the experiment of "sensitivity analysis of the param t_c" on the max-cut problem. We vary the ratio of t_c to n through the param 'ratio_tc_to_n' from 0.02 to 1.0.

parametric_sensitivity_analysis_mu.jl: The setting of the experiment of "sensitivity analysis of the param mu" on the movie recommendation problem. We vary the value of $\mu$ through the param 'mu' from 1.0 to 5.0.

Besides, we can set all the params mentioned in our papers, which are $c$, $\alpha$, $\mu$, $\beta$, $\gamma$, $\epsilon$, $\delta$ corresponding to the params 'param_c', 'param_alpha', 'mu', 'beta_scaling', 'gamma_scaling', 'epsilon', 'delta' in our code.

*************************************************************

All of the experiments are coded in Julia 1.7.2 and run on an identical configuration: a server with 2 Intel(R) Xeon(R) Platinum 8171M CPU @ 2.60Ghz CPU and 128GB of RAM. The kernel is Linux version 5.4.0-122-generic.

------------------------------------------------------------------------------------------
