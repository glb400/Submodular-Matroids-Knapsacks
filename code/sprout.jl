using Random
using JLD2
using DataFrames
using LinearAlgebra

function density_search_for_sprout(param_c::Integer,fA_value::AbstractFloat, pq::PriorityQueue, num_sol::Integer, beta_scaling::AbstractFloat, gamma_scaling::AbstractFloat, delta::AbstractFloat, f_diff, ind_add_oracle, knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing}, epsilon::AbstractFloat, opt_size_ub::Integer; verbose=true,mu::AbstractFloat=1.0)
    # println("---------------density_search_for_sprout_start--------------")
    # get the max gain 
    _, max_gain = peek(pq)

    # initialize counts of function / oracle queries
    num_fun = 0
    num_oracle = 0

    # begin by specifying upper and lower bounds for the density ratio
    lower_density_ratio = beta_scaling * max_gain * 1
    upper_density_ratio = beta_scaling * max_gain * opt_size_ub

    # initialize best solutions
    best_sol = Set([])
    best_f_val = 0
    best_density_ratio =-1

    iter = 1
    while upper_density_ratio > (1+delta) * lower_density_ratio
        # compute the new density ratio (geometric average)
        density_ratio = sqrt( lower_density_ratio * upper_density_ratio )
        # run the algorithm 
        sol, fval, num_f, num_or, knap_reject = simultaneous_greedy_alg(deepcopy(pq), num_sol, epsilon, f_diff, ind_add_oracle, knapsack_constraints, density_ratio + gamma_scaling * fA_value / param_c, opt_size_ub)

        # update the best set 
        if fval > best_f_val 
            best_f_val = fval 
            best_sol = sol
            best_density_ratio = density_ratio
        end

        # update the number of function and oracle queries 
        num_fun += num_f 
        num_oracle += num_or
        if knap_reject
            lower_density_ratio = density_ratio - (mu - 1.0) * (density_ratio - lower_density_ratio) / mu            
        else
            upper_density_ratio = density_ratio + (mu - 1.0) * (upper_density_ratio - density_ratio) / mu  
        end

        iter += 1
    end # end binary search
    return best_sol, best_f_val, num_fun, num_oracle
end


function init_sgs_params_for_sprout(num_sol::Integer, k::Integer, extendible::Bool, monotone::Bool, knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing}, epsilon::AbstractFloat)
    # test whether all sets are feasible wrt knapsack constraints
    if ignorable_knapsack(knapsack_constraints)
        run_density_search = false
        m = 0
    else
        run_density_search = true
        m = num_knapsack(knapsack_constraints)
        @assert k > 0
    end

    # compute number of solutions and beta_scaling 
    if monotone

        if num_sol == 0
            num_sol = 1
        end
        
        if extendible
            p = max(num_sol-1, k)
        else
            p = k + num_sol - 1
        end

        beta_scaling = 0.0005
    else

        if extendible
            M = max( Integer(ceil(sqrt(1 + 2*m))), k )

            # only set num_sol if it was given as 0
            if num_sol == 0
                num_sol = M + 1
            end
            p = M 
        else

            # only set num_sol if it was given as 0
            if num_sol == 0
                num_sol = Integer(floor(2 + sqrt(k + 2*m + 2)))
            end
            p = k + num_sol - 1
        end

        beta_scaling = 0.0005

        gamma_scaling = 1e-6
    end

    return num_sol, run_density_search, beta_scaling, gamma_scaling
end


# main part of SPROUT
function main_part_sprout(param_c::Integer,fA_value::AbstractFloat, gnd::Array{<:Integer}, f_diff, ind_add_oracle; num_sol::Integer=0, k::Integer=0, knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing}=nothing, extendible::Bool=True, monotone::Bool=false, epsilon::AbstractFloat=0.0, opt_size_ub::Integer=length(gnd), verbose_lvl::Integer=1,mu::AbstractFloat=1.0)
    if length(gnd)==0
        # println("------------------EMPTY-------------------")
        best_sol = Set([])
        best_f_val = 0
        num_fun = 0
        num_oracle = 0
        return best_sol, best_f_val, num_fun, num_oracle    
    end
    
    # check that inputs make sense
    @assert (num_sol > 0) || (k > 0) "At least num_sol or k need to be specified."
    @assert dimension_check(gnd, knapsack_constraints) "There are more elemnents in knapsack constraints than in the ground set."
    @assert ignorable_knapsack(knapsack_constraints) || (epsilon > 0) "Because knapsack constraints are provided, please specify a non-zero value of epsilon for running the density search."
    @assert ignorable_knapsack(knapsack_constraints) || (k > 0) "Because knapsack constraints are provided, please specify the independence system parameter k for running the density search."
    @assert 0.0 <= epsilon <= 1.0 "Epsilon needs to be set in the range [0, 1]"

    # initialize a priority queue, initialize algorithm parameters
    num_sol, run_density_search, beta_scaling, gamma_scaling  = init_sgs_params_for_sprout(num_sol, k, extendible, monotone, knapsack_constraints, epsilon)
    pq, num_fun, num_oracle = initialize_pq(gnd, f_diff, ind_add_oracle, num_sol, knapsack_constraints)

    if length(pq)==0
        # println("------------------EMPTY-------------------")
        best_sol = Set([])
        best_f_val = 0
        return best_sol, best_f_val, num_fun, num_oracle
    end

    # verbosity levels
    info_verbose = verbose_lvl >= 1
    alg_verbose = verbose_lvl >= 2

    # run the algorithms
    if run_density_search
        delta = epsilon
        best_sol, best_f_val, num_f, num_or = density_search_for_sprout(param_c, fA_value, pq, num_sol, beta_scaling, gamma_scaling, delta, f_diff, ind_add_oracle, knapsack_constraints, epsilon, opt_size_ub; verbose=alg_verbose, mu=mu)

    else
        # run the plain simultaneous greedy with density ratio = 0
        density_ratio = 0.0
        best_sol, best_f_val, num_f, num_or, _ = simultaneous_greedy_alg(pq, num_sol, epsilon, f_diff, ind_add_oracle, knapsack_constraints, density_ratio, opt_size_ub; verbose=alg_verbose)
    end

    # update the number of oracle queries
    num_fun += num_f 
    num_oracle += num_or

    return best_sol, best_f_val, num_fun, num_oracle
end