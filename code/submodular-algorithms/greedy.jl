using StatsBase

function greedy(pq::PriorityQueue, f_diff, ind_add_oracle; knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing}=nothing, density_ratio::AbstractFloat=0.0, epsilon::AbstractFloat=0.0, opt_size_ub::Integer=length(pq), verbose::Bool=false)
    best_sol, best_f_val, num_fun, num_oracle, knap_reject = simultaneous_greedy_alg(pq, 1, epsilon, f_diff, ind_add_oracle, knapsack_constraints, density_ratio, opt_size_ub, verbose=verbose)
    return best_sol, best_f_val, num_fun, num_oracle, knap_reject
end

function greedy(gnd::Array{<:Integer}, f_diff, ind_add_oracle; knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing}=nothing, density_ratio::AbstractFloat=0.0, epsilon::AbstractFloat=0.0, opt_size_ub::Integer=length(gnd), verbose::Bool=false)
    
    # initialize a priority queue, initialize algorithm parameters
    pq, num_fun, num_oracle = initialize_pq(gnd, f_diff, ind_add_oracle, 1, knapsack_constraints)

    # run greedy algorithm
    best_sol, best_f_val, num_f, num_or, knap_reject = greedy(pq, f_diff, ind_add_oracle; knapsack_constraints=knapsack_constraints, density_ratio=density_ratio, epsilon=epsilon, opt_size_ub=opt_size_ub, verbose=verbose)

    # update the number of oracle queries
    num_fun += num_f
    num_oracle += num_or

    return best_sol, best_f_val, num_fun, num_oracle, knap_reject
end

function sample_greedy(gnd::Array{<:Integer}, sample_prob::AbstractFloat, f_diff, ind_add_oracle; verbose::Bool=false)

    @assert 0 <= sample_prob <= 1

    # subsample the elements
    sample_gnd = Int64[]
    for elm in gnd 
        if rand() <= sample_prob
            push!(sample_gnd, elm)
        end
    end

    if verbose
        print("The subsampled elements are ")
        print(sample_gnd)
    end

    # run greedy with density ratio of 0.0
    best_sol, best_f_val, num_fun, num_oracle, _ = greedy(sample_gnd, f_diff, ind_add_oracle; knapsack_constraints=nothing, density_ratio=0.0, verbose=verbose)
    return best_sol, best_f_val, num_fun, num_oracle
end