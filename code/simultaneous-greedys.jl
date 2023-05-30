using CPUTime

function simultaneous_greedy_alg(pq::PriorityQueue, num_sol::Integer, epsilon::AbstractFloat, f_diff, ind_add_oracle, knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing}, density_ratio::AbstractFloat, opt_size_ub::Integer; verbose=false)

    @assert num_sol > 0
    @assert 0 <= epsilon <= 1

    # get dimensions
    n = length(pq)

    # initialize the list of k solutions, their function values
    sol_list = [Set{Int64}() for i=1:num_sol]
    sol_vals_list = [0.0 for i=1:num_sol]

    # create solution costs (knapsack) and create densities of each element
    sol_costs = init_knapsack_costs(num_sol, knapsack_constraints)

    # initialize counts of function / oracle queries / knapsack reject status
    num_fun = 0
    num_oracle = 0
    knap_reject = false

    # get the maximum gain element
    best_elm_info, max_gain = peek(pq)

    # set the threshold
    threshold = (1 - epsilon) * max_gain
    min_threshold = epsilon * max_gain / opt_size_ub
    prev_threshold = threshold

    iter = 1
    print_iter = true

    # repeat until the threshold becomes too small, or we run out of elements
    while (threshold > min_threshold) && (length(pq) > 0)

        # take out the best element / solution pair 
        (elm, sol_ind, prev_size, density), f_gain = dequeue_pair!(pq)
        add_this_elm = false

        # if prev_size = |sol|, then this element / solution pair has maximal gain -- so, add it!
        if prev_size == length(sol_list[sol_ind])

            # check that knapsack constraint is satisfied, then update (if knapsack isn't satisfied, we just go back to the pq)
            if knapsack_feasible_to_add(elm, sol_ind, sol_costs, knapsack_constraints)

                # if verbose
                #     println("\t\tWe will be adding this element $elm to solution $sol_ind")
                # end

                # update solution & solution info
                union!(sol_list[sol_ind], elm)
                sol_vals_list[sol_ind] += f_gain
                update_sol_costs!(elm, sol_ind, sol_costs, knapsack_constraints)

                # remove added element from all element / solution pairs
                keys_to_remove = [k for k in keys(pq) if k[1] == elm]
                for k in keys_to_remove
                    delete!(pq, k)
                end

                iter += 1
                print_iter = true
            end

            # if the gain is less than the threshold, update the threshold appropriately (if no  more elements, then continue)
            if (f_gain < threshold) && (length(pq) > 0)
                prev_threshold = threshold
                threshold = min( (1 - epsilon) * threshold, peek(pq)[2])
            end

        # otherwise, we need to re-evaluate f(elm | sol) and decide whether to re-enqueue
        else

            # first test if S + e remains independent (otherwise throw away element)
            num_oracle += 1
            if ind_add_oracle(elm, sol_list[sol_ind])
                
                # evaluate the new marginal gain
                num_fun += 1
                f_gain = f_diff(elm, sol_list[sol_ind])

                # next, only keep considering the element if marginal gain is above the density ratio (otherwise throw away element)
                if (f_gain > density_ratio * density)

                    # re-enqueue this element
                    prev_size = length(sol_list[sol_ind]) # this is the current set size
                    enqueue!(pq, (elm, sol_ind, prev_size, density), f_gain)

                else
                    # record if the density condition is not satisfied
                    knap_reject = true
                end
            end
        end # end decision of whether to add or re-enqueue
    end # end while -- done creating solutions

    # add the single best element to the solution list
    push!(sol_list, Set([best_elm_info[1]]))
    push!(sol_vals_list, max_gain)
    
    # return the best solution
    best_sol_ind = argmax([sol_vals_list])
    best_sol = sol_list[best_sol_ind]
    best_f_val = sol_vals_list[best_sol_ind]

    return best_sol, best_f_val, num_fun, num_oracle, knap_reject
end 

function density_search(pq::PriorityQueue, num_sol::Integer, beta_scaling::AbstractFloat, delta::AbstractFloat, f_diff, ind_add_oracle, knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing}, epsilon::AbstractFloat, opt_size_ub::Integer; verbose=true)

    # get the max gain 
    _, max_gain = peek(pq)

    # initialize counts of function / oracle queries
    num_fun = 0
    num_oracle = 0

    # begin by specifying upper and lower bounds for the density ratio
    lower_density_ratio = beta_scaling * max_gain * 1
    upper_density_ratio = beta_scaling * max_gain * opt_size_ub

    # initialize best solutions
    best_sol = nothing
    best_f_val = -Inf

    iter = 1
    while upper_density_ratio > (1+delta) * lower_density_ratio

        # compute the new density ratio (geometric average)
        density_ratio = sqrt( lower_density_ratio * upper_density_ratio )

        # run the algorithm 
        sol, fval, num_f, num_or, knap_reject = simultaneous_greedy_alg(deepcopy(pq), num_sol, epsilon, f_diff, ind_add_oracle, knapsack_constraints, density_ratio, opt_size_ub)

        # update the best set 
        if fval > best_f_val 
            best_f_val = fval 
            best_sol = sol
        end

        # update the number of function and oracle queries 
        num_fun += num_f 
        num_oracle += num_or

        # update the binary search based on knap reject status
        if knap_reject
            upper_density_ratio = density_ratio
        else
            lower_density_ratio = density_ratio
        end

        iter += 1
    end # end binary search
    
    return best_sol, best_f_val, num_fun, num_oracle
end

function init_sgs_params(num_sol::Integer, k::Integer, extendible::Bool, monotone::Bool, knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing}, epsilon::AbstractFloat)

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
        
        beta_scaling = 2 * (1 - epsilon)^2 / (p + 1 + 2*m)

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

        beta_scaling = 2 * (1 - epsilon) * (1 - 1/num_sol - epsilon) / (p + 1 + 2*m)
    end

    return num_sol, run_density_search, beta_scaling
end

function simultaneous_greedys(gnd::Array{<:Integer}, f_diff, ind_add_oracle; num_sol::Integer=0, k::Integer=0, knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing}=nothing, extendible::Bool=false, monotone::Bool=false, epsilon::AbstractFloat=0.0, opt_size_ub::Integer=length(gnd), verbose_lvl::Integer=1)
    CPUtic();

    # check that inputs make sense
    @assert (num_sol > 0) || (k > 0) "At least num_sol or k need to be specified."
    @assert dimension_check(gnd, knapsack_constraints) "There are more elemnents in knapsack constraints than in the ground set."
    @assert ignorable_knapsack(knapsack_constraints) || (epsilon > 0) "Because knapsack constraints are provided, please specify a non-zero value of epsilon for running the density search."
    @assert ignorable_knapsack(knapsack_constraints) || (k > 0) "Because knapsack constraints are provided, please specify the independence system parameter k for running the density search."
    @assert 0.0 <= epsilon <= 1.0 "Epsilon needs to be set in the range [0, 1]"

    # initialize a priority queue, initialize algorithm parameters
    num_sol, run_density_search, beta_scaling  = init_sgs_params(num_sol, k, extendible, monotone, knapsack_constraints, epsilon)
    pq, num_fun, num_oracle = initialize_pq(gnd, f_diff, ind_add_oracle, num_sol, knapsack_constraints)

    if length(pq)==0
        best_sol = nothing
        best_f_val = -Inf
        return best_sol, best_f_val, num_fun, num_oracle
    end

    # verbosity levels
    info_verbose = verbose_lvl >= 1
    alg_verbose = verbose_lvl >= 2

    # run the algorithms
    if run_density_search
        delta = epsilon
        best_sol, best_f_val, num_f, num_or = density_search(pq, num_sol, beta_scaling, delta, f_diff, ind_add_oracle, knapsack_constraints, epsilon, opt_size_ub; verbose=alg_verbose)

    else
        # run the plain simultaneous greedy with density ratio = 0
        density_ratio = 0.0
        best_sol, best_f_val, num_f, num_or, _ = simultaneous_greedy_alg(pq, num_sol, epsilon, f_diff, ind_add_oracle, knapsack_constraints, density_ratio, opt_size_ub; verbose=alg_verbose)
    end

    # update the number of oracle queries
    num_fun += num_f 
    num_oracle += num_or

    CPUtoc();

    return best_sol, best_f_val, num_fun, num_oracle
end

