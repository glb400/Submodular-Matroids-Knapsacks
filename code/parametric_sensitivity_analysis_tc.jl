using Random
using DataFrames
using DataStructures
using StatsBase
using Revise
using LinearAlgebra
using Combinatorics
using Distances
using Graphs
using CPUTime
include("./SubmodularGreedy.jl")
using .SubmodularGreedy

n = 1000
erdos_graph = erdos_renyi(n, 0.01, seed=1)

# generate the weight matrix
weight_matrix = rand(Random.seed!(1), n, n)


# step1: construct objective function
function graph_cut(graph::Graphs.AbstractGraph,sol::Set{<:Integer})
    sol = collect(sol)
    sum = 0.0
    for elm in sol
        arr = all_neighbors(graph, elm)
        arr_new=filter(x->x∉sol,arr)
        for e in arr_new
            sum += weight_matrix[elm, e]
        end
    end
    return sum
end

f_diff(elm,sol)=graph_cut(erdos_graph,Set(union(elm,sol))) - graph_cut(erdos_graph,sol)

# step2: construct constraints
# construct matroid constraints
card_limit=10

function all_matroid_feasible(sol::Set{<:Integer},cardinality_limit::Integer)
    # cardinality constraint is feasible
    if length(sol)>cardinality_limit
        return false 
    end
    return true
end

# construct marginal independence oracle and independence oracle
ind_add_oracle(elm,sol)=all_matroid_feasible(union(sol,elm),card_limit)
ind_oracle(sol)=all_matroid_feasible(sol,card_limit)

budget_param=1.0

# construct knapsack constraints
knapsack_constraints=Array{Float64,2}(undef,2,n)
budget1=100*budget_param
degree_array = degree(erdos_graph)
knapsack_constraints[1,:] = degree_array / budget1

gnd=collect(1:n)
budget2=40*budget_param
knapsack_constraints[2,:] = (gnd .% 10 .+ 1) / budget2


# run algorithms
# 1 Greedy
sol,f_val,_=greedy(gnd,f_diff,ind_add_oracle,knapsack_constraints=knapsack_constraints)

# 2 Repeated Greedy
num_sol=2
epsilon=0.25
sol,f_val,_=repeated_greedy(gnd,f_diff,ind_add_oracle,knapsack_constraints=knapsack_constraints,num_sol=num_sol,epsilon=epsilon,k=1)

# 3 Simultaneous Greedy
sol,f_val,_=simultaneous_greedys(gnd,f_diff,ind_add_oracle,knapsack_constraints=knapsack_constraints,extendible=true,num_sol=num_sol,epsilon=epsilon,k=1)

# 4 FANTOM
sol,f_val,_=fantom(gnd,f_diff,ind_add_oracle,knapsack_constraints=knapsack_constraints,epsilon=epsilon,k=1)

# 5 SPROUT++
# param_c must be 1 if SPROUT++
param_c=1
ratio_tc_to_n=0.06
tc=n*ratio_tc_to_n
mu=1.0

function knapsack_feasible(sol::Set{<:Integer},knapsack_constraints::Union{Array{<:AbstractFloat,2},Nothing})
    if knapsack_constraints===nothing 
        return true
    else
        return all(sum(knapsack_constraints[:,collect(sol)],dims=2) .<= 1)
    end
end

function sproutpp(param_c::Integer,gnd::Array{<:Integer},f_diff,ind_oracle,ind_add_oracle;num_sol::Integer=0,k::Integer=0,knapsack_constraints::Union{Array{<:AbstractFloat,2},Nothing}=nothing,extendible::Bool=false,monotone::Bool=false,epsilon::AbstractFloat=0.0,opt_size_ub::Integer=length(gnd),verbose_lvl::Integer=1,param_alpha::AbstractFloat=0.5)
    CPUtic();
    num_fun=0
    best_sol=nothing
    best_f_val=0
    # generate all combinatorial set
    gnd_combinatorials=shuffle(collect(combinations(gnd,param_c)))
    # enumerate single-element set
    fA_max = 0.0
    for elm in gnd
        fA_max = max(fA_max,f_diff(elm, Set{Int64}()))
    end
    num_fun+=length(gnd)
    cnt=0
    for gnd_combinatorial in gnd_combinatorials
        # check independent system and knapsack constraints
        if ind_oracle(Set(gnd_combinatorial))==false||knapsack_feasible(Set(gnd_combinatorial),knapsack_constraints)==false
            continue
        end
        # compute f(A)
        fA_value=0.0
        set_a=Set{Int64}()
        for elm in gnd_combinatorial
            fA_value+=f_diff(elm,set_a)
            num_fun+=1
            union!(set_a,elm)
        end
        if fA_value < param_alpha * fA_max
            continue
        end
        # println("--------------------------------$cnt iteration--------------------------------")
        if cnt>=tc
            break
        end
        cnt+=1
        # construct new objective function z(S)=f(S \cup A) - f(A)
        # obviously, the marginal oracle z_diff(e,S)=z(e|S)=f(e \cup (S \cup A)) - f(S \cup A)=f(e|(S \cup A))=f_diff(e,S \cup A)
        z_diff(elm,sol)=f_diff(elm,union(sol,gnd_combinatorial))
        num_fun+=1
        
        # update the ground set N
        # CASE:SPROUT
        # gnd_new=filter(x->x∉gnd_combinatorial&&f_diff(x,Set(gnd_combinatorial))<=(fA_value/param_c),gnd)
        # num_fun+=n-length(gnd_combinatorial)
        
        # CASE:SPROUT++
        gnd_new=filter(x->x∉gnd_combinatorial,gnd)
        # constract all matroid constraints
        card_limit_new=card_limit-param_c

        ind_add_oracle_new(elm,sol)=all_matroid_feasible(union(sol,elm),card_limit_new)
        # decrease the capacity of all knapsack-cost functions
        knapsack_constraints_new=deepcopy(knapsack_constraints)

        # knapsack constraint 1
        budget1_new=budget1-sum(degree_array[gnd_combinatorial])
        knapsack_constraints_new[1,:]= degree_array / budget1_new

        # knapsack constraint 2
        budget2_new=budget2-sum(gnd[gnd_combinatorial] .% 10 .+ 1)
        knapsack_constraints_new[2,:]= (gnd .% 10 .+ 1) / budget2_new

        # run simultaneousGreedy
        sol,f_val,oracle_calls,_=main_part_sprout(param_c,fA_value,gnd_new,z_diff,ind_add_oracle_new,knapsack_constraints=knapsack_constraints_new,extendible=true,num_sol=num_sol,epsilon=epsilon,k=1,mu=mu)
        num_fun+=oracle_calls

        # add set A at the end
        f_val+=fA_value
        if sol===nothing
            sol=Set(gnd_combinatorial)
        else
            union!(sol,gnd_combinatorial)
        end

        if f_val>best_f_val
            best_sol,best_f_val=sol,f_val
        end

        println("-------------------------")
        println("Best current value is: $best_f_val")
        println("Best current solution set is: $best_sol")
        println("-------------------------")
    end
    CPUtoc();
    return best_sol,best_f_val,num_fun
end

sol,f_val,_=sproutpp(param_c,gnd,f_diff,ind_oracle,ind_add_oracle,knapsack_constraints=knapsack_constraints,extendible=true,num_sol=num_sol,epsilon=epsilon,k=1)
