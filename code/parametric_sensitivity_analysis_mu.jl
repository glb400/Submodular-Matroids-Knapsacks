using Random
using JLD2
using DataFrames
using DataStructures
using StatsBase
using Revise
using LinearAlgebra
using Combinatorics
using Distances
Random.seed!(12)
include("./SubmodularGreedy.jl")
using .SubmodularGreedy

data_file="movie_info.jld2"
movie_info_df=load(data_file)["movie_info_df"]

n=size(movie_info_df,1)
names(movie_info_df)


# step1: construct objective function
function ij_to_ind(i::Integer,j::Integer,n::Integer)
    # ensure that i<j
    if i>j
        i,j=j,i
    end
    # map (i,j) to the correct entry in the upper triangular array
    return (i-1)*n - div(i*(i-1),2) + j
end

function get_similarity(similarity_array::Array{<:AbstractFloat},i::Integer,j::Integer,n::Integer)
    ind=ij_to_ind(i,j,n)
    return similarity_array[ind]
end

function compute_similarity_array(vec_list; sigma::AbstractFloat=1.0)
    # get dimensions
    n=length(vec_list)
    # initialize array of pair-wise similarities
    n_pairs=div(n*(n+1),2)
    similarity_array=zeros(n_pairs)
    # fill in the pair-wise similarities
    ind=1
    for i=1:n
        vi=vec_list[i]
        for j=i:n
            vj=vec_list[j]
            similarity_array[ind]= exp(-4 * euclidean(vi, vj))  
            ind+=1
        end
    end
    return similarity_array
end

similarity_array=compute_similarity_array(Array(movie_info_df[!,:vec]))

# construct marginal value oracle for efficient computation
function dispersion_diff(elm::Integer,sol::Set{<:Integer},similarity_array::Array{<:AbstractFloat},n::Integer)
    if elm in sol
        return 0.0
    else 
        # compute coverage and diversity terms
        coverage_term=sum([ get_similarity(similarity_array,i,elm,n) for i=1:n ])
        diversity_term=( 2 * sum([ get_similarity(similarity_array,i,elm,n) for i in sol]) + get_similarity(similarity_array,elm,elm,n))
        return (coverage_term - diversity_term) / n
    end
end

f_diff(elm,sol)=dispersion_diff(elm,sol,similarity_array,n)

# step2: construct constraints
# construct matroid constraints
movie_genre_df=Array(movie_info_df[!,:genre_set])

card_limit=10

genre_list=Set{String}()
for i=1:length(movie_genre_df)
    union!(genre_list,movie_genre_df[i])
end
genre_list=collect(genre_list)

genre_card_limit=2
limit_values=fill(genre_card_limit,19)

genre_limit=Dict{String,Int64}()
for i=1:length(limit_values)
    genre_limit[genre_list[i]]=limit_values[i]
end

function all_matroid_feasible(sol::Set{<:Integer},cardinality_limit::Integer,genre_limit::Dict{<:String,<:Integer},genre_list::Array{<:String},movie_genre_df::Array{<:Set{String}})
    # cardinality constraint is feasible
    if length(sol)>cardinality_limit
        return false 
    end    
    # each genre constraint is feasible
    genre_count=Dict{String,Int64}()
    for genre in genre_list
        genre_count[genre]=0
    end
    for elm in sol
        for genre in movie_genre_df[elm]
            if in(genre, genre_list)
                genre_count[genre]+=1
                if genre_count[genre]>genre_limit[genre]
                    return false
                end
            end
        end
    end
    return true
end

# construct marginal independence oracle and independence oracle
ind_add_oracle(elm,sol)=all_matroid_feasible(union(sol,elm),card_limit,genre_limit,genre_list,movie_genre_df)
ind_oracle(sol)=all_matroid_feasible(sol,card_limit,genre_limit,genre_list,movie_genre_df)

budget_param=1.0

# construct knapsack constraints
# knapsack_constraints=Array{Float64,2}(undef,3,n)
knapsack_constraints=Array{Float64,2}(undef,2,n)
# knapsack constraint 1
rating_array=Array(movie_info_df[!,:rating])

# max_rating=maximum(rating_array) # 9.7
max_rating=10
budget1=20*budget_param
knapsack_constraints[1,:]= (max_rating .-  reshape(rating_array,(1,n))) / budget1

# knapsack constraint 2
date_array=Array(movie_info_df[!,:year])

year1=1995
budget2=30*budget_param
knapsack_constraints[2,:]=abs.(year1 .- reshape(date_array,(1,n))) / budget2


# run algorithms
# 1 Greedy
gnd=collect(1:n)
sol,f_val,_=greedy(gnd,f_diff,ind_add_oracle,knapsack_constraints=knapsack_constraints)

# 2 Repeated Greedy
num_sol=2
epsilon=0.25
sol,f_val,_=repeated_greedy(gnd,f_diff,ind_add_oracle,knapsack_constraints=knapsack_constraints,num_sol=num_sol,epsilon=epsilon,k=20)

# 3 Simultaneous Greedy
sol,f_val,_=simultaneous_greedys(gnd,f_diff,ind_add_oracle,knapsack_constraints=knapsack_constraints,extendible=true,num_sol=num_sol,epsilon=epsilon,k=20)

# 4 FANTOM
sol,f_val,_=fantom(gnd,f_diff,ind_add_oracle,knapsack_constraints=knapsack_constraints,epsilon=epsilon,k=20)

# 5 SPROUT++
# set tc = 5
tc=floor(n*0.0005)

# param_c must be 1 if SPROUT++
param_c=1
param_alpha=0.5
mu=1.0

function knapsack_feasible(sol::Set{<:Integer},knapsack_constraints::Union{Array{<:AbstractFloat,2},Nothing})
    if knapsack_constraints===nothing 
        return true
    else
        return all(sum(knapsack_constraints[:,collect(sol)],dims=2) .<= 1)
    end
end

function sproutpp(param_c::Integer,gnd::Array{<:Integer},f_diff,ind_oracle,ind_add_oracle;num_sol::Integer=0,k::Integer=0,knapsack_constraints::Union{Array{<:AbstractFloat,2},Nothing}=nothing,extendible::Bool=false,monotone::Bool=false,epsilon::AbstractFloat=0.0,opt_size_ub::Integer=length(gnd),verbose_lvl::Integer=1,param_alpha::AbstractFloat=0.5)
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
        if cnt>=tc
            break
        end
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
        println("--------------------------------$cnt iteration--------------------------------")
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
        genre_limit_new=deepcopy(genre_limit)

        for elm in gnd_combinatorial
            for genre in movie_genre_df[elm]
                if in(genre, genre_list)
                    genre_limit_new[genre]=max(genre_limit_new[genre]-1,0)
                end
            end
        end
        ind_add_oracle_new(elm,sol)=all_matroid_feasible(union(sol,elm),card_limit_new,genre_limit_new,genre_list,movie_genre_df)

        # decrease the capacity of all knapsack-cost functions
        knapsack_constraints_new=deepcopy(knapsack_constraints)

        # knapsack constraint 1
        budget1_new=budget1-sum((max_rating .-  reshape(rating_array,(1,n)))[gnd_combinatorial])
        knapsack_constraints_new[1,:]= (max_rating .-  reshape(rating_array,(1,n))) / budget1_new

        # knapsack constraint 2
        budget2_new=budget2-sum(abs.(year1 .- reshape(date_array,(1,n)))[gnd_combinatorial])
        knapsack_constraints_new[2,:]=abs.(year1 .- reshape(date_array,(1,n))) / budget2_new

        # # knapsack constraint 3
        # budget3_new=budget3-sum(abs.(year2 .- reshape(date_array,(1,n)))[gnd_combinatorial])
        # knapsack_constraints_new[3,:]=abs.(year2 .- reshape(date_array,(1,n))) / budget3_new 

        # run simultaneousGreedy
        sol,f_val,oracle_calls,_=main_part_sprout(param_c,fA_value,gnd_new,z_diff,ind_add_oracle_new,knapsack_constraints=knapsack_constraints_new,extendible=true,num_sol=num_sol,epsilon=epsilon,k=20,mu=mu)
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
    return best_sol,best_f_val,num_fun
end

sol,f_val,_=sproutpp(param_c,gnd,f_diff,ind_oracle,ind_add_oracle,knapsack_constraints=knapsack_constraints,extendible=true,num_sol=num_sol,epsilon=epsilon,k=20)
