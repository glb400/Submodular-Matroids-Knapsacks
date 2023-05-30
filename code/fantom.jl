using CPUTime

function knapsack_feasible(sol::Set{<:Integer},knapsack_constraints::Union{Array{<:AbstractFloat,2},Nothing})
    if knapsack_constraints===nothing 
        return true
    else
        return all(sum(knapsack_constraints[:,collect(sol)],dims=2) .<= 1)
    end
end

function fantom(gnd::Array{<:Integer},f_diff,ind_add_oracle;knapsack_constraints::Union{Array{<:AbstractFloat,2},Nothing}=nothing,epsilon::AbstractFloat=0.5,k::Integer)
    CPUtic();
    num_fun = 0
    num_oracle  = 0

    n = length(gnd)
    vals_elements = [f_diff(elm, Set{Int64}()) for elm in gnd]
    num_fun += length(gnd)    
    M = maximum(vals_elements)
    omgea = deepcopy(gnd)

	# number of matroid & knapsack constraints
    p = k
	ell = num_knapsack(knapsack_constraints)
	gamma = (2 * p * M) / ( (p+1) * (2 * p + 1))
	R = []
	val = 1
    while val <= n
        push!(R, val * gamma)                                                                                         
		val *= (1 + epsilon)
    end
    len_R = length(R)
    println("--------------len_R: $len_R---------------")
    U = []
	for rho in R
		omega = deepcopy(gnd)
		sols, num_f, num_oracle = iterated_GDT(p,ell,f_diff,ind_add_oracle,ground=omega,knapsack_constraints=knapsack_constraints,rho=rho)
        for sol in sols
            push!(U,sol)
            current_val=0.0
            set_a=Set{Int64}()
            for elm in sol
                current_val+=f_diff(elm,set_a)
                num_f+=1
                union!(set_a,elm)
            end
        end
        num_fun += num_f
    end
    best_sol = []
    best_f_val = 0
    for sol in U
        # calculate current value
        current_val=0.0
        set_a=Set{Int64}()
        for elm in sol
            current_val+=f_diff(elm,set_a)
            num_fun+=1
            union!(set_a,elm)
        end
        
        if current_val > best_f_val
            best_f_val = current_val
            best_sol = deepcopy(sol)
        end
    end
    CPUtoc();
    return best_sol, best_f_val, num_fun, num_oracle
end


function iterated_GDT(p::Integer,ell::Integer,f_diff,ind_add_oracle;ground::Array{<:Integer},knapsack_constraints::Union{Array{<:AbstractFloat,2},Nothing}=nothing,rho::AbstractFloat)
    num_fun = 0
    num_oracle = 0
    println("--------------rho: $rho---------------")
    S_i = []
	for i in collect(1:p+1)
		S, nf, noc = GDT(p,ell,f_diff,ind_add_oracle,ground=ground,knapsack_constraints=knapsack_constraints,rho=rho)
        num_fun += nf
        num_oracle += noc
        push!(S_i,deepcopy(S))
        ground = filter(x->x∉S, ground)
    end
	return S_i, num_fun, num_oracle
end


function GDT(p::Integer,ell::Integer,f_diff,ind_add_oracle;ground::Array{<:Integer},knapsack_constraints::Union{Array{<:AbstractFloat,2},Nothing}=nothing,rho::AbstractFloat)
    num_fun = 0
    num_oracle = 0
    S = Set{Int64}()
	flag = true	
    while flag
		flag = false
		cand = -1
		cand_val = -1
		for elm in ground
            # check the feasibility
            num_oracle += 1
            if ind_add_oracle(elm, S) && knapsack_feasible(Set(union(elm, S)), knapsack_constraints)
                val = f_diff(elm, S)
                num_fun += 1
				if  val / (sum(knapsack_constraints[:,collect(elm)]) + 1e-6) >= rho
					if val > cand_val
						cand = elm
						cand_val = val
						flag = true
                    end
                end
            end
        end
		if cand != -1
            union!(S, cand)
        end
    end
	return S, num_fun, num_oracle
end




































    