module SubmodularGreedy

using DataStructures
using StatsBase
using Revise

include("helper-funs.jl")
export intersection_ind_oracle

include("greedy.jl")
export greedy

include("repeated-greedy.jl")
export repeated_greedy, deterministic_usm

include("simultaneous-greedys.jl")
export simultaneous_greedys

include("sprout.jl")
export main_part_sprout

include("fantom.jl")
export fantom

end # module
