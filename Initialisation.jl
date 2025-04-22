module Initialisation

export Patch, create_subpopulation, create_population, Population, SubPopulation

SubPopulation::Type = Vector{Tuple{Bool,Bool}} 

struct Population
    males::Matrix{SubPopulation}
    females::Matrix{SubPopulation}
end

#=
in: params
out: Population at t=0

•Population initialisation function
•Each individual has the same probability of getting each genotype
=# 
function create_population(number_of_patches::I, individuals_per_patch::I)::Population where I<:Integer
    males::Matrix{SubPopulation} = (x -> [(rand(Bool),rand(Bool)) for i=1:floor(Int,individuals_per_patch/5)]).(zeros(number_of_patches, 5))
    females::Matrix{SubPopulation} = (x -> [(rand(Bool),rand(Bool)) for i=1:floor(Int,individuals_per_patch/5)]).(zeros(number_of_patches, 5))
    return Population(males,females)
end





end