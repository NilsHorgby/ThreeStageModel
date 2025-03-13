module Initialisation

export Patch, create_subpopulation, create_population, Population, SubPopulation

SubPopulation::Type = Vector{Tuple{Bool,Bool}} 

struct Population
    males::Matrix{SubPopulation}
    females::Matrix{SubPopulation}
end

function create_population(number_of_patches::Int64, individuals_per_patch::Int64)::Population
    males::Matrix{SubPopulation} = (x -> [(rand(Bool),rand(Bool)) for i=1:floor(Int,individuals_per_patch/6)]).(zeros(number_of_patches, 5))
    females::Matrix{SubPopulation} = (x -> [(rand(Bool),rand(Bool)) for i=1:floor(Int,individuals_per_patch/6)]).(zeros(number_of_patches, 5))
    return Population(males,females)
end





end