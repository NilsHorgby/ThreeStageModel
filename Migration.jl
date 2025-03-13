module Migration
using .Main.Initialisation
using Distributions: Normal


function put_within_boundaries(new_patch::Int64, number_of_patches::Int64)
    return (1 < new_patch <= number_of_patches) ? new_patch : (1 > new_patch) ? 0 : number_of_patches
end


function put_within_boundaries(new_patch::Int64, number_of_patches::Int64)
    if new_patch <= 0
        return 1
    elseif new_patch > number_of_patches
        return number_of_patches
    else 
        return new_patch
    end
end

function discrete_distrubution(dispersal_distance::Real, location::Int64)::Int64
    return round(Int,rand(Normal(location,dispersal_distance)))
end

function migrate!(population::Vector{SubPopulation},next_population::Vector{SubPopulation}, dispersal_distance::Float64)::Vector{SubPopulation}
    number_of_patches::Int64 = length(population)
    @inbounds for location in 1:number_of_patches
        @inbounds for individual_index in 1:length(population[location]) #see if a while loop is faster
            next_location::Int64 = put_within_boundaries(discrete_distrubution(dispersal_distance,location), number_of_patches)
            push!(next_population[next_location],population[location][individual_index])
        end
    end
    return next_population
end 

function migration!(population::Population, dispersal_distance::Float64)::Population
    #newborns don't migrate
    next_males = [[Vector{Tuple{Bool,Bool}}() for i=eachindex(population.males[:,1])] for j=1:6]
    next_males = hcat(migrate!.(eachcol(population.males),next_males,dispersal_distance)...)

    next_females = [[Vector{Tuple{Bool,Bool}}() for i=eachindex(population.females[:,1])] for j=1:6]
    next_females = hcat(migrate!.(eachcol(population.females),next_females,dispersal_distance)...)

    return Population(next_males,next_females)
end


end