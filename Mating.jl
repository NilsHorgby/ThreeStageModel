module Mating

using .Main.Initialisation
using Distributions: Poisson

function poisson_gamete_production(individuals::SubPopulation, number_of_inds_in_patch::Int64, location::Int64, fittnesses_by_genotype_first_half::Vector{Float64}, number_of_patches::Int64, carring_capacity::Int64,reproductive_rate::Float64)::BitVector
    #this function produces the gametes of the next generation
    #the "middle" of the habiat is at number_of_patches ÷ 2 + 0.5
    fittnesses_by_genotype::Vector{Float64} = (location <= number_of_patches ÷ 2) ? fittnesses_by_genotype_first_half :
                                                                                    reverse(fittnesses_by_genotype_first_half)
    individual_genotypes::Vector{Int64} = sum.(individuals) 
    n_AA::Int64 = sum(individual_genotypes .== 0) 
    n_Aa::Int64 = sum(individual_genotypes .== 1) 
    n_aa::Int64 = number_of_inds_in_patch - n_AA - n_Aa
    
    mean_offspring_per_genotype::Vector{Float64} = 2*exp(reproductive_rate*(1 - number_of_inds_in_patch/carring_capacity)) .* fittnesses_by_genotype

    relized_number_of_offspring_by_genotype::Vector{Int64} = [rand(Poisson((mean_offspring_per_genotype[1]  * n_AA))),
                                                              rand(Poisson((mean_offspring_per_genotype[2]  * n_Aa))),
                                                              rand(Poisson((mean_offspring_per_genotype[3]  * n_aa)))]
    relized_number_of_offspring::Int64 = sum(relized_number_of_offspring_by_genotype)
    p::Float64 = (relized_number_of_offspring_by_genotype[1] + relized_number_of_offspring_by_genotype[2]/2)/relized_number_of_offspring
    

    return rand(relized_number_of_offspring) .> p
end 

function poisson_mating(males::SubPopulation, females::SubPopulation, location::Int64, fittnesses_by_genotype_first_half::Vector{Float64}, number_of_patches::Int64, carring_capacity::Int64,reproductive_rate::Float64)::Tuple{SubPopulation,SubPopulation}
    number_of_inds_in_patch::Int64 = length(males) + length(females)
    male_gametes::BitVector  = poisson_gamete_production(males, number_of_inds_in_patch, location, fittnesses_by_genotype_first_half, number_of_patches, carring_capacity, reproductive_rate)
    female_gametes::BitVector = poisson_gamete_production(females, number_of_inds_in_patch, location, fittnesses_by_genotype_first_half, number_of_patches, carring_capacity, reproductive_rate)
    number_offspring::Int64 = min(length(male_gametes),length(female_gametes))
    offspring_genomes::SubPopulation = tuple.(male_gametes[1:number_offspring],female_gametes[1:number_offspring])
    
    females::SubPopulation = offspring_genomes[1:floor(Int,number_offspring/2)]
    males::SubPopulation = offspring_genomes[floor(Int,number_offspring/2)+1:end]
 
    return males, females
end


function age_patch(patch::Patch, fittnesses_by_genotype_first_half::Vector{Float64}, number_of_patches::Int64, carring_capacity::Int64, reproductive_rate::Float64)
    newborn_males, newborn_females = poisson_mating(vcat(patch.adult_males...), vcat(patch.adult_females...), 1, fittnesses_by_genotype_first_half, number_of_patches, carring_capacity, reproductive_rate)
    
    return Patch(newborn_males, newborn_females,
                 patch.newborn_males, patch.newborn_females,
                 tuple(patch.juvenile_males, patch.adult_males[1], patch.adult_males[2]), tuple(patch.juvenile_females, patch.adult_females[1], patch.adult_females[2]),
                 patch.location)
end


population = create_population(120,100)


selection_cofficient = 0.3
fittnesses_by_genotype_first_half = [1, 1-selection_cofficient, 1-2*selection_cofficient]
number_of_patches = 100
carring_capacity = 100
reproductive_rate = 1.2


function global_mating(population::Population)::Population
    adult_males::Vector{SubPopulation}= vcat.(eachcol(population.males[:,3:5])...)
    adult_females::Vector{SubPopulation} = vcat.(eachcol(population.females[:,3:5])...)

    offspring::Vector{Tuple{SubPopulation, SubPopulation}} = ((adult_males, adult_females) -> poisson_mating(adult_males, adult_females, 1, fittnesses_by_genotype_first_half,number_of_patches,carring_capacity,reproductive_rate)).(adult_males, adult_females)
    males_offspring::Vector{SubPopulation} = (offspring_i -> offspring_i[1]).(offspring)
    females_offspring::Vector{SubPopulation} = (offspring_i -> offspring_i[2]).(offspring)
    
    next_males::Matrix{SubPopulation} = hcat(males_offspring, population.males[:,1:4])
    next_females::Matrix{SubPopulation} = hcat(females_offspring, population.females[:,1:4])
    return Population(next_males, next_females)
end

end