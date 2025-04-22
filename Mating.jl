#=
I comment each function before the function, and signify what line I am takning about with a superscript letter
=#

module Mating
using .Main.Initialisation
using Distributions: Poisson
using Random: bitrand

export age_population!, poisson_gamete_production
#=
in: Vector of individuals ¹
out: Vector of gametes ²

•This function produces the gametes of the next generation²
•The "middle" of the habiat is at number_of_patches ÷ 2 + 0.5³
•Assumes that each individual produces a Poisson-distuted number of gametes⁴
•The total number of gametes is a Poisson-distrubuted random variable.⁵

=#
function poisson_gamete_production(individuals::SubPopulation, #¹
                                   effective_number_of_inds_in_patch::Real,
                                   location::I,
                                   fittnesses_by_genotype_first_half::Tuple{Float64,Float64,Float64},
                                   number_of_patches::I,
                                   carring_capacity::I,
                                   reproductive_rate::Float64)::BitVector where I<:Integer#² 
    
    fittnesses_by_genotype::Tuple{Float64,Float64,Float64} = (location <= number_of_patches/2) ? fittnesses_by_genotype_first_half :
                                                                                                 reverse(fittnesses_by_genotype_first_half)#³
    individual_genotypes::Vector{I} = sum.(individuals) 
    n_AA::I = sum(individual_genotypes .== 0) 
    n_Aa::I = sum(individual_genotypes .== 1)  
    n_aa::I = length(individual_genotypes) - n_AA - n_Aa

    mean_offspring_per_genotype::Tuple{Float64,Float64,Float64} = 2*exp(reproductive_rate*(1 - effective_number_of_inds_in_patch/carring_capacity))/3 .* fittnesses_by_genotype
    mean_offspring::Float64 = sum(mean_offspring_per_genotype .* (n_AA,n_Aa,n_aa))#⁴
    relized_number_of_offspring::I = rand(Poisson(mean_offspring))#⁵
    fittness_adjusted_ps::Tuple{Float64,Float64,Float64} = fittnesses_by_genotype .* (n_AA, n_Aa, n_aa)
    p::Float64 = (fittness_adjusted_ps[1] + fittness_adjusted_ps[2]*0.5) / sum(fittness_adjusted_ps)

    return rand(relized_number_of_offspring) .> p #²
end 


#=
in: females, males¹
out: female offspring, male offspring ²

•the total number of offspring is the minimum of the number of gametes from males and females ³
•the number of male offspring is a binomially-distrubuted random variable with p = 0.5 and N = the number of offspring ⁴
•the number of female offspring is (the number of offspring - the number of male offspring)⁵
=#
function mating(adult_males::Vector{SubPopulation},#¹
                adult_females::Vector{SubPopulation},#¹
                number_of_inds_by_patch::Vector{I},
                fittnesses_by_genotype_first_half::Tuple{Float64,Float64,Float64},
                number_of_patches::I, 
                carring_capacity::I,
                reproductive_rate::Float64)::Tuple{Vector{SubPopulation},Vector{SubPopulation}} where I<:Integer#²
    curried_gamete_production(individuals::SubPopulation, number_of_inds_in_patch::I, location::I)::BitVector = poisson_gamete_production(individuals,
                                                                                                                                          number_of_inds_in_patch,
                                                                                                                                          location,
                                                                                                                                          fittnesses_by_genotype_first_half,
                                                                                                                                          number_of_patches,
                                                                                                                                          carring_capacity,
                                                                                                                                          reproductive_rate)
    #300 element vector of BitVectors with varing lenght
    male_gametes::Vector{BitVector} = curried_gamete_production.(adult_males,number_of_inds_by_patch, I.(1:300))
    female_gametes::Vector{BitVector} = curried_gamete_production.(adult_females,number_of_inds_by_patch, I.(1:300))
    #300 element vector of Ints - preallocate
    number_offspring_by_patch::Vector{I} = min.(length.(male_gametes),length.(female_gametes)) #³
    
    mitosis_one_patch(male_gametes_one_patch::BitVector,female_gametes_one_patch::BitVector,number_offspring::I)::SubPopulation = tuple.(male_gametes_one_patch[1:number_offspring],female_gametes_one_patch[1:number_offspring])
    
    offspring_by_patch::Vector{SubPopulation} = mitosis_one_patch.(male_gametes, female_gametes, number_offspring_by_patch)
    #300 element vector of Ints - preallocate
    sex_indices::Vector{I} = sum.(bitrand.(number_offspring_by_patch))

    males_offspring::Vector{SubPopulation} = ((offspring,sex_index) -> offspring[1:sex_index]).(offspring_by_patch, sex_indices)#⁴
    females_offspring::Vector{SubPopulation} = ((offspring,sex_index) -> offspring[sex_index+1:end]).(offspring_by_patch, sex_indices)#⁵

    return males_offspring, females_offspring #²
end

#=
in: The whole Population¹
out: nothing (mutates input)²

•Ages population one time step, including recruitment 
=#
function age_population!(population::Population,#¹
                         fittnesses_by_genotype_first_half::Tuple{Float64,Float64,Float64},
                         carring_capacity::I,
                         reproductive_rate::Float64,
                         number_of_patches::I)::Nothing where I<:Integer#²

    
    adult_males::Vector{SubPopulation}= vcat.(eachcol(population.males[:,3:5])...)
    adult_females::Vector{SubPopulation} = vcat.(eachcol(population.females[:,3:5])...)
    number_of_inds_by_patch::Vector{I} = sum.(eachrow(length.(population.males[:,3:5]))) .+ sum.(eachrow(length.(population.females[:,3:5])))
    population.males[:,2:5], population.females[:,2:5] = population.males[:,1:4], population.females[:,1:4]
    population.males[:,1], population.females[:,1] = mating(adult_males, adult_females, number_of_inds_by_patch, fittnesses_by_genotype_first_half, number_of_patches, carring_capacity, reproductive_rate)
    return nothing
end

end