module Migration

using .Main.Initialisation
using Distributions: Normal

export migration!

#=
in: A Patch number
out: A patch number placed between 0 and N

•Guarantees that the patch number is between 0 and N 
=#
function put_within_boundaries(new_patch::I, number_of_patches::I) where I<:Integer
    if new_patch <= 0
        return 1
    elseif new_patch > number_of_patches
        return number_of_patches
    else 
        return new_patch
    end
end

#=
in: params
out: A random integer distrubuted according to dicretized normal distrubtion
=#
function discrete_distrubution(dispersal_distance::Real, location::I)::I where I<:Integer
    return round(Int,rand(Normal(location,dispersal_distance)))
end

#=
in: one age group of a population
out: the same age group after migration

•Note, mutates input as side-effect. 

there should be a better way of doing this. Idealy we would not create a new population array. Creating a new next_population array is more difficult to avoid. 
=#
function migrate!(population::AbstractArray{SubPopulation},
                  next_population::Vector{SubPopulation},
                  dispersal_distance::Float64)::Vector{SubPopulation}
    number_of_patches::Integer = length(population)
    @inbounds for location in 1:number_of_patches
        @inbounds for individual_index in 1:length(population[location]) #see if a while loop is faster
            next_location::Integer = put_within_boundaries(discrete_distrubution(dispersal_distance,location), number_of_patches)
            push!(next_population[next_location],population[location][individual_index])
        end
    end
    #population .= next_population
    return next_population
end 



#=
in: A population
out: Population after migration

•Note, mutates input
=#
function migration!(population::Population, dispersal_distance::Float64)::Population
    #newborns don't migrate
    next_males = [[Vector{Tuple{Bool,Bool}}() for i=eachindex(population.males[:,1])] for j=1:5]

    next_males = hcat(migrate!.(eachcol(population.males),next_males,dispersal_distance)...)

    next_females = [[Vector{Tuple{Bool,Bool}}() for i=eachindex(population.females[:,1])] for j=1:5]
    next_females = hcat(migrate!.(eachcol(population.females),next_females,dispersal_distance)...)

    return Population(next_males,next_females)
end

function migration_fast!(population::Population, dispersal_distance::Float64)::Nothing #Population
    #newborns don't migrate
    next_males = [[Vector{Tuple{Bool,Bool}}() for i=eachindex(population.males[:,1])] for j=1:5]

    next_males = hcat(migrate!.(eachcol(population.males),next_males,dispersal_distance)...)

    next_females = [[Vector{Tuple{Bool,Bool}}() for i=eachindex(population.females[:,1])] for j=1:5]
    next_females = hcat(migrate!.(eachcol(population.females),next_females,dispersal_distance)...)
    population.males .= next_males
    population.females .= next_females
    return nothing
    #return Population(next_males,next_females)
end
#=
population = create_population(300,100)
@benchmark migration!(population, 1.)

BenchmarkTools.Trial: 8803 samples with 1 evaluation.
 Range (min … max):  447.980 μs …   8.464 ms  ┊ GC (min … max): 0.00% … 77.49%
 Time  (median):     507.924 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   566.649 μs ± 356.257 μs  ┊ GC (mean ± σ):  8.33% ± 11.04%

  ▄█▆▃▂▁                                                        ▁
  ██████▇▆▆▄▄▄▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▃▃▁▅▆▄▅▆▆▇▇▇▆ █
  448 μs        Histogram: log(frequency) by time       2.69 ms <

 Memory estimate: 630.44 KiB, allocs estimate: 9034.
=#


#= test
population = create_population(300,100)

test_A = collect(zip(repeat([false],50),repeat([false],50)))
test_a = collect(zip(repeat([true],50),repeat([true],50)))
test_males = migrate!([(i <= 150) ? test_A : test_a for i=1:300],[Vector{Tuple{Bool,Bool}}() for i=1:300], 1.0 )
test_males = migrate!(test_males,[Vector{Tuple{Bool,Bool}}() for i=1:300], 1.0)
using Plots
plot(sum.((x -> sum.(x)).(test_males)))

=#
end