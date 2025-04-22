cd("/home/gadus/programing-projects/julia-projects/ThreeStageModel")
include("Initialisation.jl")
include("Mating.jl")
include("Migration.jl")

using .Initialisation, .Mating, .Migration
using ProgressMeter
using Plots
using StatsBase

function calc_genotype_counts!(population, arr::Array{Int16,4}, index::Int64)::Nothing
    arr[index, :, :, 1] .= ((x,y) -> sum(sum.(x) .== 0) + sum(sum.(y) .== 0)).(population.males,population.females)
    arr[index, :, :, 2] .= ((x,y) -> sum(sum.(x) .== 1) + sum(sum.(y) .== 1)).(population.males,population.females)
    arr[index, :, :, 3] .= ((x,y) -> sum(sum.(x) .== 2) + sum(sum.(y) .== 2)).(population.males,population.females)
    return nothing
end

function run_simulation(number_of_patches::Int64,
                        individuals_per_patch::Int64,
                        carring_capacity::Int64,
                        selection_cofficient::Float64,
                        dispersal_distance::Float64,
                        reproductive_rate::Float64,
                        generations::Int64)::Array{Int16,4}

    
    genotype_counts = Array{Int16}(undef, 100 + ceil(Int,(generations-100)/10), number_of_patches, 5, 3)

    fittnesses_by_genotype_first_half::Tuple{Float64,Float64,Float64} = (1, 1-selection_cofficient, 1-2selection_cofficient)
    population = create_population(number_of_patches, individuals_per_patch)
    index = 0 
    for gen in 1:generations
        if gen <= 100 || gen%10 == 0
            index += 1
            calc_genotype_counts!(population, genotype_counts, index)
        end
        age_population!(population,
                        fittnesses_by_genotype_first_half,
                        carring_capacity,
                        reproductive_rate,
                        number_of_patches)
        population = migration!(population, dispersal_distance)
    end
    return genotype_counts
end

using JLD2
function main()
    selection_coefficients::Vector{Float64} = [0.1,0.2,0.3]
    dispersal_distances::Vector{Float64} = [2., 4., 6.]
    #cline_widths::Vector{Float64} = [10,15,20]


    carring_capacity::Int64 = 200
    reproductive_rate::Float64 = 1.1
    number_of_patches::Int64 = 300
    individuals_per_patch::Int64 = 100
    generations::Int64 = 10_000
    n_runs::Int64 = 100
   
    for selection_coefficient in selection_coefficients, dispersal_distance in dispersal_distances
        if ~isfile("results/genotype_counts_s_$(selection_coefficient)_l_$(dispersal_distance)_r_$reproductive_rate.jld2")
        #dispersal_distance::Float64 = sqrt(cline_width ^ 2 * selection_coefficient/3)

        fittnesses_by_genotype_first_half::Tuple{Float64,Float64,Float64} = (1, 1-selection_coefficient, 1-2*selection_coefficient)
        genotype_counts::Array{Int16,5} = zeros(n_runs,
                                                100 + ceil(Int,( generations-100)/10),
                                                number_of_patches,
                                                5,
                                                3)
        Base.Threads.@threads for run in 1:n_runs
            genotype_counts[run,:,:,:,:] .= run_simulation(number_of_patches,
                                                        individuals_per_patch,
                                                        carring_capacity,
                                                        selection_coefficient,
                                                        dispersal_distance,
                                                        reproductive_rate,
                                                        generations)
        end
        save_object("results/genotype_counts_s_$(selection_coefficient)_l_$(dispersal_distance)_r_$reproductive_rate.jld2", genotype_counts)
        end
    end
end
main()

#=
allele_frequencies
mean.(eachrow(allele_frequencies[end,:,:]))

fitnesses_first_half = (allele_frequencies[ :, 1,1:150] .* (1-2selection_cofficient)  .+ allele_frequencies[ :, 2,1:150] .* (1-selection_cofficient) + allele_frequencies[ :, 3,1:150]) ./ sum.(eachslice(allele_frequencies, dims=(1,3)))[1:150,:]
fitnesses_second_half = (allele_frequencies[151:end, 1,:] + allele_frequencies[151:end, 2,:] .* (1-selection_cofficient) + allele_frequencies_over_time[151:end, 3,:] .* (1-2selection_cofficient)) ./ sum.(eachslice(allele_frequencies_over_time, dims=(1,3)))[151:end,:]
fitnesses = vcat(fitnesses_first_half, fitnesses_second_half)
mean_fitnesses = mean.(eachcol(fitnesses))
                                                                  
#plot mean pop size in last 100 gens
plot(sum.(eachrow(mean(eachslice(population_sizes[(end-99):end,:,:],dims=(1))))))
#plot population size over time

t = 1:20000

plot(t ./ 2.5, mean([sum.(eachrow(population_sizes[:,i,:])) for i=1:100]), xscale=:log10,
     xlabel = "Generations (log10 scale)",
     ylabel = "Population Size",
     label = "",
     title = "Population Size Over Time\nMultistage Model")
#plot allele frequencies of last gen
plot(mean.(eachrow(allele_frequencies[end,:,:])))

using StatsBase

plot(sum.(eachrow(population_sizes[end,:,:])))

using BenchmarkTools
#Profile.Allocs.@profile sample_rate = 0.01 run_simulation(number_of_patches,individuals_per_patch,carring_capacity,selection_cofficient,dispersal_distance,reproductive_rate, 100)
#=
@benchmark run_simulation(number_of_patches,individuals_per_patch,carring_capacity,selection_cofficient,dispersal_distance,reproductive_rate, 100)

BenchmarkTools.Trial: 46 samples with 1 evaluation.
 Range (min … max):   97.078 ms … 182.021 ms  ┊ GC (min … max):  8.28% … 28.82%
 Time  (median):     103.735 ms               ┊ GC (median):     9.28%
 Time  (mean ± σ):   109.966 ms ±  15.493 ms  ┊ GC (mean ± σ):  10.75% ±  5.25%

   █▄▃ ▁                                                         
  ▄███▇█▄▄▄▁▁▁▆▆▇▆▆▄▁▁▁▄▁▁▁▄▁▁▁▁▁▁▁▁▁▁▁▁▁▄▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▄ ▁
  97.1 ms          Histogram: frequency by time          182 ms <

 Memory estimate: 196.21 MiB, allocs estimate: 2569898.
=#
=#