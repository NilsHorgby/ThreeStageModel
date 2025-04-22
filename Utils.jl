module Utils
export predicted_allele_frequency
function predicted_allele_frequency(x::Float64, selection_cofficient::Float64, dispersal_distance::Float64)::Float64
    if x >= 0
        return -1/2 + 3/2 * tanh(sqrt(selection_cofficient/(2*dispersal_distance^2))*x+ atanh(sqrt(2/3)))^2
    else
        return 3/2 - 3/2 * tanh(-sqrt(selection_cofficient/(2*dispersal_distance^2))*x + atanh(sqrt(2/3)))^2
    end
end
#using .Main.Initialisation
#using PrettyTables
#=
export summary

function Base.summary(patch::Patch)::Nothing
    patch_info::Array{Union{String,Int64}} = ["Newborns"              length(patch.newborn_males)     length(patch.newborn_females);
                                              "Juveniles"             length(patch.juvenile_males)    length(patch.juvenile_females);
                                              "First Year Adults"     length(patch.adult_males[1])    length(patch.adult_females[1]);
                                              "Second Year Adults"    length(patch.adult_males[2])    length(patch.adult_females[2]);
                                              "Third Year Adults"     length(patch.adult_males[3])    length(patch.adult_females[3])]

    patch_info = vcat(patch_info,  reshape(vcat("Total", sum.(eachcol(patch_info[:,2:3]))), (1,3)))
    print(pretty_table(patch_info, header=["","Males", "Females"]))
end


function Base.summary(population::Population, axis::Int64=1, patch::Union{Int64,UndefInitializer}=undef)::Nothing
    if patch != undef
        let females::Vector{Int64} = length.(population.females[patch,:]),
        males::Vector{Int64} = length.(population.males[patch,:])
        col_names::Vector{String} = ["Newborns","Juveniles","First Year Adults","Second Year Adults","Third Year Adults"]
        patch_info::Array{Union{String,Int64}} = hcat(col_names, males, females)
        patch_info = vcat(patch_info,  reshape(vcat("Total", sum.(eachcol(patch_info[:,2:3]))), (1,3)))
        print(pretty_table(patch_info, header=["","Males", "Females"]))
        end
    else
        if axis == 1
            let females::Vector{Int64} = sum.(eachcol(length.(population.females))),
            males::Vector{Int64} = sum.(eachcol(length.(population.males)))
            col_names::Vector{String} = ["Newborns","Juveniles","First Year Adults","Second Year Adults","Third Year Adults"]
            patch_info::Array{Union{String,Int64}} = hcat(col_names, males, females)
            patch_info = vcat(patch_info,  reshape(vcat("Total", sum.(eachcol(patch_info[:,2:3]))), (1,3)))
            print(pretty_table(patch_info, header=["","Males", "Females"]))
            end
        elseif axis == 2
            population::Population = create_population(100,100)
            females::Vector{Int64} = sum.(eachrow(length.(population.females)))
            males::Vector{Int64} = sum.(eachrow(length.(population.males)))
            print(pretty_table(hcat(males, females), header=["Males", "Females"]))
        end
    end
end
=#
end