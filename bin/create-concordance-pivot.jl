#!/usr/bin/env julia

# loading necessary packages
using CSV, DataFrames

# locate input files based on command line arguments
const cycle_table = islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]
const sample_manifest = islink(ARGS[2]) ? readlink(ARGS[2]) : ARGS[2]

# create a vector of detailed report files
const results_files = filter(x -> occursin("_Detailed_Results", x), readdir())

# find available cycles
available_cycles = Int64[]
for i in results_files
    cycle = parse(Int64, replace(split(split(i, ".")[1], "_")[end], r"cycle" => ""))
    push!(available_cycles, cycle)
end

# next data prep in a function, as per advice at https://docs.julialang.org/en/v1/manual/performance-tips/#kernel-functions
function prep_cycle_data(cycle_df::DataFrame, results_files::AbstractArray)
    
    # Read in the list of SNPs
    SNPs = CSV.read(results_files[1], DataFrame, skipto=17, header=16)
    snp_split = split.(SNPs.ID, '-')
    snp_list = sort(unique(getindex.(snp_split, 2)))

    # make an empty data frame that will subsequently be filled for each SNP
    concordance_table = DataFrame(SNP=snp_list, Assay="missing", Cycle=missing, Comments=missing)

    # determine optimal cycles
    sort!(cycle_df, [:SNP])
    concordance_table[:Cycle] = cycle_df[:cycle]
    concordance_table = filter(row -> row[:Cycle] != 0, concordance_table)
    cycle_df = filter(row -> row[:cycle] != 0, cycle_df)

    return concordance_table

end

function reorient_with_manifest(concordance_table::DataFrame, sample_manifest::String)

    # read in and format the sample manifest
    samples = CSV.read(sample_manifest, DataFrame, skipto=3, header=2)
    samples = dropmissing!(samples, :Relationship)

    # parse out the relationship "groups"
    groups = replace.(samples.Relationship, "Dam " => "")
    groups = replace.(groups, "Progeny " => "")
    group_int = parse.(Int64, groups)
    samples[:, :group] = group_int

    # go through each group and add them to new vectors that will eventually become
    # new column headers for the concordance pivot
    ordered_samples = []
    ordered_relations = []
    for group in unique(samples.group)

        # go group by group
        sub = filter(row -> group == row[:group], samples)

        # first add the dams
        append!(ordered_samples, filter(row -> occursin("Dam ", row[:Relationship]), sub)[Symbol("Animal ID")])
        append!(ordered_relations, filter(row -> occursin("Dam ", row[:Relationship]), sub)[:Relationship])

        # then add the progeny
        append!(ordered_samples, filter(row -> occursin("Progeny ", row[:Relationship]), sub)[Symbol("Animal ID")])
        append!(ordered_relations, filter(row -> occursin("Progeny ", row[:Relationship]), sub)[:Relationship])

    end
    
    # reorder columns in the concordance pivot
    ordered_cols = vcat(names(concordance_table)[1:4], ordered_samples)
    concordance_table = concordance_table[:, ordered_cols]

    # create new top row
    top_row = vcat("Animal ID", missing, missing, missing, ordered_samples)

    # create new second row and rejigger data frame
    row_two = vcat(names(concordance_table)[1:4], ordered_relations)
    for i in range(1,length(names(concordance_table)))
        old = names(concordance_table)[i]
        new = row_two[i]
        rename!(concordance_table, Symbol(old) => Symbol(new))
    end
    top_row = reshape(top_row, 1, length(top_row))
    row_two = reshape(row_two, 1, length(row_two))
    final_table = convert(Array, concordance_table)
    final_table = vcat(top_row, row_two, final_table)

    # make the matrix a data frame and write out
    final_table = DataFrame(final_table)

    return final_table

end

# nest loops in a function to improve performance
function generate_pivot(cycle_table::String, sample_manifest::String, results_files::AbstractArray)

    # Create an empty dictionary to store the results dataframes
    df_dict = Dict()

    # Read in the results files and store them in the list
    for file in results_files
        i = parse(Int64, String(split(split(file, ".")[1], "_cycle")[2]))
        df_dict[i] = CSV.read(file, DataFrame, skipto=17, header=16)
    end

    # read optimal cycles CSV into memory
    cycle_df = CSV.read(cycle_table, DataFrame)

    # prep data for loop
    concordance_table = prep_cycle_data(cycle_df, results_files)

    # Loop through each SNP to collate calls
    for i in range(1,nrow(concordance_table))

        # define SNP for this iteration
        snp = concordance_table[i,:SNP]
        
        # determine optimal cycle for that SNP
        optimal_cycle = cycle_df[i,:cycle]

        if !(optimal_cycle in available_cycles)
            continue
        end
    
        # Read in the correct BioMark output file based on that cycle
        data = df_dict[optimal_cycle]
        col_split = split.(data.ID, '-')
        data.sample = getindex.(col_split, 1)
        data.SNP = getindex.(col_split, 2)

        # make table of SNPs
        snp_table = filter(row -> occursin(snp, row[:SNP]), data)

        # Record SNP assay code in the concordance pivot
        concordance_table[i,:Assay] = unique(snp_table[:Assay])[1]

        if snp == concordance_table[:SNP][1]
            for animal in unique(sort(snp_table[:Name]))
                call = filter(row -> occursin(animal, row[:Name]), snp_table)[:Final][1]
                concordance_table[:, Symbol(animal)] = "No Call"
                concordance_table[i, Symbol(animal)] = call
            end
        else
            for animal in unique(sort(snp_table[:Name]))
                call = filter(row -> occursin(animal, row[:Name]), snp_table)[:Final][1]
                concordance_table[i, Symbol(animal)] = call
            end
        end
    end

    # change assay names to SNP names
    concordance_table.SNP = replace.(concordance_table.SNP, "A" => "SNP")

    # reorient with sample manifest
    final_table = reorient_with_manifest(concordance_table, sample_manifest)

    # export to CSV-formatted spreadsheet
    CSV.write("concordance-pivot.csv", final_table, writeheader=false)

    # return final_table

end

# Generate the pivot using the function we just defined
generate_pivot(cycle_table, sample_manifest, results_files)
