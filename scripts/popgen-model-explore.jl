using DataFrames
using CategoricalArrays
using CairoMakie
using AlgebraOfGraphics

include("../src/popgen-three-locus.jl")
include("../src/popgen-turing-vm.jl")



# ============================================================
# Plotting deterministic model trajectories of different
# fitness effects at L2 and L3
# ============================================================
# - "is there a scenario where double het VI-FC can increase?"
# where the demoted pair doesn
# - yes, but (requires?) a scenario where either VF and IC are the best, or 
# VC and IF are the best. Since we expect both IC and IF to be unfit,
# we should expect VI-FC to decrease in the data
# - well, VI-FC actually is overall stable, but this could be from short generation time or 
# from IC being recessive (for example, VI-CC also appear stable)

function deterministic_trajectory(
    selection,
    dominance,
    x0;
    rAB = 0.05,
    rBC = 0.05,
    generations = 20,
    scenario = "scenario"
)
    W = build_fitness_matrix(
        selection,
        dominance;
        checks = false
    )
    model = build_model(
        W;
        rAB = rAB,
        rBC = rBC,
        checks = false
    )
    simulation = simulate(
        model,
        x0,
        generations;
        checks = true
    )

    records = DataFrame(
        scenario = String[],
        generation = Int[],
        genotype = String[],
        proportion = Float64[]
    )

    for generation in 0:generations
        probabilities = vm_two_locus_genotype_probs(
            simulation.trajectory[:, generation + 1]
        )
        for genotype_index in eachindex(L2_L3_GENOTYPES)
            push!(
                records,
                (
                    scenario,
                    generation,
                    L2_L3_GENOTYPES[genotype_index],
                    probabilities[genotype_index]
                )
            )
        end
    end

    return records, W
end

# Haplotype order: VVF, VVC, VIF, VIC, LVF, LVC, LIF, LIC.
# A uniform initial distribution makes the initial two-locus genotype
# distribution easy to interpret and keeps scenarios comparable.
x0 = fill(1 / N_HAPLOTYPES, N_HAPLOTYPES)

# Selection is interpreted through W[i,j] = 1 - h[i]s[i] - h[j]s[j].
# For s > 0:
#   h =  1: the haplotype is selected against in homozygotes and heterozygotes
#   h =  0: no heterozygote effect
#   h = -1: overdominant effect relative to the selected haplotype
scenarios = [
    # (
    #     name = "VIF/VIC, LIF/LIC, h = 1",
    #     selection = [0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 0.8, 0.8],
    #     dominance = [0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0]
    # ),
    (
        name = "VIF/VIC, LIF/LIC, h = 0",
        selection = [0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 0.8, 0.8],
        dominance = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    ),
    # (
    #     name = "VIF/VIC, LIF/LIC, h = -1",
    #     selection = [0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 0.8, 0.8],
    #     dominance = [0.0, 0.0, -1.0, -1.0, 0.0, 0.0, -1.0, -1.0]
    # ),
    (
        name = "VC and IF fatal, IC nuetral",
        selection = [0.0, 1, 1, 0, 0.0, 1, 1, 0],
        dominance = [0.0, 0.0, 0.1, -1.0, 0.0, 0.0, 0.1, -1.0]
    ),
    (
        name = "VC and IF fatal, IC overdominant",
        selection = [0.0, 0.0, 1, 0.5, 0.0, 0.0, 1, 0.6],
        dominance = [0.0, 0.0, 0.1, -1, 0.0, 0.0, 0.1, -1]
    ),
    (
        name = "VC and IF fatal, IC additive",
        selection = [0.0, 0.0, 1, 0.5, 0.0, 0.0, 1, 0.6],
        dominance = [0.0, 0.0, 0.1, 0.5, 0.0, 0.0, 0.1, 0.5]
    )
]

build_fitness_matrix(scenarios[3].selection, scenarios[3].dominance; checks=false)

trajectory_tables = DataFrame[]
fitness_matrices = Dict{String, Matrix{Float64}}()
for scenario in scenarios
    trajectory, W = deterministic_trajectory(
        scenario.selection,
        scenario.dominance,
        x0;
        scenario = scenario.name,
        generations = 20
    )
    push!(trajectory_tables, trajectory)
    fitness_matrices[scenario.name] = W
end

trajectory_data = reduce(vcat, trajectory_tables)
trajectory_data.genotype = categorical(
    trajectory_data.genotype;
    levels = L2_L3_GENOTYPES,
    ordered = true
)

spec = data(trajectory_data) *
    mapping(
        :generation,
        :proportion;
        color = :genotype,
        col = :genotype,
        row = :scenario
    ) *
    visual(Lines)

figure = draw(
    spec;
    figure = (; size = (1400, 1100)),
    axis = (; limits = (nothing, (0, 1))),
    facet = (; linkxaxes = :all, linkyaxes = :all)
)
