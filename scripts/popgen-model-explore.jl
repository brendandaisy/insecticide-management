using DataFrames
using CategoricalArrays
using CairoMakie
using AlgebraOfGraphics

include("../src/popgen-three-locus.jl")
include("../src/popgen-turing-vm.jl")

function two_locus_trajectory(
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

    for generation in 1:generations
        probabilities = vm_two_locus_genotype_probs(
            simulation[:, generation]
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

# ============================================================
# Exploring the effect of the recombination parameters
# ============================================================
# outcomes are certainly similar, but not definitively exchangeable. Proceeding
# to a proper practical identifiability analysis

selection_neutral = fill(0., 8)
dominance_neutral = fill(1., 8)
Wneutral = build_fitness_matrix(selection_neutral, dominance_neutral)

selection_guess = [0, 0.05, 0.95, 0.8, 0.2, 0.15, 0.99, 0.85]
dominance_guess = [0, 0.5, 0.95, 0.1, 0.5, 0.5, 0.95, 0.3]
Wguess = build_fitness_matrix(selection_guess, dominance_guess)

xinit = [1, 1, 1, 3, 1, 2, 1, 10]
xinit = xinit ./ sum(xinit)

scenario1 = two_locus_trajectory(
    selection_neutral, dominance_neutral, xinit;
    rAB=0.1, rBC=0.9,
    scenario="W neutral, rAB low, rBC high"
)

scenario2 = two_locus_trajectory(
    selection_neutral, dominance_neutral, xinit;
    rAB=0.9, rBC=0.1,
    scenario="W neutral, rAB high, rBC low"
)

scenario3 = two_locus_trajectory(
    selection_guess, dominance_guess, xinit;
    rAB=0.1, rBC=0.9,
    scenario="W guess, rAB low, rBC high"
)

scenario4 = two_locus_trajectory(
    selection_guess, dominance_guess, xinit;
    rAB=0.9, rBC=0.1,
    scenario="W guess, rAB high, rBC low"
)

trajectories = reduce(
    vcat,
    first.([scenario1, scenario2, scenario3, scenario4])
)

spec = data(trajectories) *
    mapping(:generation, :proportion; color=:genotype, layout=:scenario) *
    visual(ScatterLines)

draw(spec, scales(Color=(;palette=:tab10)))

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
