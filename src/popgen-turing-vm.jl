using Distributions
using Turing

const L2_L3_GENOTYPES = [
    "VV-FF",
    "VV-FC",
    "VV-CC",
    "VI-FF",
    "VI-FC",
    "VI-CC",
    "II-FF",
    "II-FC",
    "II-CC"
]

const HAPLOTYPES_1016_1534 = ["VF", "VC", "IF", "IC"]

"""
Convert nine two-locus diploid genotype counts into expected counts
of the four 1016/1534 haplotypes: VF, VC, IF, and IC.

The phase of VI-FC is unobserved. Its count is split equally between
the VF/IC and VC/IF phases.
"""
function genotype_counts_to_haplotype_counts(counts)
    length(counts) == length(L2_L3_GENOTYPES) ||
        throw(DimensionMismatch("counts must have length 9."))

    VV_FF, VV_FC, VV_CC,
    VI_FF, VI_FC, VI_CC,
    II_FF, II_FC, II_CC = counts

    return [
        2 * VV_FF + VV_FC + VI_FF + 0.5 * VI_FC,
        2 * VV_CC + VV_FC + VI_CC + 0.5 * VI_FC,
        2 * II_FF + VI_FF + II_FC + 0.5 * VI_FC,
        2 * II_CC + VI_CC + II_FC + 0.5 * VI_FC
    ]
end

"""
Construct eight three-locus pseudo-counts from four two-locus counts.

The 410 locus is unobserved, so each 1016/1534 haplotype is split
equally across the V and L backgrounds.
"""
function haplotype_counts_to_three_locus_counts(counts)
    two_locus_counts = genotype_counts_to_haplotype_counts(counts)
    return vcat(two_locus_counts ./ 2, two_locus_counts ./ 2)
end

@model function vm_fitness_priors(
    n_sites,
    ::Type{T}=Float64
) where {T}
    selection = Matrix{T}(undef, n_sites, N_HAPLOTYPES)
    dominance_raw = Matrix{T}(undef, n_sites, N_HAPLOTYPES)
    dominance = Matrix{T}(undef, n_sites, N_HAPLOTYPES)

    selection_phi = Vector{T}(undef, N_HAPLOTYPES)
    selection_lambda = Vector{T}(undef, N_HAPLOTYPES)
    dominance_phi = Vector{T}(undef, N_HAPLOTYPES)
    dominance_lambda = Vector{T}(undef, N_HAPLOTYPES)

    for haplotype in 1:N_HAPLOTYPES
        selection_phi[haplotype] ~ Beta(1, 1)
        selection_lambda[haplotype] ~ Pareto(1, 1.5)

        dominance_phi[haplotype] ~ Beta(1, 1)
        dominance_lambda[haplotype] ~ Pareto(1, 1.5)
    end

    for site in 1:n_sites
        for haplotype in 1:N_HAPLOTYPES
            selection[site, haplotype] ~ Beta(
                selection_lambda[haplotype] * selection_phi[haplotype],
                selection_lambda[haplotype] * (
                    1 - selection_phi[haplotype]
                )
            )

            dominance_raw[site, haplotype] ~ Beta(
                dominance_lambda[haplotype] * dominance_phi[haplotype],
                dominance_lambda[haplotype] * (
                    1 - dominance_phi[haplotype]
                )
            )

            dominance[site, haplotype] =
                2 * dominance_raw[site, haplotype] - 1
        end
    end

    return (selection = selection, dominance = dominance)
end

"""
Convert the four 1016/1534 joint haplotype frequencies into
nine diploid genotype probabilities.
"""
function haplotype23_to_genotype_probs(x23)
    length(x23) == 4 || throw(DimensionMismatch("x23 must have length 4."))

    VF, VC, IF, IC = x23

    probabilities = [
        VF^2,
        2 * VF * VC,
        VC^2,
        2 * VF * IF,
        2 * VF * IC + 2 * VC * IF,
        2 * VC * IC,
        IF^2,
        2 * IF * IC,
        IC^2
    ]

    return probabilities
end

function haplotype23_to_genotype_probs!(probabilities, x23)
    length(probabilities) == length(L2_L3_GENOTYPES) ||
        throw(DimensionMismatch("probabilities must have length 9."))
    length(x23) == 4 || throw(DimensionMismatch("x23 must have length 4."))

    VF, VC, IF, IC = x23

    probabilities[1] = VF^2
    probabilities[2] = 2 * VF * VC
    probabilities[3] = VC^2
    probabilities[4] = 2 * VF * IF
    probabilities[5] = 2 * VF * IC + 2 * VC * IF
    probabilities[6] = 2 * VC * IC
    probabilities[7] = IF^2
    probabilities[8] = 2 * IF * IC
    probabilities[9] = IC^2

    return probabilities
end

"""
Map one three-locus haplotype distribution to the observed two-locus
1016/1534 genotype probabilities.
"""
function vm_two_locus_genotype_probs(x)
    return haplotype23_to_genotype_probs(
        get_1016_1534_frequencies(x)
    )
end

function vm_two_locus_genotype_probs!(probabilities, x)
    VF = x[1] + x[5]
    VC = x[2] + x[6]
    IF = x[3] + x[7]
    IC = x[4] + x[8]

    probabilities[1] = VF^2
    probabilities[2] = 2 * VF * VC
    probabilities[3] = VC^2
    probabilities[4] = 2 * VF * IF
    probabilities[5] = 2 * VF * IC + 2 * VC * IF
    probabilities[6] = 2 * VC * IC
    probabilities[7] = IF^2
    probabilities[8] = 2 * IF * IC
    probabilities[9] = IC^2

    return probabilities
end

"""
Turing model for the two-locus VM observations.

The latent population is three-locus, while observations contain only
loci 1016 and 1534. Fitness and recombination parameters are shared by
replicates within each site, while recombination parameters are globally shared across sites; 
each site/replicate has its own initial three-locus haplotype distribution.
"""
@model function vm_two_locus_model(
    counts,
    sample_sizes,
    n_generations,
    group_site,
    n_sites,
    n_groups,
    ::Type{T}=Float64
) where {T}
    # n_genotypes, n_generation_columns, n_groups_in_counts = size(counts)
    # n_groups_in_counts == n_groups ||
    #     throw(DimensionMismatch("counts has an unexpected group dimension."))

    rAB ~ Beta(1, 3)
    rBC ~ Beta(1, 3)
    delta ~ truncated(Normal(0, 3), 0, Inf)

    fitness_params ~ to_submodel(vm_fitness_priors(n_sites), false)
    selection = fitness_params.selection
    dominance = fitness_params.dominance
    initial_haplotypes = Vector{Vector{T}}(undef, n_groups)
    for group in 1:n_groups
        initial_haplotype_counts = haplotype_counts_to_three_locus_counts(
            counts[:, 1, group]
        )
        initial_haplotypes[group] ~ Dirichlet(
            T.(initial_haplotype_counts) .+ delta
        )
    end

    site_models = [
        build_model(
            build_fitness_matrix(
                selection[site, :],
                dominance[site, :];
                checks = false
            );
            rAB = rAB,
            rBC = rBC,
            checks = false
        )
        for site in 1:n_sites
    ]

    predicted = similar(counts, T)

    for group in 1:n_groups
        site = group_site[group]
        popgen_model = site_models[site]

        x = copy(initial_haplotypes[group])
        x_next = similar(x)
        probabilities = similar(x, length(L2_L3_GENOTYPES))

        for generation in 1:n_generations[group]
            vm_two_locus_genotype_probs!(probabilities, x)
            predicted[:, generation, group] .= probabilities

            counts[:, generation, group] ~ Multinomial(
                sample_sizes[generation, group],
                probabilities
            )

            if generation < n_generations[group]
                step!(x_next, x, popgen_model; checks = false)
                x, x_next = x_next, x
            end
        end
    end

    return predicted
end

"""
Turing model for the two-locus VM observations.

The latent population is three-locus, while observations contain only
loci 1016 and 1534. 

If multiple sites are given, fitness and recombination parameters are shared by
replicates within each site, while recombination parameters are globally shared across sites; 
each site/replicate has its own initial three-locus haplotype distribution.
"""
@model function vm_two_locus_model(
    counts,
    sample_sizes,
    generations,
    ::Type{T}=Float64
) where {T}
    # n_genotypes, n_generation_columns, n_groups_in_counts = size(counts)
    # n_groups_in_counts == n_groups ||
    #     throw(DimensionMismatch("counts has an unexpected group dimension."))

    rAB ~ Beta(1, 3)
    rBC ~ Beta(1, 3)
    delta ~ truncated(Normal(0, 3), 0, Inf)

    selection ~ filldist(Beta(1, 1), N_HAPLOTYPES-1)
    pushfirst!(selection, 0)

    dominance ~ filldist(Beta(1, 1), N_HAPLOTYPES-1)
    dominance = 2 * dominance_raw[site, haplotype] - 1
    pushfirst!(dominance, 0)

    for haplotype in 1:N_HAPLOTYPES
        selection_phi[haplotype] ~ Beta(1, 1)
        selection_lambda[haplotype] ~ Pareto(1, 1.5)

        dominance_phi[haplotype] ~ Beta(1, 1)
        dominance_lambda[haplotype] ~ Pareto(1, 1.5)
    end

    for site in 1:n_sites
        for haplotype in 1:N_HAPLOTYPES
            selection[site, haplotype] ~ Beta(
                selection_lambda[haplotype] * selection_phi[haplotype],
                selection_lambda[haplotype] * (
                    1 - selection_phi[haplotype]
                )
            )

            dominance_raw[site, haplotype] ~ Beta(
                dominance_lambda[haplotype] * dominance_phi[haplotype],
                dominance_lambda[haplotype] * (
                    1 - dominance_phi[haplotype]
                )
            )

            
        end
    end

    fitness_params ~ to_submodel(vm_fitness_priors(n_sites), false)
    selection = fitness_params.selection
    dominance = fitness_params.dominance
    initial_haplotypes = Vector{Vector{T}}(undef, n_groups)
    for group in 1:n_groups
        initial_haplotype_counts = haplotype_counts_to_three_locus_counts(
            counts[:, 1, group]
        )
        initial_haplotypes[group] ~ Dirichlet(
            T.(initial_haplotype_counts) .+ delta
        )
    end

    site_models = [
        build_model(
            build_fitness_matrix(
                selection[site, :],
                dominance[site, :];
                checks = false
            );
            rAB = rAB,
            rBC = rBC,
            checks = false
        )
        for site in 1:n_sites
    ]

    predicted = similar(counts, T)

    for group in 1:n_groups
        site = group_site[group]
        popgen_model = site_models[site]

        x = copy(initial_haplotypes[group])
        x_next = similar(x)
        probabilities = similar(x, length(L2_L3_GENOTYPES))

        for generation in 1:n_generations[group]
            vm_two_locus_genotype_probs!(probabilities, x)
            predicted[:, generation, group] .= probabilities

            counts[:, generation, group] ~ Multinomial(
                sample_sizes[generation, group],
                probabilities
            )

            if generation < n_generations[group]
                step!(x_next, x, popgen_model; checks = false)
                x, x_next = x_next, x
            end
        end
    end

    return predicted
end
