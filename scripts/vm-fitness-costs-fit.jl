using Distributions
using InvertedIndices
using CSV, DataFrames, DataFramesMeta
using CairoMakie, AlgebraOfGraphics
using Turing, FlexiChains
using Statistics
using Serialization
# using HDF5, MCMCChainsStorage

include("../src/popgen-three-locus.jl")
include("../src/popgen-turing-vm.jl")

vm_two_loci = CSV.read("data-proc/vm-geno-counts-1016-1534.csv", DataFrame)

function get_generations(data)
	grouped = unique(select(data, [:site, :rep, :generation]))
	sort!(grouped, [:site, :rep, :generation])

	group_generations = [
		Int.(group.generation)
		for group in groupby(grouped, [:site, :rep]; sort=false)
	]

	all_generations = Int.(grouped.generation)

	return (
		values = all_generations,
		by_group = group_generations,
		n = maximum.(group_generations)
	)
end

function prepare_vm_two_loci(data)
	required_columns = [:site, :rep, :generation, :genotype, :N, :count]
	all(column -> column in propertynames(data), required_columns) ||
		throw(ArgumentError("VM data is missing a required column."))

	genotype_index = Dict(
		genotype => index
		for (index, genotype) in enumerate(L2_L3_GENOTYPES)
	)

	site_levels = sort!(unique(data.site))
	site_index = Dict(site => index for (index, site) in enumerate(site_levels))

	group_keys = unique(select(data, [:site, :rep]))
	sort!(group_keys, [:site, :rep])
	group_index = Dict(
		(row.site, row.rep) => index
		for (index, row) in enumerate(eachrow(group_keys))
	)

	grouped = unique(select(data, [:site, :rep, :generation]))
	sort!(grouped, [:site, :rep, :generation])
	generation_data = get_generations(data)
	max_generation = maximum(generation_data.n)

	counts = zeros(Int, length(L2_L3_GENOTYPES), max_generation, nrow(group_keys))
	sample_sizes = zeros(Int, max_generation, nrow(group_keys))

	for key in eachrow(grouped)
		group = group_index[(key.site, key.rep)]
		generation = Int(key.generation)
		rows = data[
			(data.site .== key.site) .&
			(data.rep .== key.rep) .&
			(data.generation .== key.generation),
			:
		]

		nrow(rows) == length(L2_L3_GENOTYPES) ||
			throw(ArgumentError(
				"Each VM site/replicate/generation must have exactly nine genotype rows."
			))

		seen = falses(length(L2_L3_GENOTYPES))
		for row in eachrow(rows)
			index = get(genotype_index, String(row.genotype), 0)
			index == 0 && throw(ArgumentError(
				"Unexpected two-locus genotype: $(row.genotype)."
			))
			seen[index] && throw(ArgumentError(
				"Duplicate genotype $(row.genotype) in VM data."
			))
			seen[index] = true
			counts[index, generation, group] = Int(row.count)
			sample_sizes[generation, group] = Int(row.N)
		end

		all(seen) || throw(ArgumentError(
			"A VM observation is missing at least one genotype category."
		))
		sum(counts[:, generation, group]) == sample_sizes[generation, group] ||
			throw(ArgumentError(
				"VM genotype counts do not sum to N for " *
				"$(key.site), replicate $(key.rep), generation $(key.generation)."
			))
	end

	return (
		counts = counts,
		sample_sizes = sample_sizes,
		n_generations = generation_data.n,
		group_site = [site_index[row.site] for row in eachrow(group_keys)],
		n_sites = length(site_levels),
		n_groups = nrow(group_keys),
		site_levels = site_levels,
		group_keys = group_keys
	)
end

vm_data = prepare_vm_two_loci(vm_two_loci)

vm_model = vm_two_locus_model(
	vm_data.counts,
	vm_data.sample_sizes,
	vm_data.n_generations,
	vm_data.group_site,
	vm_data.n_sites,
	vm_data.n_groups
)

# ============================================================
# Prior predictive check
# ============================================================

function plot_predictions(predictions, site_rep_group)
    site = vm_data.group_keys[!, :site][site_rep_group]
    replicate = vm_data.group_keys[!, :rep][site_rep_group]
    n_generations = vm_data.n_generations[site_rep_group]

    predictive_draws = [
        predictions[iteration, 1][
            :,
            1:n_generations,
            site_rep_group
        ]
        for iteration in axes(predictions, 1)
    ]

    predictive_mean = reduce(
        .+,
        predictive_draws
    ) ./ length(predictive_draws)

    predictive_lower = similar(predictive_mean)
    predictive_upper = similar(predictive_mean)
    for genotype_index in axes(predictive_mean, 1)
        for generation_index in axes(predictive_mean, 2)
            values = [
                draw[genotype_index, generation_index]
                for draw in predictive_draws
            ]
            predictive_lower[genotype_index, generation_index] =
                quantile(values, 0.025)
            predictive_upper[genotype_index, generation_index] =
                quantile(values, 0.975)
        end
    end

    observed_counts = vm_data.counts[
        :,
        1:n_generations,
        site_rep_group
    ]
    observed_sizes = vm_data.sample_sizes[
        1:n_generations,
        site_rep_group
    ]
    observed_proportions = observed_counts ./ reshape(observed_sizes, 1, :)

    predictive_summary = DataFrame(
        generation = repeat(
            1:n_generations,
            inner = length(L2_L3_GENOTYPES)
        ),
        genotype = repeat(
            L2_L3_GENOTYPES,
            outer = n_generations
        ),
        mean = vec(predictive_mean),
        lower = vec(predictive_lower),
        upper = vec(predictive_upper)
    )

    observed_summary = DataFrame(
        generation = repeat(
            1:n_generations,
            inner = length(L2_L3_GENOTYPES)
        ),
        genotype = repeat(
            L2_L3_GENOTYPES,
            outer = n_generations
        ),
        proportion = vec(observed_proportions)
    )

    prior_spec = data(predictive_summary) *
        (mapping(
            :generation,
            :lower,
            :upper,
            color = :genotype,
            col = :genotype
        ) * visual(Band; alpha = 0.25) +
        mapping(
            :generation,
            :mean,
            color = :genotype,
            col = :genotype
        ) * visual(Lines; linewidth = 2))

    observed_spec = data(observed_summary) *
        mapping(
            :generation,
            :proportion,
            color = :genotype,
            col = :genotype
        ) * visual(Scatter; color = :black, markersize = 7)

    draw(
        prior_spec + observed_spec;
        figure = (
            title = "predictions for site $site and replicate $replicate",
            size = (1200, 500)
        ),
        axis = (; limits = (nothing, (0, 1))),
        facet = (; linkxaxes = :all, linkyaxes = :all)
    )
end

ch_prior = sample(
	vm_model,
	Prior(),
	5000;
	progress=false
)

prior_predictions = returned(vm_model, ch_prior)

plot_predictions(6)

# save("figs/vm-prior-predictive.png", prior_predictive_figure.figure)

# ============================================================
# Posterior fit (two locus data only)
# ============================================================

ch_fit = sample(
	vm_model,
	NUTS(1000, 0.75),
	MCMCThreads(),
    500,
    8
)

summarystats(ch_fit)

using Serialization

# Save the FlexiChain object to a file
serialize("fits/vm-two-locus-fit-9-23.jls", ch_fit)

post_sim = returned(vm_model, ch_fit)
plot_predictions(post_sim, 3)

# ============================================================
# Summarizing parameter inferences
# ============================================================

# Diagnostics
using PairPlots

pairplot(
    ch_fit,
    [@varname(rAB), @varname(rBC), @varname(delta)];
    pool_chains=false,
    divergences=:numerical_error
)

# Recombination parameters
ch_fit_recom = ch_fit[[@varname(rAB), @varname(rBC)]]

# TODO 9/23 these looks horrible. hopefully they get better...
Makie.plot(ch_fit_recom)

df_fit_recom = stack(DataFrame(ch_fit_recom), [:rAB, :rBC])
df_prior_recom = stack(DataFrame(ch_prior[[@varname(rAB), @varname(rBC)]]), [:rAB, :rBC])

spec_post = data(df_fit_recom) *
    mapping(:value; color=:variable) *
    AlgebraOfGraphics.density(;datalimits=(0, 1))

spec_prior = data(df_prior_recom) *
    mapping(:value; strokecolor=:variable) *
    visual(Density; boundary=(0, 1), color=:transparent, strokewidth=2)
    # AlgebraOfGraphics.density(;datalimits=(0, 1)) *
    # subvisual(Band, color = :transparent)

draw(spec_post + spec_prior)

# Fitness parameters
ch_fit_fitness = ch_fit[[
    @varname(selection_phi), @varname(selection_lambda),
    @varname(dominance_phi), @varname(dominance_lambda),
]]

Makie.plot(ch_fit_fitness)

df_fit_fitness = FlexiChains.transform_values(
    ch_fit_fitness,
    [@varname(selection_phi), @varname(selection_lambda)] =>
        ((x, y) -> x ./ (x .+ y)) =>
        @varname(mean_selection),
    # TODO confirm that transforming the hierachical mean in this way is reasonable
    [@varname(dominance_phi), @varname(dominance_lambda)] =>
        ((x, y) -> 2(x ./ (x .+ y)) .- 1) =>
        @varname(mean_dominance),
) |> DataFrame

df_fit_fitness = stack(df_fit_fitness, r"mean")
select!(df_fit_fitness, :variable, :value)
@rtransform!(df_fit_fitness, :parameter=:variable[6:end-3])

df_fit_fitness_summary = @chain df_fit_fitness begin
    groupby([:variable, :parameter])
    @combine(
        :mean=mean(:value), 
        :lower=quantile(:value, 0.025), 
        :upper=quantile(:value, 0.975)
    )
end

spec_post = data(df_fit_fitness_summary) *
    (mapping(:mean, :variable; layout=:parameter) *
    visual(Scatter) +
    mapping(:variable, :lower, :upper; layout=:parameter) *
    visual(Rangebars; direction=:x))

draw(spec_post; facet=(;linkxaxes=:none))