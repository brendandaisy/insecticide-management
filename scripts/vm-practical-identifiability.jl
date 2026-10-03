using Distributions
using InvertedIndices
using CSV, DataFrames, DataFramesMeta
using CairoMakie, AlgebraOfGraphics
using Turing, FlexiChains
using Statistics
using Random

using Optim
using LikelihoodProfiler, OrdinaryDiffEqTsit5

include("../src/popgen-three-locus.jl")
include("../src/popgen-turing-vm.jl")
include("../src/popgen-data-proc.jl")

function multistart_mle(
    model,
    fixed_parameters=NamedTuple();
    seed::Integer,
    maxiters::Integer,
    n_starts::Integer=1
)
    n_starts >= 1 || throw(ArgumentError("n_starts must be positive."))
    conditioned_model = model | fixed_parameters
    best_candidate = nothing
    messages = String[]

    for start in 1:n_starts
        try
            # initial_params = Turing.InitFromParams(
            #     profile_initial_params(anchor, fixed_parameters)
            # )

            fit = maximum_likelihood(Xoshiro(seed+start-1), conditioned_model; maxiters)
            result = fit.optim_result
            success = string(result.retcode) == "Success"
            candidate = (
                fit = fit,
                objective = result.objective,
                success = success,
                message = string(result.retcode)
            )

            if isfinite(result.objective) &&
               (best_candidate === nothing || result.objective < best_candidate.objective)
                best_candidate = candidate
            end

            if !success
                push!(
                    messages,
                    "start $start: $(result.retcode), objective=$(result.objective)"
                )
            end
        catch exception
            push!(messages, sprint(showerror, exception))
        end
    end

    if best_candidate !== nothing
        return merge(
            best_candidate,
            (message=join(messages, "\n"),)
        )
    end

    # if all starts didn't succeed, return a failed fit result
    return (
        fit = nothing,
        objective = Inf,
        success = false,
        message = join(messages, "\n")
    )
end

function vm_expected_genotype_probs(parameters, counts, generations)
    selection = vcat(0.0, parameters.selection_raw)
    dominance = vcat(0.0, 2 .* parameters.dominance_raw .- 1)
    fitness = build_fitness_matrix(selection, dominance; checks=false)
    population_model = build_model(
        fitness;
        rAB=parameters.rAB,
        rBC=parameters.rBC,
        checks=false
    )

    initial_haplotype_counts = haplotype_counts_to_three_locus_counts(
        counts[:, 1]
    ) .+ parameters.delta
    initial_haplotypes = initial_haplotype_counts ./ sum(initial_haplotype_counts)
    trajectory = simulate(
        population_model,
        initial_haplotypes,
        generations;
        checks=false
    )

    expected = zeros(Float64, length(L2_L3_GENOTYPES), generations)
    for generation in 1:generations
        vm_two_locus_genotype_probs!(
            view(expected, :, generation),
            view(trajectory, :, generation)
        )
    end
    return expected
end

# Rebuild full 7-element raw vectors, putting fixed haplotypes back at 1.
function vm_full_parameters(params)
    raw(name, i) = i in fixed_raw_indices ? 1.0 :
        params[name == :selection_raw ? @varname(selection_raw[i]) : @varname(dominance_raw[i])]
    return (
        rAB=params[@varname(rAB)],
        rBC=params[@varname(rBC)],
        delta=params[@varname(delta)],
        selection_raw=[raw(:selection_raw, i) for i in 1:(N_HAPLOTYPES - 1)],
        dominance_raw=[raw(:dominance_raw, i) for i in 1:(N_HAPLOTYPES - 1)]
    )
end

# `adjustments` overrides MLE parameters (natural scale), e.g. (; rAB=0.5).
function vm_mle_plot_data(
    mle,
    counts,
    sample_sizes,
    generations;
    adjustments::NamedTuple=(;)
)
    parameters = vm_full_parameters(mle.params)
    series_parameters = ["MLE expected" => parameters]
    if !isempty(adjustments)
        label = "Adjusted: " *
            join(("$k=$(round.(v; digits=3))" for (k, v) in pairs(adjustments)), ", ")
        push!(series_parameters, label => merge(parameters, adjustments))
    end

    rows = NamedTuple[]
    for generation in 1:generations, genotype_index in eachindex(L2_L3_GENOTYPES)
        push!(
            rows,
            (
                generation=generation,
                genotype=L2_L3_GENOTYPES[genotype_index],
                proportion=counts[genotype_index, generation] /
                    sample_sizes[generation],
                series="Observed"
            )
        )
    end

    for (label, series_params) in series_parameters
        expected = vm_expected_genotype_probs(series_params, counts, generations)
        for generation in 1:generations, genotype_index in eachindex(L2_L3_GENOTYPES)
            push!(
                rows,
                (
                    generation=generation,
                    genotype=L2_L3_GENOTYPES[genotype_index],
                    proportion=expected[genotype_index, generation],
                    series=label
                )
            )
        end
    end

    return DataFrame(rows)
end

function plot_vm_mle_comparison(comparison_data)
    figure = Figure(size=(1200, 900))

    for (genotype_index, genotype) in enumerate(L2_L3_GENOTYPES)
        panel_data = filter(:genotype => ==(genotype), comparison_data)
        expected_data = filter(:series => !=("Observed"), panel_data)
        observed_data = filter(:series => ==("Observed"), panel_data)

        specification =
            data(expected_data) *
                mapping(:generation, :proportion; color=:series => presorted) *
                visual(Lines; linewidth=2) +
            data(observed_data) *
                mapping(:generation, :proportion) *
                visual(Scatter; color=:black, markersize=8)

        row = div(genotype_index - 1, 3) + 1
        column = mod1(genotype_index, 3)
        grid = draw!(
            figure[row, column],
            specification;
            axis=(;
                title=genotype,
                xlabel="generation",
                ylabel="proportion",
                limits=(nothing, (0, 1))
            )
        )
        genotype_index == 1 && legend!(figure[:, 4], grid)
    end

    Label(
        figure[0, :],
        "MLE dynamics and observed genotypes: site $vm_group_site, replicate $vm_group_rep";
        fontsize=18
    )
    return figure
end

vm_two_loci = CSV.read("data-proc/vm-geno-counts-1016-1534.csv", DataFrame)

vm_data = prepare_vm_two_loci(vm_two_loci)

vm_group = 2
vm_group_site = vm_data.group_keys[!, :site][vm_group]
vm_group_rep = vm_data.group_keys[!, :rep][vm_group]

group_counts = vm_data.counts[:, :, vm_group]
group_samples = vm_data.sample_sizes[:, vm_group]
group_generations = vm_data.n_generations[vm_group]

# Fixed haplotypes have selection_raw = dominance_raw = 1; raw index = haplotype index - 1.
fixed_haplotypes = ["VIF", "LIF"]
fixed_raw_indices = [
    findfirst(==(haplotype), HAPLOTYPE_NAMES) - 1 for haplotype in fixed_haplotypes
]
free_raw_indices = setdiff(1:(N_HAPLOTYPES - 1), fixed_raw_indices)
fixed_values = Dict(
    vcat(
        [@varname(selection_raw[i]) => 1.0 for i in fixed_raw_indices],
        [@varname(dominance_raw[i]) => 1.0 for i in fixed_raw_indices]
    )
)

group_model = vm_two_locus_model(
	group_counts,
    group_samples,
    group_generations
) | fixed_values

mle = multistart_mle(group_model; seed=115, maxiters=20000, n_starts=10);
mle.fit.params

vm_mle_plot_data(
    mle.fit,
    group_counts,
    group_samples,
    group_generations
) |> plot_vm_mle_comparison

profile_objective(theta, _) =
    -2 * Turing.LogDensityProblems.logdensity(mle.fit.ldf, theta)

profile_optfunction = OptimizationFunction(
    profile_objective,
    AutoForwardDiff()
)
profile_optproblem = OptimizationProblem(
    profile_optfunction,
    copy(mle.fit.optim_result.u)
)

varname_ranges = mle.fit.ldf._varname_ranges.data
# Free parameters are ordered rAB, rBC, delta, free selection_raw, free dominance_raw.
length(mle.fit.optim_result.u) == 3 + 2 * length(free_raw_indices) ||
    error("Unexpected number of free parameters.")
selection_indices = 3 .+ eachindex(free_raw_indices)
profile_indices = vcat(
    first(varname_ranges.rAB.range),
    first(varname_ranges.rBC.range),
    selection_indices
)
profile_names = vcat(
    :rAB,
    :rBC,
    [Symbol("selection ", HAPLOTYPE_NAMES[i + 1]) for i in free_raw_indices]
)
profile_lower = min.(-40.0, mle.fit.optim_result.u[profile_indices] .- 10)
profile_upper = max.(40.0, mle.fit.optim_result.u[profile_indices] .+ 10)
profile_target = ParameterTarget(
    profile_indices,
    profile_lower,
    profile_upper,
    profile_names
)
profile_problem = ProfileLikelihoodProblem(
    profile_optproblem,
    copy(mle.fit.optim_result.u),
    profile_target;
    conf_level=0.95
)

optimization_profiler = OptimizationProfiler(optimizer = LBFGS(), stepper = AdaptiveStep())
optimization_profile = solve(
    profile_problem,
    optimization_profiler;
    maxiters=3000,
    verbose=false
)

retcodes(optimization_profile)

integration_profiler = IntegrationProfiler(
    integrator=Tsit5(),
    integrator_opts=(;dtmax=0.1),
    matrix_type=:hessian
)
integration_profile = solve(
    profile_problem,
    integration_profiler;
    maxiters=3000,
    verbose=false
)

retcodes(integration_profile)

function profile_plot_data(solution, profile_index, parameter_name)
    parameter_index = profile_indices[profile_index]
    profile = DataFrame(solution[profile_index])
    linked_values = profile[!, Symbol("x", parameter_index)]
    mle_value = mle.fit.optim_result.u[parameter_index]
    profile[!, :linked_value] = linked_values
    profile[!, :delta_deviance] = profile.objective .-
        2 * mle.fit.optim_result.objective
    profile[!, :branch] = ifelse.(linked_values .<= mle_value, "lower", "upper")
    profile[!, :parameter] .= string(parameter_name)
    profile[!, :mle] .= mle_value
    return profile[!, [:parameter, :linked_value, :delta_deviance, :branch, :mle]]
end

function plot_profiles(profile)
    parameters = unique(profile.parameter)
    threshold_data = DataFrame(
        parameter=parameters,
        threshold=fill(profile_problem.threshold, length(parameters))
    )
    mle_data = unique(profile[!, [:parameter, :mle]])

    specification = data(profile) *
        mapping(
            :linked_value => "transformed parameter value",
            :delta_deviance => "profile likelihood-ratio statistic";
            color=:branch,
            layout=:parameter => presorted
        ) * visual(ScatterLines; markersize=5) +
        data(threshold_data) *
        mapping(:threshold; layout=:parameter => presorted) *
        visual(HLines; color=:black, linestyle=:dash) +
        data(mle_data) *
        mapping(:mle; layout=:parameter => presorted) *
        visual(VLines; color=:gray, linestyle=:dot)

    draw(
        specification;
        figure=(; size=(1400, 900), title="site $vm_group_site, replicate $vm_group_rep"),
        facet=(; linkxaxes=:none, linkyaxes=:none),
    )
end

profile_data = reduce(
    vcat,
    [
        profile_plot_data(integration_profile, i, name)
        for (i, name) in enumerate(profile_names)
    ]
)

plot_profiles(profile_data)

#=
9/29

from profiles of recombination paramters, the MLE of both is 0 (probably not really but 
these clearly produce a reasonable fit).

with this in mind, when r_AB≈0, we find increasing r_BC->0.5 does not change the likelihood ratio,
clearly a degenerate surface. I'm somewhat suprised that with the crossover definitions, changing
r_BC should change the model dynamics, and theoretically the unambigous genotypes should give info
about r_BC. 

TODO simply making a heatmap of the loglikelihood for both params with the others at an MLE may give you
the visual you want
=#

mle = multistart_mle(group_model; seed=115, maxiters=20000, n_starts=10);
mle.fit.params

vm_mle_plot_data(
    mle.fit,
    group_counts,
    group_samples,
    group_generations
    # adjustments=(;rAB=0.99)
) |> plot_vm_mle_comparison

mle_other_params = [
    k => v
    for (k, v) in pairs(mle.fit.params) 
    if k ∉ (@varname(rAB), @varname(rBC))
]

recomb_grid = allcombinations(
    DataFrame,
    rAB=collect(range(1e-4, 1-1e-4, 10)),
    rBC=collect(range(1e-4, 1-1e-4, 10))
)

recomb_loglik = @rtransform(
    recomb_grid,
    :loglik=loglikelihood(
        group_model, 
        VarNamedTuple(vcat(
            mle_other_params,
            @varname(rAB) => :rAB, @varname(rBC) => :rBC
        ))
    )
)

spec = data(recomb_loglik) *
    mapping(:rAB, :rBC, :loglik) *
    visual(Heatmap)

draw(spec)