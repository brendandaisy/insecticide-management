using Printf
using Plots

# ============================================================
# LOAD THE THREE-LOCUS MODEL
# ============================================================

# Load all constants, structures, and functions from the main
# three-locus model script.
#
# Using @__DIR__ ensures that Julia looks in the directory
# containing this experiment script, regardless of the directory
# from which Julia was launched.
include(
    joinpath(
        @__DIR__,
        "popgen-three-locus_parametric_type.jl"
    )
)

# ============================================================
# FREQUENCY-PLOTTING FUNCTION
# ============================================================

"""
Plot the frequencies of all eight haplotypes over generations.

Inputs
------
trajectory :
    A (generations + 1) × 8 matrix returned by simulate.

experiment_name :
    A descriptive name used in the plot title.

Returns
-------
frequency_plot :
    A Plots.jl plot object.
"""
function plot_frequency_trajectory(
    trajectory::AbstractMatrix;
    experiment_name::AbstractString = "Three-Locus Simulation"
)
    size(trajectory, 2) == N_HAPLOTYPES ||
        throw(DimensionMismatch(
            "The trajectory must have 8 columns. " *
            "Its size is $(size(trajectory))."
        ))

    # Row 1 contains generation 0.
    #
    # Therefore, a trajectory with n rows corresponds to
    # generations 0 through n - 1.
    generation_values = 0:(size(trajectory, 2) - 1)

    # Construct an initially empty plot.
    frequency_plot = plot(
        xlabel = "Generation",
        ylabel = "Haplotype frequency",
        title = experiment_name,
        ylim = (0.0, 1.0),
        legend = :outerright,
        linewidth = 2,
        grid = true,
        size = (950, 600)
    )

    # Add one trajectory for each of the eight haplotypes.
    for haplotype_index in 1:N_HAPLOTYPES
        plot!(
            frequency_plot,
            generation_values,
            trajectory[haplotype_index, :],
            label = HAPLOTYPE_NAMES[haplotype_index],
            linewidth = 2
        )
    end

    return frequency_plot
end

# ============================================================
# MEAN-FITNESS PLOTTING FUNCTION
# ============================================================

"""
Plot mean population fitness over generations.

The mean-fitness vector returned by simulate contains one value
for every transition from generation g to generation g + 1.
"""
function plot_mean_fitness(
    mean_fitness::AbstractVector;
    experiment_name::AbstractString = "Three-Locus Simulation"
)
    # mean_fitness[1] corresponds to the transition from
    # generation 0 to generation 1.
    generation_values = 0:(length(mean_fitness) - 1)

    mean_fitness_plot = plot(
        generation_values,
        mean_fitness,
        xlabel = "Generation",
        ylabel = "Mean fitness",
        title = "$experiment_name: Mean Fitness",
        label = "Mean fitness",
        linewidth = 2,
        grid = true,
        size = (850, 500)
    )

    return mean_fitness_plot
end

# ============================================================
# REUSABLE EXPERIMENT FUNCTION
# ============================================================

"""
Construct and run one three-locus simulation experiment.

Inputs
------
experiment_name :
    Name used for printed output and saved files.

s, h :
    Selection and dominance parameter vectors.

rAB, rBC :
    Recombination fractions for intervals A-B and B-C.

x0 :
    Initial haplotype-frequency distribution.

generations :
    Number of generations to simulate.

save_plots :
    If true, save the plots as PNG files.

output_directory :
    Directory in which plots are saved.

Returns
-------
A named tuple containing:

    W
        Constructed fitness matrix.

    model
        Constructed ThreeLocusModel.

    result
        Simulation output containing trajectory and mean fitness.

    frequency_plot
        Plot of all haplotype frequencies.

    mean_fitness_plot
        Plot of mean fitness.
"""
function run_experiment(;
    experiment_name::AbstractString,
    s::AbstractVector{<:Real},
    h::AbstractVector{<:Real},
    rAB::Real,
    rBC::Real,
    x0::AbstractVector{<:Real},
    generations::Integer,
    save_plots::Bool = true,
    output_directory::AbstractString = joinpath(
        @__DIR__,
        "simulation_output"
    )
)
    println()
    println("============================================================")
    println("RUNNING EXPERIMENT: $experiment_name")
    println("============================================================")

    # Construct the fitness matrix from s and h.
    #
    # checks=false is used here because the supplied main script
    # refers to check_constructed_fitness_matrix, but that validation
    # function is not currently present in the supplied file.
    #
    # The general fitness-matrix validation is performed immediately
    # afterward using check_fitness_matrix.
    W = build_fitness_matrix(
        s,
        h;
        checks = false
    )

    check_fitness_matrix(W)

    # Construct the complete selection-recombination model.
    model = build_model(
        W;
        rAB = rAB,
        rBC = rBC,
        checks = true
    )

    # Run the simulation.
    result = simulate(
        model,
        x0,
        generations;
        checks = true
    )

    # Construct the frequency plot.
    frequency_plot = plot_frequency_trajectory(
        result.trajectory;
        experiment_name = "$experiment_name: Haplotype Frequencies"
    )

    # Construct the mean-fitness plot.
    mean_fitness_plot = plot_mean_fitness(
        result.mean_fitness;
        experiment_name = experiment_name
    )

    # Print the parameter values used in the experiment.
    println()
    println("Recombination parameters:")
    @printf("rAB = %.6f\n", rAB)
    @printf("rBC = %.6f\n", rBC)

    println()
    println("Initial haplotype distribution:")
    print_distribution(x0)

    println()
    println("Final haplotype distribution after $generations generations:")
    print_distribution(result.trajectory[:, end])

    println()
    @printf(
        "Final mean fitness: %.10f\n",
        result.mean_fitness[end]
    )

    # Save the plots if requested.
    if save_plots
        mkpath(output_directory)

        # Replace spaces with underscores to create simple filenames.
        file_prefix = replace(
            lowercase(experiment_name),
            " " => "_"
        )

        frequency_filename = joinpath(
            output_directory,
            "$(file_prefix)_frequencies.png"
        )

        mean_fitness_filename = joinpath(
            output_directory,
            "$(file_prefix)_mean_fitness.png"
        )

        savefig(
            frequency_plot,
            frequency_filename
        )

        savefig(
            mean_fitness_plot,
            mean_fitness_filename
        )

        println()
        println("Frequency plot saved to:")
        println(frequency_filename)

        println()
        println("Mean-fitness plot saved to:")
        println(mean_fitness_filename)
    end

    return (
        W = W,
        model = model,
        result = result,
        frequency_plot = frequency_plot,
        mean_fitness_plot = mean_fitness_plot
    )
end

# ============================================================
# EXPERIMENT 1: PARAMETER VALUES
# ============================================================

# Selection parameters.
s = Float64[
    0.0,    # s[1]: VVF -- fixed reference value
    0.05,   # s[2]: VVC
    0.95,   # s[3]: VIF
    0.8,    # s[4]: VIC
    0.2,    # s[5]: LVF
    0.15,   # s[6]: LVC
    0.99,   # s[7]: LIF
    0.85    # s[8]: LIC
]

# Dominance parameters.
h = Float64[
    0.0,    # h[1]: VVF -- fixed reference value
    0.5,    # h[2]: VVC
    0.95,   # h[3]: VIF
    0.1,    # h[4]: VIC
    0.5,    # h[5]: LVF
    0.5,    # h[6]: LVC
    0.95,   # h[7]: LIF
    0.3     # h[8]: LIC
]

# Adjacent-locus recombination fractions.
#
# Change these values to investigate different levels of
# recombination between loci A-B and B-C.
rAB = 0.10
rBC = 0.10

# Initial haplotype-frequency distribution.
#
# These eight values must be nonnegative and must sum to 1.
#
# This example begins with all eight haplotypes at equal frequency.
x0 = Float64[
    0.125,   # VVF
    0.125,   # VVC
    0.125,   # VIF
    0.125,   # VIC
    0.125,   # LVF
    0.125,   # LVC
    0.125,   # LIF
    0.125    # LIC
]

# Number of generations to simulate.
generations = 100

# ============================================================
# RUN EXPERIMENT 1
# ============================================================

experiment_1 = run_experiment(
    experiment_name = "Baseline Experiment",
    s = s,
    h = h,
    rAB = rAB,
    rBC = rBC,
    x0 = x0,
    generations = generations,
    save_plots = true
)

# Display both plots in Julia or the VS Code plot pane.
display(experiment_1.frequency_plot)
display(experiment_1.mean_fitness_plot)

# ============================================================
# ACCESSING THE RESULTS
# ============================================================

# Complete frequency trajectory:
trajectory = experiment_1.result.trajectory

# Mean fitness at each transition:
mean_fitness = experiment_1.result.mean_fitness

# Final haplotype-frequency distribution:
final_distribution = trajectory[:, end]

# Constructed fitness matrix:
W = experiment_1.W

# ============================================================
# HAPLOTYPE-FREQUENCY PLOT
# ============================================================

"""
Plot the frequencies of all eight haplotypes over generations.

The trajectory matrix has:
    rows    = generations 0, 1, ..., generations
    columns = the eight haplotypes
"""
function plot_haplotype_frequencies(
    trajectory::AbstractMatrix;
    title_text::AbstractString = "Haplotype Frequencies Over Generations"
)
    size(trajectory, 2) == N_HAPLOTYPES ||
        throw(DimensionMismatch(
            "The trajectory must have 8 columns. " *
            "Its size is $(size(trajectory))."
        ))

    # The first row of trajectory is generation 0.
    generations = 0:(size(trajectory, 1) - 1)

    haplotype_plot = plot(
        xlabel = "Generation",
        ylabel = "Haplotype frequency",
        title = title_text,
        ylim = (0.0, 1.0),
        legend = :outerright,
        grid = true,
        size = (950, 600)
    )

    # Add one line for each haplotype.
    for i in 1:N_HAPLOTYPES
        plot!(
            haplotype_plot,
            generations,
            trajectory[i, :],
            label = HAPLOTYPE_NAMES[i],
            linewidth = 2
        )
    end

    return haplotype_plot
end

# ============================================================
# LOCUS-1 ALLELE-FREQUENCY PLOT
# ============================================================

"""
Plot the frequencies of alleles V and L at locus 1 over generations.

The plot also displays V + L. This total should remain equal
to 1 at every generation.
"""
function plot_locus1_allele_frequencies(
    trajectory::AbstractMatrix;
    title_text::AbstractString =
        "Allele Frequencies at Locus 1",
    checks::Bool = true,
    atol::Real = 1e-12
)
    allele_frequencies = locus1_allele_frequencies(
        trajectory;
        checks = checks,
        atol = atol
    )

    generations = 0:(size(trajectory, 1) - 1)

    allele_plot = plot(
        generations,
        allele_frequencies.V_frequency,
        label = "V",
        xlabel = "Generation",
        ylabel = "Allele frequency",
        title = title_text,
        linewidth = 3,
        ylim = (0.0, 1.05),
        grid = true,
        legend = :right,
        size = (850, 550)
    )

    plot!(
        allele_plot,
        generations,
        allele_frequencies.L_frequency,
        label = "L",
        linewidth = 3
    )

    # Plot the sum of the two allele frequencies.
    #
    # This line should remain at 1 for every generation.
    plot!(
        allele_plot,
        generations,
        allele_frequencies.total_frequency,
        label = "V + L",
        linewidth = 2,
        linestyle = :dash
    )

    # Report the largest numerical departure from 1.
    maximum_sum_error = maximum(
        abs.(
            allele_frequencies.total_frequency .-
            one(eltype(allele_frequencies.total_frequency))
        )
    )

    @printf(
        "Maximum error in V + L = 1: %.3e\n",
        maximum_sum_error
    )

    return allele_plot
end

locus1_plot = plot_locus1_allele_frequencies(
    experiment_1.result.trajectory;
    title_text = "Baseline Experiment: Locus 1"
)

display(locus1_plot)

# ============================================================
# JOINT MARGINAL-FREQUENCY PLOT FOR LOCI 2 AND 3
# ============================================================

"""
Plot the joint marginal frequencies VF, VC, IF, and IC over
generations.

The allele at locus 1 is ignored.

The plot also displays:

    VF + VC + IF + IC

which should remain equal to 1 at every generation.
"""
function plot_loci23_joint_frequencies(
    trajectory::AbstractMatrix;
    title_text::AbstractString =
        "Joint Marginal Frequencies of Loci 2 and 3",
    checks::Bool = true,
    atol::Real = 1e-12
)
    joint_frequencies = loci23_joint_frequencies(
        trajectory;
        checks = checks,
        atol = atol
    )

    # Row 1 corresponds to generation 0.
    generations = 0:(size(trajectory, 1) - 1)

    joint_plot = plot(
        generations,
        joint_frequencies.VF_frequency,
        label = "VF",
        xlabel = "Generation",
        ylabel = "Joint marginal frequency",
        title = title_text,
        linewidth = 3,
        ylim = (0.0, 1.05),
        grid = true,
        legend = :outerright,
        size = (950, 600)
    )

    plot!(
        joint_plot,
        generations,
        joint_frequencies.VC_frequency,
        label = "VC",
        linewidth = 3
    )

    plot!(
        joint_plot,
        generations,
        joint_frequencies.IF_frequency,
        label = "IF",
        linewidth = 3
    )

    plot!(
        joint_plot,
        generations,
        joint_frequencies.IC_frequency,
        label = "IC",
        linewidth = 3
    )

    # Plot the sum of all four joint marginal frequencies.
    #
    # This dashed line should remain at 1 for every generation.
    plot!(
        joint_plot,
        generations,
        joint_frequencies.total_frequency,
        label = "VF + VC + IF + IC",
        linewidth = 2,
        linestyle = :dash
    )

    # Calculate the largest numerical departure from 1.
    maximum_sum_error = maximum(
        abs.(
            joint_frequencies.total_frequency .-
            one(eltype(joint_frequencies.total_frequency))
        )
    )

    @printf(
        "Maximum error in VF + VC + IF + IC = 1: %.3e\n",
        maximum_sum_error
    )

    return joint_plot
end

loci23_plot = plot_loci23_joint_frequencies(
    experiment_1.result.trajectory;
    title_text =
        "Baseline Experiment: Joint Frequencies of Loci 2 and 3",
    checks = true
)

display(loci23_plot)