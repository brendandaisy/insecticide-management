using Printf

# ============================================================
# FIXED THREE-LOCUS HAPLOTYPE INFORMATION
# ============================================================

# There are 2^3 = 8 possible three-locus haplotypes.
const N_HAPLOTYPES = 8

# Julia indices:
#   1     2     3     4     5     6     7     8
#
# Mathematical indices:
#   0     1     2     3     4     5     6     7
const HAPLOTYPE_NAMES = (
    "VVF",
    "VVC",
    "VIF",
    "VIC",
    "LVF",
    "LVC",
    "LIF",
    "LIC"
)

# Binary representation of each haplotype.
#
# The three entries correspond to loci A, B, and C.
# Lowercase allele = 0
# Uppercase allele = 1
const HAPLOTYPE_BITS = (
    (0, 0, 0),   # VVF
    (0, 0, 1),   # VVC
    (0, 1, 0),   # VIF
    (0, 1, 1),   # VIC
    (1, 0, 0),   # LVF
    (1, 0, 1),   # LVC
    (1, 1, 0),   # LIF
    (1, 1, 1)    # LIC
)

# ============================================================
# FITNESS MATRIX PARAMETERIZATIONS
# ============================================================

"""
Construct the symmetric 8 × 8 fitness matrix W from s and h.

The numerical element type is inferred from s and h. Therefore,
this function supports Float64, Float32, BigFloat, and
ForwardDiff.Dual numbers used by Turing.
"""
function build_fitness_matrix(
    s::AbstractVector{<:Real},
    h::AbstractVector{<:Real};
    checks::Bool = true,
    atol::Real = 1e-12
)
    length(s) == N_HAPLOTYPES ||
        throw(DimensionMismatch(
            "s must have length $N_HAPLOTYPES. Got $(length(s))."
        ))

    length(h) == N_HAPLOTYPES ||
        throw(DimensionMismatch(
            "h must have length $N_HAPLOTYPES. Got $(length(h))."
        ))

    # Find a common concrete numerical type.
    #
    # Examples:
    #   Float64 and Float64        -> Float64
    #   Float64 and BigFloat       -> BigFloat
    #   Float64 and Dual           -> Dual
    T = promote_type(eltype(s), eltype(h))

    # Convert both vectors to the common type.
    s_typed = T.(s)
    h_typed = T.(h)

    # Allocate W using T rather than Float64.
    W = zeros(T, N_HAPLOTYPES, N_HAPLOTYPES)

    zero_T = zero(T)
    one_T  = one(T)

    # Diagonal entries.
    for i in 1:N_HAPLOTYPES
        W[i, i] = one_T - s_typed[i]
    end

    # Off-diagonal entries.
    for i in 1:(N_HAPLOTYPES - 1)
        for j in (i + 1):N_HAPLOTYPES
            fitness = max(
                zero_T,
                one_T -
                h_typed[i] * s_typed[i] -
                h_typed[j] * s_typed[j]
            )

            W[i, j] = fitness
            W[j, i] = fitness
        end
    end

    if checks
        check_constructed_fitness_matrix(
            W,
            s_typed,
            h_typed;
            atol = atol
        )
    end

    return W
end

# ============================================================
# MODEL CONTAINER
# ============================================================

"""
Container for a fixed three-locus selection-recombination model.

Fields
------
W : 8 × 8 relative-fitness matrix

P : 8 × 8 × 8 recombination tensor, where
        P[k, i, j]
    is the probability that the diploid i/j produces gamete k.

K : 8 × 8 × 8 precomputed combined tensor,
        K[k, i, j] = W[i, j] * P[k, i, j]

The tensors P and K are built only once and reused at every generation.
"""
struct ThreeLocusModel{T<:Real}
    W::Matrix{T}
    P::Array{T, 3}
    K::Array{T, 3}
end

# ============================================================
# VALIDATION FUNCTIONS
# ============================================================

"""
Check four recombination-pattern probabilities.

Requirements:
    q0 + q1 + q2 + q12 = 1
    each q value is finite and nonnegative
"""
function check_recombination_parameters(
    q0::Real,
    q1::Real,
    q2::Real,
    q12::Real;
    atol::Real = 1e-12
)
    # promote returns values with one common concrete type
    q0_t, q1_t, q2_t, q12_t = promote(q0, q1, q2, q12)

    q = (q0_t, q1_t, q2_t, q12_t)

    all(isfinite, q) ||
        throw(ArgumentError(
            "All recombination probabilities must be finite."
        ))

    all(value -> value >= zero(value), q) ||
        throw(ArgumentError(
            "All recombination probabilities must be nonnegative. Got $q."
        ))

    q_sum = sum(q)

    isapprox(
        q_sum,
        one(q_sum);
        atol = atol,
        rtol = 0
    ) ||
        throw(ArgumentError(
            "Recombination probabilities must sum to 1. " *
            "Their sum is $q_sum."
        ))

    return true
end


"""
Check an 8 × 8 relative-fitness matrix.
"""
function check_fitness_matrix(
    W::AbstractMatrix{<:Real};
    atol::Real = 1e-12
)
    size(W) == (N_HAPLOTYPES, N_HAPLOTYPES) ||
        throw(DimensionMismatch(
            "W must be an 8 × 8 matrix. Its size is $(size(W))."
        ))

    all(isfinite, W) ||
        throw(ArgumentError(
            "Every entry of W must be finite."
        ))

    minimum(W) >= 0.0 ||
        throw(ArgumentError(
            "Relative fitness values cannot be negative."
        ))

    maximum(W) > 0.0 ||
        throw(ArgumentError(
            "At least one genotype must have positive fitness."
        ))

    isapprox(W, transpose(W); atol = atol, rtol = 0.0) ||
        throw(ArgumentError(
            "W must be symmetric because W[i,j] = W[j,i]."
        ))

    return true
end


"""
Check a haplotype-frequency distribution.

Requirements:
    length(x) = 8
    all entries are finite
    all entries are nonnegative up to numerical tolerance
    sum(x) = 1
"""
function check_distribution(
    x::AbstractVector{<:Real};
    atol::Real = 1e-12
)
    length(x) == N_HAPLOTYPES ||
        throw(DimensionMismatch(
            "The haplotype distribution must have length 8. " *
            "Its length is $(length(x))."
        ))

    all(isfinite, x) ||
        throw(ArgumentError(
            "Every haplotype frequency must be finite."
        ))

    minimum(x) >= -atol ||
        throw(ArgumentError(
            "Haplotype frequencies cannot be negative. " *
            "The smallest value is $(minimum(x))."
        ))

    isapprox(sum(x), 1.0; atol = atol, rtol = 0.0) ||
        throw(ArgumentError(
            "Haplotype frequencies must sum to 1. " *
            "Their sum is $(sum(x))."
        ))

    return true
end

# ============================================================
# RECOMBINATION PARAMETER CONSTRUCTION
# ============================================================

"""
Convert adjacent-locus recombination fractions into the four
three-locus recombination-pattern probabilities.

This function assumes independent recombination in intervals A-B and B-C.
"""
function q_from_independent_recombination(
    rAB::Real,
    rBC::Real;
    checks::Bool = true,
    atol::Real = 1e-12
)
    rAB_t, rBC_t = promote(rAB, rBC)

    zero_T = zero(rAB_t)
    one_T  = one(rAB_t)

    all(isfinite, (rAB_t, rBC_t)) ||
        throw(ArgumentError(
            "Recombination fractions must be finite."
        ))

    zero_T <= rAB_t <= one_T ||
        throw(ArgumentError(
            "rAB must be between 0 and 1. Got rAB = $rAB_t."
        ))

    zero_T <= rBC_t <= one_T ||
        throw(ArgumentError(
            "rBC must be between 0 and 1. Got rBC = $rBC_t."
        ))

    q0  = (one_T - rAB_t) * (one_T - rBC_t)
    q1  = rAB_t * (one_T - rBC_t)
    q2  = (one_T - rAB_t) * rBC_t
    q12 = rAB_t * rBC_t

    if checks
        check_recombination_parameters(
            q0,
            q1,
            q2,
            q12;
            atol = atol
        )
    end

    return (
        q0  = q0,
        q1  = q1,
        q2  = q2,
        q12 = q12
    )
end

# ============================================================
# RECOMBINATION LOOKUP TENSOR
# ============================================================

"""
Convert a binary haplotype such as (1,0,1) into its Julia index.

Examples:
    (0,0,0) -> 1, corresponding to VVF
    (1,0,1) -> 6, corresponding to LVC
    (1,1,1) -> 8, corresponding to LIC
"""
function haplotype_to_index(h)
    return 1 + 4 * h[1] + 2 * h[2] + h[3]
end


"""
Check the complete recombination tensor P.
"""
function check_recombination_kernel(
    P::AbstractArray{<:Real, 3};
    atol::Real = 1e-12
)
    size(P) == (
        N_HAPLOTYPES,
        N_HAPLOTYPES,
        N_HAPLOTYPES
    ) ||
        throw(DimensionMismatch(
            "P must have size 8 × 8 × 8. Its size is $(size(P))."
        ))

    all(isfinite, P) ||
        throw(ArgumentError(
            "Every entry of P must be finite."
        ))

    minimum(P) >= -atol ||
        throw(ArgumentError(
            "P contains a negative probability."
        ))

    # For each parental diplotype i/j, gamete probabilities must sum to 1.
    for i in 1:N_HAPLOTYPES
        for j in 1:N_HAPLOTYPES
            probability_sum = sum(@view P[:, i, j])

            isapprox(
                probability_sum,
                one(probability_sum);
                atol = atol,
                rtol = 0.0
            ) ||
                throw(ArgumentError(
                    "Gamete probabilities for parental pair " *
                    "($i,$j) sum to $probability_sum rather than 1."
                ))
        end
    end

    # Parental order must not matter:
    # P[k,i,j] = P[k,j,i].
    for k in 1:N_HAPLOTYPES
        for i in 1:N_HAPLOTYPES
            for j in 1:N_HAPLOTYPES
                isapprox(
                    P[k, i, j],
                    P[k, j, i];
                    atol = atol,
                    rtol = 0.0
                ) ||
                    throw(ArgumentError(
                        "P is not symmetric in parental indices " *
                        "for k=$k, i=$i, j=$j."
                    ))
            end
        end
    end

    return true
end


"""
Construct P[k,i,j] = probability that diploid i/j produces gamete k.

The tensor is constructed once and used as a lookup table during simulation.
"""
function build_recombination_kernel(
    q0::Real,
    q1::Real,
    q2::Real,
    q12::Real;
    checks::Bool = true,
    atol::Real = 1e-12
)
    # Convert all recombination probabilities to one common
    # concrete numerical type.
    q0_t, q1_t, q2_t, q12_t = promote(
        q0,
        q1,
        q2,
        q12
    )

    # Ensure that the numerical type supports division.
    #
    # Examples:
    #   Float64            -> Float64
    #   Float32            -> Float32
    #   BigFloat           -> BigFloat
    #   ForwardDiff.Dual   -> ForwardDiff.Dual
    T = typeof(q0_t / 2)

    q0_t  = convert(T, q0_t)
    q1_t  = convert(T, q1_t)
    q2_t  = convert(T, q2_t)
    q12_t = convert(T, q12_t)

    if checks
        check_recombination_parameters(
            q0_t,
            q1_t,
            q2_t,
            q12_t;
            atol = atol
        )
    end

    # First index: resulting gamete k
    # Second index: parental haplotype i
    # Third index: parental haplotype j
    P = zeros(
        T,
        N_HAPLOTYPES,
        N_HAPLOTYPES,
        N_HAPLOTYPES
    )

    for i in 1:N_HAPLOTYPES
        for j in 1:N_HAPLOTYPES

            # Binary parental haplotypes.
            u = HAPLOTYPE_BITS[i]
            v = HAPLOTYPE_BITS[j]

            # ------------------------------------------------
            # 1. No recombination
            #
            # Gametes:
            #   (u1,u2,u3)
            #   (v1,v2,v3)
            # ------------------------------------------------

            P[haplotype_to_index(u), i, j] += q0_t / 2
            P[haplotype_to_index(v), i, j] += q0_t / 2

            # ------------------------------------------------
            # 2. Recombination in interval A-B only
            #
            # Gametes:
            #   (u1,v2,v3)
            #   (v1,u2,u3)
            # ------------------------------------------------

            gamete_1 = (u[1], v[2], v[3])
            gamete_2 = (v[1], u[2], u[3])

            P[haplotype_to_index(gamete_1), i, j] += q1_t / 2
            P[haplotype_to_index(gamete_2), i, j] += q1_t / 2

            # ------------------------------------------------
            # 3. Recombination in interval B-C only
            #
            # Gametes:
            #   (u1,u2,v3)
            #   (v1,v2,u3)
            # ------------------------------------------------

            gamete_1 = (u[1], u[2], v[3])
            gamete_2 = (v[1], v[2], u[3])

            P[haplotype_to_index(gamete_1), i, j] += q2_t / 2
            P[haplotype_to_index(gamete_2), i, j] += q2_t / 2

            # ------------------------------------------------
            # 4. Recombination in both intervals
            #
            # Gametes:
            #   (u1,v2,u3)
            #   (v1,u2,v3)
            # ------------------------------------------------

            gamete_1 = (u[1], v[2], u[3])
            gamete_2 = (v[1], u[2], v[3])

            P[haplotype_to_index(gamete_1), i, j] += q12_t / 2
            P[haplotype_to_index(gamete_2), i, j] += q12_t / 2
        end
    end

    if checks
        check_recombination_kernel(P; atol = atol)
    end

    return P
end

# ============================================================
# MODEL CONSTRUCTION
# ============================================================

"""
Check the combined tensor K[k,i,j] = W[i,j] * P[k,i,j].

For each i,j:
    sum over k of K[k,i,j] must equal W[i,j]
because sum over k of P[k,i,j] equals 1.
"""
function check_combined_kernel(
    K::AbstractArray{<:Real, 3},
    W::AbstractMatrix{<:Real};
    atol::Real = 1e-12
)
    size(K) == (
        N_HAPLOTYPES,
        N_HAPLOTYPES,
        N_HAPLOTYPES
    ) ||
        throw(DimensionMismatch(
            "K must have size 8 × 8 × 8."
        ))

    all(isfinite, K) ||
        throw(ArgumentError(
            "Every entry of K must be finite."
        ))

    minimum(K) >= -atol ||
        throw(ArgumentError(
            "K contains a negative value."
        ))

    for i in 1:N_HAPLOTYPES
        for j in 1:N_HAPLOTYPES
            kernel_sum = sum(@view K[:, i, j])

            isapprox(
                kernel_sum,
                W[i, j];
                atol = atol,
                rtol = 0.0
            ) ||
                throw(ArgumentError(
                    "For parental pair ($i,$j), " *
                    "sum(K[:,i,j]) = $kernel_sum, " *
                    "but W[i,j] = $(W[i,j])."
                ))
        end
    end

    return true
end


"""
Construct a model directly from q0, q1, q2, and q12.
"""
function build_model_from_q(
    W::AbstractMatrix{<:Real},
    q0::Real,
    q1::Real,
    q2::Real,
    q12::Real;
    checks::Bool = true,
    atol::Real = 1e-12
)
    # Find one common concrete numerical type for the fitness
    # matrix and all four recombination probabilities.
    T0 = promote_type(
        eltype(W),
        typeof(q0),
        typeof(q1),
        typeof(q2),
        typeof(q12)
    )

    # Ensure that the numerical type supports division.
    T = typeof(one(T0) / 2)

    # Store a concrete copy of the fitness matrix using the
    # promoted numerical type.
    W_typed = Matrix{T}(W)

    q0_t  = convert(T, q0)
    q1_t  = convert(T, q1)
    q2_t  = convert(T, q2)
    q12_t = convert(T, q12)

    if checks
        check_fitness_matrix(W_typed; atol = atol)
    end

    # Build the recombination lookup tensor once.
    P = build_recombination_kernel(
        q0_t,
        q1_t,
        q2_t,
        q12_t;
        checks = checks,
        atol = atol
    )

    # Reshape W to 1 × 8 × 8 so it broadcasts over the k index.
    #
    # K[k,i,j] = P[k,i,j] * W[i,j].
    K = P .* reshape(
        W_typed,
        1,
        N_HAPLOTYPES,
        N_HAPLOTYPES
    )

    if checks
        check_combined_kernel(
            K,
            W_typed;
            atol = atol
        )
    end

    return ThreeLocusModel(
        W_typed,
        P,
        K
    )
end


"""
Convenience constructor using independent recombination fractions rAB and rBC.
"""
function build_model(
    W::AbstractMatrix{<:Real};
    rAB::Real,
    rBC::Real,
    checks::Bool = true,
    atol::Real = 1e-12
)
    q = q_from_independent_recombination(
        rAB,
        rBC;
        checks = checks,
        atol = atol
    )

    return build_model_from_q(
        W,
        q.q0,
        q.q1,
        q.q2,
        q.q12;
        checks = checks,
        atol = atol
    )
end

# ============================================================
# ONE-GENERATION UPDATE
# ============================================================

"""
Advance the model by one generation.

Inputs
------
x_next : preallocated vector that will receive the next distribution

x      : current haplotype-frequency distribution

model  : fixed ThreeLocusModel

Returns
-------
mean_fitness : normalization constant for the generation

The function modifies x_next in place.
"""
function step!(
    x_next::AbstractVector{<:Real},
    x::AbstractVector{<:Real},
    model::ThreeLocusModel;
    checks::Bool = false,
    atol::Real = 1e-12
)
    length(x_next) == N_HAPLOTYPES ||
        throw(DimensionMismatch(
            "x_next must have length 8."
        ))

    if checks
        check_distribution(x; atol = atol)
    end

    # Erase values left over from the previous generation.
    fill!(x_next, zero(eltype(x_next)))

    # For each ordered parental haplotype pair i/j:
    for i in 1:N_HAPLOTYPES
        for j in 1:N_HAPLOTYPES

            # Random-mating frequency of ordered zygote i/j.
            zygote_frequency = x[i] * x[j]

            # Skip calculations when the parental pair is absent.
            if iszero(zygote_frequency)
                continue
            end

            # Distribute the selected zygote contribution among
            # the eight possible gametes.
            for k in 1:N_HAPLOTYPES
                x_next[k] += (
                    model.K[k, i, j] *
                    zygote_frequency
                )
            end
        end
    end

    # Before normalization, sum(x_next) equals mean fitness.
    mean_fitness = sum(x_next)

    # This condition must always be checked because division by
    # zero would make the next generation undefined.
    if !isfinite(mean_fitness) ||
       mean_fitness <= zero(mean_fitness)

        throw(DomainError(
            mean_fitness,
            "Mean fitness must be finite and strictly positive."
        ))
    end

    # Normalize the new gamete distribution.
    x_next ./= mean_fitness

    if checks
        check_distribution(x_next; atol = atol)
    end

    return mean_fitness
end

# ============================================================
# MULTI-GENERATION SIMULATION
# ============================================================

"""
Simulate a complete trajectory.

Returns a named tuple with:

    trajectory
        A (generations + 1) × 8 matrix.
        Row 1 is generation 0.
        Row g+1 is generation g.

    mean_fitness
        A vector containing mean fitness for each transition.
"""
function simulate(
    model::ThreeLocusModel,
    x0,
    generations::Integer;
    checks::Bool = true,
    atol::Real = 1e-12
)
    generations >= 0 ||
        throw(ArgumentError(
            "The number of generations cannot be negative."
        ))

    # Collect the initial input before determining the common
    # numerical type.
    x_raw = collect(x0)

    # Find a common concrete numerical type for the initial
    # distribution and the model tensors.
    #
    # If the model contains ForwardDiff.Dual values, the
    # trajectory will also contain ForwardDiff.Dual values.
    T0 = promote_type(
        eltype(x_raw),
        eltype(model.K)
    )

    # Ensure that the numerical type supports division.
    T = typeof(one(T0) / 2)

    # Convert the initial input into a vector with the common type.
    x_current = T.(x_raw)

    if checks
        check_distribution(x_current; atol = atol)
    end

    # Allocate storage for the complete frequency trajectory.
    trajectory = zeros(
        T,
        generations + 1,
        N_HAPLOTYPES
    )

    # Store generation 0.
    trajectory[1, :] .= x_current

    # Mean fitness is defined for each transition t -> t+1.
    mean_fitness = zeros(T, generations)

    # Reusable next-generation vector.
    x_next = similar(x_current)

    for generation in 1:generations

        mean_fitness[generation] = step!(
            x_next,
            x_current,
            model;
            checks = checks,
            atol = atol
        )

        # Store the newly calculated distribution.
        trajectory[generation + 1, :] .= x_next

        # Exchange the roles of the two vectors without allocating
        # new vectors.
        x_current, x_next = x_next, x_current
    end

    return (
        trajectory = trajectory,
        mean_fitness = mean_fitness
    )
end


"""
Simulate without saving the complete trajectory.

Use this function when only the final distribution is needed.
"""
function simulate_final(
    model::ThreeLocusModel,
    x0,
    generations::Integer;
    checks::Bool = true,
    atol::Real = 1e-12
)
    generations >= 0 ||
        throw(ArgumentError(
            "The number of generations cannot be negative."
        ))

    # Collect the initial input before determining the common
    # numerical type.
    x_raw = collect(x0)

    # Find a common concrete numerical type for the initial
    # distribution and the model tensors.
    T0 = promote_type(
        eltype(x_raw),
        eltype(model.K)
    )

    # Ensure that the numerical type supports division.
    T = typeof(one(T0) / 2)

    x_current = T.(x_raw)

    if checks
        check_distribution(x_current; atol = atol)
    end

    x_next = similar(x_current)

    for _ in 1:generations
        step!(
            x_next,
            x_current,
            model;
            checks = checks,
            atol = atol
        )

        x_current, x_next = x_next, x_current
    end

    return x_current
end

# ============================================================
# OUTPUT UTILITY
# ============================================================

"""
Print one haplotype-frequency distribution with labels.
"""
function print_distribution(x)
    check_distribution(x)

    println("Math index   Haplotype   Frequency")
    println("----------------------------------")

    for julia_index in 1:N_HAPLOTYPES
        mathematical_index = julia_index - 1

        @printf(
            "%5d        %-3s       %.10f\n",
            mathematical_index,
            HAPLOTYPE_NAMES[julia_index],
            x[julia_index]
        )
    end

    return nothing
end

# ============================================================
# LOCUS-1 ALLELE FREQUENCIES
# ============================================================
function get_410_frequencies(x)
    # Sum the frequencies of all haplotypes carrying V
    # at the first locus.
    V_frequency = sum(x[1:4])

    # Sum the frequencies of all haplotypes carrying L
    # at the first locus.
    L_frequency = sum(x[5:8])

    [V_frequency, L_frequency]
end

# ============================================================
# JOINT MARGINAL FREQUENCIES FOR LOCI 2 AND 3
# ============================================================

"""
Calculate the joint marginal frequencies of loci 2 and 3.

The allele at locus 1 is ignored. Therefore:

    VF includes:
        VVF and LVF

    VC includes:
        VVC and LVC

    IF includes:
        VIF and LIF

    IC includes:
        VIC and LIC

At every generation, the four joint marginal frequencies must
add to 1:

    VF + VC + IF + IC = 1
"""
function get_1016_1534_frequencies(x)
    # VF appears in haplotypes VVF and LVF.
    VF_frequency = x[1] .+ x[5]

    # VC appears in haplotypes VVC and LVC.
    VC_frequency = x[2] .+ x[6]

    # IF appears in haplotypes VIF and LIF.
    IF_frequency = x[3] .+ x[7]

    # IC appears in haplotypes VIC and LIC.
    IC_frequency = x[4] .+ x[8]

    return [VF_frequency, VC_frequency, IF_frequency, IC_frequency]
end

