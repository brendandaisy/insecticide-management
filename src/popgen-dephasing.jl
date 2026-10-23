using Distributions
using Random
using SpecialFunctions: loggamma

# Requires popgen-three-locus.jl (HAPLOTYPE_NAMES, N_HAPLOTYPES) to be included first.

# Per-locus genotype classes ordered (susceptible homozygote, heterozygote, resistant homozygote).
const LOCUS_CLASSES = (("VV", "VL", "LL"), ("VV", "VI", "II"), ("FF", "FC", "CC"))

const GENOTYPES_3L = [
    join((a, b, c), '-')
    for a in LOCUS_CLASSES[1] for b in LOCUS_CLASSES[2] for c in LOCUS_CLASSES[3]
]

const HAPLOTYPE_INDEX = Dict(name => i for (i, name) in enumerate(HAPLOTYPE_NAMES))

"""
All unordered haplotype pairs `(i, j)`, `i <= j`, that are compatible with an
unphased three-locus genotype string such as "VL-VI-FC".
"""
function phase_pairs(genotype)
    loci = split(genotype, '-')
    het = [i for i in 1:3 if loci[i][1] != loci[i][2]]
    pairs = Set{Tuple{Int,Int}}()
    for mask in 0:(2^length(het) - 1)
        h1 = [loci[i][1] for i in 1:3]
        h2 = [loci[i][2] for i in 1:3]
        for (bit, locus) in enumerate(het)
            if (mask >> (bit - 1)) & 1 == 1
                h1[locus], h2[locus] = h2[locus], h1[locus]
            end
        end
        push!(pairs, minmax(HAPLOTYPE_INDEX[String(h1)], HAPLOTYPE_INDEX[String(h2)]))
    end
    return sort!(collect(pairs))
end

const PHASE_TABLE = [phase_pairs(g) for g in GENOTYPES_3L]

"Counts over `GENOTYPES_3L` from genotype strings and their counts."
function genotype_count_vector(genotypes, counts)
    v = zeros(Int, length(GENOTYPES_3L))
    for (g, c) in zip(genotypes, counts)
        i = findfirst(==(String(g)), GENOTYPES_3L)
        isnothing(i) && throw(ArgumentError("Unexpected three-locus genotype: $g."))
        v[i] += c
    end
    return v
end

# ============================================================
# PER-LOCUS HWE
# ============================================================

"Counts of (susceptible homozygote, heterozygote, resistant homozygote) at one locus."
function locus_genotype_counts(counts, locus)
    n = zeros(Int, 3)
    for (g, c) in zip(GENOTYPES_3L, counts)
        n[findfirst(==(split(g, '-')[locus]), LOCUS_CLASSES[locus])] += c
    end
    return n
end

"Exact HWE test for a biallelic locus (Wigginton et al. 2005, enumerating all heterozygote counts)."
function hwe_exact_pvalue(n_AA, n_AB, n_BB)
    N = n_AA + n_AB + n_BB
    n_A = 2n_AA + n_AB
    n_B = 2n_BB + n_AB
    logp(h) = loggamma(N + 1) - loggamma((n_A - h) ÷ 2 + 1) - loggamma(h + 1) -
        loggamma((n_B - h) ÷ 2 + 1) + h * log(2) +
        loggamma(n_A + 1) + loggamma(n_B + 1) - loggamma(2N + 1)
    lps = logp.((n_A % 2):2:min(n_A, n_B))
    return min(1.0, sum(exp(lp) for lp in lps if lp <= logp(n_AB) + 1e-9))
end

"Observed/expected heterozygosity, F-hat and exact p-value for one locus from `(n_AA, n_AB, n_BB)`."
function hwe_summary(n)
    N = sum(n)
    p_resistant = (2n[3] + n[2]) / (2N)
    H_obs = n[2] / N
    H_exp = 2p_resistant * (1 - p_resistant)
    return (;
        N, p_resistant, H_obs, H_exp,
        F_hat=H_exp > 0 ? 1 - H_obs / H_exp : NaN,
        p_exact=hwe_exact_pvalue(n...),
    )
end

"Fisher's method for combining independent p-values."
fisher_combine(p) = ccdf(Chisq(2length(p)), -2sum(log.(max.(p, floatmin()))))

# ============================================================
# HAPLOTYPE MODEL WITH INBREEDING F
# ============================================================

# F is the probability that the two haplotypes are identical by descent.
@inline function pair_prob(x, F, i, j)
    return i == j ? (1 - F) * x[i]^2 + F * x[i] : 2 * (1 - F) * x[i] * x[j]
end

"Probabilities of the 27 genotypes in `GENOTYPES_3L` under haplotype frequencies `x` and inbreeding `F`."
function genotype_probs(x, F)
    return [sum(pair_prob(x, F, i, j) for (i, j) in pairs) for pairs in PHASE_TABLE]
end

function loglik(x, F, counts)
    ll = 0.0
    for (pairs, c) in zip(PHASE_TABLE, counts)
        c == 0 && continue
        ll += c * log(sum(pair_prob(x, F, i, j) for (i, j) in pairs))
    end
    return ll
end

"""
EM for haplotype frequencies `x` (and inbreeding `F`) from unphased genotype counts.

With a latent identical-by-descent indicator, both updates are closed form:
a homozygous pair contributes 1 haplotype draw if IBD and 2 otherwise, and
`F` is the expected IBD fraction. Pass `F=value` to hold `F` fixed.
"""
function dephase_em(counts; F=nothing, starts=10, tol=1e-9, maxiter=5000, rng=Random.default_rng())
    free_F = isnothing(F)
    N = sum(counts)
    best = nothing

    for _ in 1:starts
        x = rand(rng, Dirichlet(ones(N_HAPLOTYPES)))
        F_cur = free_F ? 0.5rand(rng) : Float64(F)
        ll_old = -Inf
        ll = -Inf
        converged = false

        for _ in 1:maxiter
            draws = zeros(N_HAPLOTYPES)
            ibd = 0.0
            for (pairs, c) in zip(PHASE_TABLE, counts)
                c == 0 && continue
                probs = [pair_prob(x, F_cur, i, j) for (i, j) in pairs]
                total = sum(probs)
                for (k, (i, j)) in enumerate(pairs)
                    w = c * probs[k] / total
                    if i == j
                        denom = (1 - F_cur) * x[i] + F_cur
                        q = denom > 0 ? F_cur / denom : 0.0
                        draws[i] += w * (2 - q)
                        ibd += w * q
                    else
                        draws[i] += w
                        draws[j] += w
                    end
                end
            end
            x = draws ./ sum(draws)
            free_F && (F_cur = clamp(ibd / N, 0.0, 1 - 1e-8))
            ll = loglik(x, F_cur, counts)
            if abs(ll - ll_old) < tol
                converged = true
                break
            end
            ll_old = ll
        end

        if isnothing(best) || ll > best.loglik
            best = (; x, F=F_cur, loglik=ll, converged)
        end
    end
    return best
end

"Profile-likelihood interval for `F` on a grid (chi-square(1) cutoff)."
function profile_F_interval(counts, ll_max; grid=0:0.025:0.975, starts=3, cutoff=3.84, rng=Random.default_rng())
    lls = [dephase_em(counts; F=Float64(f), starts, rng).loglik for f in grid]
    ll_ref = max(ll_max, maximum(lls))
    keep = [f for (f, ll) in zip(grid, lls) if 2(ll_ref - ll) <= cutoff]
    return (lo=minimum(keep), hi=maximum(keep))
end

"""
Test F = 0 against F free by likelihood ratio. Reports the boundary-corrected
asymptotic p-value (0.5 * chi-square(1)) and a parametric bootstrap p-value.
"""
function inbreeding_lrt(counts; B=200, starts=10, boot_starts=5, rng=Random.default_rng())
    fit0 = dephase_em(counts; F=0.0, starts, rng)
    fit1 = dephase_em(counts; starts, rng)
    stat = max(2(fit1.loglik - fit0.loglik), 0.0)

    probs = genotype_probs(fit0.x, 0.0)
    N = sum(counts)
    boot = map(1:B) do _
        c = rand(rng, Multinomial(N, probs))
        f0 = dephase_em(c; F=0.0, starts=boot_starts, rng)
        f1 = dephase_em(c; starts=boot_starts, rng)
        max(2(f1.loglik - f0.loglik), 0.0)
    end

    return (;
        fit0, fit1, stat,
        p_asymptotic=stat > 0 ? 0.5 * ccdf(Chisq(1), stat) : 1.0,
        p_bootstrap=(1 + count(>=(stat - 1e-9), boot)) / (B + 1),
        F_ci=profile_F_interval(counts, fit1.loglik; rng),
    )
end
