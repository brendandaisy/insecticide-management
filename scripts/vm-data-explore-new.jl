using Distributions
using InvertedIndices
using CSV, DataFrames, DataFramesMeta
using CairoMakie, AlgebraOfGraphics

include("../src/popgen-three-locus.jl")
include("../src/popgen-helpers.jl")
include("../src/popgen-dephasing.jl")

vm_two_loci = CSV.read(
    "data-proc/vm-geno-counts-1016-1534.csv",
    DataFrame
)

vm_three_loci = CSV.read(
    "data-proc/vm-geno-counts-F1-F8.csv", 
    DataFrame
)

# plot F1 and F8 three locus data on average across sites

# by genotype
vm3_by_genotype = @by vm_three_loci [:generation, :genotype] begin
    :prop=mean(:prop)
end

spec = data(vm3_by_genotype) *
    mapping(:generation, :prop; layout=:genotype, color=:genotype) *
    visual(ScatterLines)

draw(spec; figure=(;size=(600, 600)))

# by haplotype: per site/rep frequencies, then averaged across sites and reps
vm3_by_haplotype = @chain vm_three_loci begin
    @groupby [:site, :rep, :generation]
    @combine begin
        :haplotype = collect(HAPLOTYPE_NAMES)
        :freq = genotype_counts_to_haplotype_freq(:genotype, :count)
    end
    @by([:generation, :haplotype],
        :mean_freq = mean(:freq),
        :freq_lo = quantile(:freq, 0.025),
        :freq_hi = quantile(:freq, 0.975)
    )
end

spec = data(vm3_by_haplotype) *
    # (mapping(:generation, :freq_lo, :freq_hi; color=:haplotype) *
    # visual(Errorbars) +
    mapping(:generation, :mean_freq; color=:haplotype) *
    visual(ScatterLines)

draw(spec, scales(;Color=(;palette=:tab10)); figure=(;size=(500, 400)))

# by allele
vm3_by_allele = @chain vm_three_loci begin
    @groupby [:site, :rep, :generation]
    @combine @astable begin
        L, I, C = genotype_counts_to_allele_freq(:genotype, :count)
        :freq_L = L
        :freq_I = I
        :freq_C = C
    end
    stack(r"freq"; value_name=:frequency)
end

vm3_by_allele_summ = @by vm3_by_allele [:generation, :variable] begin
    :frequency=mean(:frequency)
end

spec = (data(vm3_by_allele) *
    visual(Scatter; alpha=0.3) +
    data(vm3_by_allele_summ) *
    visual(ScatterLines; markersize=12)) *
    mapping(
        :generation, 
        :frequency; 
        color=:variable, col=:variable => presorted
    )

draw(spec, scales(;Color=(;palette=:tab10)); figure=(;size=(600, 300)))

# ============================================================
# F1 DEPHASING: HWE TESTS AND EM WITH INBREEDING
# ============================================================

f1 = @rsubset vm_three_loci :generation == 1
sites = sort(unique(f1.site))

# Scheme A: everything pooled. Scheme B: each site separately
counts_by_scheme = Pair{String,Vector{Int}}[
    "pooled" => genotype_count_vector(f1.genotype, f1.count);
    [s => genotype_count_vector(f1.genotype[f1.site .== s], f1.count[f1.site .== s])
        for s in sites]
]

# --- per-locus HWE ---
loci_names = ["410 (L)", "1016 (I)", "1534 (C)"]

hwe_table = DataFrame(
    (; scheme, locus=loci_names[l], hwe_summary(locus_genotype_counts(counts, l))...)
    for (scheme, counts) in counts_by_scheme for l in 1:3
)

hwe_table

hwe_fisher = @chain hwe_table begin
    @rsubset :scheme != "pooled"
    @by :locus :p_fisher = fisher_combine(:p_exact)
end

# --- EM fits, F = 0 vs F free, with bootstrap LRT ---
lrt_results = [
    scheme => inbreeding_lrt(counts; B=200, rng=Xoshiro(i))
    for (i, (scheme, counts)) in enumerate(counts_by_scheme)
]

lrt_table = DataFrame(
    (;
        scheme,
        N=sum(counts),
        F_hat=r.fit1.F,
        F_lo=r.F_ci.lo,
        F_hi=r.F_ci.hi,
        ll_F0=r.fit0.loglik,
        ll_Ffree=r.fit1.loglik,
        lrt_stat=r.stat,
        p_asymptotic=r.p_asymptotic,
        p_bootstrap=r.p_bootstrap,
        converged=r.fit0.converged && r.fit1.converged,
    )
    for ((scheme, r), (_, counts)) in zip(lrt_results, counts_by_scheme)
)

# --- pooled vs per-site model comparison (7 free haplotype frequencies + F per fit) ---
function pooled_vs_site_comparison(lrt_results, counts_by_scheme)
    pooled = last(first(lrt_results))
    site_fits = last.(lrt_results[2:end])
    n_total = sum(sum(last(c)) for c in counts_by_scheme[2:end])
    n_sites = length(site_fits)

    map([(:F0, 7, :fit0), (:Ffree, 8, :fit1)]) do (model, k, field)
        ll_pooled = getfield(pooled, field).loglik
        ll_sites = sum(getfield(f, field).loglik for f in site_fits)
        df = k * (n_sites - 1)
        (;
            model,
            ll_pooled, ll_sites,
            lrt_stat=2(ll_sites - ll_pooled),
            df,
            p_asymptotic=ccdf(Chisq(df), 2(ll_sites - ll_pooled)),
            AIC_pooled=2k - 2ll_pooled,
            AIC_sites=2k * n_sites - 2ll_sites,
            BIC_pooled=k * log(n_total) - 2ll_pooled,
            BIC_sites=k * n_sites * log(n_total) - 2ll_sites,
        )
    end |> DataFrame
end

scheme_comparison = pooled_vs_site_comparison(lrt_results, counts_by_scheme)

# --- haplotype frequencies: equal-split baseline vs EM ---
dephased = DataFrame(
    (; scheme, haplotype=HAPLOTYPE_NAMES[h], method, freq)
    for ((scheme, counts), (_, r)) in zip(counts_by_scheme, lrt_results)
    for (method, x) in [
        "equal split" => let nz = counts .> 0
            genotype_counts_to_haplotype_freq(GENOTYPES_3L[nz], counts[nz])
        end,
        "EM, F = 0" => r.fit0.x,
        "EM, F free" => r.fit1.x,
    ]
    for (h, freq) in enumerate(x)
)

spec = data(dephased) *
    mapping(
        :haplotype => presorted, :freq;
        color=:method, dodge=:method, layout=:scheme => presorted
    ) *
    visual(BarPlot)

draw(spec; figure=(;size=(1000, 800)), facet=(;linkxaxes=:none))

spec = data(lrt_table) * (
    mapping(:scheme => presorted, :F_hat) * visual(Scatter) +
    mapping(:scheme => presorted, :F_lo, :F_hi) * visual(Rangebars)
)

draw(spec; figure=(;size=(600, 350)))

# --- checks ---
@assert all(r -> isapprox(sum(r.fit1.x), 1; atol=1e-8), last.(lrt_results))
@assert all(r -> r.fit1.loglik >= r.fit0.loglik - 1e-6, last.(lrt_results))

# EM allele frequencies should match direct allele frequencies from genotypes
# Holds exactly only for F = 0; with F > 0 the haplotype marginals shift.
let counts = last(first(counts_by_scheme)), x = last(first(lrt_results)).fit0.x
    nz = counts .> 0
    direct = genotype_counts_to_allele_freq(GENOTYPES_3L[nz], counts[nz])
    from_em = (
        L=sum(x[h] for h in 5:8),
        I=sum(x[h] for h in (3, 4, 7, 8)),
        C=sum(x[h] for h in (2, 4, 6, 8)),
    )
    @assert all(isapprox.(values(direct), values(from_em); atol=1e-6))
end

# Parameter recovery on simulated data (known x and F)
let rng = Xoshiro(7), x_true = rand(rng, Dirichlet(fill(0.5, N_HAPLOTYPES)))
    for F_true in (0.0, 0.3)
        sim = rand(rng, Multinomial(5000, genotype_probs(x_true, F_true)))
        fit = dephase_em(sim; starts=20, rng)
        @info "recovery" F_true F_hat=fit.F max_abs_err_x=maximum(abs.(fit.x .- x_true))
    end
end

