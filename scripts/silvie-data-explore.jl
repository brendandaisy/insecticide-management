using Distributions
using InvertedIndices
using CSV, DataFrames, DataFramesMeta
using CairoMakie, AlgebraOfGraphics

set_aog_theme!()

rr0 = CSV.read("data-raw/Resistance.Reversal_ALL.kdr_2025-07-24.csv", DataFrame)

rr = @select rr0 begin
    :rep=:Replicate
    :generation=:Generation
    :box=:Box
    :sex=:Sex
    :gen_410=replace.(:gen_410, "S" => "V", "R" => "L")
    :gen_1016=replace.(:gen_1016, "S" => "V", "R" => "I")
    :gen_1534
end

# won't be considering 1534 since problem with assay
@transform!(rr, @byrow :genotype="$(:gen_410)-$(:gen_1016)")

rr_geno_counts = @by(rr, [:rep, :generation, :genotype], :count=length(:genotype))

CSV.write("data-proc/rr-two-locus-geno-counts.csv", rr_geno_counts)

spec = data(rr_geno_counts) *
    mapping(:generation, :count, color=:genotype, col=:rep) *
    visual(ScatterLines)

draw(spec, scales(Color=(;palette=:tab10)))

# view trends of genotypes at the two individual loci
rr_loc_counts = @chain rr begin
    stack([:gen_410, :gen_1016, :gen_1534])
    @transform(:locus=replace.(:variable, "gen_" => ""), :genotype=:value)
    groupby([:rep, :generation, :locus, :genotype])
    @combine(:count=length(:genotype))
end

spec = data(@subset(rr_loc_counts, :locus .!= "1534")) *
    mapping(:generation, :count, color=:genotype, row=:rep, col=:locus) *
    visual(ScatterLines)

draw(spec)

# fit a fitness costs model with two loci and recombination
using Turing
using FlexiChains
using CategoricalArrays
import IterTools: product

include("../src/popgen-two-locus.jl")
include("../src/popgen-turing-two-locus.jl")

unique(rr_geno_counts.genotype)

mono_pairs = product(["VV", "VL", "LL"], ["VV", "VI", "II"])
genotypes = vec(collect("$h1-$h2" for (h1, h2) in mono_pairs))

rr_geno_counts = @chain rr_geno_counts begin
    @transform!(:genotype=CategoricalArray(:genotype; levels=genotypes))
    fillcombinations([:rep, :generation, :genotype]; fill=0)
    sort([:generation, :rep, :genotype])
end

# gmiss_gen = Vector{Union{Missing, Int}}(missing, 9)
# counts = [missing for gen in 1:maximum(rr_geno_counts.generation), rep in 1:3]
# counts = Matrix{Union{Missing, Vector{Int}}}(missing, maximum(rr_geno_counts.generation), 3)
# n_obs = fill(50, size(counts))
# for (key, sdf) in pairs(groupby(rr_geno_counts, [:rep, :generation]))
#     counts[key.generation, key.rep] = sdf.count
#     n_obs[key.generation, key.rep] = sum(skipmissing(sdf.count))
# end
counts = Matrix{Vector{Int}}(undef, 5, 3)
n_obs = Matrix{Int}(undef, size(counts))
for (key, sdf) in pairs(groupby(rr_geno_counts, [:rep, :generation]))
    counts[key.generation .÷ 2, key.rep] = sdf.count
    n_obs[key.generation .÷ 2, key.rep] = sum(sdf.count)
end

mod_rr = popgen_model_rr(counts, n_obs, 10, 2:2:10)

ch_fit = sample(mod_rr, NUTS(1000, 0.75; init_ϵ=0.05), MCMCThreads(), 500, 2)

using HDF5, MCMCChainsStorage

h5open("fits/rr-two-locus-fit-7-22.h5", "w") do f
  	write(f, ch_fit)
end

ch_fit = h5open("fits/rr-two-locus-fit-7-22.h5", "r") do f
    from_mcmcchains(read(f, MCMCChains.Chains))
end

fitness_coefs = @chain DataFrame(FlexiChains.FlexiChain(ch_fit)) begin
    @select(Between("h₁[1]", "s₃[3]"))
    stack()
    @transform(
        :param=string.(getindex.(:variable, 1)),
        :locus=replace.(string.(getindex.(:variable, 2)), "₁" => "LV", "₂" => "VI", "₃" => "LI"),
        :rep=string.(getindex.(:variable, 6))
    )
end

spec = data(fitness_coefs) *
    mapping(:value, color=:rep, col=:locus => presorted, row=:param) *
    AlgebraOfGraphics.density()

draw(spec; facet=(;linkxaxes=:none, linkyaxes=:minimal))

save("figs/silvie-two-locus-fitness.pdf", draw(spec; facet=(;linkxaxes=:none, linkyaxes=:minimal)))

mod_rr_pred = popgen_model_rr([missing for gen in 1:size(counts, 1), rep in 1:size(counts, 2)], n_obs, 10, 2:2:10)
ch_yhat = predict(mod_rr_pred, ch_fit)

yhat = @chain DataFrame(ch_yhat) begin
    stack(Not(:iteration, :chain))
    @rtransform(
        :generation=2parse(Int, match(r"\[(\d)", :variable).captures[1]),
        :rep=parse(Int, match(r"(\d)\]\[", :variable).captures[1]),
        :genotype=genotypes[parse(Int, match(r"\[(\d+)\]", :variable).captures[1])]
    )
    groupby([:rep, :generation, :genotype])
    @combine(
        :mean=mean(:value),
        :lo=quantile(:value, 0.025),
        :hi=quantile(:value, 0.975)
    )
    @transform(:genotype=CategoricalArray(:genotype; levels=genotypes))
    # combine(Not(:iterati|on, :chain) .=> [mean x->quantile(x, 0.025)])
end

spec_pred = data(yhat) *
    (mapping(:generation, :lo, :hi, col=:genotype, row=:rep) *
    visual(Band; alpha=0.5) +
    mapping(:generation, :mean, col=:genotype, row=:rep) *
    visual(Lines; alpha=0.5))

spec_dat = data(rr_geno_pred) *
    mapping(:generation, :count, col=:genotype, row=:rep) *
    visual(ScatterLines)

draw(spec_pred + spec_dat)

### TODO kinda wish I knew why this didn't work, but whatever
rr_geno_counts_resrt = sort(rr_geno_counts, [:generation, :rep, :genotype])
geno_count_ret = returned(mod_rr, ch_fit)

ghat_stack = stack(x->vcat(x...), vec(geno_count_ret))

vcat(vec(geno_count_ret)[1]...)
rr_geno_pred = @chain rr_geno_counts_resrt begin
    @transform(
        :g_mean=mean.(eachrow(ghat_stack)),
        :g_lo=quantile.(eachrow(ghat_stack), 0.025),
        :g_hi=quantile.(eachrow(ghat_stack), 0.975)
    )
    groupby([:rep, :generation])
    @transform(:prop=:count ./ sum(:count))
end

spec_pred = data(rr_geno_pred) *
    (mapping(:generation, :g_lo, :g_hi, col=:genotype, row=:rep) *
    visual(Band; alpha=0.5) +
    mapping(:generation, :g_mean, col=:genotype, row=:rep) *
    visual(Lines; alpha=0.5))

spec_dat = data(rr_geno_pred) *
    mapping(:generation, :prop, col=:genotype, row=:rep) *
    visual(ScatterLines)

draw(spec_pred + spec_dat)