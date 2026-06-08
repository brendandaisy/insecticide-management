using InvertedIndices
using LinearAlgebra
using Turing
using HDF5, MCMCChainsStorage
using Distributions
using Random
using StatsFuns
using SpecialFunctions
using CategoricalArrays
using CSV, DataFrames, DataFramesMeta
using CairoMakie, AlgebraOfGraphics

set_aog_theme!()

include("src/Quasinomials.jl")
include("src/helpers.jl")

vm = CSV.read("data-proc/vera-maloof.csv", DataFrame);

vm_long = @chain vm begin
    select(Not(Cols(r"IF"), Between("VF", "IC")))
    stack(Between(:VFVF, :ICIC); variable_name="genotype", value_name="count")
    @rtransform(:rep = "$(:site)$(:rep)")
    @groupby([:site, :rep, :generation])
    @transform(
        :N = sum(:count),
        :genotype = CategoricalArray(:genotype)
    )
end

Y = Matrix{Int}[]
n_obs = Vector{Int}[]
for vmt in groupby(vm_long, :generation)
    df = @select(vmt, :rep, :genotype, :count)
    df = unstack(df, :rep, :count)
    dm = Matrix{Int}(df[:, Not(1)])
    push!(Y, dm)
    push!(n_obs, vec(sum(dm, dims=1)))
end

reps = unique(vm_long.rep)

inits = norm_floor.(eachcol(Y[1]), 0) |> stack

### Huibjen lab RR data
rr = CSV.read("data-raw/Resistance.Reversal_ALL.kdr_2025-07-24.csv", DataFrame)

rr_long = @chain rr begin
    select(Not(1))
    @rtransform(
        :gen_410=replace(:gen_410, "S" => "V", "R" => "L"),
        :gen_1016=replace(:gen_1016, "S" => "V", "R" => "I"),
    )
    @rselect(
        :site="RR",
        :generation=:Generation,
        :rep="RR" * string(:Replicate),
        :genotype="$(:gen_410)-$(:gen_1016)"
    )
    @groupby All()
    combine(nrow => :count)
end
# TODO how to get to haplotype?

# original haplotypes: VV, LV, LI (aka no VI!)
geno2haplo = Dict(
    "LL-II" => "LI/LI",
    "LL-VI" => "LV/LI",
    "LL-VV" => "LV/LV",
    "VL-II" => "VI/LI", #### pretty common :(
    "VL-VI" => ["VV/LI", "LV/VI"], # OR LV/VI
    "VL-VV" => "VV/LV",
    "VV-II" => "VI/VI", #### very rare
    "VV-VI" => "VV/VI", #### rare
    "VV-VV" => "VV/VV"
)

for g in rr_long[!, :genotype]
    println(geno2haplo[g])
end

for rr_t in groupby(rr_long, :generation)
    rrtwide = unstack(rr_t, :rep, :count; fill=0)
    println(Matrix{Int}(@select(rrtwide, Cols(r"RR\d+"))))
    # push!(ysomething, Matrix{Int}(@select(rrtwide, Cols(r"RR\\d+"))))
end

rr_viz = data(rr_long) *
    mapping(:generation, :count; color=:rep, layout=:genotype) *
    visual(ScatterLines)

draw(rr_viz)




