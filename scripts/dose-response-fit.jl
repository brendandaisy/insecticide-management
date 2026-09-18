using InvertedIndices
using LinearAlgebra
using Distributions
using Turing
using HDF5, MCMCChainsStorage
# using StatsFuns
# using SpecialFunctions
using CSV, DataFrames, DataFramesMeta
using CategoricalArrays
using CairoMakie, AlgebraOfGraphics

include("../src/dose-response.jl")

set_aog_theme!()

dose0 = CSV.read("data-raw/Silvie data - dose-response insecticide mortality data.csv", DataFrame; normalizenames=true);

# ╔═╡ 5b7bd0fb-2c39-4614-b5b4-b195fc078397
dose = @chain dose0 begin
	@subset(:dose_per_mosq .> 0)
	@rtransform(:genotype="$(:locus3)-$(:locus1)-$(:locus2)")
end

viz = data(dose) * 
	mapping(:dose_per_mosq => log, :mortality => (x->quantile(Normal(), x)); color=:genotype) * 
	visual(Scatter)

draw(viz; axis=(limits=(nothing, (-2, 2)),))

# process data in required format for Turing model
loci = @select(dose, :locus3, :locus1, :locus2) |> Matrix

U = occursin.(r"[LIC]", loci)
V = occursin.(r"(\w)(?!\1)(\w)", loci)

Y = dose.dead
N = dose.total
X = log.(dose.dose_per_mosq)

dr_model = dose_response(U, V, X, N)

# TODO continue incorporating dose-response-model.jl using new conditioning format and
# switch to BetaBinomial