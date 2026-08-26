using Distributions
using InvertedIndices
using CSV, DataFrames, DataFramesMeta
using CairoMakie, AlgebraOfGraphics
using Turing
using HDF5, MCMCChainsStorage

include("../src/popgen-three-locus.jl")
include("../src/popgen-turing.jl")

vm410 = CSV.read("data-proc/vm-geno-counts-410.csv", DataFrame)

