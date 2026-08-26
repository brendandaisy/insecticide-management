using Distributions
using InvertedIndices
using CSV, DataFrames, DataFramesMeta
using CairoMakie, AlgebraOfGraphics

vm01 = CSV.read("data-proc/dose-response.csv", DataFrame)

vm02 = XLSX.readtable()