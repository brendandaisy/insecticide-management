using Distributions
using InvertedIndices
using CSV, DataFrames, DataFramesMeta
using XLSX
using CairoMakie, AlgebraOfGraphics

column_labels = ["site", "generation", "rep", "N", "VV-FF", "VV-FC", "VV-CC", "IV-FF", "IV-FC", "IV-CC", "II-FF", "II-FC", "II-CC"]
vma0 = XLSX.readtable(
    "data-raw/Vera-Maloof tableS2.xlsx", 1, "A:M"; 
    header=false, column_labels
) |> DataFrame

unique(vma0[!, :site])

filldown(v) = accumulate((x, y) -> coalesce(y, x), v, init = v[1])

vma = vma0[Not(vcat(1:2, 53:55, 109:112, 166:169, 223:226, 280:283)),:]
vma = @chain vma begin
    @transform(:site=filldown(:site), :generation=filldown(:generation))
    @subset(:rep .!= "Total") # implicitely removes rows with :rep === missing
    transform(
        :generation => x->parse.(Int, last.(x)),
        Between(:rep, "II-CC") .=> ByRow(Int);
        renamecols=false
    )
    stack(Between("VV-FF", "II-CC"); variable_name="genotype", value_name="count")
    @by([:site, :rep, :generation], :N=sum(:count), :count, :prop=:count ./ sum(:count))
end

CSV.write("data-proc/vm-geno-counts-1016-1534.csv", vma)

vma_agg = @chain vma begin
    groupby([:generation, :genotype])
    @combine(:prop=sum(:count) / sum(:N))
end

spec = data(vma_agg) *
    mapping(:generation, :prop, color=:genotype, layout=:genotype) *
    visual(ScatterLines)

draw(spec)

### now the 410 data
column_labels = ["id", "Dz.1", "Dz.8", "Tap.1", "Tap.8", "Mer3.1", "Mer3.8", "Mer2.1", "Mer2.8", "Acp.1", "Acp.8", "Ac.1", "Ac.8", "Mer1.1", "Mer1.8", "Co.1", "Co.8"]
vmb0 = XLSX.readtable("data-raw/Old DNA PCR results(in).xlsx"; column_labels, missing_strings="N/A") |> DataFrame

vmb = @chain vmb0[Not(49:53), :] begin
    stack(Not(:id))
    @subset(.!ismissing.(:value))
    @select(
        :site=first.(split.(:variable, ".")),
        :generation=parse.(Int, last.(split.(:variable, "."))),
        :genotype=reverse.(replace.(:value, "G" => "V", "T" => "L"))
    )
    @by([:site, :generation, :genotype], :count=length(:genotype))
    fillcombinations([:site, :generation, :genotype]; fill=0)
    groupby([:site, :generation])
    @transform(:N=sum(:count), :prop=:count ./ sum(:count))
end

CSV.write("data-proc/vm-geno-counts-410.csv", vmb)

spec = data(vmb) *
    mapping(
        :generation => nonnumeric, 
        :prop, 
        color=:site, 
        col=:genotype => sorter(["VV", "VL", "LL"])
    ) *
    visual(ScatterLines)

set_aog_theme!()
draw(spec, scales(Color=(;palette=:tab10)); figure=(;size=(620, 350)))