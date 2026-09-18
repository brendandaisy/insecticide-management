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
    @by([:site, :rep, :generation], :genotype, :N=sum(:count), :count, :prop=:count ./ sum(:count))
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

##~~~~

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

# ============================================================
# updated F1 and F8 data for all three loci!
# ============================================================

# ============================================================
# Input / output
# ============================================================

input_file = "data-raw/2026-08-26_F1_F8_genotypes.xlsx"
input_cols = ["mosquito_id", "Dz.1", "Dz.8", "Tap.1", "Tap.8", "Mer3.1", "Mer3.8", "Mer2.1", "Mer2.8", "Acp.1", "Acp.8", "Ac.1", "Ac.8", "Mer1.1", "Mer1.8", "Co.1", "Co.8"]

# Sheets/loci, in the order they should appear in the
# combined genotype
loci = ["410", "1016", "1534"]


# ============================================================
# Genotype recoding
# ============================================================

genotype_maps = Dict(
    "410" => Dict(
        "TT" => "LL",
        "TG" => "VL",
        "GG" => "VV"
    ),

    "1016" => Dict(
        "AA" => "II",
        "AG" => "VI",
        "GG" => "VV"
    ),

    "1534" => Dict(
        "GG" => "CC",
        "TG" => "FC",
        "TT" => "FF"
    )
)


# ============================================================
# Read a locus sheet and convert to long format
# ============================================================

function read_locus(sheet, locus_name)
    raw = DataFrame(XLSX.readtable(input_file, sheet; column_labels=input_cols))

    # # First column is mosquito ID
    # rename!(raw, names(raw)[1] => :mosquito_id)

    # Remove summary/total rows at the bottom.
    # Actual mosquito IDs are numeric.
    raw = filter(row -> row.mosquito_id isa Number, raw)

    raw.mosquito_id = Int.(raw.mosquito_id)

    # Convert site x generation columns to long format
    long = stack(
        raw,
        Not(:mosquito_id),
        variable_name = :site_generation,
        value_name = :genotype
    )

    # Extract site and generation from column names

    # long.site = String[]
    # long.generation = Int[]

    # for x in long.site_generation

    #     parts = split(String(x), r"\s+")

    #     # Everything except final F1/F8 is the site name
    #     site = join(parts[1:end-1], " ")

    #     # Convert F1 -> 1, F8 -> 8, etc.
    #     generation = parse(
    #         Int,
    #         replace(parts[end], "F" => "")
    #     )

    #     push!(long.site, site)
    #     push!(long.generation, generation)
    # end

    @select!(
        long,
        :mosquito_id,
        :site=first.(split.(:site_generation, ".")),
        :generation=parse.(Int, last.(split.(:site_generation, "."))),
        :genotype
    )

    # Standardize missing/blank genotype calls
    long.genotype = map(long.genotype) do x

        if ismissing(x)
            missing
        elseif isempty(strip(String(x)))
            missing
        else
            strip(String(x))
        end

    end

    # Recode the locus-specific genotype
    genotype_map = genotype_maps[locus_name]

    long.genotype = map(long.genotype) do x

        if ismissing(x)
            missing
        elseif haskey(genotype_map, x)
            genotype_map[x]
        else
            error(
                "Unexpected genotype '$x' found at locus $locus_name"
            )
        end

    end

    # Rename genotype according to locus
    rename!(
        long,
        :genotype => Symbol("locus_", locus_name)
    )

    return long
end


# ============================================================
# Read all three loci
# ============================================================

locus_data = Dict(
    locus => read_locus(locus, locus)
    for locus in loci
)


# ============================================================
# Join the three loci by mosquito, site, and generation
# ============================================================

combined = select(
    locus_data[loci[1]],
    :mosquito_id,
    :site,
    :generation,
    Symbol("locus_", loci[1])
)

for locus in loci[2:end]

    tmp = select(
        locus_data[locus],
        :mosquito_id,
        :site,
        :generation,
        Symbol("locus_", locus)
    )

    combined = leftjoin(
        combined,
        tmp,
        on = [:mosquito_id, :site, :generation]
    )

end


# ============================================================
# Keep only mosquitoes with complete calls at all three loci
# ============================================================

locus_columns = Symbol.("locus_", loci)

combined = filter(
    row -> all(
        !ismissing(row[c])
        for c in locus_columns
    ),
    combined
)


# ============================================================
# Construct the three-locus genotype
#
# Example:
#
#   410  = TT -> LL
#   1016 = AG -> VI
#   1534 = TG -> FC
#
# becomes:
#
#   LL-VI-FC
# ============================================================

combined.genotype = [
    join(
        String[row[c] for c in locus_columns],
        "-"
    )
    for row in eachrow(combined)
]


# ============================================================
# Calculate N
#
# N is the total number of mosquitoes with complete genotype
# calls at all three loci within each site x generation.
# ============================================================

N_by_site_generation = combine(
    groupby(
        combined,
        [:site, :generation]
    ),
    nrow => :N
)


# ============================================================
# Count each three-locus genotype
# ============================================================

result = combine(
    groupby(
        combined,
        [:site, :generation, :genotype]
    ),
    nrow => :count
)


# ============================================================
# Add N to each genotype group
# ============================================================

result = leftjoin(
    result,
    N_by_site_generation,
    on = [:site, :generation]
)


# ============================================================
# Arrange columns
# ============================================================

select!(
    result,
    :site,
    :generation,
    :genotype,
    :count,
    :N
)


# ============================================================
# Sort output
# ============================================================

sort!(
    result,
    [:site, :generation, :genotype]
)

# ============================================================
# Plot results
# ============================================================
spec = data(result) *
        mapping(
            :generation, 
            (:count, :N) => ((x, y) -> x ./ y) =>  "proportion"; 
            color=:genotype, col=:site
        ) *
        visual(ScatterLines)

draw(spec, scales(Color=(;palette=:tab20)))

# ============================================================
# Write CSV
# ============================================================

CSV.write("data-proc/vm-geno-counts-new.csv", result)
