function get_generations(data)
	grouped = unique(select(data, [:site, :rep, :generation]))
	sort!(grouped, [:site, :rep, :generation])

	group_generations = [
		Int.(group.generation)
		for group in groupby(grouped, [:site, :rep]; sort=false)
	]

	all_generations = Int.(grouped.generation)

	return (
		values = all_generations,
		by_group = group_generations,
		n = maximum.(group_generations)
	)
end

function prepare_vm_two_loci(data)
	required_columns = [:site, :rep, :generation, :genotype, :N, :count]
	all(column -> column in propertynames(data), required_columns) ||
		throw(ArgumentError("VM data is missing a required column."))

	genotype_index = Dict(
		genotype => index
		for (index, genotype) in enumerate(L2_L3_GENOTYPES)
	)

	site_levels = sort!(unique(data.site))
	site_index = Dict(site => index for (index, site) in enumerate(site_levels))

	group_keys = unique(select(data, [:site, :rep]))
	sort!(group_keys, [:site, :rep])
	group_index = Dict(
		(row.site, row.rep) => index
		for (index, row) in enumerate(eachrow(group_keys))
	)

	grouped = unique(select(data, [:site, :rep, :generation]))
	sort!(grouped, [:site, :rep, :generation])
	generation_data = get_generations(data)
	max_generation = maximum(generation_data.n)

	counts = zeros(Int, length(L2_L3_GENOTYPES), max_generation, nrow(group_keys))
	sample_sizes = zeros(Int, max_generation, nrow(group_keys))

	for key in eachrow(grouped)
		group = group_index[(key.site, key.rep)]
		generation = Int(key.generation)
		rows = data[
			(data.site .== key.site) .&
			(data.rep .== key.rep) .&
			(data.generation .== key.generation),
			:
		]

		nrow(rows) == length(L2_L3_GENOTYPES) ||
			throw(ArgumentError(
				"Each VM site/replicate/generation must have exactly nine genotype rows."
			))

		seen = falses(length(L2_L3_GENOTYPES))
		for row in eachrow(rows)
			index = get(genotype_index, String(row.genotype), 0)
			index == 0 && throw(ArgumentError(
				"Unexpected two-locus genotype: $(row.genotype)."
			))
			seen[index] && throw(ArgumentError(
				"Duplicate genotype $(row.genotype) in VM data."
			))
			seen[index] = true
			counts[index, generation, group] = Int(row.count)
			sample_sizes[generation, group] = Int(row.N)
		end

		all(seen) || throw(ArgumentError(
			"A VM observation is missing at least one genotype category."
		))
		sum(counts[:, generation, group]) == sample_sizes[generation, group] ||
			throw(ArgumentError(
				"VM genotype counts do not sum to N for " *
				"$(key.site), replicate $(key.rep), generation $(key.generation)."
			))
	end

	return (
		counts = counts,
		sample_sizes = sample_sizes,
		n_generations = generation_data.n,
		group_site = [site_index[row.site] for row in eachrow(group_keys)],
		n_sites = length(site_levels),
		n_groups = nrow(group_keys),
		site_levels = site_levels,
		group_keys = group_keys
	)
end