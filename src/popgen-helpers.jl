

"""
Convert three-locus diploid genotype counts (strings such as "VL-VI-FC",
ordered 410-1016-1534) into frequencies of the eight three-locus haplotypes,
ordered as `HAPLOTYPE_NAMES`.

Unphased genotypes with h >= 2 heterozygous loci are split equally across
their 2^(h-1) phase resolutions.
"""
function genotype_counts_to_haplotype_freq(genotypes, counts; checks=true)
    checks && length(genotypes) != length(counts) &&
        throw(DimensionMismatch("genotypes and counts must have the same length."))

    haplotype_index = Dict(
        name => index for (index, name) in enumerate(HAPLOTYPE_NAMES)
    )
    totals = zeros(Float64, N_HAPLOTYPES)

    for (genotype, count) in zip(genotypes, counts)
        loci = split(genotype, '-')
        checks && !(length(loci) == 3 && all(locus -> length(locus) == 2, loci)) &&
            throw(ArgumentError("Unexpected three-locus genotype: $genotype."))

        heterozygous = [i for i in 1:3 if loci[i][1] != loci[i][2]]
        n_phases = 2^max(length(heterozygous) - 1, 0)

        for phase in 0:(n_phases - 1)
            haplotype_1 = [loci[i][1] for i in 1:3]
            haplotype_2 = [loci[i][2] for i in 1:3]

            # The first heterozygous locus keeps its orientation; the rest are swapped by bit.
            for (bit, locus) in enumerate(heterozygous[2:end])
                if (phase >> (bit - 1)) & 1 == 1
                    haplotype_1[locus], haplotype_2[locus] =
                        haplotype_2[locus], haplotype_1[locus]
                end
            end

            for haplotype in (haplotype_1, haplotype_2)
                index = get(haplotype_index, String(haplotype), 0)
                checks && index == 0 && throw(
                    ArgumentError("Unexpected alleles in genotype: $genotype.")
                )
                totals[index] += count / n_phases
            end
        end
    end

    return totals ./ sum(totals)
end

"""
Convert three-locus diploid genotype counts (strings such as "VL-VI-FC",
ordered 410-1016-1534) into resistant allele frequencies, returned as
`(L=, I=, C=)` for 410, 1016 and 1534 respectively.
"""
function genotype_counts_to_allele_freq(genotypes, counts; checks=true)
    checks && length(genotypes) != length(counts) &&
        throw(DimensionMismatch("genotypes and counts must have the same length."))

    resistant = ('L', 'I', 'C')
    resistant_copies = zeros(Float64, 3)

    for (genotype, count) in zip(genotypes, counts)
        loci = split(genotype, '-')
        checks && !(length(loci) == 3 && all(locus -> length(locus) == 2, loci)) &&
            throw(ArgumentError("Unexpected three-locus genotype: $genotype."))

        for locus in 1:3
            resistant_copies[locus] += count * count_alleles(loci[locus], resistant[locus])
        end
    end

    frequencies = resistant_copies ./ (2 * sum(counts))
    return (L=frequencies[1], I=frequencies[2], C=frequencies[3])
end

count_alleles(locus, allele) = count(==(allele), locus)

"""
Convert nine two-locus diploid genotype counts into expected counts
of the four 1016/1534 haplotypes: VF, VC, IF, and IC.

The phase of VI-FC is unobserved. Its count is split equally between
the VF/IC and VC/IF phases.
"""
function L2_L3_genotype_counts_to_haplotype_counts(counts)
    length(counts) == length(L2_L3_GENOTYPES) ||
        throw(DimensionMismatch("counts must have length 9."))

    VV_FF, VV_FC, VV_CC,
    VI_FF, VI_FC, VI_CC,
    II_FF, II_FC, II_CC = counts

    return [
        2 * VV_FF + VV_FC + VI_FF + 0.5 * VI_FC,
        2 * VV_CC + VV_FC + VI_CC + 0.5 * VI_FC,
        2 * II_FF + VI_FF + II_FC + 0.5 * VI_FC,
        2 * II_CC + VI_CC + II_FC + 0.5 * VI_FC
    ]
end

"""
Construct eight three-locus pseudo-counts from four two-locus counts.

The 410 locus is unobserved, so each 1016/1534 haplotype is split
equally across the V and L backgrounds.
"""
function haplotype_counts_to_three_locus_counts(counts)
    error("Not sure where I need this yet")
    # two_locus_counts = L2_L3_genotype_counts_to_haplotype_counts(counts)
    # return vcat(two_locus_counts ./ 2, two_locus_counts ./ 2)
end