using Distributions, Turing

function haplo2geno_410(x1)
    g = [x1[1]^2, x1[1]x1[2], x1[2]^2]
    g ./ sum(g)
end

# assumes haplottypes are in the ordering:
# "VF"
# "VC"
# "IF"
# "IC"
# assumes genotypes are in the ordering:
# "VV-FF"
# "VV-FC"
# "VV-CC"
# "VI-FF"
# "VI-FC"
# "VI-CC"
# "II-FF"
# "II-FC"
# "II-CC"
function haplo2geno_1016_1534(x23)
    g = [
        x23[1]^2, x23[1]x23[2], x23[2]^2, x23[1]x23[3], x23[1]x23[4]+x23[2]x23[3], 
        x23[2]x23[4], x23[3]^2, x23[3]x23[4], x23[4]^2
    ]
    g ./ sum(g)
end

# TODO test me please!
# TODO also, consider moving to different script
function initial_haplotype_probs(Y1, Y23)
    num_rep = size(Y23, 3)
    probs = zeros(8, num_rep)

    VF_ind = [1, 5]
    VC_ind = [2, 6]
    IF_ind = [3, 7]
    IC_ind = [4, 8]

    for rep in 1:num_rep
        counts1 = view(Y1, :, 1, rep)
        probs_rep = view(probs, :, rep)
        probs_rep[1:4] += counts1[1]
        # probs_rep += counts1[2] not necessary yeah?
        probs_rep[5:8] += counts1[3]

        counts23 = view(Y23, :, 1, rep)
        probs_rep[VF_ind] += counts23[1]
        probs_rep[VF_ind] += counts23[2] / 2
        probs_rep[VC_ind] += counts23[2] / 2
        probs_rep[VC_ind] += counts23[3]
        probs_rep[VF_ind] += counts23[4] / 2
        probs_rep[IF_ind] += counts23[4] / 2
        probs_rep[VC_ind] += counts23[6] / 2
        probs_rep[IC_ind] += counts23[6] / 2
        probs_rep[IF_ind] += counts23[7]
        probs_rep[IF_ind] += counts23[8] / 2
        probs_rep[IC_ind] += counts23[8] / 2
        probs_rep[IC_ind] += counts23[9]

        probs_rep += sum(probs_rep)
    end

    return probs_rep
end

@model function popgen_model_vm(Y1, Y23, N1, N23, ::Type{T}=Float64) where {T}
    num_rep = size(Y23, 3)
    num_gen = size(Y23, 2)

    # priors
    r12 ~ Beta(1, 3)
    r23 ~ Beta(1, 3)

    γ ~ Beta(1, 1)
    λ ~ Beta(1, 1)

	h ~ filldist(TriangularDist(-1, 1, γ), num_rep)
	s ~ filldist(TriangularDist(-1, 1, λ), num_rep)
    pushfirst!.(h, 0)
    pushfirst!.(s, 0)
    Ws = build_fitness_matrix.(s, h; checks=false)

    # initial conditions
    haplos_init = to_submodel(initial_haplotypes(Y1, Y23), false)

    for rep in 1:num_rep
        sim = build_model(Ws[rep]; rAB = r12, rBC = r23, checks = false)
        x = view(haplos_init, :, :, rep)
        xnext = similar(h)
        for gen in 1:num_gen
            if gen > 1
                p1, p23 = get_marginal_frequencies(x)
                Y1[:, gen, rep] .~ Multinomial(N1[gen, rep], p1)
                Y23[:, gen, rep] .~ Multinomial(N23[gen, rep], haplo2geno_1016_1534(p23))
            end
            step!(xnext, x, sim; checks=false)
            x = xnext
        end
    end   
end