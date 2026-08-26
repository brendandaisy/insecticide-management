using Distributions, Turing

# assumes genotypes are in the ordering:
#  "VV-VV"
#  "VL-VV"
#  "LL-VV"
#  "VV-VI"
#  "VL-VI"
#  "LL-VI"
#  "VV-II"
#  "VL-II"
#  "LL-II"
function haplo2geno_rr(h)
    g = [
        h[1]^2, h[1]h[2], h[2]^2, h[1]h[3], h[1]h[4]+h[2]h[3], h[2]h[4], h[3]^2, h[3]h[4], h[4]^2
    ]
    g ./ sum(g)
end

@model function fitness_priors_rr(num_reps)
	h₁ ~ filldist(Uniform(-1, 1), num_reps)
	h₂ ~ filldist(Uniform(-1, 1), num_reps)
    h₃ ~ filldist(Uniform(-1, 1), num_reps)

    s₁ ~ filldist(Uniform(0, 1), num_reps)
    s₂ ~ filldist(Uniform(0, 1), num_reps)
    s₃ ~ filldist(Uniform(0, 1), num_reps)

    map(1:num_reps) do i
        fitness_weights([0, h₁[i], h₂[i], h₃[i]], [0, s₁[i], s₂[i], s₃[i]])
    end
end

# TODO almost definitely want to make inner vector of matrices static vectors
@model function popgen_model_rr(Y, N, num_gen, obs_gen, ::Type{T}=Float64) where {T}
    num_reps = size(Y, 2)

    r ~ Truncated(Normal(0, 0.1), 0, 1)
    Ws ~ to_submodel(fitness_priors_rr(num_reps), false)

    init_haplos = Vector{Vector{T}}(undef, num_reps)

    haplo_freq = Matrix{Vector{T}}(undef, num_gen, num_reps)
    geno_probs = Matrix{Vector{T}}(undef, size(Y))

	for rep in 1:num_reps
        modf = x -> update_two_locus(x, Ws[rep], r)

        init_haplos[rep] ~ Dirichlet(fill(1, 3))
        X₀ = zeros(T, 4)
        X₀[[1, 2, 4]] = init_haplos[rep]

        @views begin
            haplo_freq[:, rep] = iter_popgen(modf, X₀, num_gen)            
            geno_probs[:, rep] = haplo2geno_rr.(haplo_freq[obs_gen, rep])
        end

        for obs_t in axes(Y, 1)
            # @show (rep, t)
            # println(length(geno_probs[obs_t, rep]))
            try
                Y[obs_t, rep] ~ Multinomial(N[obs_t, rep], geno_probs[obs_t, rep])
            catch e
                @show (rep, obs_t)
                println(Ws[rep])
            end
            
        end
    end
    return geno_probs
    # Y ~ arraydist(Multinomial.(N, geno_probs))
end