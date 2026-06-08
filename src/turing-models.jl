@model function fitness_priors(num_reps)
	φ₁ ~ Beta(1, 1); λ₁ ~ Pareto(0.5, 1.5)
	φ₂ ~ Beta(1, 1); λ₂ ~ Pareto(0.5, 1.5)

	h₁ ~ filldist(Beta(λ₁*φ₁, λ₁*(1-φ₁)), num_reps)
	h₂ ~ filldist(Beta(λ₂*φ₂, λ₂*(1-φ₂)), num_reps)
	
	μ₁ ~ Gamma(1, 1); ν₁ ~ truncated(Cauchy(0, 2.5); lower=num_reps)
	μ₂ ~ Gamma(1, 1); ν₂ ~ truncated(Cauchy(0, 2.5); lower=num_reps)
	
	s₁ ~ filldist(truncated(Normal(μ₁, ν₁), lower=0), num_reps)
	s₂ ~ filldist(truncated(Normal(μ₂, ν₂), lower=0), num_reps)

	return fitness_weights.(h₁, h₂, s₁, s₂) # broadcase over all reps
end

@model function fitness_priors_old()
	φ₁ ~ Beta(1.5, 5); λ₁ ~ Pareto(1, 2)
	φ₂ ~ Beta(1.5, 3); λ₂ ~ Pareto(1, 2)
	φ₃ ~ Beta(5, 1.5); λ₃ ~ Pareto(1, 1.5)
	φ₄ ~ Beta(3, 1.5); λ₄ ~ Pareto(1, 1.5)

	h₁ ~ filldist(Beta(λ₁*φ₁, λ₁*(1-φ₁)), 8)
	h₂ ~ filldist(Beta(λ₂*φ₂, λ₂*(1-φ₂)), 8)
	s₁ ~ filldist(Beta(λ₃*φ₃, λ₃*(1-φ₃)), 8)
	s₂ ~ filldist(Beta(λ₄*φ₄, λ₄*(1-φ₄)), 8)

	return fitness_weights_old.(h₁, h₂, s₁, s₂)
end

@model function ir_model(y, A₀, C, Tₘ, N, M, ::Type{T}=Float64) where {T}
	δ ~ Beta(1, 100)
	ω ~ to_submodel(fitness_priors(8), false)
	ωrep = repeat(ω, inner=3)
    ν ~ Gamma(1, 1)

	A = [Matrix{T}(undef, 6, num_rep) for t=1:num_gen]
	p = Matrix{T}(undef, 6, num_rep)

	A[1] = A₀ .+ δ
	for t in 2:Tₘ
		for rep in 1:M
			h = C * (ωrep[rep] .* A[t-1][:,rep])
			h /= sum(h)
			p[:,rep] = [h[1]^2, h[2]^2, h[3]^2, 2h[1]*h[2], 2h[1]*h[3], 2h[2]*h[3]]
			A[t][:,rep] = p[:,rep] / sum(p[:,rep])
		end
		y[t] ~ arraydist([QuasiMultinomial(N[t][r], A[t][:,r], ν) for r=1:M])
	end
end

@model function fitness_priors_comb(Mvm, Mrr)
	M = Mvm + Mrr

	# dominance param for 410 - RR only
	φ₁ ~ Beta(1, 1); λ₁ ~ Pareto(0.5, 1.5)
	h₁ ~ filldist(Beta(λ₁*φ₁, λ₁*(1-φ₁)), Mrr)

	# dominance for 1016 - both VM and RR
	φ₂ ~ Beta(1, 1); λ₂ ~ Pareto(0.5, 1.5)
	h₂ ~ filldist(Beta(λ₂*φ₂, λ₂*(1-φ₂)), M)

	# dominance for 1536 - 
	φ₂ ~ Beta(1, 1); λ₂ ~ Pareto(0.5, 1.5)
	h₂ ~ filldist(Beta(λ₂*φ₂, λ₂*(1-φ₂)), M)
	
	μ₁ ~ Gamma(1, 1); ν₁ ~ truncated(Cauchy(0, 2.5); lower=M)
	μ₂ ~ Gamma(1, 1); ν₂ ~ truncated(Cauchy(0, 2.5); lower=M)
	
	s₁ ~ filldist(truncated(Normal(μ₁, ν₁), lower=0), M)
	s₂ ~ filldist(truncated(Normal(μ₂, ν₂), lower=0), M)

	return fitness_weights.(h₁, h₂, s₁, s₂) # broadcase over all reps
end

# TODO this will only apply to the VM data assuming we need to add recombination to the RR data
function proj_mat_vm(g, ω, r=fill(1, 6))
	F♀ = r .* ω
	F♂ = diagm(ω)
	# q = (C * F♂ * g) / sum(C * F♂ * g)
	q = (C * F♂ * g)
	A = [q[1]*F♀[1] 0 0 0.5*q[1]*F♀[4] 0.5*q[1]*F♀[5] 0;
	0 q[2]*F♀[2] 0 0.5*q[2]*F♀[4] 0  0.5*q[2]*F♀[6];
	0 0 q[3]*F♀[3] 0 0.5*q[3]*F♀[5] 0.5*q[3]*F♀[6];
	q[2]*F♀[1] q[1]*F♀[2] 0 0.5*(q[1]+q[2])*F♀[4] 0.5*q[2]*F♀[5] 0.5*q[1]*F♀[6];
	q[3]*F♀[1] 0 q[1]*F♀[3] 0.5*q[3]*F♀[4] 0.5*(q[1]+q[3])*F♀[5] 0.5*q[1]*F♀[6];
	0 q[3]*F♀[2] q[2]*F♀[3] 0.5*q[3]*F♀[4] 0.5*q[2]*F♀[5] 0.5*(q[2]+q[3])*F♀[6]]
	A / sum(A*g)
end

@model function ir_model_comb(y, A₀, C, Tₘ, N, M, ::Type{T}=Float64) where {T}
	δ ~ Beta(1, 100)
	ω ~ to_submodel(fitness_priors(), false)
	ωrep = repeat(ω, inner=3)
    ν ~ Gamma(1, 1)

	A = [Matrix{T}(undef, 6, M) for t=1:num_gen]
	p = Matrix{T}(undef, 6, M)

	A[1] = A₀ .+ δ
	for t in 2:Tₘ
		for rep in 1:M
			h = C * (ωrep[rep] .* A[t-1][:,rep])
			h /= sum(h)
			p[:,rep] = [h[1]^2, h[2]^2, h[3]^2, 2h[1]*h[2], 2h[1]*h[3], 2h[2]*h[3]]
			A[t][:,rep] = p[:,rep] / sum(p[:,rep])
		end
		y[t] ~ arraydist([QuasiMultinomial(N[t][r], A[t][:,r], ν) for r=1:M])
	end
end