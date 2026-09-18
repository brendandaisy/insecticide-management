using Distributions
using Turing

function latent_risk(i, θ, U, V, X)
	α₀, α, β₀, γ = θ
	allele_risk = @. α*γ^(V[i, :])  * U[i, :] 
	α₀ + β₀*X[i] + sum(allele_risk)
end

# ╔═╡ acfdaa61-1bfd-4a39-a0d5-3a8211a3b9da
@model function dose_response(U, V, X, N)
	α₀ ~ Normal(0, 2)
	α ~ filldist(truncated(Normal(0, 2), upper=0), 3)
	β₀ ~ truncated(Normal(0, 3); lower=0)
	# β ~ filldist(Normal(0, 1), 3)
	ν ~ Gamma(1, 1)
	γ ~ filldist(Beta(1, 1), 3)
	
	# σ_g ~ truncated(Cauchy(0, 2.5); lower=0)
	# g_raw ~ filldist(Normal(0, 1), 6)
	# g = g_raw .* σ_g

	η = map(i->latent_risk(i, (α₀, α, β₀, γ), U, V, X), 1:length(y))
	r = max.(1e-4, min.(0.9999, cdf.(Normal(), η)))
	y ~ arraydist(QuasiBinomial.(N, r, ν))

	# return r
end