struct QuasiBinomial{T<:Real} <: DiscreteUnivariateDistribution
    n::Int
    p::T
    β::T
    M::Int
end

function QuasiBinomial(n::Integer, p::T, β::T, M::Int) where {T<:Real}
    @assert n >= 1
    @assert β >= 0
    @assert p >= 0 && p <= 1 "p is not a probability."

    return QuasiBinomial{T}(n, p, β, M)
end

QuasiBinomial(n::Int, p::T, β::T) where {T <: Real} = QuasiBinomial{T}(n, p, β, 50)

Distributions.params(d::QuasiBinomial) = (d.n, d.p, d.β)
Distributions.mean(d::QuasiBinomial) = d.n * d.p

function Distributions.insupport(d::QuasiBinomial, x::T) where T<:Real
    if !isinteger(x) || x < 0 || x > d.n
        return false
    end
    return true
end

function Distributions.logpdf(d::QuasiBinomial, x::T) where T<:Real
    p = d.p
    n = d.n
	β = d.β
    S = eltype(p)
    R = promote_type(T, S)
    insupport(d, x) || return -R(Inf)
    s = R(loggamma(n+1) - loggamma(x+1) - loggamma(n-x+1))
	s -= R(xlogy(n-1, 1+β*n))
	s += R(log(p))
	s += R(log(1-p))
	s += R(xlogy(x-1, p+β*x))
	s += R(xlogy(n-x-1, 1 - p + β*(n-x)))
	return s
end

function acc_rate_qb_bb(x, a₁, a₂, n)
	num = gamma(a₁+n) * gamma(a₂+1) * (a₁+x)^(x-1) * (a₂+n-x)^(n-x-1)
	denom = gamma(a₁+x) * gamma(a₂+n-x) * (a₁+n)^(n-1)
	num / denom
end

function accept_reject_qb(rng::AbstractRNG, a₁::T, a₂::T, n, M) where {T<:Real}
    swap = false
    if a₂ < a₁
        tmp = a₁
        a₁ = a₂
        a₂ = tmp
        swap = true
    end
    bb = BetaBinomial(n, a₁, a₂)
    qb = QuasiBinomial(n, a₁/(a₁+a₂), 1/(a₁+a₂))
    x = rand(rng, bb)
    u = rand(rng, T)
    # while u > acc_rate_qb_bb(x, a₁, a₂, n)
    while u > pdf(qb, x) / (M*pdf(bb, x))
        x = rand(rng, bb)
        u = rand(rng, T)
    end
    swap ? n - x : x
end

function Distributions.rand(rng::AbstractRNG, d::QuasiBinomial{T}) where {T <: Real}
    p = d.p
    n = d.n
    β = d.β
    a₁ = p / β
    a₂ = (1 - p) / β
    accept_reject_qb(rng, a₁, a₂, n, d.M)
end

struct QuasiMultinomial{T<:Real, TV<:AbstractVector{T}} <: DiscreteMultivariateDistribution
    n::Int
    p::TV
    β::T
    QuasiMultinomial{T, TV}(n::Int, p::TV, β::T) where {T <: Real, TV <: AbstractVector{T}} = new{T, TV}(n, p, β)
end

function QuasiMultinomial(n::Integer, p::AbstractVector{T}, β::T) where {T<:Real}
    @assert n >= 1
    @assert β >= 0
    @assert isprobvec(p) "p is not a probability vector."

    return QuasiMultinomial{T, typeof(p)}(n, p, β)
end

Distributions.length(d::QuasiMultinomial) = length(d.p)
Distributions.params(d::QuasiMultinomial) = (d.n, d.p, d.β)

function Distributions.insupport(d::QuasiMultinomial, x::AbstractVector{T}) where T<:Real
    k = length(d)
    length(x) == k || return false
    s = 0.0
    for i = 1:k
        xi = x[i]
        if !(isinteger(xi) && xi >= 0)
            return false
        end
        s += xi
    end
    return s == d.n
end

function Distributions.logpdf(d::QuasiMultinomial, x::AbstractVector{T}) where T<:Real
    p = d.p
    n = d.n
	β = d.β
    S = eltype(p)
    R = promote_type(T, S)
    insupport(d, x) || return -R(Inf)
    s = R(loggamma(n + 1))
	s -= R(xlogy(n-1, 1+β*n))
    for i in eachindex(p)
        xi = x[i]
        p_i = p[i]
        s -= R(loggamma(R(xi) + 1))
		s += log(p_i)
        s += xlogy(xi-1, p_i + β*xi)
    end
    return s
end

function Distributions.rand(rng::AbstractRNG, d::QuasiMultinomial)
	p = d.p
    n = d.n
	β = d.β
	k = length(d)

	x = zeros(typeof(n), k)
    z = zero(eltype(p))
    rp = oftype(z + z, 1) # remaining total probability (widens type if needed)
    i = 0
    km1 = k - 1

    while i < km1 && n > 0
        i += 1
        p_i = p[i]
        if p_i < rp
            xi = rand(rng, QuasiBinomial(n, p_i / rp, β))
            x[i] = xi
            n -= xi
            rp -= p_i
        else
            # In this case, we don't even have to sample
            # from Binomial. Just assign remaining counts
            # to xi.

            x[i] = n
            n = 0
        end
    end

    if i == km1
        x[k] = n
	end
    # else  # n must have been zero
    #     z = zero(eltype(x))
    #     for j = i+1 : k
    #         x[j] = z
    #     end
    # end

    return x
end