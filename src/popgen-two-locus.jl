using LinearAlgebra
import IterTools: iterated

function fitness_weights(h, s)
    W = diagm(1 .- s)
    for i in 1:(size(W, 1)-1)
        for j in (i+1):size(W, 2)
            W[i, j] = max(0, 1 - h[i]s[i] - h[j]s[j])
            W[j, i] = W[i, j]
        end
    end
    return W
end

function marginal_fitness(x, W)
    [sum(wj .* x) for wj in eachrow(W)]
end

linkage_disequilibrium(x) = x[1]*x[4] - x[2]*x[3]

linkage_disequilibrium(x, W) = W[1, 4]*x[1]*x[4] - W[2, 3]*x[2]*x[3]

function mean_fitness(x, W)
    sum(W .* (x * x'))
end

function iter_popgen(modf, xinit, num_gen)
    Iterators.take(iterated(modf, xinit), num_gen) |> collect
end

function update_two_locus(x, W, r)
    wi = marginal_fitness(x, W)
    rD = r * linkage_disequilibrium(x, W)
    W̄ = mean_fitness(x, W)
    @. (x * wi + [-1, 1, 1, -1] * rD) / W̄
end

# α = 0.2
# β = 0.1
# γ = -0.1
# δ = -0.2
# W = [1-δ 1-β 1-γ 1;
#     1-β 1-α 1 1-γ;
#     1-γ 1 1-α 1-β;
#     1 1-γ 1-β 1-δ]

# marginal_fitness(xinit, W)

# xinit = [0.4, 0.4, 0.1, 0.1]
# res = iter_popgen(x->update_two_locus(x, W, 0.05), xinit, 20)

# resdf = @chain stack(res)' begin
#     DataFrame(["AB", "Ab", "aB", "ab"])
#     @transform(:t=eachindex(res))
#     stack(Not(:t); variable_name=:haplotype)
# end

# spec = data(resdf) *
#     mapping(:t, :value, color=:haplotype) *
#     visual(ScatterLines)

# draw(spec)