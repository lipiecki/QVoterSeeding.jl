"""
    seeding(network::Union{AbstractString, SimpleGraph}, strategy::Symbol, zealots::Bool, q::Integer, budgets::AbstractVector{<:AbstractFloat}, rng::AbstractRNG = Random.default_rng(); kwargs...)
Run network seeding experiments with `q`-voter dynamics with for a specified seeding `strategy` and seeding `budgets`. If `zealots = true`, seeds will act as zealots.

If `network` is an `AbstractString`, it is passed to [`loadnetwork`](@ref), along with the optional `rng` argument.

Available options for `strategy`:
- `:cc` - seeding nodes with the highest **c**omplex **c**entrality (adapted to q-voter dynamics)
- `:hd` - seeding nodes with the **h**ighest **d**egree centrality
- `:pr` - seeding nodes with the highest **p**age**r**ank centrality
- `:onehop` - seeding nodes using **one**-**hop**
- `:random` - seeding **random** nodes

Return a vector of final concentrations of active vertices for each seeding budget.    
"""
function seeding(network::AbstractString, strategy::Symbol, zealots::Bool, q::Integer, budgets::AbstractVector{<:AbstractFloat}, rng::AbstractRNG = Random.default_rng(); kwargs...)
    g = loadnetwork(network, rng)
    return seeding(g, strategy, zealots, q, budgets; kwargs...)
end

function seeding(g::SimpleGraph, strategy::Symbol, zealots::Bool, q::Integer, budgets::AbstractVector{<:AbstractFloat}; kwargs...)
    strategy ∈ [:hd, :pr, :cc, :onehop, :random] || error("unknown strategy")
    if strategy == :hd
        rank = sortperm(degree(g), rev = true)
    elseif strategy == :pr
        rank = sortperm(pagerank(g), rev = true)
    elseif strategy == :cc
        rank = sortperm(competitive_complexcentrality(g, q), rev = true)
    end
    results = zeros(length(budgets))
    active = falses(nv(g))
    flexible = trues(nv(g))
    for s in eachindex(budgets)
        fill!(active, false)
        fill!(flexible, true)
        seedsize = Int(ceil(budgets[s]*nv(g)))
        if strategy == :onehop
            active[onehop(g, seedsize)] .= true
        elseif strategy == :random
            active[sample(vertices(g), seedsize, replace = false)] .= true
        else
            active[rank[1:seedsize]] .= true
        end
        if zealots 
            for i in eachindex(active)
                active[i] && (flexible[i] = false)
            end
        end
        results[s] = qvoter(g, active, q; flexible = flexible, kwargs...)
    end
    return results
end
