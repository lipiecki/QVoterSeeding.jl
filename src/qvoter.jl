"""
    qvoter(g::SimpleGraph{<:Integer}, active::BitVector, q::Integer; ε::AbstractFloat = 0.0, flexible::BitVector = trues(nv(g)), replace::Bool = false, warmup::Integer = 10^5, averaging::Integer = 10^3)
Run a Monte Carlo simulation of the [q-voter model](https://doi.org/10.1103/PhysRevE.80.041129).
### Arguments
- `g`: Graph representing the structure of the agent system.
- `active`: BitVector representing states of the agents (the vector is mutated during simulation).
- `q`: Number of neighbors in the influence group.
### Keyword arguments
- `ε = 0.0`: Probability of flipping state when q-panel is not unanimous.
- `flexible = trues(nv(g))`: BitVector represeting flexibility of the agents, i.e., whether an agent changes its state during simulation.
- `replace = false`: If true, cosntruct the q-panel by drawing neighbors with replacement, otherwise draw without replacement.
- `warmup = 10^5`: Number of Monte Carlo steps calculated during the simulation warm up.
- `averaging = 10^3`: Number of Monte Carlo steps over which the concentration is averaged.

Return the final concentration of active agents. If the system reaches an absorbing state, the final concentration is returned immediately.
"""
function qvoter(g::SimpleGraph{<:Integer}, active::BitVector, q::Integer; ε::AbstractFloat = 0.0, flexible::BitVector = trues(nv(g)), replace::Bool = false, warmup::Integer = 10^5, averaging::Integer = 10^3)
    c̄ = 0.0
    c = sum(active)
    is_noise = ε > 0.0
    flexible_vertices = filter(v -> flexible[v], vertices(g))
    qpanel = Vector{Int}(undef, q)
    for t in 1:warmup+averaging
        for _ in 1:nv(g)
            target = rand(flexible_vertices)
            if !replace && degree(g, target) <= q
                if active[target]
                    if !any(@view active[neighbors(g, target)]) || (is_noise && rand() < ε)
                        active[target] = false
                        c -= 1
                    end
                elseif all(@view active[neighbors(g, target)]) || (is_noise && rand() < ε)
                    active[target] = true
                    c += 1
                end
            else
                sample!(neighbors(g, target), qpanel, replace = replace)
                if active[target]
                    if !any(@view active[qpanel]) || (is_noise && rand() < ε)
                        active[target] = false
                        c -= 1
                    end
                elseif all(@view active[qpanel]) || (is_noise && rand() < ε)
                    active[target] = true
                    c += 1
                end
            end
        end
        if c == nv(g) || c == 0
            return c/nv(g)
        end
        if t > warmup
            c̄ += c/nv(g)
        end
    end
    return c̄/averaging
end
