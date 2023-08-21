const StateOrSystem = Union{SimulationState, MosiSystem}

"""
    NeighborList(state_or_system, cutoff; [initial_capacity=5N])

Create a neighbor list for the given system or state.

You can specify and the initial capacity of the neighbor list, which is set to `5N` by default,
where N is the number of particles in the system.
"""
struct NeighborList
    cutoff²::Float64
    N::Int
    end_indices::Vector{Int}
    neighbors::Vector{Int}
end

function NeighborList(s::StateOrSystem, cutoff::Real; capacity::Integer = 5 * natoms(s))
    if s isa SimulationState
        s = system(s)
    end
    N = natoms(s)
    end_indices = zeros(Int, N - 1)
    neighbors = Vector{Int}(undef, capacity)
    nl = NeighborList(Float64(cutoff^2), N, end_indices, neighbors)
    update!(nl, s)
    nl
end

"""
    update!(neighbor_list, state_or_system)

Update the neighbor list for the given system or state.

!!! note
    Users should guarantee that the number of particles in the system is not changed.
"""
update!(nl::NeighborList, s::StateOrSystem) = update!(nl, positions(s), distance_function(s))
function update!(nl::NeighborList, rs::AbstractVector, dist::F) where {F}
    N = nl.N
    cutoff² = nl.cutoff²
    end_indices = nl.end_indices
    neighbors = nl.neighbors
    n = length(neighbors)
    index = 0
    # @assert N == length(rs)
    for i in 1:N-1
        @inbounds rᵢ = rs[i]
        for j in i+1:N
            @inbounds rⱼ = rs[j]
            d = dist(rᵢ, rⱼ)
            norm_sqr(d) > cutoff² && continue
            index += 1
            if n < index
                n = Base._nextpow2(index)
                resize!(neighbors, n)
            end
            @inbounds neighbors[index] = j
        end
        @inbounds end_indices[i] = index
    end
    nl
end

Base.length(nl::NeighborList) = nl.end_indices[end]
Base.eltype(::NeighborList) = Pair{Int, Int}
function Base.iterate(nl::NeighborList, state::Pair{Int, Int} = (1, 0))
    (i, ii) = state
    N = nl.N
    ii += 1
    # only need to loop until N - 1
    while i < N
        end_index = @inbounds nl.end_indices[i]
        ii ≤ end_index && return ((i, @inbounds nl.neighbors[ii]), (i, ii))
        i += 1
    end
    nothing
end

atom_pairs(nl::NeighborList) = nl
