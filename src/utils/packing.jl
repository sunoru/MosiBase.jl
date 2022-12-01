# Support arbitrary dimensions
struct SimplePack{T}
    ms::T
    spaces::T
end

Base.eltype(::Type{SimplePack{T}}) where T = T
Base.IteratorSize(::Type{<:SimplePack}) = Base.IsInfinite()

function Base.iterate(it::SimplePack, state=zeros(Int, length(it.ms)))
    isnothing(state) && return nothing
    ndim = length(state)
    v = @inbounds [
        (state[i] - (it.ms[i] - 1) / 2) * it.spaces[i]
        for i in 1:ndim
    ]
    i = 1
    @inbounds while i ≤ ndim
        state[i] += 1
        state[i] < it.ms[i] && return v, state
        state[i] = 0
        i += 1
    end
    return v, state
end

function simple_pack(box::T, N::Integer) where T
    ndim = length(box)
    d = (prod(box) / N) ^ (1 / ndim)
    ms = ceil.(box / d)
    spaces = box ./ ms
    Iterators.take(SimplePack(ms, spaces), N)
end

function ceil_root(n::Real, k::Integer)
    m = n ^ (1 / k)
    p = round(m)
    if p ^ k < n
        p + 1
    else
        p
    end
end
# Actually not requiring the box to be a cube.
function fcc_pack(box::Vector3, N::Integer)
    m = ceil_root(N / 4, 3)
    A, B, C = box
    sA, sB, sC = box / m
    Iterators.take(Iterators.flatten((
        Vector3(i * sA - A / 2, j * sB - B / 2, k * sC - C / 2),
        Vector3((i + 0.5) * sA - A / 2, (j + 0.5) * sB - B / 2, k * sC - C / 2),
        Vector3(i * sA - A / 2, (j + 0.5) * sB - B / 2, (k + 0.5) * sC - C / 2),
        Vector3((i + 0.5) * sA - A / 2, j * sB - B / 2, (k + 0.5) * sC - C / 2)
    ) for i in 0:m-1 for j in 0:m-1 for k in 0:m-1), N)
end

fcc_pack(L::Real, N::Integer) = fcc_pack(Vector3(L, L, L), N)

