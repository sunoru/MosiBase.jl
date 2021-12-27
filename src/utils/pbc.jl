# Return the vector in the box and the number of periods moved.
function period_check(vec::T, box::T) where T <: MosiVector
    t = round.(vec ./ box)
    vec -= t .* box
    vec, t
end

function update_periods!(vecs::AbstractVector{T}, periods::AbstractVector{T}, box::T) where T <: MosiVector
    N = length(vecs)
    @inbounds @simd for i in 1:N
        vec, t = period_check(vecs[i], box)
        vecs[i] = vec
        periods[i] += t
    end
    vecs, periods
end

default_distance(a::T, b::T) where T <: MosiVector = a - b
pbc_distance(pbox::T) where T <: MosiVector = (a::T, b::T) -> period_check(a - b, pbox)[1]

original_vector(vec::T, period::T, box::T) where T <: MosiVector = vec + period .* box

function original_vectors(vecs::AbstractVector{T}, periods::AbstractVector{T}, box::T) where T <: MosiVector
    N = length(vecs)
    if length(periods) ≡ 0
        return copy(vecs)
    end
    @inbounds [
        original_vector(vecs[i], periods[i], box)
        for i = 1:N
    ]
end