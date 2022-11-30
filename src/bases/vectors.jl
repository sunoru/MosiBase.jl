using Base: ReinterpretArray

import StaticArrays
using StaticArrays: SVector

const Vector2 = SVector{2, Float64}
const Vector3 = SVector{3, Float64}
const MosiVector = Union{Vector2, Vector3}

const Vector2s = Vector{Vector2}
const Vector3s = Vector{Vector3}
const AbstractVector2s = AbstractVector{Vector2}
const AbstractVector3s = AbstractVector{Vector3}

is_3d(x) = false
is_3d(::Type{Vector3}) = true
is_3d(::Type{<:AbstractVector3s}) = true

flatten(x::Vector{T}) where T <: MosiVector = reinterpret(Float64, x)
unflatten(::Type{T}, x::Vector{Float64}) where T <: MosiVector = reinterpret(T, x)
unflatten(x::ReinterpretArray{Float64, 1, T, Vector{T}}) where T <: MosiVector = x.parent

to_matrix(v::Vector2s) = reshape(flatten(v), 2, length(v))
to_matrix(v::Vector3s) = reshape(flatten(v), 3, length(v))
