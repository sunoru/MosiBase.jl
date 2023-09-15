using Base: ReinterpretArray

using StaticArrays: SVector

const Vector2 = SVector{2, Float64}
const Vector3 = SVector{3, Float64}
const MosiVector = Union{Vector2, Vector3}

const Matrix2 = SMatrix{2, 2, Float64, 4}
const Matrix3 = SMatrix{3, 3, Float64, 9}
const MosiMatrix = Union{Matrix2, Matrix3}

const Vector2s = Vector{Vector2}
const Vector3s = Vector{Vector3}
const AbstractVector2s = AbstractVector{Vector2}
const AbstractVector3s = AbstractVector{Vector3}

is_3d(x) = false
is_3d(::Type{Vector3}) = true
is_3d(::Type{<:AbstractVector3s}) = true

flatten(x::AbstractVector{<:Real}) = x
flatten(x::AbstractVector{<:MosiVector}) = reinterpret(Float64, x)
unflatten(::Type{T}, x::Vector{Float64}) where {T <: MosiVector} = reinterpret(T, x)
unflatten(x::ReinterpretArray{Float64, 1, T, Vector{T}}) where {T <: MosiVector} = x.parent

to_matrix(v::Vector2s) = reshape(flatten(v), 2, length(v))
to_matrix(v::Vector3s) = reshape(flatten(v), 3, length(v))

zero_similar(a::AbstractArray) = fill!(similar(a), zero(eltype(a)))

# Not a good practice but we need these
Base.zero(::Type{Adjoint{Float64, T}}) where T <: MosiVector = adjoint(zero(T))
# https://github.com/JuliaSparse/SparseArrays.jl/pull/433
function Base.convert(::Type{T}, x::Int) where T <: MosiVector
    @assert x == 0
    zero(T)
end
Base.:(-)(b::Bool, A::MosiMatrix) = UniformScaling(b) - A
