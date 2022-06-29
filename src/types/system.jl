import Base: copy

import LinearAlgebra: norm, normalize
import Statistics: mean

abstract type MosiSystem{T <: MosiVector} end
abstract type Molecule{T} <: MosiSystem{T} end

natoms(s::MosiSystem) = length(positions(s))
positions(::MosiSystem) = error("Unimplemented")
velocities(::MosiSystem{T}) where {T} = T[]
periods(::MosiSystem{T}) where {T} = T[]
pbc_box(::MosiSystem{T}) where {T} = zero(T)
has_pbc(s::MosiSystem{T}) where {T} = pbc_box(s) ≢ zero(T)
original_positions(s::MosiSystem) =
    if has_pbc(s)
        original_vectors(positions(s), periods(s), pbc_box(s))
    else
        positions(s)
    end
distance_function(s::MosiSystem) =
    if has_pbc(s)
        pbc_distance(pbc_box(s))
    else
        default_distance
    end
update_periods!(s::MosiSystem) =
    if has_pbc(s)
        update_periods!(positions(s), periods(s), pbc_box(s))
        s
    else
        s
    end

struct MolecularSystem{T <: MosiVector, AT <: AbstractVector{T}} <: MosiSystem{T}
    positions::AT
    velocities::AT
    periods::AT
    box::T
end

function MolecularSystem(rs, vs, box = nothing)
    N = length(rs)
    @assert N === length(vs)
    T = eltype(rs)
    periods = if isnothing(box)
        box = zero(T)
        T[]
    else
        periods = zeros(T, N)
        update_periods!(rs, periods, box)
        periods
    end
    ms = MolecularSystem(rs, vs, periods, box)
end

positions(s::MolecularSystem) = s.positions
velocities(s::MolecularSystem) = s.velocities
periods(s::MolecularSystem) = s.periods
pbc_box(s::MolecularSystem) = s.box

struct ConfigurationSystem{T <: MosiVector, AT1 <: AbstractVector{T}, AT2 <: AbstractVector{T}} <:
       MosiSystem{T}
    positions::AT1
    periods::AT2
    box::T
end

function ConfigurationSystem(rs; box = nothing, update_periods = true)
    T = eltype(rs)
    ps, box = if box ≡ nothing || box ≡ zero(T)
        T[], zero(T)
    else
        zeros(T, length(rs)), box
    end
    conf = ConfigurationSystem(rs, ps, box)
    if update_periods
        update_periods!(conf)
    end
    conf
end
function ConfigurationSystem((rs, ps)::Tuple; box, update_periods = true)
    conf = ConfigurationSystem(rs, ps, box)
    if update_periods
        update_periods!(conf)
    end
    conf
end
ConfigurationSystem(s::MosiSystem; box = pbc_box(s)) = ConfigurationSystem(positions(s), periods(s), box)
copy(s::ConfigurationSystem) = ConfigurationSystem(copy(positions(s)), copy(periods(s)), pbc_box(s))

positions(s::ConfigurationSystem) = s.positions
velocities(::ConfigurationSystem{T}) where {T} = T[]
periods(s::ConfigurationSystem) = s.periods
pbc_box(s::ConfigurationSystem) = s.box

center_of_mass(s::MosiSystem, args...) = center_of_mass(original_positions(s), args...)
center_of_mass(rs::AbstractVector{<:MosiVector}) = mean(rs)
function center_of_mass(rs::AbstractVector{<:MosiVector}, mass_func::Function)
    s = 0.0
    M = 0.0
    @inbounds @simd for i in 1:length(rs)
        m = mass_func(i)
        s += rs[i] * m
        M += m
    end
    s / M
end
