import Base: copy

import LinearAlgebra: norm, normalize
import Statistics: mean

abstract type System{T <: MosiVector} end
abstract type Molecule{T} <: System{T} end

natoms(s::System) = length(positions(s))
positions(::System) = error("Unimplemented")
velocities(::System{T}) where T = T[]
periods(::System{T}) where T = T[]
box(::System{T}) where T = zero(T)
has_pbc(s::System{T}) where T = box(s) ≢ zero(T)
original_positions(s::System) = if has_pbc(s)
    original_vectors(
        positions(s), periods(s), box(s)
    )
else
    positions(s)
end
distance_function(s::System) = if has_pbc(s)
    pbc_distance(box(s))
else
    default_distance
end
update_periods!(s::System) = if has_pbc(s)
    update_periods!(positions(s), periods(s), box(s))
    s
else
    nothing
end

struct MolecularSystem{T <: MosiVector, AT <: AbstractVector{T}} <: System{T}
    positions::AT
    velocities::AT
    periods::AT
    box::T
end

function MolecularSystem(
    rs, vs, box = nothing
)
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
box(s::MolecularSystem) = s.box

struct ConfigurationSystem{T <: MosiVector, AT <: AbstractVector{T}} <: System{T}
    positions::AT
    periods::AT
    box::T
end

function ConfigurationSystem(rs; box = nothing)
    T = eltype(rs)
    @assert box ≡ nothing || box ≡ zero(T)
    ConfigurationSystem(rs, [], zero(T))
end
ConfigurationSystem((rs, ps)::Tuple; box) = update_periods!(ConfigurationSystem(rs, ps, box))
ConfigurationSystem(s::System; box = box(s)) = ConfigurationSystem(positions(s), periods(s), box)
copy(s::ConfigurationSystem) = ConfigurationSystem(copy(positions(s)), copy(periods(s)), box(s))

positions(s::ConfigurationSystem) = s.positions
velocities(::ConfigurationSystem{T}) where T = T[]
periods(s::ConfigurationSystem) = s.periods
box(s::ConfigurationSystem) = s.box
