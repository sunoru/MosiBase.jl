import Base: copy

import LinearAlgebra: norm, normalize
import Statistics: mean

abstract type MosiSystem{T <: MosiVector} end
abstract type Molecule{T} <: MosiSystem{T} end

natoms(s::MosiSystem) = length(positions(s))
positions(::MosiSystem) = error("Unimplemented")
velocities(::MosiSystem{T}) where T = T[]
periods(::MosiSystem{T}) where T = T[]
box(::MosiSystem{T}) where T = zero(T)
has_pbc(s::MosiSystem{T}) where T = box(s) ≢ zero(T)
original_positions(s::MosiSystem) = if has_pbc(s)
    original_vectors(
        positions(s), periods(s), box(s)
    )
else
    positions(s)
end
distance_function(s::MosiSystem) = if has_pbc(s)
    pbc_distance(box(s))
else
    default_distance
end
update_periods!(s::MosiSystem) = if has_pbc(s)
    update_periods!(positions(s), periods(s), box(s))
    s
else
    nothing
end

struct MolecularSystem{T <: MosiVector, AT <: AbstractVector{T}} <: MosiSystem{T}
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

struct ConfigurationSystem{T <: MosiVector, AT <: AbstractVector{T}} <: MosiSystem{T}
    positions::AT
    periods::AT
    box::T
end

function ConfigurationSystem(rs; box = nothing)
    T = eltype(rs)
    ps, box = if box ≡ nothing || box ≡ zero(T)
        T[], zero(T)
    else
        zeros(T, length(rs)), box
    end
    update_periods!(ConfigurationSystem(rs, ps, box))
end
ConfigurationSystem((rs, ps)::Tuple; box) = update_periods!(ConfigurationSystem(rs, ps, box))
ConfigurationSystem(s::MosiSystem; box = box(s)) = ConfigurationSystem(positions(s), periods(s), box)
copy(s::ConfigurationSystem) = ConfigurationSystem(copy(positions(s)), copy(periods(s)), box(s))

positions(s::ConfigurationSystem) = s.positions
velocities(::ConfigurationSystem{T}) where T = T[]
periods(s::ConfigurationSystem) = s.periods
box(s::ConfigurationSystem) = s.box
