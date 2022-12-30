abstract type MosiModel{T <: MosiVector} end

name(::T) where T <: MosiModel = string(T)
vectype(::Type{<:MosiModel{T}}) where T = T
vectype(::MosiModel{T}) where T = T
is_3d(model::Union{MosiModel, Type{<:MosiModel}}) = is_3d(vectype(model))
natoms(model::MosiModel) = model.N
constraints(::MosiModel, rs) = Float64[]
constraint_gradients(::MosiModel{T}, rs) where T = Vector{T}[]
constraint_gradients(::MosiModel{T}, rs, i) where T = T[]

potential_energy_function(::MosiModel, rs) = error("Unimplemented")
force_function(::MosiModel, rs; inplace=similar(rs)) = (error("Unimplemented"); inplace)
force_function(::MosiModel, rs, i) = error("Unimplemented")
function potential_energy_gradients(model::MosiModel, rs; inplace=similar(rs))
    force_function(model, rs; inplace)
    inplace .= -inplace
end
potential_energy_gradients(model::MosiModel, rs, i) = -force_function(model, rs, i)
mass(::MosiModel, i) = 1.0

pbc_box(model::MosiModel) = model.box
has_pbc(model::MosiModel{T}) where T = pbc_box(model) ≢ zero(T)
distance_function(model_or_system::Union{MosiModel, MosiSystem}) = if has_pbc(model_or_system)
    pbc_distance(pbc_box(model_or_system))
else
    default_distance
end

generate_initial(::MosiModel, ST::Type{<:MosiSystem}=MolecularSystem; rng::AbstractRNG=GLOBAL_RNG) = error("Unimplemented")

struct UnknownModel <: MosiModel{MosiVector} end
name(::UnknownModel) = "UnknownModel"
