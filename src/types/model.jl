abstract type MosiModel{T <: MosiVector} end

name(::MosiModel) = error("Unimplemented")
is_3d(::MosiModel) = T ≡ Vector3
natoms(model::MosiModel) = model.N
constraints(::MosiModel, rs) = Float64[]
constraint_gradients(::MosiModel{T}, rs) where T = Vector{T}[]
constraint_gradients(::MosiModel{T}, rs, i) where T = T[]

potential_energy_function(::MosiModel, rs) = error("Unimplemented")
force_function(::MosiModel, rs) = error("Unimplemented")
force_function(::MosiModel, rs, i) = error("Unimplemented")
potential_energy_gradients(model::MosiModel, rs) = -force_function(model, rs)
potential_energy_gradients(model::MosiModel, rs, i) = -force_function(model, rs, i)

has_pbc(model::MosiModel) = false
pbc_box(model::MosiModel) = model.box
distance_function(model_or_system::Union{MosiModel, MosiSystem}) = if has_pbc(model_or_system)
    pbc_distance(pbc_box(model_or_system))
else
    default_distance
end

struct UnknownModel <: MosiModel{MosiVector} end
name(::UnknownModel) = "UnknownModel"
