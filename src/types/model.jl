abstract type MosiModel end

name(::MosiModel) = error("Unimplemented")
natoms(model::MosiModel) = model.N
constraints(::MosiModel) = error("Unimplemented")
constraint_gradients(::MosiModel) = error("Unimplemented")

potential_energy_function(::MosiModel) = error("Unimplemented")
force_function(::MosiModel) = error("Unimplemented")
potential_energy_gradients(model::MosiModel) = (rs) -> -force_function(model)(rs)

has_pbc(model::MosiModel) = false
box(model::MosiModel) = model.box
distance_function(model_or_system::Union{MosiModel, System}) = if has_pbc(model_or_system)
    pbc_distance(box(model_or_system))
else
    default_distance
end

struct UnknownModelType <: MosiModel end
const UnknownModel = UnknownModelType()
name(::UnknownModelType) = "UnknownModel"
