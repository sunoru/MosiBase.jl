"""
    potential_energy(system, model)

Compute the potential energy of the system with the given model.
"""
potential_energy(s::MosiSystem, model::MosiModel) = potential_energy_function(model, positions(s))
"""
    kinetic_energy(system, model)

Compute the kinetic energy of the system with the given model.
"""
kinetic_energy(s::MosiSystem, model::MosiModel) = error("Unimplemented")
kinetic_energy(::ConfigurationSystem, ::MosiModel) = 0.0
mechanical_energy(s::MosiSystem, model::MosiModel) = potential_energy(s, model) + kinetic_energy(s, model)

degree_of_freedom(s::MosiSystem, model::MosiModel) = 0
_temperature(K::Float64, dof::Int) = K * 2 / dof
function temperature(s::MosiSystem, model::MosiModel)
    dof = degree_of_freedom(s, model)
    T = _temperature(kinetic_energy(s, model), dof)
end
