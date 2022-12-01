potential_energy(s::MosiSystem, model::MosiModel) = potential_energy_function(model, positions(s))
kinetic_energy(s::MosiSystem, model::MosiModel) = error("Unimplemented")
kinetic_energy(::ConfigurationSystem, ::MosiModel) = 0.0
mechanical_energy(s::MosiSystem, model::MosiModel) = potential_energy(s, model) + kinetic_energy(s, model)
