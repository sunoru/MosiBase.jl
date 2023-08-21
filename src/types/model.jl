abstract type MosiModel{T <: MosiVector} end

name(::T) where {T <: MosiModel} = string(T)
vectype(::Type{<:MosiModel{T}}) where {T} = T
vectype(::MosiModel{T}) where {T} = T
is_3d(model::Union{MosiModel, Type{<:MosiModel}}) = is_3d(vectype(model))
natoms(model::MosiModel) = model.N
constraints(::MosiModel, rs) = Float64[]
constraint_gradients(::MosiModel{T}, rs) where {T} = Vector{T}[]
constraint_gradients(::MosiModel{T}, rs, i) where {T} = T[]

interaction_pairs(model::MosiModel) =
    let N = natoms(model)
        ((i, j) for i in 1:N-1 for j in i+1:N)
    end
# Must guarantee that i < j ≤ length(rs)
potential_energy_pair(::MosiModel, rs, i, j) = error("Unimplemented")
force_pair(::MosiModel, rs, i, j) = error("Unimplemented")

potential_energy_function(model::MosiModel, rs; neighbor_list = interaction_pairs(model)) =
    sum(potential_energy_pair(model, rs, i, j) for (i, j) in neighbor_list)
function force_function(::MosiModel, rs; inplace = similar(rs), neighbor_list = interaction_pairs(model))
    for (i, j) in neighbor_list
        fij = force_pair(model, rs, i, j)
        inplace[i] += fij
        inplace[j] -= fij
    end
    inplace
end
force_function(::MosiModel, rs, i) = error("Unimplemented")
function potential_energy_gradients(model::MosiModel, rs; inplace = similar(rs))
    force_function(model, rs; inplace)
    inplace .= -inplace
end
potential_energy_gradients(model::MosiModel, rs, i) = -force_function(model, rs, i)
mass(::MosiModel, i) = 1.0

pbc_box(model::MosiModel) = model.box
has_pbc(model::MosiModel{T}) where {T} = pbc_box(model) ≢ zero(T)
distance_function(model_or_system::Union{MosiModel, MosiSystem}) =
    if has_pbc(model_or_system)
        pbc_distance(pbc_box(model_or_system))
    else
        default_distance
    end

generate_initial(::MosiModel, ST::Type{<:MosiSystem}; rng::AbstractRNG = GLOBAL_RNG) = error("Unimplemented")
# a default implementation for MolecularSystem
function generate_initial(
    model::MosiModel,
    ::Type{MolecularSystem} = MolecularSystem;
    rng::AbstractRNG = GLOBAL_RNG,
)
    conf = generate_initial(model, ConfigurationSystem; rng)
    box = pbc_box(model)
    rs = positions(conf)
    vs = randn(rng, Vector3, length(rs))
    # Remove center-of-mass momentum
    vs_cm = mean(vs)
    vs = Vector3[v - vs_cm for v in vs]
    MolecularSystem(rs, vs, box)
end

struct UnknownModel <: MosiModel{MosiVector} end
name(::UnknownModel) = "UnknownModel"
