module MosimoBase

using Reexport: @reexport
@reexport using LinearAlgebra, Statistics, StaticArrays, Random, DataFrames, JLD
using MacroTools: @capture
@reexport using Random: AbstractRNG, GLOBAL_RNG
@reexport using MLStyle

using LinearAlgebra: norm_sqr
export norm_sqr

export Nullable, ∞
include("./bases/bases.jl")

export Vector2, Vector3, MosiVector,
    Vector2s, Vector3s,
    AbstractVector2s, AbstractVector3s,
    is_3d, flatten, unflatten,
    to_matrix
include("./bases/vectors.jl")

export read_vector_all,
    read_vectors, read_vectors_all
include("./utils/io.jl")

export MosiRNG,
    make_seed, new_rng, restore_rng, rng_state
include("./utils/rng.jl")

export period_check, update_periods!,
    default_distance, pbc_distance, pbc_box,
    original_vector, original_vectors
include("./utils/pbc.jl")

export simple_pack, fcc_pack
include("./utils/packing.jl")

export polar_to_cartesian, cartesian_to_polar,
    random_2d_direction,
    random_point_on_circle, random_point_on_disk
include("./utils/2d_math.jl")

export spherical_to_cartesian, cartesian_to_spherical,
    ρ_hat, θ_hat, φ_hat,
    ∂Ω_∂ρ, ∂Ω_∂θ, ∂Ω_∂φ
include("./utils/spherical_coordinates.jl")

export random_3d_direction, random_direction_on_plane,
    rotate_around, project_onto, decompose_vector
include("./utils/3d_math.jl")

export orthogonalize, projection_matrix
include("./utils/projection_matrix.jl")

export MosiSystem, Molecule, MolecularSystem, ConfigurationSystem,
    natoms, positions, velocities, periods, has_pbc, original_positions,
    distance_function, update_periods!,
    center_of_mass
include("./types/system.jl")

export MosiModel, UnknownModel,
    name, vectype, constraints, constraint_gradients,
    potential_energy_function, force_function, potential_energy_gradients,
    mass,
    generate_initial
include("./types/model.jl")

export potential_energy, kinetic_energy, mechanical_energy
include("./types/energy.jl")

export SimulationState, SimulationTape, SimulationSetup, SimulationResult, SimulationError,
    SimpleTape,
    system, times, mosi_model, get_configuration_func, init_state, tape, observables
include("./types/simulation.jl")

export TapeFiles, MultiFileMemoryMapTape,
    has_vs, update_tape, configuration
include("./types/tape_files.jl")

export Ensemble, CanonicalEnsemble, MicrocanonicalEnsemble, PotentialEnergyLandscapeEnsemble
include("./types/ensembles.jl")

end
