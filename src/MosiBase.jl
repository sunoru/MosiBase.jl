module MosiBase

using LinearAlgebra
import Base: @pure

export Nullable
include("./bases/bases.jl")

export Vector2, Vector3, MosiVector,
    Vector2s, Vector3s,
    AbstractVector2s, AbstractVector3s,
    flatten, unflatten,
    to_matrix, ndims
include("./bases/vectors.jl")

export read_vector_all,
    read_vectors, read_vectors_all
include("./utils/io.jl")

export Xoshiro256StarStar,
    make_seed, new_rng, restore_rng, rng_state
include("./utils/rng.jl")

export period_check, update_periods!,
    default_distance, pbc_distance, pbc_box,
    original_vector, original_vectors
include("./utils/pbc.jl")

export polar_to_cartesian, cartesian_to_polar,
    random_2d_direction
include("./utils/2d_math.jl")

export spherical_to_cartesian, cartesian_to_spherical,
    random_3d_direction, random_direction_plane
include("./utils/3d_math.jl")

export orthogonalize, projection_matrix
include("./utils/projection_matrix.jl")

export MosiSystem, Molecule, MolecularSystem, ConfigurationSystem,
    natoms, positions, velocities, periods, has_pbc, original_positions,
    distance_function, update_periods!
include("./types/system.jl")

export MosiModel, UnknownModel,
    name, is_3d, vectype, constraints, constraint_gradients,
    potential_energy_function, force_function, potential_energy_gradients
include("./types/model.jl")

export potential_energy
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
