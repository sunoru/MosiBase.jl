abstract type SimulationState end

system(::SimulationState) = error("Unimplemented")
Base.time(::SimulationState) = error("Unimplemented")
positions(state::SimulationState) = positions(system(state))
velocities(state::SimulationState) = velocities(system(state))
periods(state::SimulationState) = periods(system(state))

abstract type SimulationTape{TM <: MosiModel} end

struct SimpleTape{T, TM <: MosiModel{T}} <: SimulationTape{TM}
    times::Vector{Float64}
    positions::Vector{Vector{T}}
    velocities::Vector{Vector{T}}
    periods::Vector{Vector{T}}
    model::TM
end
SimulationTape(ts, rs) = SimpleTape(ts, rs, [], [], UnknownModel())
SimulationTape(ts, rs, vs, ps) = SimpleTape(ts, rs, vs, ps, UnknownModel())
has_vs(tape::SimulationTape) = length(tape.velocities) > 0
mosi_model(tape::SimulationTape) = tape.model
natoms(tape::SimulationTape{UnknownModel}) = length(tape.positions[1])
natoms(tape::SimulationTape) = natoms(tape.model)
Base.length(tape::SimulationTape) = length(tape.times)
times(tape::SimulationTape) = tape.times
positions(tape::SimulationTape) = tape.positions
velocities(tape::SimulationTape) = tape.velocities
periods(tape::SimulationTape) = tape.periods

times(tape::SimulationTape, i) = tape.times[i]
positions(tape::SimulationTape, i) = tape.positions[i]
velocities(tape::SimulationTape, i) = tape.velocities[i]
periods(tape::SimulationTape{T}, i) where T = length(tape.periods) ≡ 0 ? T[] : tape.periods[i]

function get_configuration_func(
    tape::SimulationTape;
    atom_range = 1:natoms(tape)
)
    rss = positions(tape)
    pss = periods(tape)
    model = mosi_model(tape)
    box = pbc_box(model)
    (i::Int) -> ConfigurationSystem(
        view(rss[i], atom_range),
        view(pss[i], atom_range),
        box
    )
end

function system(tape::SimulationTape, i::Integer)
    model = mosi_model(tape)
    has_vs(mosi_model(tape)) || return get_configuration_func(tape)(i)
    rs = positions(tape, i)
    vs = velocities(tape, i)
    ps = periods(tape, i)
    box = pbc_box(model, i)
    MolecularSystem(rs, vs, ps, box)
end

abstract type SimulationSetup end
mosi_model(::SimulationSetup) = UnknownModel()
init_state(::SimulationSetup) = error("Unimplemented")

abstract type SimulationResult end
tape(result::SimulationResult) = result.tape
observables(result::SimulationResult) = result.observables

struct SimulationError <: Exception
    msg::String
end

Base.showerror(io::IO, e::SimulationError) = print(io, e.msg)
