abstract type SimulationState end

system(::SimulationState) = error("Unimplemented")
Base.time(::SimulationState) = error("Unimplemented")
positions(state::SimulationState) = positions(system(state))
velocities(state::SimulationState) = velocities(system(state))
periods(state::SimulationState) = periods(system(state))

abstract type SimulationTape{T <: MosiVector} end

struct SimpleTape{T} <: SimulationTape{T}
    times::Vector{Float64}
    positions::Vector{Vector{T}}
    velocities::Vector{Vector{T}}
    periods::Vector{Vector{T}}
end
SimulationTape(ts, rs) = SimpleTape(ts, rs, [], [])
SimulationTape(ts, rs, vs, ps) = SimpleTape(ts, rs, vs, ps)
natoms(tape::SimulationTape) = length(tape.positions[1])
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
    tape::SimulationTape,
    model::MosiModel;
    atom_range = 1:natoms(tape)
)
    rss = positions(tape)
    pss = periods(tape)
    tbox = box(model)
    (i::Int) -> ConfigurationSystem(
        view(rss[i], atom_range),
        view(pss[i], atom_range),
        tbox
    )
end

abstract type SimulationSetup end
model(::SimulationSetup) = UnknownModel()
init_state(::SimulationSetup) = error("Unimplemented")

abstract type SimulationResult end
tape(result::SimulationResult) = result.tape
observables(result::SimulationResult) = result.observables

struct SimulationError <: Exception
    msg::String
end

Base.showerror(io::IO, e::SimulationError) = print(io, e.msg)


