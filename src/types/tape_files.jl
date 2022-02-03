struct TapeFiles{T1 <: Nullable{IOStream}, T2 <: Nullable{IOStream}}
    ts::IOStream
    rs::IOStream
    vs::T1
    ps::T2
end

function TapeFiles(
    ts_file,
    rs_file,
    vs_file = nothing,
    ps_file = nothing
)
    pos = position(rs_file)
    @assert isnothing(vs_file) || pos ≡ position(vs_file)
    @assert isnothing(ps_file) || pos ≡ position(ps_file)
    T1 = isnothing(vs_file) ? Nothing : IOStream
    T2 = isnothing(ps_file) ? Nothing : IOStream
    TapeFiles{T1, T2}(ts_file, rs_file, vs_file, ps_file)
end

function TapeFiles(
    ts_filename::AbstractString,
    rs_filename::AbstractString,
    vs_filename::Nullable{AbstractString} = nothing,
    ps_filename::Nullable{AbstractString} = nothing;
    read = true
)
    files = [
        if isnothing(filename) || read && !isfile(filename)
            nothing
        else
            open(filename, read ? "r" : "a+")
        end for filename in (
            ts_filename, rs_filename,
            vs_filename, ps_filename
        )
    ]
    TapeFiles(files...)
end

TapeFiles(
    datadir::AbstractString;
    filename_prefix = "tape",
    has_vs = true,
    has_pbc = false,
    read = true
) = TapeFiles(
    joinpath(datadir, "$filename_prefix-ts.dat"),
    joinpath(datadir, "$filename_prefix-rs.dat"),
    has_vs ? joinpath(datadir, "$filename_prefix-vs.dat") : nothing,
    has_pbc ? joinpath(datadir, "$filename_prefix-ps.dat") : nothing;
    read = read
)

has_vs(tf::TapeFiles) = !isnothing(tf.vs)
has_pbc(tf::TapeFiles) = !isnothing(tf.ps)

function read_tape(tf::TapeFiles, ::Type{T}, N::Integer) where T <: MosiVector
    times = read_vector_all(tf.ts, Float64)
    rs = read_vectors_all(tf.rs, T, N)
    vs = has_vs(tf) ? read_vectors_all(tf.vs, T, N) : Vector{T}[]
    periods = has_pbc(tf) ? read_vector3s_all(tf.ps, T, N) : Vector{T}[]
    SimulationTape(times, rs, vs, periods)
end

function update_tape(tf::TapeFiles, s::SimulationState)
    write(tf.ts, time(s))
    write(tf.rs, positions(s))
    if has_vs(tf)
        write(tf.vs, velocities(s))
    end
    if has_pbc(tf)
        write(tf.ps, periods(s))
    end
end
all_files(tf::TapeFiles) = let a = has_vs(tf), b = has_pbc(tf)
    a ? b ?
        (tf.ts, tf.rs, tf.vs, tf.ps) :
        (tf.ts, tf.rs, tf.vs) :
    b ?
        (tf.ts, tf.rs, tf.ps) :
        (tf.ts, tf.rs)
end
Base.close(tf::TapeFiles) = close.(all_files(tf))
Base.flush(tf::TapeFiles) = flush.(all_files(tf))
Base.seekend(tf::TapeFiles) = seekend.(all_files(tf))
Base.seekstart(tf::TapeFiles) = seekstart.(all_files(tf))

struct MultiFileMemoryMapTape{
    TM <: MosiModel, T1 <: Nullable{IOStream}, T2 <: Nullable{IOStream}
} <: SimulationTape{TM}
    tape_files::TapeFiles{T1, T2}
    len::Int
    model::TM
end

MultiFileMemoryMapTape(tape_files::TapeFiles{T1, T2}, len, model; is3d=true) where {T1, T2} =
    MultiFileMemoryMapTape{is3d ? Vector3 : Vector2, T1, T2}(tape_files, len, model)

function MultiFileMemoryMapTape(
    datadir::AbstractString, model;
    filename_prefix = "tape",
    has_vs = true,
    has_pbc = false
)
    tape_files = TapeFiles(datadir; filename_prefix = filename_prefix, has_vs = has_vs, has_pbc = has_pbc)
    len = filesize(joinpath(datadir, "$filename_prefix-ts.dat")) ÷ sizeof(Float64)
    MultiFileMemoryMapTape(tape_files, len, model)
end

function fetch_data!(file::IOStream, ::Type{T}, len::Integer, i::Integer) where T
    record_size = sizeof(T) * len
    offset = record_size * (i - 1)
    seek(file, offset)
    vec = Vector{T}(undef, len)
    read!(file, vec)
    vec
end

has_pbc(tape::MultiFileMemoryMapTape) = has_pbc(tape.tape_files)
Base.length(tape::MultiFileMemoryMapTape) = tape.len
natoms(tape::MultiFileMemoryMapTape) = natoms(tape.model)
Base.time(tape::MultiFileMemoryMapTape, i) = fetch_data!(tape.tape_files.ts, Float64, 1, i)[1]
times(tape::MultiFileMemoryMapTape) = read_vector_all(tape.tape_files.ts, Float64)
positions(tape::MultiFileMemoryMapTape{TM}, i) where TM = fetch_data!(
    tape.tape_files.rs, vectype(TM),
    tape.N, i
)
positions(tape::MultiFileMemoryMapTape{TM}) where TM = read_vectors_all(tape.tape_files.rs, vectype(TM), tape.N)
velocities(tape::MultiFileMemoryMapTape{TM}, i) where TM = fetch_data!(
    tape.tape_files.vs, vectype(TM),
    tape.N, i
)
velocities(tape::MultiFileMemoryMapTape{TM}) where TM = if has_vs(tape)
    read_vectors_all(tape.tape_files.vs, vectype(TM), tape.N)
else
    Vector{vectype(TM)}[]
end
periods(tape::MultiFileMemoryMapTape{TM}, i) where TM = if has_pbc(tape)
    fetch_data!(
        tape.tape_files.ps, vectype(TM),
        tape.N, i
    )
else
    T[]
end
periods(tape::MultiFileMemoryMapTape{TM}) where TM = if has_pbc(tape)
    read_vectors_all(tape.tape_files.ps, vectype(TM), tape.N)
else
    Vector{vectype(TM)}[]
end
configuration(tape::MultiFileMemoryMapTape, i) = ConfigurationSystem(
    positions(tape, i),
    periods(tape, i),
    pbc_box(tape.model)
)

function get_configuration_func(
    tape::MultiFileMemoryMapTape;
    atom_range = 1:natoms(tape),
    buffer_size = 256
)
    model = tape.model
    T = vectype(model)
    pbc = has_pbc(tape)
    box = pbc_box(model)
    N = natoms(tape)
    a, b = first(atom_range), last(atom_range)
    @assert a > 0 && b ≤ N
    s1 = sizeof(T)
    rsize = s1 * N
    a, b = (a - 1) * s1, b * s1
    buffer_len = buffer_size * rsize
    buffer_rs = zeros(UInt8, buffer_len)
    buffer_ps = pbc ? zeros(UInt8, buffer_len) : nothing
    range = 1:buffer_size
    f_rs = tape.tape_files.rs
    f_ps = tape.tape_files.ps
    seekstart(f_rs)
    readbytes!(f_rs, buffer_rs; all=false)
    if pbc
        seekstart(f_ps)
        readbytes!(f_ps, buffer_ps; all=false)
    end
    (i::Int) -> begin
        if i ∉ range
            t = (i - 1) ÷ buffer_size
            range = t * buffer_size + 1 : (t + 1) * buffer_size
            seek(f_rs, rsize * t * buffer_size)
            readbytes!(f_rs, buffer_rs; all=false)
            if pbc
                seek(f_ps, rsize * t * buffer_size)
                readbytes!(f_ps, buffer_ps; all=false)
            end
        end
        j = (i - 1) % buffer_size
        i_start = j * rsize
        rs = reinterpret(T, view(buffer_rs, i_start + 1 + a:i_start + b))
        if !pbc
            ConfigurationSystem(rs)
        else
            ps = reinterpret(T, view(buffer_ps, i_start + 1 + a:i_start + b))
            ConfigurationSystem(rs, ps, box)
        end
    end
end
