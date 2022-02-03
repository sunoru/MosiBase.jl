function read_vector_all(io::IOStream, ::Type{T}) where T
    t = position(io)
    seekend(io)
    len = position(io) ÷ sizeof(T)
    seekstart(io)
    vec = read!(io, Vector{T}(undef, len))
    seek(io, t)
    vec
end

read_vectors(io::IOStream, ::Type{T}, N::Integer) where T <: MosiVector =
    read!(io, Vector{T}(undef, N))

function read_vectors_all(io::IOStream, ::Type{T}, N::Integer) where T <: MosiVector
    seekend(io)
    len = position(io) ÷ (sizeof(T) * N)
    seekstart(io)
    [read_vectors(io, T, N) for _ in 1:len]
end
