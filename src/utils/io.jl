function read_vector_all(io::IOStream, ::Type{T}) where T
    t = position(io)
    seekend(io)
    len = position(io) ÷ sizeof(T)
    seekstart(io)
    vec = read!(io, Vector{T}(undef, len))
    seek(io, t)
    vec
end
read_vector_all(filename::AbstractString, T::Type) = open(filename) do fi
    read_vector_all(fi, T)
end

read_vectors(io::IOStream, ::Type{T}, N::Integer) where T <: MosiVector =
    read!(io, Vector{T}(undef, N))

function read_vectors_all(io::IOStream, ::Type{T}, N::Integer) where T <: MosiVector
    seekend(io)
    len = position(io) ÷ (sizeof(T) * N)
    seekstart(io)
    [read_vectors(io, T, N) for _ in 1:len]
end
read_vectors_all(filename::AbstractString, T::Type, N::Integer) = open(filename) do fi
    read_vectors_all(fi, T, N)
end
