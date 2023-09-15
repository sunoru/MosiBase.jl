function orthogonalize(x, bases::AbstractVector)
    F = 1.0
    x_new = copy(x)
    @inbounds for j in 1:length(bases)
        cji = bases[j] ⋅ x
        F -= cji^2
        x_new -= cji * bases[j]
    end
    x_new /= √F
    x_new
end

# N_vector × N_bases
function orthogonalize(x, bases::AbstractMatrix)
    F = 1.0
    x_new = copy(x)
    N_bases = size(bases, 2)
    @inbounds for j in 1:N_bases
        base = @view(bases[:, j])
        cji = base ⋅ x
        F -= cji^2
        x_new -= cji * base
    end
    x_new /= √F
    x_new
end

function normalize_columns!(A::AbstractMatrix)
    N = size(A, 2)
    @inbounds for i in 1:N
        normalize!(@view(A[:, i]))
    end
    A
end
normalize_columns(A::AbstractMatrix) = normalize_columns!(copy(A))

@doc raw"""All Constraints. Fixed lengths, potential energy, etc.

```math
P = 1 - \sum_{i}^{N_C}|C_i'\rangle\langle C_i'|
```
"""
function projection_matrix(∇Cs::AbstractVector; normed = false)
    N_C = length(∇Cs)
    ∇C_unit_new = if normed
        copy(∇Cs)
    else
        normalize.(∇Cs)
    end
    P = I
    @inbounds for i in 1:N_C
        tmp = ∇C_unit_new[i] = orthogonalize(∇C_unit_new[i], ∇C_unit_new[1:i-1])
        P -= tmp * tmp'
    end
    P
end

function projection_matrix(∇Cs::AbstractMatrix; normed = false)
    N_C = size(∇Cs, 2)
    ∇C_unit_new = if normed
        copy(∇Cs)
    else
        normalize_columns(∇Cs)
    end
    P = I
    @inbounds for i in 1:N_C
        tmp = ∇C_unit_new[:, i] = orthogonalize(@view(∇C_unit_new[:, i]), @view(∇C_unit_new[:, 1:i-1]))
        P -= tmp * tmp'
    end
    P
end