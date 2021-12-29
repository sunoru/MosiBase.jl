function orthogonalize(x, bases)
    F = 1.0
    x_new = copy(x)
    @inbounds for j = 1:length(bases)
        cji = bases[j] ⋅ x
        F -= cji ^ 2
        x_new -= cji * bases[j]
    end
    x_new /= √F
    x_new
end

@doc raw"""All Constraints. Fixed lengths, potential energy, etc.

```math
P = 1 - \sum_{i}^{N_C}|C_i'\rangle\langle C_i'|
```
"""
function projection_matrix(∇Cs)
    N_C = length(∇Cs)
    ∇C_unit_new = normalize.(∇Cs)
    P = I
    @inbounds for i = 1:N_C
        ∇C_unit_new[i] = orthogonalize(∇C_unit_new[i], ∇C_unit_new[1:i-1])
        tmp = flatten(∇C_unit_new[i])
        P -= tmp * tmp'
    end
    P
end
