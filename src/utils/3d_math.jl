using LinearAlgebra: ×, norm, norm_sqr, normalize
using Random: AbstractRNG, GLOBAL_RNG

function spherical_to_cartesian(ρ, θ, ϕ)
    x = ρ * sin(θ) * cos(ϕ)
    y = ρ * sin(θ) * sin(ϕ)
    z = ρ * cos(θ)
    Vector3(x, y, z)
end
spherical_to_cartesian(spherical::Vector3) = spherical_to_cartesian(spherical[1], spherical[2], spherical[3])

function cartesian_to_spherical(x, y, z)
    ρ = √(x ^ 2 + y ^ 2 + z ^ 2)
    θ = acos(z / ρ)
    ϕ = atan(y, x)
    Vector3(ρ, θ, ϕ)
end
cartesian_to_spherical(xyz::Vector3) = cartesian_to_spherical(xyz[1], xyz[2], xyz[3])

function random_3d_direction(rng::AbstractRNG = GLOBAL_RNG)
    ϕ = (2 * rand(rng) - 1) * π
    θ = acos(1 - 2 * rand(rng))
    spherical_to_cartesian(1.0, θ, ϕ)
end

function random_direction_plane(rng::AbstractRNG, normal_vector::Vector3)
    a, b, c = normal_vector
    u = normalize(Vector3(b - c, -a + c, a - b))
    v = normalize(normal_vector × u)
    ϕ = 2 * π * rand(rng)
    w = cos(ϕ) * u + sin(ϕ) * v
end

random_direction_plane(normal_vector::Vector3) = random_direction_plane(GLOBAL_RNG, normal_vector)
