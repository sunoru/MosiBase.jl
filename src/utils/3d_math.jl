function rotate_around(v::Vector3, k::Vector3, φ::Real)
    k = normalize(k)
    v_parallel = (v ⋅ k) * k
    v_perpendicular = v - v_parallel
    ω = k × v_perpendicular
    v_parallel + v_perpendicular * cos(φ) + ω * sin(φ)
end

@_with_default_rng function random_3d_direction()
    θ = 2π * rand(rng)
    z = 2 * rand(rng) - 1
    r = √(1 - z ^ 2)
    Vector3(r * cos(θ), r * sin(θ), z)
end

"""
    random_direction_on_plane([rng], normal_vector)
    random_direction_on_plane([rng], u, v)

Generate a random direction on a plane defined by a normal vector or two vectors on the plane.
"""
random_direction_on_plane(rng::AbstractRNG, u::Vector3, v::Vector3) = let φ = 2π * rand(rng)
    w = cos(φ) * u + sin(φ) * v
end
@_with_default_rng function random_direction_on_plane(normal_vector::Vector3)
    a, b, c = normal_vector
    u = normalize(Vector3(b - c, -a + c, a - b))
    v = normalize(normal_vector × u)
    random_direction_on_plane(rng, u, v)
end
