using LinearAlgebra: norm_sqr

function polar_to_cartesian(ρ, θ)
    x = ρ * cos(θ)
    y = ρ * sin(θ)
    Vector2(x, y)
end
polar_to_cartesian(polar::Vector2) = polar_to_cartesian(polar[1], polar[2])

function cartesian_to_polar(x, y)
    ρ = √(x ^ 2 + y ^ 2)
    θ = acos(x / ρ)
    Vector2(ρ, θ)
end
cartesian_to_polar(xy::Vector3) = cartesian_to_polar(xy[1], xy[2])

@_with_default_rng function random_point_on_circle(radius::Real)
    θ = 2π * rand(rng)
    polar_to_cartesian(radius, θ)
end

@_with_default_rng random_2d_direction() = random_point_on_circle(rng, 1.0)

@_with_default_rng function random_point_on_disk(radius::Real)
    ρ = √(rand(rng) * radius ^ 2)
    random_point_on_circle(rng, ρ)
end
