import LinearAlgebra: ×, norm, norm_sqr, normalize
import Random: AbstractRNG

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

function random_2d_direction(rng::AbstractRNG)
    θ = 2 * π * rand(rng)
    polar_to_cartesian(1.0, θ)
end
