function spherical_to_cartesian(ρ, θ, φ)
    x = ρ * sin(θ) * cos(φ)
    y = ρ * sin(θ) * sin(φ)
    z = ρ * cos(θ)
    Vector3(x, y, z)
end
spherical_to_cartesian(spherical::Vector3) = spherical_to_cartesian(spherical[1], spherical[2], spherical[3])

function cartesian_to_spherical(x, y, z)
    ρ = √(x ^ 2 + y ^ 2 + z ^ 2)
    θ = acos(z / ρ)
    φ = atan(y, x)
    Vector3(ρ, θ, φ)
end
cartesian_to_spherical(xyz::Vector3) = cartesian_to_spherical(xyz[1], xyz[2], xyz[3])

ρ_hat(θ, φ) = Vector3(sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ))
ρ_hat(Ω::Vector3) = normalize(Ω)
θ_hat(θ, φ) = Vector3(cos(θ) * cos(φ), cos(θ) * sin(φ), -sin(θ))
φ_hat(θ, φ) = Vector3(-sin(φ), cos(φ), 0)
∂Ω_∂ρ(θ, φ) = 0.0
∂Ω_∂θ(θ, φ) = Vector3(cos(θ) * cos(φ), cos(θ) * sin(φ), -sin(θ))
∂Ω_∂φ(θ, φ) = Vector3(-sin(θ) * sin(φ), sin(θ) * cos(φ), 0)

for func in (:θ_hat, :φ_hat, :∂Ω_∂ρ, :∂Ω_∂θ, :∂Ω_∂φ)
    @eval $func(Ω::Vector3) = let (_, θ, φ) = cartesian_to_spherical(Ω)
        $func(θ, φ)
    end
end
