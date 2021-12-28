abstract type Ensemble end

struct CanonicalEnsemble <: Ensemble end
struct MicrocanonicalEnsemble <: Ensemble end
struct PotentialEnergyLandscapeEnsemble <: Ensemble
    E_L::Float64
end
