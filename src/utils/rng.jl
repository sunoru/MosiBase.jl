import RandomNumbers: gen_seed
import RandomNumbers.Xorshifts: Xoshiro256StarStar

make_seed(seed::Nullable{Integer} = nothing) = seed === nothing ? gen_seed(UInt64) : (seed % UInt64)
new_rng(seed::UInt64) = Xoshiro256StarStar(seed)

function restore_rng(state::NTuple{4, UInt64})
    rng = Xoshiro256StarStar(state)
    rng.x, rng.y, rng.z, rng.w = state
    rng
end

rng_state(rng::Xoshiro256StarStar) = [rng.x, rng.y, rng.z, rng.w]
