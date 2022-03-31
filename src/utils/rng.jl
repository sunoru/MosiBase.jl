import Random
using Random: Xoshiro

const MosiRNG = Xoshiro

new_rng(seed) = MosiRNG(seed)

restore_rng(state::NTuple{4, UInt64}) = MosiRNG(state...)

rng_state(rng::MosiRNG) = (rng.s0, rng.s1, rng.s2, rng.s3)

make_seed() = UInt64(Random.make_seed()[1])
make_seed(n::Integer) = UInt64(Random.make_seed(n)[1])
