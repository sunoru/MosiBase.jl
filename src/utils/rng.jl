using Random: make_seed, Xoshiro

const MosiRNG = Xoshiro

new_rng(seed) = MosiRNG(seed)

restore_rng(state::NTuple{4, UInt64}) = MosiRNG(state...)

rng_state(rng::MosiRNG) = (rng.s0, rng.s1, rng.s2, rng.s3)
