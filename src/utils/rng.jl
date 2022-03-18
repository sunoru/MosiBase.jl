using Random: make_seed, Xoshiro

new_rng(seed) = Xoshiro(seed)

restore_rng(state::NTuple{4, UInt64}) = Xoshiro(state...)

rng_state(rng::Xoshiro) = (rng.s0, rng.s1, rng.s2, rng.s3)
