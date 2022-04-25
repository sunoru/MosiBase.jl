using Random: Xoshiro

const MosiRNG = Xoshiro

new_rng(seed) = MosiRNG(seed)

restore_rng(state::NTuple{4, UInt64}) = MosiRNG(state...)

rng_state(rng::MosiRNG) = (rng.s0, rng.s1, rng.s2, rng.s3)

make_seed() = UInt64(Random.make_seed()[1])
make_seed(n::Integer) = UInt64(Random.make_seed(n)[1])

macro _with_default_rng(func_def)
    @capture(func_def, function f_(__) _ end) ||
        @capture(func_def, f_(__) = _) ||
            error("Could not parse function definition")
    insert!(func_def.args[1].args, 2, :(rng::AbstractRNG))
    esc(quote
        $func_def
        $f(args...; kwargs...) = $f(GLOBAL_RNG, args...; kwargs...)
    end)
end