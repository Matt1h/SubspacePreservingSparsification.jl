using Random
using SSA: separated_min_max

v = rand(MersenneTwister(0), Float32, 20)
@test separated_min_max(v, 0.2f0) == (0.009370327f0, 0.14352024f0, 0.23794222f0, 0.9680481f0)
@test separated_min_max([1, 2, 3], 2) == (1, 2, 2, 3)
