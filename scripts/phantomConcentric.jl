using Serialization

include("scenes.jl")

scene = concentric(
    [34.0, 2.5, 50.0], 10^6,
    shifts=[[-10.0, 0.0, 0.0], [0.0, 0.0, 0.0], [10.0, 0.0, 0.0]]
)

open(f -> serialize(f, scene), "concentric.jls", "w")
