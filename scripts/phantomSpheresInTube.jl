using Serialization

include("scenes.jl")

scene = spheres([12.0, 12.0, 36.0], 10^4, 0.0025)

open(f -> serialize(f, scene), "spheresInTube.jls", "w")
