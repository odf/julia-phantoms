using Serialization

include("scenes.jl")

scene = spheres([12.0, 12.0, 12.0], 4, 0.0075)

open(f -> serialize(f, scene), "spheresInTubeSimple.jls", "w")
