using Serialization

include("scenes.jl")

tile, isAdditive = concentric([5.0, 2.5, 5.0], 20000)

scene = Array{SceneItem{Float64}, 1}()

for shiftX in -15.0 : 5.0 : 15.0
    for shiftZ in -15.0 : 5.0 : 15.0
        s = shift(Vec3(shiftX, 0.0, shiftZ))

        for item in tile[2:end]
            push!(scene, SceneItem(s(item.object), item.value))
        end
    end
end

push!(scene, SceneItem(Cuboid(35.0, 2.5, 35.0), 0.05))

open(f -> serialize(f, (scene, isAdditive)), "tiled.jls", "w")
