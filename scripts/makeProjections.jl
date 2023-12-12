using Serialization
using AnalyticPhantoms

scenePath = ARGS[1]
dataDir = ARGS[2]

(scene, additive) = open(deserialize, scenePath)

geoms = open(readProjectionGeometries, "$(dataDir)/proj.dat")

pixelSize = 0.278
camSizePx = 1504

open("$(dataDir)/info.dat") do f
    for line in readlines(f)
        key, val = map(strip, split(line, ':'))
        if key == "pixelSize"
            pixelSize = parse(Float64, val)
        elseif key == "camSizePx"
            camSizePx = parse(Float64, val)
        end
    end
end


projectAndWrite(
    scene,
    geoms,
    pixelSize,
    camSizePx,
    outdir=dataDir,
    additive=additive,
    oversamplingFactor=2,
    parallel=true
)
