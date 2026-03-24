using Serialization

@everywhere include("../AnalyticPhantoms.jl")
using ..AnalyticPhantoms

scenePath = ARGS[1]
dataDir = ARGS[2]
binning = parse(Int64, ARGS[3])

(scene, additive) = open(deserialize, scenePath)

geoms = open(readProjectionGeometries, "$(dataDir)/proj.dat")

pixelSize = 0.139 * binning
camSizePx = div(3008, binning)

print("pixelSize is $(pixelSize)\n")
print("camSizePx is $(camSizePx)\n")

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
