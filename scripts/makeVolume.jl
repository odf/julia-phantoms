using Serialization
using AnalyticPhantoms

scenePath = ARGS[1]
voxelSize = parse(Float64, ARGS[2])
outputShape = map(s -> parse(Float64, s), split(ARGS[3], ','))

(scene, additive) = open(deserialize, scenePath)

outdir = "vol_raw"
mkpath(outdir)

rasterizeAndWrite(
    scene,
    outputShape,
    voxelSize,
    additive=additive,
    outdir=outdir,
    oversamplingFactor=2,
    parallel=true
)
