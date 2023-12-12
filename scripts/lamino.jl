using AnalyticPhantoms
using ArgParse
using Distributed
using Random

@everywhere include("$($(@__DIR__))/scenes.jl")


function parseCommandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--binFactor", "-b"
            help = "Bin factor"
            arg_type = Int64
            default = 4
        "--cylindrical", "-c"
            help = "use a cylindrical (helical) trajectory"
            action = :store_true
        "--coldRun", "-C"
            help = "only computes metadata, but no projections"
            action = :store_true
        "--heightByDiameter", "-H"
            help = "Volume height as a multiple of volume diameter"
            arg_type = Float64
            default = 0.5
        "--magnification", "-m"
            help = "the desired magnification factor"
            arg_type = Float64
            default = 50.0
        "--nrElements", "-n"
            help = "the number of scene elements"
            arg_type = Int64
            default = 1000000
        "--outdir", "-o"
            help = "the output directory"
            default = "-"
        "--projectionDensity", "-p"
            help = "projection density per area relative to unwrapped SFT"
            arg_type = Float64
            default = 1.0
        "--oversampling", "-r"
            help = "oversampling factor"
            arg_type = Int64
            default = 2
        "--sourceDistance", "-s"
            help = "the desired source distance in millimetres"
            arg_type = Float64
            default = 5.0
        "--rasterize", "-z"
            help = "produce a volume image instead of a projection set"
            action = :store_true
    end

    return parse_args(s)
end


args = parseCommandline()

coldRun = args["coldRun"]
outdir = if (args["outdir"] == "-") "concentric" else args["outdir"] end



# -- Experiment geometry settings

volDepthMm = 5.0
volWidthMm = 7.0 * volDepthMm
volHeightMm = volWidthMm * args["heightByDiameter"]

sourceSeparationMm = 0.5

detectorSizeMm = 460.0
pixelSizeMm = 0.15 * args["binFactor"]
maxCamLengthMm = 835.0



# -- Experiment geometry calculations

helical = args["cylindrical"]

minSourceDistMm = sourceSeparationMm +
    if helical
        sqrt(volWidthMm^2 + volDepthMm^2) / 2
    else
        volDepthMm / 2
    end

sourceDistanceMm = max(args["sourceDistance"], minSourceDistMm)
cameraLengthMm = min(args["magnification"] * sourceDistanceMm, maxCamLengthMm)

magnification = cameraLengthMm / sourceDistanceMm
voxelSizeMm = pixelSizeMm / magnification



# -- Scene computation

volumeShape = [ volWidthMm, volDepthMm, volHeightMm ]

if coldRun
    scene = [ SceneItem(Sphere(Vec3(0.0, 0.0, 0.0), 1.0), 0.1) ]
    additive = true
else
    Random.seed!(135796531)
    scene, additive = concentric(volumeShape, args["nrElements"])
end


# -- Output generation

mkpath(outdir)

scanWidthMm = 0.4 * volWidthMm * (50.0 / magnification)
scanHeightMm = volHeightMm * (50.0 / magnification)

if args["rasterize"]
    outputShape = [ scanWidthMm, volDepthMm, scanHeightMm ] .+ 1.0

    rasterizeAndWrite(
        scene,
        outputShape,
        voxelSizeMm,
        additive=additive,
        outdir=outdir,
        oversamplingFactor=args["oversampling"],
        coldRun=coldRun
    )
else
    geom = BeamGeometry(sourceDistanceMm, cameraLengthMm, detectorSizeMm)

    trajectory =
        if helical
            helicalTrajectory(
                geom, volHeightMm, voxelSizeMm,
                numProjFactor=args["projectionDensity"]
            )
        else
            planeGridTrajectory(
                geom, scanHeightMm, scanWidthMm, voxelSizeMm,
                numProjFactor=args["projectionDensity"]
            )
        end

    projectAndWrite(
        scene,
        geom,
        trajectory,
        voxelSizeMm,
        additive=additive,
        outdir=outdir,
        oversamplingFactor=args["oversampling"],
        coldRun=coldRun
    )
end
