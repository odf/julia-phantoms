using AnalyticPhantoms
using ArgParse


function parseCommandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--binFactor", "-b"
            help = "bin factor"
            arg_type = Int64
            default = 2
        "--objectDepth", "-d"
            help = "scanned object depth in mm"
            arg_type = Float64
            default = 2.5
        "--scanHeight", "-H"
            help = "scanned region height in mm"
            arg_type = Float64
            default = 10.0
        "--magnification", "-m"
            help = "magnification factor"
            arg_type = Float64
            default = 40.0
        "--outdir", "-o"
            help = "output directory"
            default = "trajectory"
        "--projectionDensity", "-p"
            help = "projection density per area relative to unwrapped SFT"
            arg_type = Float64
            default = 1.0
        "--sourceDistance", "-s"
            help = "source distance in mm"
            arg_type = Float64
            default = 5.0
        "--scanWidth", "-w"
            help = "scanned region width in mm"
            arg_type = Float64
            default = 10.0
    end

    return parse_args(s)
end


args = parseCommandline()

binFactor = args["binFactor"]
objectDepthMm = args["objectDepth"]
scanHeightMm = args["scanHeight"]
magnification = args["magnification"]
projDensity = args["projectionDensity"]
sourceDistanceMm = args["sourceDistance"]
scanWidthMm = args["scanWidth"]

cameraLengthMm = sourceDistanceMm * magnification
pixelSizeMm = 0.139 * binFactor
voxelSizeMm = pixelSizeMm / magnification
camSizePx = 3008 / binFactor
camSizeMm = camSizePx * pixelSizeMm

geom = BeamGeometry(sourceDistanceMm, cameraLengthMm, camSizeMm)

trajectory = planeGridTrajectory(
    geom, scanHeightMm, scanWidthMm, voxelSizeMm,
    objectDepth=objectDepthMm,
    numProjFactor=projDensity
)

projGeoms = makeProjectionGeometries(geom, trajectory)


outdir = args["outdir"]
mkdir(outdir)

open("$(outdir)/geom.dat", "w") do f
    writeConeBeamGeometry(f, geom, trajectory)
end

open("$(outdir)/proj.dat", "w") do f
    writeProjectionGeometries(f, projGeoms)
end

open("$(outdir)/info.dat", "w") do f
    write(f, "pixelSize: $(pixelSizeMm)\n")
    write(f, "camSizePx: $(camSizePx)\n")
end
