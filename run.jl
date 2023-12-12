import Distributed


function linearResample(seriesIn::Vector{T}, nOut::Int) where T
    nIn = length(seriesIn)
    scale = (nIn - 1) / (nOut - 1)

    function pluck(i)
        k0, a = divrem(scale * (i - 1) + 1, 1)
        k = Int(k0)

        if k >= nIn
            seriesIn[nIn]
        else
            (1.0 - a) * seriesIn[k] + a * seriesIn[k+1]
        end
    end

    return map(pluck, 1:nOut)
end


function projectAndWrite(
    sceneAt::Union{Function, Vector{SceneItem{T}}},
    projGeoms::Vector{ProjectionGeometry{T}},
    pixelSize::T,
    camSizePx::Int;

    outdir::String=".",
    additive=false,
    oversamplingFactor::Int=1,
    parallel::Bool=true
) where T <: Real

    if typeof(sceneAt) == Vector{SceneItem{T}}
        scene = sceneAt
        sceneAt = (i, n) -> scene
    end

    open("$(outdir)/info.dat", "a") do f
        write(f, "pixelSize: $(pixelSize)\n")
        write(f, "camSizePx: $(camSizePx)\n")
    end

    open("$(outdir)/proj.dat", "w") do f
        writeProjectionGeometries(f, projGeoms)
    end

    writeImage = function(img, i)
        name = "$(outdir)/proj-$(@sprintf("%05d", i)).raw"
        open(name, "w") do f
            write(f, map(Float32, img))
        end
    end

    if parallel
        Distributed.@sync for (i, p) in enumerate(Distributed.workers())
            Distributed.@spawnat p projectScene(
                sceneAt,
                projGeoms,
                pixelSize,
                camSizePx,
                writeImage,
                oversamplingFactor,
                additive,
                i,
                Distributed.nworkers())
        end
    else
        projectScene(
            sceneAt,
            projGeoms,
            pixelSize,
            camSizePx,
            writeImage,
            oversamplingFactor,
            additive)
    end
end


function projectAndWrite(
    sceneAt::Union{Function, Vector{SceneItem{T}}},
    geom::BeamGeometry{T},
    trajectory::Trajectory{T},
    voxelSize::T;

    outdir::String=".",
    additive=false,
    oversamplingFactor::Int=1,
    misalignment::Misalignment{T}=makeMisalignment(0.0),
    motionVectors::Vector{Vec3{T}}=[Vec3{T}(0.0, 0.0, 0.0)],
    projectionGeometries::Union{Nothing, Vector{ProjectionGeometry{T}}}=nothing,
    pushOut::Function=identity,
    coldRun::Bool=false,
    parallel::Bool=true
) where T <: Real

    if typeof(sceneAt) == Vector{SceneItem{T}}
        scene = sceneAt
        sceneAt = (i, n) -> scene
    end

    if projectionGeometries == nothing
        motionVectors = linearResample(motionVectors, trajectory.numProjTotal)
        projGeoms = map(pushOut, makeProjectionGeometries(
            geom,
            trajectory,
            misalignment=misalignment,
            motionVectors=motionVectors
        ))
    else
        projGeoms = projectionGeometries
    end

    pixelSize = voxelSize * geom.cameraLength / geom.sourceDistance
    camSizePx = 2 * Int(ceil(0.5 * geom.cameraWidth / pixelSize)) + 1

    open("$(outdir)/geom.dat", "w") do f
        writeConeBeamGeometry(f, geom, trajectory)
    end

    open("$(outdir)/geom-actual.dat", "w") do f
        writeConeBeamGeometry(f, geom, trajectory, misalignment=misalignment)
    end

    open("$(outdir)/info.dat", "a") do f
        write(f, "pixelSize: $(pixelSize)\n")
        write(f, "camSizePx: $(camSizePx)\n")
    end

    open("$(outdir)/proj.dat", "w") do f
        writeProjectionGeometries(f, projGeoms)
    end

    if coldRun
        return
    end

    writeImage = function(img, i)
        name = "$(outdir)/proj-$(@sprintf("%05d", i)).raw"
        open(name, "w") do f
            write(f, map(Float32, img))
        end
    end

    if parallel
        Distributed.@sync for (i, p) in enumerate(Distributed.workers())
            Distributed.@spawnat p projectScene(
                sceneAt,
                projGeoms,
                pixelSize,
                camSizePx,
                writeImage,
                oversamplingFactor,
                additive,
                i,
                Distributed.nworkers())
        end
    else
        projectScene(
            sceneAt,
            projGeoms,
            pixelSize,
            camSizePx,
            writeImage,
            oversamplingFactor,
            additive)
    end
end


function rasterizeAndWrite(
    scene::Vector{SceneItem{T}},
    volumeShape::Vector{T},
    voxelSize::T;

    outdir::String=".",
    additive=false,
    oversamplingFactor::Int=1,
    coldRun::Bool=false,
    parallel::Bool=true
) where T <: Real

    sx, sy, sz = map(s -> 2 * Int(ceil(0.5 * s)) + 1, volumeShape / voxelSize)
    cx, cy, cz = div.((sx, sy, sz), 2)

    open("$(outdir)/info.dat", "a") do f
        write(f, "voxelSize: $(voxelSize)\n")
        write(f, "sizeXYZ: $(sx) $(sy) $(sz)\n")
    end

    if coldRun
        return
    end

    xSlice = Array{T}(undef, sy, sz)
    ySlice = Array{T}(undef, sx, sz)

    writeSlice = function(img, axis, i)
        name = "$(outdir)/slice$(axis)-$(@sprintf("%05d", i)).raw"
        open(name, "w") do f
            write(f, map(Float32, img))
        end
    end

    writeImage = function(img, i)
        name = "$(outdir)/volz-$(@sprintf("%05d", i)).raw"
        open(name, "w") do f
            write(f, map(Float32, img))
        end

        xSlice[:, i] = img[cx, :]
        ySlice[:, i] = img[:, cy]

        if i == cz
            writeSlice(img, "Z", i - 1)
        end
    end

    if parallel
        Distributed.@sync for (i, p) in enumerate(Distributed.workers())
            Distributed.@spawnat p rasterizeScene(
                scene,
                volumeShape,
                voxelSize,
                writeImage,
                oversamplingFactor,
                additive,
                128,
                i,
                Distributed.nworkers())
        end
    else
        rasterizeScene(
            scene,
            volumeShape,
            voxelSize,
            writeImage,
            oversamplingFactor,
            additive)
    end

    writeSlice(xSlice, "X", cx - 1)
    writeSlice(ySlice, "Y", cy - 1)
end


# -- Alternate methods for backward compatibility.

function projectAndWrite(
    sceneAt::Union{Function, Vector{SceneItem{T}}},
    supportRadius::T,
    targetHeight::T,
    voxelSize::T;

    outdir::String=".",
    additive=false,
    oversamplingFactor::Int=1,
    magnification::T=50.0,
    sourceDistMin::T=0.0,
    cameraWidth::T=0.0,
    misalignment::Misalignment{T}=makeMisalignment(0.0),
    motionVectors::Vector{Vec3{T}}=[Vec3{T}(0.0, 0.0, 0.0)],
    numScanPasses::Int=1,
    numProjFactor::T=1.0,
    isSpaceFilling::Bool=true,
    isDoubleHelix::Bool=false,
    trajectory::Union{Nothing, HelicalTrajectory{T}}=nothing,
    projectionGeometries::Union{Nothing, Vector{ProjectionGeometry{T}}}=nothing,
    pushOut::Function=identity,
    coldRun::Bool=false,
    parallel::Bool=true
) where T <: Real

    if trajectory == nothing
        geom = computeGeometry(
            supportRadius,
            magnification=magnification,
            sourceDistMin=sourceDistMin,
            cameraWidth=cameraWidth
        )

        trajectory = helicalTrajectory(
            geom,
            targetHeight,
            voxelSize,
            numScanPasses=numScanPasses,
            numProjFactor=numProjFactor,
            isSpaceFilling=isSpaceFilling,
            isDoubleHelix=isDoubleHelix
        )
    end

    projectionAndWrite(
        sceneAt,
        geom,
        trajectory,
        voxelSize,
        outdir=outdir,
        additive=additive,
        oversamplingFactor=oversamplingFactor,
        misalignment=misalignment,
        motionVectors=motionVectors,
        projectionGeometries=projectionGeometries,
        pushOut=pushOut,
        coldRun=coldRun,
        parallel=parallel
    )
end


function rasterizeAndWrite(
    scene::Vector{SceneItem{T}},
    supportRadius::T,
    targetHeight::T,
    voxelSize::T;

    outdir::String=".",
    additive=false,
    oversamplingFactor::Int=1,
    coldRun::Bool=false,
    parallel::Bool=true
) where T <: Real

    volumeShape = [ 2 * supportRadius, 2 * supportRadius, targetHeight ]
    rasterizeAndWrite(
        scene,
        volumeShape,
        voxelSize,
        outdir=outdir,
        additive=additive,
        oversamplingFactor=oversamplingFactor,
        coldRun=coldRun,
        parallel=parallel
    )
end
