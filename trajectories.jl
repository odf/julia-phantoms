abstract type Trajectory{T <: Real} end


struct HelicalTrajectory{T <: Real} <: Trajectory{T}
    isDoubleHelix::Bool
    numProjTotal::Int
    numProjPerPass::Int
    numProjPerRevolution::Int
    dZ::T
    dTheta::T
end


struct PlaneGridTrajectory{T <: Real} <: Trajectory{T}
    isTwoSided::Bool
    numProjTotal::Int
    numberOfRows::Int
    dZ::T
    dX::T
end


struct StadiumTrajectory{T <: Real} <: Trajectory{T}
    numProjTotal::Int
    numProjPerRound::Int
    width::T
    radius::T
    dZ::T
    dT::T
end


struct BeamGeometry{T <: Real}
    sourceDistance::T
    cameraLength::T
    cameraWidth::T
end


struct Misalignment{T <: Real}
    sourceDistanceOffset::T
    cameraLengthOffset::T
    originOffsetHorizontal::T
    originOffsetVertical::T
    cameraTiltPhi::T
    cameraTiltTheta::T
    cameraTiltPsi::T
end


function makeMisalignment(
    sourceDistanceOffset::T,
    cameraLengthOffset::T=0.0,
    originOffsetHorizontal::T=0.0,
    originOffsetVertical::T=0.0,
    cameraTiltPhi::T=0.0,
    cameraTiltTheta::T=0.0,
    cameraTiltPsi::T=0.0
) where T
    return Misalignment(
        sourceDistanceOffset,
        cameraLengthOffset,
        originOffsetHorizontal,
        originOffsetVertical,
        cameraTiltPhi,
        cameraTiltTheta,
        cameraTiltPsi
    )
end


function writeConeBeamGeometryCommon(
    f::IOStream,
    geom::BeamGeometry{T},
    misalignment::Misalignment{T}=makeMisalignment(0.0)
) where T
    m = misalignment
    sd = geom.sourceDistance + m.sourceDistanceOffset
    cl = geom.cameraLength + m.cameraLengthOffset

    write(f, "specimen_distance: $(sd)\n")
    write(f, "camera_length: $(cl)\n")
    write(f, "origin_offset_horizontal: $(m.originOffsetHorizontal)\n")
    write(f, "origin_offset_vertical: $(m.originOffsetVertical)\n")
    write(f, "detector_tilt_phi: $(m.cameraTiltPhi)\n")
    write(f, "detector_tilt_theta: $(m.cameraTiltTheta)\n")
    write(f, "detector_tilt_psi: $(m.cameraTiltPsi)\n")
    write(f, "beam_hardening_corr_exp: 1\n")
    write(f, "phase_retrieval_reg_parm: 0\n")
end


function writeConeBeamGeometry(
    f::IOStream,
    geom::BeamGeometry{T},
    traj::HelicalTrajectory{T};
    misalignment::Misalignment{T}=makeMisalignment(0.0)
) where T
    num360sPerRev = traj.numProjPerRevolution * traj.dTheta / (2 * pi)
    writeConeBeamGeometryCommon(f, geom, misalignment)

    if traj.isDoubleHelix
        write(f, "trajectory: DOUBLE_HELIX\n")
    else
        write(f, "trajectory: HELIX\n")
    end

    write(f, "projections_per_revolution: $(traj.numProjPerRevolution)\n")
    write(f, "num360s_per_revolution: $(num360sPerRev)\n")
    write(f, "roi: 0\n")
    write(f, "pitch: $(traj.dZ * traj.numProjPerRevolution)\n")
    write(f, "angle_offset: 0\n")
end


function writeConeBeamGeometry(
    f::IOStream,
    geom::BeamGeometry{T},
    traj::PlaneGridTrajectory{T};
    misalignment::Misalignment{T}=makeMisalignment(0.0)
) where T
    writeConeBeamGeometryCommon(f, geom, misalignment)

    # TODO use correct trajectory types once Mango supports them
    if traj.isTwoSided
        write(f, "trajectory: DOUBLE_HELIX\n")
    else
        write(f, "trajectory: HELIX\n")
    end

    write(f, "projections_per_revolution: $(traj.numberOfRows)\n")
    write(f, "num360s_per_revolution: 1\n")
    write(f, "roi: 0\n")
    write(f, "pitch: $(traj.dZ)\n")
    write(f, "angle_offset: 0\n")
end


function writeConeBeamGeometry(
    f::IOStream,
    geom::BeamGeometry{T},
    traj::StadiumTrajectory{T};
    misalignment::Misalignment{T}=makeMisalignment(0.0)
) where T
    tPerRound = 2 * (traj.width + pi * geom.sourceDistance)
    num360sPerRound = traj.numProjPerRound * traj.dT / tPerRound

    writeConeBeamGeometryCommon(f, geom, misalignment)

    # TODO use correct trajectory types once Mango supports them
    write(f, "trajectory: HELIX\n")

    write(f, "projections_per_revolution: $(traj.numProjPerRound)\n")
    write(f, "num360s_per_revolution: $(num360sPerRound)\n")
    write(f, "roi: 0\n")
    write(f, "pitch: $(traj.dZ * traj.numProjPerRound)\n")
    write(f, "angle_offset: 0\n")
end


function idealPitch(sourceDist::T, virtualCamWidth::T)::T where T <: Real
    a = virtualCamWidth / 2.0
    r = sourceDist
    return a / (1.0 + a^2 / r^2) / (atan(a, r) / (2.0 * pi) + 0.25)
end


function projPerRowForSpaceFilling(rowWidth::T, dZ::T)::T where T <: Real
    cos60 = sqrt(3.0) / 2.0
    goldenRatio = (sqrt(5.0) + 1.0) / 2.0

    return floor(sqrt(cos60 * rowWidth / dZ)) + 1.0 / goldenRatio
end


function projPer360ForSpaceFilling(sourceDist::T, dZ::T)::T where T <: Real
    return projPerRowForSpaceFilling(2.0 * pi * sourceDist, dZ)
end


function computeGeometry(
    supportRadius::T;
    coneAngle::T=pi/3,
    magnification::T=50.0,
    sourceDistMin::T=0.0,
    cameraWidth::T=0.0
) where T <: Real

    sourceDist = max(sourceDistMin, supportRadius / sin(coneAngle / 2.0))
    virtualCameraWidth = (2.0 * sourceDist * supportRadius /
                          sqrt(sourceDist^2 - supportRadius^2))

    cameraLength = sourceDist * magnification
    if cameraWidth == 0.0
        cameraWidth = virtualCameraWidth * magnification
    end

    return BeamGeometry(sourceDist, cameraLength, cameraWidth)
end


function helicalTrajectory(
    geom::BeamGeometry{T},
    targetHeight::T,
    voxelSize::T;
    numScanPasses::Int=1,
    numProjFactor::T=1.0,
    isSpaceFilling::Bool=true,
    isDoubleHelix::Bool=false
) where T <: Real

    sourceDist = geom.sourceDistance
    magnification = geom.cameraLength / sourceDist
    halfWidth = 0.5 * geom.cameraWidth / magnification
    supportRadius = halfWidth * sourceDist / sqrt(halfWidth^2 + sourceDist^2)

    pitchPerRev = idealPitch(sourceDist, 2 * halfWidth)
    numProjPerRev = Int(floor(
        3.0 * numProjFactor * supportRadius / voxelSize / numScanPasses))

    dZ = pitchPerRev / numProjPerRev

    numProjPer360 = if isSpaceFilling
        projPer360ForSpaceFilling(sourceDist, dZ)
    else
        numProjPerRev
    end

    dTheta = 2.0 * pi / numProjPer360

    numHelices = if isDoubleHelix 2 else 1 end
    numProjPerPass = Int(ceil((targetHeight + pitchPerRev) / dZ)) * numHelices
    numProjTotal = numProjPerPass * numScanPasses

    # -- Build the result
    return HelicalTrajectory(
        isDoubleHelix, numProjTotal, numProjPerPass, numProjPerRev, dZ, dTheta
    )
end


function planeGridTrajectory(
    geom::BeamGeometry{T},
    scanHeight::T,
    scanWidth::T,
    voxelSize::T;
    objectDepth::T=geom.sourceDist,
    numProjFactor::T=1.0,
    isTwoSided::Bool=true
) where T <: Real

    sd = geom.sourceDistance
    cl = geom.cameraLength
    mag = cl / sd

    camSizePx = geom.cameraWidth / mag / voxelSize
    voxelCoverage = 0.5 * numProjFactor * camSizePx * objectDepth / sd

    dX = dZ = camSizePx / sqrt(voxelCoverage) * voxelSize

    numberOfRows = Int(ceil(scanHeight / dZ))
    numberOfColums = Int(ceil(scanWidth / dX))
    numberOfSides = if isTwoSided 2 else 1 end

    numProjTotal = numberOfRows * numberOfColums * numberOfSides

    return PlaneGridTrajectory(isTwoSided, numProjTotal, numberOfRows, dZ, dX)
end


function stadiumTrajectory(
    geom::BeamGeometry{T},
    targetHeight::T,
    rectangleWidth::T,
    voxelSize::T,
    numProjFactor::T=1.0,
    isSpaceFilling::Bool=true
) where T <: Real

    sourceDist = geom.sourceDistance
    magnification = geom.cameraLength / sourceDist
    pitchPerRound = idealPitch(sourceDist, geom.cameraWidth / magnification)

    distPerRound = 2.0 * (pi * sourceDist + rectangleWidth)
    numProjPerRound = Int(floor(distPerRound * numProjFactor / voxelSize))

    dZ = pitchPerRound / numProjPerRound

    numProjPer360 = if isSpaceFilling
        projPerRowForSpaceFilling(distPerRound, dZ)
    else
        numProjPerRound
    end

    dT = distPerRound / numProjPer360

    numProjTotal = Int(ceil((targetHeight + pitchPerRound) / dZ))

    # -- Build the result
    return StadiumTrajectory(
        numProjTotal, numProjPerRound, rectangleWidth, sourceDist, dZ, dT
    )
end
