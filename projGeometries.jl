using Printf

struct ProjectionGeometry{T <: Real}
    sourcePos::Vec3{T}
    detectorPos::Vec3{T}
    wStep::Vec3{T}
    hStep::Vec3{T}
end


shift(projGeoms::Vector{ProjectionGeometry{T}}, v::Vec3{T}) where T = map(
    pg -> ProjectionGeometry(
        pg.sourcePos + v, pg.detectorPos + v, pg.wStep, pg.hStep
    ),
    projGeoms
)


function writeProjectionGeometries(
    f::IOStream,
    projGeoms::Vector{ProjectionGeometry{T}}
) where T
    write(f,
          "# source_x,   source_y,   source_z,   ",
          "  detector_x,   detector_y,   detector_z,   ",
          "   wstep_x,    wstep_y,    wstep_z,   ",
          "   hstep_x,    hstep_y,    hstep_z\n")

    for pg in projGeoms
        (s, d, w, h) = (pg.sourcePos, pg.detectorPos, pg.wStep, pg.hStep)

        write(f,
              (@sprintf "%10.7f, %10.7f, %10.7f,   " s.x s.y s.z),
              (@sprintf "%12.7f, %12.7f, %12.7f,   " d.x d.y d.z),
              (@sprintf "%10.7f, %10.7f, %10.7f,   " w.x w.y w.z),
              (@sprintf "%10.7f, %10.7f, %10.7f\n" h.x h.y h.z))
    end
end


function readProjectionGeometries(f::IOStream)
    geoms::Vector{ProjectionGeometry{Float64}} = []

    for line in map(strip, readlines(f))
        if line[1] != '#'
            fields = map(s -> parse(Float64, s), split(line, ','))

            sourcePos = Vec3(fields[1], fields[2], fields[3])
            detectorPos = Vec3(fields[4], fields[5], fields[6])
            wStep = Vec3(fields[7], fields[8], fields[9])
            hStep = Vec3(fields[10], fields[11], fields[12])

            item = ProjectionGeometry(sourcePos, detectorPos, wStep, hStep)
            push!(geoms, item)
        end
    end

    return geoms
end


function misalignmentMatrix(misalignment::Misalignment{T}) where T
    lShift = 0.0
    wShift = misalignment.originOffsetHorizontal
    hShift = misalignment.originOffsetVertical

    phi = misalignment.cameraTiltPhi * pi / 180.0
    sinPhi, cosPhi = sin(phi), cos(phi)

    theta = misalignment.cameraTiltTheta * pi / 180.0
    sinTheta, cosTheta = sin(theta), cos(theta)

    psi = misalignment.cameraTiltPsi * pi / 180.0
    sinPsi, cosPsi = sin(psi), cos(psi)

    shift = [
        1 0 0 lShift
        0 1 0 wShift
        0 0 1 hShift
        0 0 0      1
    ]

    rotPhi = [
        1       0      0 0
        0  cosPhi sinPhi 0
        0 -sinPhi cosPhi 0
        0       0      0 1
    ]

    rotTheta = [
         cosTheta sinTheta 0 0
        -sinTheta cosTheta 0 0
                0        0 1 0
                0        0 0 1
    ]

    rotPsi = [
         cosPsi 0 sinPsi 0
              0 1      0 0
        -sinPsi 0 cosPsi 0
              0 0      0 1
    ]

    return shift * rotPhi * rotTheta * rotPsi
end


function coefficients(n::Int)
    a = [[0, 0], [1, 0], [2, 0]]
    b = [[0, 0], [0, 1], [0, 2]]
    f = 1
    s = [[0, 0]]

    while length(s) < n
        f = f / 3
        s = [u * f + v * f + w for u in b for v in a for w in s]
    end

    return hcat(s[1:n]...)
end


function rotateVector(vec::Vector{T}, thetas::Vector{T}) where T <: Real
    sines = sin.(thetas)
    cosines = cos.(thetas)

    return Vec3.(
        cosines * vec[1] - sines * vec[2],
        sines * vec[1] + cosines * vec[2],
        vec[3]
    )
end


function convertProjectionGeometries(
    cl::T,
    sd::T,
    thetas::Vector{T},
    zs::Vector{T},
    xs::Vector{T},
    misalignment::Misalignment{T},
    motionVectors::Vector{Vec3{T}}
) where T <: Real
    mat = inv(misalignmentMatrix(misalignment))

    sourcePos = rotateVector(
        [-sd; 0.0; 0.0; 0.0],
        thetas
    ) + Vec3.(xs, 0.0, zs)

    detectorPos = rotateVector(
        mat * [0.0; 0.0; 0.0; 1.0] + [cl - sd; 0.0; 0.0; 0.0],
        thetas
    ) + Vec3.(xs, 0.0, zs)

    wSteps = rotateVector(
        mat * [0.0; 1.0; 0.0; 0.0],
        thetas
    )

    hSteps = rotateVector(
        mat * [0.0; 0.0; 1.0; 0.0],
        thetas
    )

    dSteps = normalize.(cross.(wSteps, hSteps))
    bases = tuple.(wSteps, hSteps, dSteps)

    detectorPos = detectorPos .- (bases .* motionVectors)

    return ProjectionGeometry.(sourcePos, detectorPos, wSteps, hSteps)
end


function makeProjectionGeometries(
    geom::BeamGeometry{T},
    traj::HelicalTrajectory{T};
    misalignment::Misalignment{T}=makeMisalignment(0.0),
    motionVectors::Vector{Vec3{T}}=[Vec3{T}(0.0, 0.0, 0.0)],
    recursive::Bool=false
) where T <: Real
    sd = geom.sourceDistance + misalignment.sourceDistanceOffset
    cl = geom.cameraLength + misalignment.cameraLengthOffset
    dZ = traj.dZ
    dTheta = traj.dTheta

    indices = 0 : traj.numProjTotal - 1

    nrPasses = div(traj.numProjTotal, traj.numProjPerPass)
    numHelices = if traj.isDoubleHelix 2 else 1 end
    projPerHelix = div(traj.numProjPerPass, numHelices)

    pass = div.(indices, traj.numProjPerPass)
    idxInPass = mod.(indices, traj.numProjPerPass)
    helix = div.(idxInPass, projPerHelix)
    idxInHelix = (
        mod.(idxInPass, projPerHelix) .* (-1).^helix
        + helix * (projPerHelix - 1)
    )

    tau = 2 * pi
    n = ceil(tau / dTheta)
    u = [dTheta, dZ]
    v = [(n * dTheta) % tau, n * dZ]

    if recursive
        c = coefficients(nrPasses)
        s = hcat(u + v, v) * c
    else
        c = hcat([(i - 1) * [0.5545497, 0.308517] for i in 1:nrPasses]...)
        s = hcat(v, u) * c
    end

    sZ = [s[2, i + 1] for i in pass]
    sTheta = [s[1, i + 1] for i in pass]

    dZPass = mod.(sZ, dZ)
    dThetaPass = mod.(sTheta - dTheta * div.(sZ, dZ), tau)

    dThetaHelix = helix * pi

    zs = idxInHelix * dZ + dZPass
    zs = zs .- (maximum(zs) - minimum(zs)) / 2.0

    thetas = mod.(idxInHelix * dTheta + dThetaPass + dThetaHelix, tau)

    xs = zs * 0.0

    return convertProjectionGeometries(
        cl, sd, thetas, zs, xs, misalignment, motionVectors
    )
end


function makeProjectionGeometries(
    geom::BeamGeometry{T},
    traj::PlaneGridTrajectory{T};
    misalignment::Misalignment{T}=makeMisalignment(0.0),
    motionVectors::Vector{Vec3{T}}=[Vec3{T}(0.0, 0.0, 0.0)]
) where T <: Real

    sd = geom.sourceDistance + misalignment.sourceDistanceOffset
    cl = geom.cameraLength + misalignment.cameraLengthOffset
    dZ = traj.dZ
    dX = traj.dX

    numSides = if traj.isTwoSided 2 else 1 end
    projPerSide = div(traj.numProjTotal, numSides)

    numRows = traj.numberOfRows
    numCols = projPerSide / numRows

    indices = 0 : traj.numProjTotal - 1
    side = div.(indices, projPerSide)
    idxOnSide = mod.(indices, projPerSide)

    flipped = (-1).^side
    thetas = flipped * (pi / 2)

    col = div.(idxOnSide, numCols)
    row = mod.(idxOnSide, numCols)

    height = dZ * (numCols - 1)
    width = dX * (numRows - 1)

    zs = row .* dZ .- height / 2.0
    xs = col .* dX .- width / 2.0

    return convertProjectionGeometries(
        cl, sd, thetas, zs, xs .* flipped, misalignment, motionVectors
    )
end


function makeProjectionGeometries(
    geom::BeamGeometry{T},
    traj::StadiumTrajectory{T};
    misalignment::Misalignment{T}=makeMisalignment(0.0),
    motionVectors::Vector{Vec3{T}}=[Vec3{T}(0.0, 0.0, 0.0)]
) where T <: Real
    sd = geom.sourceDistance + misalignment.sourceDistanceOffset
    cl = geom.cameraLength + misalignment.cameraLengthOffset

    width = traj.width
    r = geom.sourceDistance
    arclen = pi * r

    indices = collect(0 : traj.numProjTotal - 1)
    tInRound = mod.(indices .* traj.dT, 2.0 * (arclen + width))

    inFirstStretch = tInRound .< width
    inFirstBend = width .<= tInRound .< (width + arclen)
    inSecondStretch = (width + arclen) .<= tInRound .< (2 * width + arclen)
    inSecondBend = (2 * width + arclen) .<= tInRound

    zs = indices * traj.dZ
    zs = zs .- (maximum(zs) - minimum(zs)) / 2.0

    xs =
        (tInRound .- width / 2) .* inFirstStretch .+
        width / 2 .* inFirstBend .+
        (width / 2 .- (tInRound .- width .- arclen)) .* inSecondStretch .+
        -width / 2 .* inSecondBend

    thetas =
        -pi / 2 .* inFirstStretch .+
        ((width .- tInRound) ./ r .- pi / 2) .* inFirstBend .+
        pi / 2 .* inSecondStretch .+
        ((2 * width .- arclen .- tInRound) ./ r .+ pi / 2) .* inSecondBend

    return convertProjectionGeometries(
        cl, sd, thetas, zs, xs, misalignment, motionVectors
    )
end


function pushOutFunction(
    pushFactor,
    maxCamLenFactor::T,
    maxSourcePosY::T
) where T <: Real

    function convert(pg::ProjectionGeometry{T})
        srcIn, camIn, wStep, hStep =
            pg.sourcePos, pg.detectorPos, pg.wStep, pg.hStep

        sd = norm(Vec3(srcIn.x, srcIn.y, 0.0))
        cl = norm(Vec3(camIn.x - srcIn.x, camIn.y - srcIn.y, 0.0))

        if srcIn.y * pushFactor(srcIn) > maxSourcePosY
            m = (Vec3(-1.0,  0.0, 0.0),
                 Vec3(0.0, -1.0, 0.0),
                 Vec3(0.0,  0.0, 1.0))
            srcIn = m * srcIn
            camIn = m * camIn
            wStep = m * wStep
            hStep = m * hStep
        end

        fSrc = pushFactor(srcIn)
        if fSrc <= maxCamLenFactor
            fCam = fSrc
        else
            fCam = (cl * maxCamLenFactor - sd * fSrc) / (cl - sd)
        end

        srcOut = Vec3(fSrc * srcIn.x, fSrc * srcIn.y, srcIn.z)
        camOut = Vec3(fCam * camIn.x, fCam * camIn.y, camIn.z)

        return ProjectionGeometry(srcOut, camOut, wStep, hStep)
    end

    return convert
end


function pushFactorCircular(
    obstacleCenterY::T, obstacleRadius::T
) where T <: Real

    function pushFactor(src)
        sd = norm(Vec3(src.x, src.y, 0.0))
        return max(1.0, (obstacleCenterY * src.y / sd + obstacleRadius) / sd)
    end

    return pushFactor
end


function pushFactorPolygonal(
    obstacleCorners::Vector{Tuple{T, T}}, obstacleOffset::T
) where T <: Real

    function pushFactor(src)
        sd = norm(Vec3(src.x, src.y, 0.0))
        srcDir = src / sd
        dists = map(p -> dot(Vec3(p[1], p[2], 0.0), srcDir), obstacleCorners)
        return max(1.0, (maximum(dists) - obstacleOffset) / sd)
    end

    return pushFactor
end
