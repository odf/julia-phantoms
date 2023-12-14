struct SceneItem{T <: Real}
    object::SceneObject{T}
    value::T
end


function projection(geom::ProjectionGeometry{T}, pixelSize::T) where T <: Real
    wStep = pixelSize * geom.wStep
    hStep = pixelSize * geom.hStep
    n = normalize(cross(wStep, hStep))
    m = inv((wStep, hStep, n))
    v = geom.detectorPos - geom.sourcePos
    s = geom.sourcePos
    nv = dot(n, v)

    function proj(p::Vec3{T})
        p = p - s
        return m * (nv / dot(n, p) * p - v)
    end

    return proj
end


function scanRanges(
    geom::ProjectionGeometry{T},
    pixelSize::T,
    object::SceneObject{T},
    xlo::Int,
    xhi::Int,
    ylo::Int,
    yhi::Int
) where T <: Real

    projectPoint = projection(geom, pixelSize)

    xmin = ymin = Inf
    xmax = ymax = -Inf

    corners::Vector{Vec3{T}} = hullCorners(object)

    for c in corners
        p = projectPoint(c)
        xmin = min(xmin, p.x)
        xmax = max(xmax, p.x)
        ymin = min(ymin, p.y)
        ymax = max(ymax, p.y)
    end

    xrange = max(xlo, Int(floor(xmin))) : min(xhi, Int(ceil(xmax)))
    yrange = max(ylo, Int(floor(ymin))) : min(yhi, Int(ceil(ymax)))

    return xrange, yrange
end


function scanRanges(
    geom::ProjectionGeometry{T},
    pixelSize::T,
    object::Sphere{T},
    xlo::Int,
    xhi::Int,
    ylo::Int,
    yhi::Int
) where T <: Real

    src = geom.sourcePos
    cam = geom.detectorPos

    camNormal = normalize(cross(geom.wStep, geom.hStep))
    p = object.center - src
    center = src + dot(camNormal, cam - src) / dot(camNormal, p) * p
    toSource = src - center
    distanceToSource = norm(toSource)

    sin_a = object.radius / norm(p)
    cos_a = sqrt(max(0.0, 1.0 - sin_a^2))

    sin_b = abs(dot(camNormal, toSource)) / distanceToSource
    cos_b = -sqrt(max(0.0, 1.0 - sin_b^2))

    maxDist = distanceToSource * abs(sin_a / (sin_a * cos_b + cos_a * sin_b))

    c = center - cam
    cx = dot(c, normalize(geom.wStep))
    cy = dot(c, normalize(geom.hStep))

    xmin = (cx - maxDist) / pixelSize
    xmax = (cx + maxDist) / pixelSize
    ymin = (cy - maxDist) / pixelSize
    ymax = (cy + maxDist) / pixelSize

    xrange = max(xlo, Int(floor(xmin))) : min(xhi, Int(ceil(xmax)))
    yrange = max(ylo, Int(floor(ymin))) : min(yhi, Int(ceil(ymax)))

   return xrange, yrange
end


function projectObject!(
    geom::ProjectionGeometry{T},
    pixelSize::T,
    object::SceneObject{T},
    value::T,
    buffer::Array{Buffer{T}, 2}
) where T <: Real

    nx, ny = size(buffer)
    hx, hy = div(nx, 2), div(ny, 2)
    xrange, yrange = scanRanges(
        geom, pixelSize, object, -hx, -hx + nx - 1, -hy, -hy + nx - 1
    )

    wStep = pixelSize * geom.wStep
    hStep = pixelSize * geom.hStep

    for y in yrange
        for x in xrange
            p = geom.detectorPos + wStep * T(x) + hStep * T(y)
            s = buffer[x + hx + 1, y + hy + 1]
            pushIntersections!(object, geom.sourcePos, p,
                               (near, far) -> update!(s, near, far, value))
        end
    end
end


function projectObject!(
    geom::ProjectionGeometry{T},
    pixelSize::T,
    object::Cylinder{T},
    value::T,
    buffer::Array{Buffer{T}, 2}
) where T <: Real

    nx, ny = size(buffer)
    hx, hy = div(nx, 2), div(ny, 2)
    xrange, yrange = scanRanges(
        geom, pixelSize, object, -hx, -hx + nx - 1, -hy, -hy + nx - 1
    )

    wStep = pixelSize * geom.wStep
    hStep = pixelSize * geom.hStep

    rayOrigin = geom.sourcePos
    c = rayOrigin.x^2 + rayOrigin.y^2 - object.radius^2
    bottom = -object.height / 2.0 - rayOrigin.z
    top = object.height / 2.0 - rayOrigin.z

    for y in yrange
        for x in xrange
            rayTarget = geom.detectorPos + wStep * T(x) + hStep * T(y)
            rayDir = normalize(rayTarget - rayOrigin)

            a = rayDir.x^2 + rayDir.y^2
            b = rayOrigin.x * rayDir.x + rayOrigin.y * rayDir.y
            discriminant = b^2 - a * c

            s = buffer[x + hx + 1, y + hy + 1]

            if a == 0
                if c < 0
                    if rayDir.z > 0
                        update!(s, bottom, top, value);
                    else
                        update!(s, -top, -bottom, value);
                    end
                end
            elseif discriminant > 0
                d = sqrt(discriminant)

                if rayDir.z == 0
                    if bottom <= 0 <= top
                        update!(s, -b - d, -b + d, value)
                    end
                else
                    t1 = bottom / rayDir.z
                    t2 = top / rayDir.z
                    near = max((-b - d) / a, min(t1, t2))
                    far = min((-b + d) / a, max(t1, t2))

                    if far > near
                        update!(s, near, far, value)
                    end
                end
            end
        end
    end
end


function projectObject!(
    geom::ProjectionGeometry{T},
    pixelSize::T,
    object::Sphere{T},
    value::T,
    buffer::Array{Buffer{T}, 2}
) where T <: Real

    nx, ny = size(buffer)
    hx, hy = div(nx, 2), div(ny, 2)
    xrange, yrange = scanRanges(
        geom, pixelSize, object, -hx, -hx + nx - 1, -hy, -hy + nx - 1
    )

    wStep = pixelSize * geom.wStep
    hStep = pixelSize * geom.hStep

    relativeOrigin = geom.sourcePos - object.center
    c = dot(relativeOrigin, relativeOrigin) - object.radius^2
    limit = sqrt(max(0.0, c))

    for y in yrange
        for x in xrange
            rayTarget = geom.detectorPos + wStep * T(x) + hStep * T(y)
            rayDir = normalize(rayTarget - geom.sourcePos)
            b = dot(relativeOrigin, rayDir)

            if abs(b) > limit
                d = sqrt(b^2 - c)
                s = buffer[x + hx + 1, y + hy + 1]
                update!(s, -b - d, -b + d, value)
            end
        end
    end
end


function projectObject!(
    geom::ProjectionGeometry{T},
    pixelSize::T,
    object::Shifted{T},
    value::T,
    buffer::Array{Buffer{T}, 2}
) where T <: Real

    shiftedGeom = ProjectionGeometry(
        geom.sourcePos - object.shift,
        geom.detectorPos - object.shift,
        geom.wStep,
        geom.hStep
    )

    projectObject!(shiftedGeom, pixelSize, object.object, value, buffer)
end


function projectScene(
    sceneAt::Function,
    projGeoms::Vector{ProjectionGeometry{T}},
    pixelSize::T,
    camSizePx::Int,
    saveImg::Function,
    oversamplingFactor::Int=1,
    additive::Bool=false,
    batchNr::Int=1,
    nrBatches::Int=1
) where T <: Real

    pixelSizeScan = pixelSize / oversamplingFactor
    camSizePxScan = (camSizePx + 1) * oversamplingFactor - 1

    buffer = Array{Buffer{T}}(undef, camSizePxScan, camSizePxScan)

    for idx in LinearIndices(buffer)
        buffer[idx] = if additive
            AdditiveValueBuffer{T}()
        else
            ValueBuffer{T}()
        end
    end

    for i in batchNr : nrBatches : size(projGeoms)[1]
        scene = sceneAt(i, length(projGeoms))
        geom = projGeoms[i]

        if length(scene) > 0
            foreach(clear!, buffer)

            for k in length(scene) : -1 : 1
                item = scene[k]
                projectObject!(geom, pixelSizeScan,
                               item.object, item.value, buffer)
            end
            img = downsample(map(get, buffer), oversamplingFactor)
            saveImg(img, i)
        end
    end
end
