function goodRange(
    values::Vector{T}, low::Int, high::Int
)::UnitRange{Int} where T <: Real

    l = max(Int(floor(min(values...))), low)
    h = min(Int(ceil(max(values...))), high)
    return l : h
end


function rasterizeObject!(
    voxelSize::T,
    object,
    value::T,
    buffer::Array{Buffer{T}, 2}
) where T <: Real

    nx, ny = size(buffer)
    hx, hy = div(nx, 2), div(ny, 2)

    corners = hullCorners(object)
    xrange = goodRange(map(p -> p.x / voxelSize, corners), -hx, nx - hx - 1)
    yrange = goodRange(map(p -> p.y / voxelSize, corners), -hy, ny - hy - 1)

    wStep, hStep = voxelSize * [Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0)]
    downZ = Vec3(0.0, 0.0, 1.0)

    for y in yrange
        for x in xrange
            p = wStep * T(x) + hStep * T(y)
            s = buffer[x + hx + 1, y + hy + 1]
            pushIntersections!(object, p, p + downZ,
                               (near, far) -> update!(s, near, far, value))
        end
    end
end


function rasterize(
    buffer::Buffer{T},
    voxelSize::T,
    nrVoxelsTotal::Int,
    nrVoxelsRequested::Int,
    offset::Int
) where T <: Real

    h = div(nrVoxelsTotal, 2) - offset

    out = zeros(T, nrVoxelsRequested)

    for (near, far, value) in get(buffer)
        lo = clamp(near / voxelSize + h + 0.5, 0, nrVoxelsRequested)
        hi = clamp(far / voxelSize + h + 0.5, 0, nrVoxelsRequested)

        for i in max(1, Int(ceil(lo))) : min(nrVoxelsRequested, Int(ceil(hi)))
            out[i] += value * (min(i, hi) - max(i - 1, lo))
        end
    end

    return out
end


function rasterizeScene(
    scene::Vector{SceneItem{T}},
    volumeShape::Vector{T},
    voxelSize::T,
    saveImg::Function,
    oversamplingFactor::Int=1,
    additive::Bool=false,
    slabSize::Int=128,
    batchNr::Int=1,
    nrBatches::Int=1
) where T <: Real

    shapeVx = map(s -> 2 * Int(ceil(0.5 * s)) + 1, volumeShape / voxelSize)

    rawVoxelSz = voxelSize / oversamplingFactor
    rawShapeVx = map(n -> (n + 1) * oversamplingFactor - 1, shapeVx[1:2])

    buffer = Array{Buffer{T}}(undef, rawShapeVx[1], rawShapeVx[2])

    for idx in LinearIndices(buffer)
        buffer[idx] = if additive
            AdditiveIntervalBuffer{T}()
        else
            IntervalBuffer{T}()
        end
    end

    for item in scene[length(scene):-1:1]
        rasterizeObject!(rawVoxelSz, item.object, item.value, buffer)
    end

    n = shapeVx[3]
    start = Int(floor(n / nrBatches * (batchNr - 1)))
    stop = Int(floor(n / nrBatches * batchNr))

    for offset in start : slabSize : stop
        m = min(slabSize, stop - offset)

        lines = map(b -> rasterize(b, voxelSize, n, m, offset), buffer)

        slab = Array{T}(undef, rawShapeVx[1], rawShapeVx[2], m)
        for x in 1 : rawShapeVx[1]
            for y in 1 : rawShapeVx[2]
                slab[x, y, :] = lines[x, y]
            end
        end

        for z in 1 : m
            saveImg(downsample(slab[:, :, z], oversamplingFactor), z + offset)
        end
    end
end
