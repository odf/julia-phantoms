using AnalyticPhantoms


set(val, obj) = SceneItem(obj, val)


function centralObstacleDistance(initialRadius, finalRadius, center, nrSteps)
    a = (finalRadius / initialRadius) ^ (1 / nrSteps)
    m = div(nrSteps, 10)

    return (pos, i) -> norm(pos - center) - initialRadius * a^(i - i % m)
end


function gridDimensions(volumeShape, numberOfCells)
    cellSize = (prod(volumeShape) / numberOfCells) ^ (1 / size(volumeShape, 1))
    return Int.(floor.(volumeShape ./ cellSize))
end


function signedDistanceField(xhi, yhi, zhi, pos)
    sdf = fill(Inf, (xhi, yhi, zhi))

    for x in 1 : xhi
        for y in 1 : yhi
            for z in 1 : zhi
                p = pos(x, y, z)
                sdf[x, y, z] = min(
                    norm(p - pos(1, y, z)), norm(p - pos(xhi, y, z)),
                    norm(p - pos(x, 1, z)), norm(p - pos(x, yhi, z)),
                    norm(p - pos(x, y, 1)), norm(p - pos(x, y, zhi))
                )
            end
        end
    end

    return sdf
end


function getDist(sdf, x, y, z)
    xhi, yhi, zhi = size(sdf)

    if 1 <= x <= xhi && 1 <= y <= yhi && 1 <= z <= zhi
        sdf[x, y, z]
    else
        0
    end
end


function insertSphere(
    sdf::AbstractArray{Float64, 3}, radius, cx, cy, cz, pos; reverse=false
)
    sign = reverse ? -1 : 1
    xhi, yhi, zhi = size(sdf)

    sdf[cx, cy, cz] = -sign * radius
    queue = [(cx, cy, cz)]

    while length(queue) > 0
        x0, y0, z0 = pop!(queue)

        for x in max(1, x0 - 1) : min(xhi, x0 + 1)
            for y in max(1, y0 - 1) : min(yhi, y0 + 1)
                for z in max(1, z0 - 1) : min(zhi, z0 + 1)
                    d = norm(pos(x, y, z) - pos(cx, cy, cz)) - radius
                    if d < sign * sdf[x, y, z]
                        sdf[x, y, z] = sign * d
                        push!(queue, (x, y, z))
                    end
                end
            end
        end
    end
end


insertSphereReverse(sdf, radius, cx, cy, cz, pos) =
    insertSphere(sdf, radius, cx, cy, cz, pos, reverse=true)


function insertSphericalShell(
    sdf::AbstractArray{Float64, 3}, innerRadius, outerRadius, cx, cy, cz, pos
)
    insertSphere(sdf, outerRadius, cx, cy, cz, pos)

    xhi, yhi, zhi = size(sdf)

    sdf[cx, cy, cz] = innerRadius
    queue = [(cx, cy, cz)]

    while length(queue) > 0
        x0, y0, z0 = pop!(queue)

        for x in max(1, x0 - 1) : min(xhi, x0 + 1)
            for y in max(1, y0 - 1) : min(yhi, y0 + 1)
                for z in max(1, z0 - 1) : min(zhi, z0 + 1)
                    d = norm(pos(x, y, z) - pos(cx, cy, cz)) - innerRadius
                    if d < 0 && sdf[x, y, z] < 0
                        sdf[x, y, z] = -d
                        push!(queue, (x, y, z))
                    end
                end
            end
        end
    end
end


function findPlacement(cx, cy, cz, maxDist, dist)
    cd = dist(cx, cy, cz)

    while true
        bestx, besty, bestz, bestd = cx, cy, cz, cd

        for x in cx - 1 : cx + 1
            for y in cy - 1 : cy + 1
                for z in cz - 1 : cz + 1
                    d = dist(x, y, z)
                    if d > bestd
                        bestx, besty, bestz, bestd = x, y, z, d
                    end
                end
            end
        end

        if (bestx, besty, bestz) == (cx, cy, cz) || bestd >= maxDist
            break
        else
            cx, cy, cz, cd = bestx, besty, bestz, bestd
        end
    end

    return cx, cy, cz, min(maxDist, cd)
end


function transform(center, radius)
    (x, y, z) -> Vec3((((x, y, z) .- center) ./ (center .- 1) .* radius)...)
end


function makeBoundingCylinder(
    sdf::AbstractArray{Float64, 3}, radius, cx, cy, pos
)
    xhi, yhi, zhi = size(sdf)

    for x in 1 : xhi
        for y in 1 : yhi
            d = radius - norm(pos(x, y, 1) - pos(cx, cy, 1))
            for z in 1 : zhi
                if d < sdf[x, y, z]
                    sdf[x, y, z] = d
                end
            end
        end
    end
end


function makeSpheres(sceneRadius, height, nrElements, maxElemRadius, extraDist)
    xhi, yhi, zhi = shape = (399, 399, 399)
    xmid, ymid, zmid = center = div.(shape .+ 1, 2)

    pos = transform(center, (sceneRadius, sceneRadius, height / 2.0))

    sdf = signedDistanceField(xhi, yhi, zhi, pos)
    makeBoundingCylinder(sdf, sceneRadius, xmid, ymid, pos)

    spheres = []

    for i in 1 : nrElements
        cx, cy, cz, radius = findPlacement(
            rand(1 : xhi), rand(1 : yhi), rand(1 : zhi), maxElemRadius,
            (x, y, z) -> min(getDist(sdf, x, y, z), extraDist(pos(x, y, z), i))
        )

        if radius > 0
            insertSphere(sdf, radius, cx, cy, cz, pos)
            push!(spheres, (pos(cx, cy, cz), radius))
        end
    end

    return spheres
end


function spheres(volShape, nrElements, baseSize, extraDist=((pos,i) -> Inf))
    fillValues = [-0.05, -0.02, -0.01, -0.005, 0.005, 0.01, 0.02, 0.05]

    height = volShape[3]
    radius = 0.5 * min(volShape[1], volShape[2])
    sampleRadius = 0.79 * radius

    maxRadius = 40 * baseSize * minimum(volShape)

    scene = Array{SceneItem{Float64}, 1}()

    for (c, r) in makeSpheres(
        sampleRadius, height, nrElements, maxRadius, extraDist
    )
        push!(scene, set(rand(fillValues), Sphere(c, r)))
    end

    push!(scene, set(0.05, Cylinder(sampleRadius, height)))
    push!(scene, set(0.08, Cylinder(0.85 * radius, 1.1 * height)))
    push!(scene, set(-0.08, Cylinder(0.8 * radius, 1.1 * height)))

    return scene, true
end


function board(volShape, nrElements; maxRadius=Inf, extraDist=((pos,i) -> Inf))
    width, depth, height = volShape

    xhi, yhi, zhi = shape = gridDimensions((width, depth, height), 1e8)
    center = div.(shape .+ 1, 2)
    pos = transform(center, (width / 2.0, depth / 2.0, height / 2.0))
    sdf = signedDistanceField(xhi, yhi, zhi, pos)

    scene = Array{SceneItem{Float64}, 1}()
    push!(scene, set(0.05, Cuboid(width, depth, height)))

    for i in 1 : nrElements
        cx, cy, cz, radius = findPlacement(
            rand(1 : xhi), rand(1 : yhi), rand(1 : zhi), maxRadius,
            (x, y, z) -> min(getDist(sdf, x, y, z), extraDist(pos(x, y, z), i))
        )

        if radius > 0
            insertSphere(sdf, radius, cx, cy, cz, pos)

            candidates = if radius > depth / 20.0
                [-0.05, -0.01, 0.0, 0.01]
            elseif radius > depth / 100.0
                [-0.01, 0.0, 0.01]
            else
                [-0.01, 0.0, 0.01, 0.05]
            end

            value = rand(candidates)

            if value != 0
                push!(scene, set(value, Sphere(pos(cx, cy, cz), radius)))
            end
        end
    end

    return scene, true
end


function pcb(volShape, nrElements; maxRadius=Inf)
    width, depth, height = volShape

    xhi, yhi, zhi = shape = gridDimensions((width, depth, height), 1e8)
    center = div.(shape .+ 1, 2)
    pos = transform(center, (width / 2.0, depth / 2.0, height / 2.0))
    sdf1 = signedDistanceField(xhi, yhi, zhi, pos)
    sdf2 = signedDistanceField(xhi, yhi, zhi, pos)

    scene = Array{SceneItem{Float64}, 1}()
    push!(scene, set(0.05, Cuboid(width, depth, height)))

    for i in 1 : nrElements
        mode = rand([1, 2])

        if mode == 1
            cx, cy, cz, radius = findPlacement(
                rand(1 : xhi), rand(1 : yhi), rand(1 : zhi), maxRadius,
                (x, y, z) -> getDist(sdf1, x, y, z)
            )

            value = rand([-0.05, 0.0, 0.05])

            insertSphere(sdf1, radius, cx, cy, cz, pos)
            if value > 0
                insertSphere(sdf2, radius, cx, cy, cz, pos)
                insertSphereReverse(sdf2, radius * 2 / 3, cx, cy, cz, pos)
            end

            if value != 0
                push!(scene, set(value, Sphere(pos(cx, cy, cz), radius)))
            end
        else
            cx, cy, cz, radius = findPlacement(
                rand(1 : xhi), rand(1 : yhi), rand(1 : zhi), maxRadius / 6,
                (x, y, z) -> -getDist(sdf2, x, y, z)
            )

            insertSphereReverse(sdf2, radius, cx, cy, cz, pos)

            value = rand([-0.1, -0.05, 0.0])
            if value != 0
                push!(scene, set(value, Sphere(pos(cx, cy, cz), radius)))
            end
        end
    end

    return scene, true
end


function concentric(volShape, nrElements; shifts=[[0.0, 0.0, 0.0]])
    width, depth, height = volShape

    xhi, yhi, zhi = shape = gridDimensions((width, depth, height), 1e8)
    center = div.(shape .+ 1, 2)
    pos = transform(center, (width / 2.0, depth / 2.0, height / 2.0))
    sdf = signedDistanceField(xhi, yhi, zhi, pos)

    scene = Array{SceneItem{Float64}, 1}()

    box = Cuboid(width, depth, height)
    push!(scene, set(0.05, box))

    maxRadius = min(width, depth, height) / 2

    r = 0.0125
    radii = []
    while r < maxRadius
        push!(radii, r)
        r *= sqrt(2)
    end

    for shift in shifts
        p = pos(center...) + Vec3(shift...)
        cx, cy, cz = center .+ Int.(round.(shift ./ volShape .* shape))

        for r in reverse(radii)
            push!(scene, set(-0.05, Sphere(p, r)))
            push!(scene, set(0.05, Sphere(p, 0.99 * r)))

            insertSphericalShell(sdf, 0.99 * r, r, cx, cy, cz, pos)
        end
    end

    sdf = min.(sdf, signedDistanceField(xhi, yhi, zhi, pos))

    for i in 1 : nrElements
        cx, cy, cz, radius = findPlacement(
            rand(1 : xhi), rand(1 : yhi), rand(1 : zhi), Inf,
            (x, y, z) -> getDist(sdf, x, y, z)
        )

        if radius > 0
            insertSphere(sdf, radius, cx, cy, cz, pos)

            candidates = if radius > 0.2
                [-0.05, -0.01, 0.0, 0.01]
            elseif radius > 0.1
                [-0.01, 0.0, 0.01]
            else
                [-0.01, 0.0, 0.01, 0.05]
            end

            value = rand(candidates)

            if value != 0
                push!(scene, set(value, Sphere(pos(cx, cy, cz), radius)))
            end
        end
    end

    return scene, true
end


sceneBuilders = Dict(
    "spheres" => spheres,
    "board" => board,
    "pcb" => pcb,
    "concentric" => concentric
)
