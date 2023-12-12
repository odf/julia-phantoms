struct Sphere{T <: Real} <: SceneObject{T}
    center::Vec3{T}
    radius::T
end


function pushIntersections!(
    s::Sphere{T},
    rayOrigin::Vec3{T},
    rayTarget::Vec3{T},
    accumulate::Function
) where T <: Real
    rayDir = normalize(rayTarget - rayOrigin)
    relativeOrigin = rayOrigin - s.center

    b = dot(relativeOrigin, rayDir)
    c = dot(relativeOrigin, relativeOrigin) - s.radius^2

    discriminant = b^2 - c

    if discriminant > 0.0
        d = sqrt(discriminant)
        accumulate(-b - d, -b + d)
    end
end


function hullCorners(s::Sphere{T})::Vector{Vec3{T}} where T <: Real
    r = s.radius
    c = s.center

    return [
        Vec3(c.x - r, c.y - r, c.z - r),
        Vec3(c.x - r, c.y - r, c.z + r),
        Vec3(c.x - r, c.y + r, c.z - r),
        Vec3(c.x - r, c.y + r, c.z + r),
        Vec3(c.x + r, c.y - r, c.z - r),
        Vec3(c.x + r, c.y - r, c.z + r),
        Vec3(c.x + r, c.y + r, c.z - r),
        Vec3(c.x + r, c.y + r, c.z + r)
    ]
end


struct Cylinder{T <: Real} <: SceneObject{T}
    radius::T
    height::T
end


function pushIntersections!(
    cyl::Cylinder{T},
    rayOrigin::Vec3{T},
    rayTarget::Vec3{T},
    accumulate::Function
) where T <: Real
    rayDir = normalize(rayTarget - rayOrigin)

    a = rayDir.x^2 + rayDir.y^2
    b = 2.0 * (rayOrigin.x * rayDir.x + rayOrigin.y * rayDir.y)
    c = rayOrigin.x^2 + rayOrigin.y^2 - cyl.radius^2

    discriminant = if a == 0.0 (-c) else b^2 - 4.0 * a * c end

    if discriminant > 0.0
        near = (-b - sqrt(discriminant)) / (2.0 * a)
        far = (-b + sqrt(discriminant)) / (2.0 * a)

        if rayDir.z == 0.0
            if !(-cyl.height <= 2.0 * rayOrigin.z <= cyl.height)
                far = near
            end
        else
            t1 = (-cyl.height / 2.0 - rayOrigin.z) / rayDir.z
            t2 = (cyl.height / 2.0 - rayOrigin.z) / rayDir.z
            near = max(near, min(t1, t2))
            far = min(far, max(t1, t2))
        end

        if far > near
            accumulate(near, far)
        end
    end
end


function hullCorners(cyl::Cylinder{T})::Vector{Vec3{T}} where T <: Real
    r = cyl.radius
    h = cyl.height / 2

    return [
        Vec3(-r, -r, -h),
        Vec3(-r, -r,  h),
        Vec3(-r,  r, -h),
        Vec3(-r,  r,  h),
        Vec3( r, -r, -h),
        Vec3( r, -r,  h),
        Vec3( r,  r, -h),
        Vec3( r,  r,  h)
    ]
end


struct Cuboid{T <: Real} <: SceneObject{T}
    xsize::T
    ysize::T
    zsize::T
end


function pushIntersections!(
    cub::Cuboid{T},
    rayOrigin::Vec3{T},
    rayTarget::Vec3{T},
    accumulate::Function
) where T <: Real
    rayDir = normalize(rayTarget - rayOrigin)

    near = -Inf
    far = Inf

    if rayDir.x == 0.0
        if 2 * abs(rayOrigin.x) > cub.xsize
            return
        end
    else
        t1 = (-cub.xsize / 2.0 - rayOrigin.x) / rayDir.x
        t2 = (cub.xsize / 2.0 - rayOrigin.x) / rayDir.x
        near = max(near, min(t1, t2))
        far = min(far, max(t1, t2))

        if near >= far
            return
        end
    end

    if rayDir.y == 0.0
        if 2 * abs(rayOrigin.y) > cub.ysize
            return
        end
    else
        t1 = (-cub.ysize / 2.0 - rayOrigin.y) / rayDir.y
        t2 = (cub.ysize / 2.0 - rayOrigin.y) / rayDir.y
        near = max(near, min(t1, t2))
        far = min(far, max(t1, t2))

        if near >= far
            return
        end
    end

    if rayDir.z == 0.0
        if 2 * abs(rayOrigin.z) > cub.zsize
            return
        end
    else
        t1 = (-cub.zsize / 2.0 - rayOrigin.z) / rayDir.z
        t2 = (cub.zsize / 2.0 - rayOrigin.z) / rayDir.z
        near = max(near, min(t1, t2))
        far = min(far, max(t1, t2))

        if near >= far
            return
        end
    end

    accumulate(near, far)
end


function hullCorners(cub::Cuboid{T})::Vector{Vec3{T}} where T <: Real
    a = cub.xsize / 2
    b = cub.ysize / 2
    c = cub.zsize / 2

    return [
        Vec3(-a, -b, -c),
        Vec3(-a, -b,  c),
        Vec3(-a,  b, -c),
        Vec3(-a,  b,  c),
        Vec3( a, -b, -c),
        Vec3( a, -b,  c),
        Vec3( a,  b, -c),
        Vec3( a,  b,  c)
    ]
end
