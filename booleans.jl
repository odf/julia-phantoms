mutable struct Range{T}
    near::T
    far::T
end


struct And{T} <: SceneObject{T}
    objA::SceneObject{T}
    objB::SceneObject{T}
end


function pushIntersections!(
    and::And{T},
    rayOrigin::Vec3{T},
    rayTarget::Vec3{T},
    accumulate::Function
) where T
    listA = Vector{Range{T}}()
    pushIntersections!(and.objA, rayOrigin, rayTarget,
                       (near, far) -> push!(listA, Range(near, far)))

    listB = Vector{Range{T}}()
    pushIntersections!(and.objB, rayOrigin, rayTarget,
                       (near, far) -> push!(listB, Range(near, far)))

    a() = listA[1]
    b() = listB[1]

    while !isempty(listA) && !isempty(listB)
        if a().far <= b().near
            popfirst!(listA)
        elseif b().far <= a().near
            popfirst!(listB)
        elseif a().far < b().far
            accumulate(max(a().near, b().near), a().far)
            b().near = a().far
            popfirst!(listA)
        else
            accumulate(max(a().near, b().near), b().far)
            a().near = b().far
            popfirst!(listB)
        end
    end
end


function hullCorners(and::And{T})::Vector{Vec3{T}} where T
    return hullCorners(and.objA)
end


struct AndNot{T} <: SceneObject{T}
    objA::SceneObject{T}
    objB::SceneObject{T}
end


function pushIntersections!(
    andnot::AndNot{T},
    rayOrigin::Vec3{T},
    rayTarget::Vec3{T},
    accumulate::Function
) where T
    listA = Vector{Range{T}}()
    pushIntersections!(andnot.objA, rayOrigin, rayTarget,
                       (near, far) -> push!(listA, Range(near, far)))

    listB = Vector{Range{T}}()
    pushIntersections!(andnot.objB, rayOrigin, rayTarget,
                       (near, far) -> push!(listB, Range(near, far)))

    a() = listA[1]
    b() = listB[1]

    while !isempty(listA) && !isempty(listB)
        if a().far <= b().near
            accumulate(a().near, a().far)
            popfirst!(listA)
        elseif b().far <= a().near
            popfirst!(listB)
        else
            if b().near > a().near
                accumulate(a().near, b().near)
            end

            if a().far <= b().far
                popfirst!(listA)
            else
                a().near = b().far
                popfirst!(listB)
            end
        end
    end

    for range in listA
        accumulate(range.near, range.far)
    end
end


function hullCorners(andnot::AndNot{T})::Vector{Vec3{T}} where T
    return hullCorners(andnot.objA)
end


struct Or{T} <: SceneObject{T}
    objA::SceneObject{T}
    objB::SceneObject{T}
end


function pushIntersections!(
    or::Or{T},
    rayOrigin::Vec3{T},
    rayTarget::Vec3{T},
    accumulate::Function
) where T
    listA = Vector{Range{T}}()
    pushIntersections!(or.objA, rayOrigin, rayTarget,
                       (near, far) -> push!(listA, Range(near, far)))

    listB = Vector{Range{T}}()
    pushIntersections!(or.objB, rayOrigin, rayTarget,
                       (near, far) -> push!(listB, Range(near, far)))

    a() = listA[1]
    b() = listB[1]

    while !isempty(listA) && !isempty(listB)
        if a().far < b().near
            accumulate(a().near, a().far)
            popfirst!(listA)
        elseif b().far < a().near
            accumulate(b().near, b().far)
            popfirst!(listB)
        elseif a().far < b().far
            b().near = min(a().near, b().near)
            popfirst!(listA)
        else
            a().near = min(a().near, b().near)
            popfirst!(listB)
        end
    end

    for range in listA
        accumulate(range.near, range.far)
    end
    for range in listB
        accumulate(range.near, range.far)
    end
end


function hullCorners(or::Or{T})::Vector{Vec3{T}} where T
    return vcat(hullCorners(or.objA), hullCorners(or.objB))
end


struct Xor{T} <: SceneObject{T}
    objA::SceneObject{T}
    objB::SceneObject{T}
end


function pushIntersections!(
    xor::Xor{T},
    rayOrigin::Vec3{T},
    rayTarget::Vec3{T},
    accumulate::Function
) where T
    listA = Vector{Range{T}}()
    pushIntersections!(xor.objA, rayOrigin, rayTarget,
                       (near, far) -> push!(listA, Range(near, far)))

    listB = Vector{Range{T}}()
    pushIntersections!(xor.objB, rayOrigin, rayTarget,
                       (near, far) -> push!(listB, Range(near, far)))

    a() = listA[1]
    b() = listB[1]

    while !isempty(listA) && !isempty(listB)
        if a().far < b().near
            accumulate(a().near, a().far)
            popfirst!(listA)
        elseif b().far < a().near
            accumulate(b().near, b().far)
            popfirst!(listB)
        else
            if a().near < b().near
                accumulate(a().near, b().near)
            elseif b().near < a().near
                accumulate(b().near, a().near)
            end
            if a().far < b().far
                b().near = a().far
                popfirst!(listA)
            else
                a().near = b().far
                popfirst!(listB)
            end
        end
    end

    for range in listA
        accumulate(range.near, range.far)
    end
    for range in listB
        accumulate(range.near, range.far)
    end
end


function hullCorners(xor::Xor{T})::Vector{Vec3{T}} where T
    return vcat(hullCorners(xor.objA), hullCorners(xor.objB))
end
