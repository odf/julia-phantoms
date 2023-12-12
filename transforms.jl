# -- The underlying Transformed and Shifted data types and their operations

struct Transformed{T <: Real} <: SceneObject{T}
    mat::Mat3{T}
    invmat::Mat3{T}
    shift::Vec3{T}
    object::SceneObject{T}
    Transformed{T}(mat, shift, obj) where T <: Real =
        new(mat, inv(mat), shift, obj)
end


function pushIntersections!(
    trans::Transformed{T},
    rayOrigin::Vec3{T},
    rayTarget::Vec3{T},
    accumulate::Function
) where T <: Real
    origin = trans.invmat * (rayOrigin - trans.shift)
    target = trans.invmat * (rayTarget - trans.shift)
    scale = norm(rayTarget - rayOrigin) / norm(target - origin)

    pushIntersections!(trans.object, origin, target,
                       (near, far) -> accumulate(near * scale, far * scale))
end


function hullCorners(trans::Transformed{T})::Vector{Vec3{T}} where T <: Real
    return map(p -> trans.mat * p + trans.shift, hullCorners(trans.object))
end


struct Shifted{T <: Real} <: SceneObject{T}
    shift::Vec3{T}
    object::SceneObject{T}
end


function pushIntersections!(
    trans::Shifted{T},
    rayOrigin::Vec3{T},
    rayTarget::Vec3{T},
    accumulate::Function
) where T <: Real
    origin = rayOrigin - trans.shift
    target = rayTarget - trans.shift

    pushIntersections!(trans.object, origin, target, accumulate)
end


function hullCorners(trans::Shifted{T})::Vector{Vec3{T}} where T <: Real
    return map(p -> p + trans.shift, hullCorners(trans.object))
end


# -- Optimizing overloaded factory functions on top of those types

function _transformed(
    mat::Mat3{T}, shift::Vec3{T}, obj::SceneObject{T}
) where T <: Real
    return Transformed{T}(mat, shift, obj)
end


function _transformed(
    mat::Mat3{T}, shift::Vec3{T}, base::Transformed{T}
) where T <: Real
    return _transformed(mat * base.mat, shift + mat * base.shift, base.object)
end


function _transformed(
    mat::Mat3{T}, shift::Vec3{T}, base::Shifted{T}
) where T <: Real
    return _transformed(mat, shift + mat * base.shift, base.object)
end


function _transformed(
    mat::Mat3{T}, shift::Vec3{T}, bool::And{T}
) where T <: Real
    return And(_transformed(mat, shift, bool.objA),
               _transformed(mat, shift, bool.objB))
end


function _transformed(
    mat::Mat3{T}, shift::Vec3{T}, bool::AndNot{T}
) where T <: Real
    return AndNot(_transformed(mat, shift, bool.objA),
                  _transformed(mat, shift, bool.objB))
end


function _transformed(
    mat::Mat3{T}, shift::Vec3{T}, bool::Or{T}
) where T <: Real
    return Or(_transformed(mat, shift, bool.objA),
              _transformed(mat, shift, bool.objB))
end


function _transformed(
    mat::Mat3{T}, shift::Vec3{T}, bool::Xor{T}
) where T <: Real
    return Xor(_transformed(mat, shift, bool.objA),
               _transformed(mat, shift, bool.objB))
end


function _shifted(shift::Vec3{T}, obj::SceneObject{T}) where T <: Real
    return Shifted{T}(shift, obj)
end


function _shifted(shift::Vec3{T}, obj::Sphere{T}) where T <: Real
    return Sphere(obj.center + shift, obj.radius)
end


function _shifted(shift::Vec3{T}, base::Transformed{T}) where T <: Real
    return _transformed(base.mat, shift + base.shift, base.object)
end


function _shifted(shift::Vec3{T}, base::Shifted{T}) where T <: Real
    return _shifted(shift + base.shift, base.object)
end


function _shifted(shift::Vec3{T}, bool::And{T}) where T <: Real
    return And(_shifted(shift, bool.objA), _shifted(shift, bool.objB))
end


function _shifted(shift::Vec3{T}, bool::AndNot{T}) where T <: Real
    return AndNot(_shifted(shift, bool.objA), _shifted(shift, bool.objB))
end


function _shifted(shift::Vec3{T}, bool::Or{T}) where T <: Real
    return Or(_shifted(shift, bool.objA), _shifted(shift, bool.objB))
end


function _shifted(shift::Vec3{T}, bool::Xor{T}) where T <: Real
    return Xor(_shifted(shift, bool.objA), _shifted(shift, bool.objB))
end


# -- Some curried high-level functions to form the API

function transform(mat::Mat3{T}, shift::Vec3{T}) where T <: Real
    return object -> _transformed(mat, shift, object)
end


function transform(mat::Mat3{T}) where T <: Real
    return transform(mat, Vec3(0.0, 0.0, 0.0))
end


function shift(vector::Vec3{T}) where T <: Real
    return object -> _shifted(vector, object)
end


function scale(xscale::T, yscale::T, zscale::T) where T <: Real
    mat = (Vec3(xscale, 0.0, 0.0),
           Vec3(0.0, yscale, 0.0),
           Vec3(0.0, 0.0, zscale))

    return transform(mat)
end


function rotate(axis::Vec3{T}, angle::T) where T <: Real
    u = normalize(axis)

    mat = (cos(angle) * (Vec3(1.0, 0.0, 0.0),
                         Vec3(0.0, 1.0, 0.0),
                         Vec3(0.0, 0.0, 1.0))
           +
           sin(angle) * (Vec3( 0.0, -u.z,  u.y),
                         Vec3( u.z,  0.0, -u.x),
                         Vec3(-u.y,  u.x,  0.0))
           +
           (1.0 - cos(angle)) * (Vec3(u.x*u.x, u.x * u.y, u.x * u.z),
                                 Vec3(u.y*u.x, u.y * u.y, u.y * u.z),
                                 Vec3(u.z*u.x, u.z * u.y, u.z * u.z))
           )

    return transform(mat)
end
