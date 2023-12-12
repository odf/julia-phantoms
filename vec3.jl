import Base: +, -, *, /, inv
using LinearAlgebra


struct Vec3{T <: Real}
    x::T
    y::T
    z::T
end

(-w::Vec3{T}) where T = Vec3(-w.x, -w.y, -w.z)

(v::Vec3{T} + w::Vec3{T}) where T = Vec3(v.x + w.x, v.y + w.y, v.z + w.z)
(v::Vec3{T} - w::Vec3{T}) where T = Vec3(v.x - w.x, v.y - w.y, v.z - w.z)

(s::T * v::Vec3{T}) where T = Vec3(s * v.x, s * v.y, s * v.z)
(v::Vec3{T} * s::T) where T = s * v
(v::Vec3{T} / s::T) where T = (1.0 / s) * v


dot(v::Vec3{T}, w::Vec3{T}) where T = v.x * w.x + v.y * w.y + v.z * w.z

cross(v::Vec3{T}, w::Vec3{T}) where T =
    Vec3(v.y * w.z - v.z * w.y,
         v.z * w.x - v.x * w.z,
         v.x * w.y - v.y * w.x)

norm(v::Vec3{T}) where T = sqrt(dot(v, v))

normalize(v::Vec3{T}) where T = v / norm(v)


Mat3{T} = Tuple{Vec3{T}, Vec3{T}, Vec3{T}}


function inv(b::Mat3{T})::Mat3{T} where T
    det = dot(b[1], cross(b[2], b[3]))

    c1 = cross(b[2], b[3]) / det
    c2 = cross(b[3], b[1]) / det
    c3 = cross(b[1], b[2]) / det

    return (
        Vec3(c1.x, c2.x, c3.x),
        Vec3(c1.y, c2.y, c3.y),
        Vec3(c1.z, c2.z, c3.z)
    )
end


(a::Mat3{T} + b::Mat3{T}) where T = (a[1] + b[1], a[2] + b[2], a[3] + b[3])

(b::Mat3{T} * v::Vec3{T}) where T = v.x * b[1] + v.y * b[2] + v.z * b[3]
(a::Mat3{T} * b::Mat3{T}) where T = (a * b[1], a * b[2], a * b[3])
(s::T * b::Mat3{T}) where T = (s * b[1], s * b[2], s * b[3])
