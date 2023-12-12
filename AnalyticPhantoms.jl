module AnalyticPhantoms

abstract type SceneObject{T <: Real} end

include("vec3.jl")
include("rayBuffer.jl")
include("primitives.jl")
include("booleans.jl")
include("transforms.jl")
include("trajectories.jl")
include("projGeometries.jl")
include("downsample.jl")
include("project.jl")
include("rasterize.jl")
include("run.jl")


export Vec3, norm
export SceneItem, SceneObject, Sphere, Cylinder, Cuboid
export transform, shift, scale, rotate, And, Or, Xor, AndNot
export BeamGeometry, computeGeometry, writeConeBeamGeometry
export Trajectory, HelicalTrajectory, PlaneGridTrajectory, StadiumTrajectory
export helicalTrajectory, planeGridTrajectory, stadiumTrajectory
export makeMisalignment, convertProjectionGeometries, makeProjectionGeometries
export writeProjectionGeometries, readProjectionGeometries
export pushOutFunction, pushFactorCircular, pushFactorPolygonal
export projectScene, rasterizeScene, projectAndWrite, rasterizeAndWrite


end
