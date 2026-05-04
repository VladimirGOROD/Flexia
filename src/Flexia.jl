module Flexia
# using Pkg; Pkg.activate(".")
using StaticArrays
using Makie
using LinearAlgebra
# export MBSystem
# export Body2D

export Body2D
export FixedJoint, HingeJoint, TorsionalSpring, TrajectoryJoint
export MBSystem2D

export set_position_on_first_body!, set_position_on_second_body!, set_position!, set_rotation!
export add!, assemble!, get_body_position_dofs, get_body_velocity_dofs, number_of_dofs, last_body_dof, last_lm_dof, get_boundary_points, get_lms, get_spring_moment 
export cros!, static_solver!

export animate
export Marker
export test_func


test_func() = 1

abstract type AbstractBody2D end
abstract type AbstractJoint2D end
abstract type AbstractMarker2D end

include("solvers.jl")
include("system.jl")
include("bodies.jl")
include("joints.jl")
include("visualize.jl")
include("markers.jl")

function getdofs(sys::MBSystem2D, body::Body2D)
    if (body.index == -1)
        return 0:0
    elseif (body.index == 1)
        return 1:sys.bodiesdofs[1]
    else
        return (sys.bodiesdofs[body.index-1]+1):sys.bodiesdofs[body.index]
    end
end


end



"""
# end # module Flexia

"""