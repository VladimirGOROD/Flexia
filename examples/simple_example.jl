using Pkg; Pkg.activate("./examples")
using Flexia
using GLMakie
using ForwardDiff
using JSON
using LinearAlgebra

using StaticArrays

const g = 9.81

bd1 = Body2D(1, 1)
bd2 = Body2D(1, 1)

bd2.forces[2] = (x, t) -> -bd2.mass * g

jnt1 = FixedJoint(bd1)
jnt2 = HingeJoint(bd1, bd2)

set_position_on_second_body!(jnt2, SA[-1., 0])

sys = MBSystem2D()

add!(sys, bd1)
add!(sys, bd2)

add!(sys, jnt1)
add!(sys, jnt2)

if (!assemble!(sys))
    println("Assembling failed!")
end

func = sys.rhs

jacoby = (x) -> ForwardDiff.jacobian(func, x)

bd1_x_ind, bd1_y_ind, _ = get_body_position_dofs(sys, bd1)
bd2_x_ind, bd2_y_ind, bd2_t_ind = get_body_position_dofs(sys, bd2)

initial = zeros(number_of_dofs(sys))

func(initial)
jacoby(initial)

mass = get_mass_matrix(sys)

time_start = 0
time_end = 10
time_step = 1

time_span = 0:0.001:10

sol2 = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))
cros!(sol2, initial, mass, func, jacoby, step(time_span))
animate(sys, sol2, time_span, "out/simple_bar.mp4"; framerate = floor(Int64, 1.0 / step(time_span)), limits = (-2,13, -2, 13))
