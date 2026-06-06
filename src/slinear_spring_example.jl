using Pkg; Pkg.activate("./examples")
using Flexia
using ForwardDiff
using GLMakie
using StaticArrays
# using CSV
# using DataFrames
const g = 9.81

ground = Body2D(1000, 1000, length = 10)

bd1 = Body2D(1, 1, length= 0.1)

fx_jnt = FixedJoint(ground)

sl_jnt = SliderJoint(ground, bd1)

set_position_on_first_body!(sl_jnt, SA[0., 0.1])
set_position_on_second_body!(sl_jnt, SA[0., 0.])
set_direction_on_first_body!(sl_jnt, SA[1., 0.])
set_direction_on_second_body!(sl_jnt, SA[1., 0.])

spring = LinearSpring(sl_jnt; stiffness = 1.)


sys = MBSystem2D()

add!(sys, ground)
add!(sys, bd1)
add!(sys, fx_jnt)
add!(sys, sl_jnt)
add!(sys, spring)

assemble!(sys)

func = sys.rhs

jacoby = (x) -> ForwardDiff.jacobian(func, x)

initial = zeros(number_of_dofs(sys))


bd1_x_ind, bd1_y_ind, bd1_t_ind = get_body_position_dofs(sys, bd1)

initial[bd1_x_ind] = 0.1
initial[bd1_y_ind] = 0.1

mass = get_mass_matrix(sys)
t_end = 10
time_span = LinRange(0,t_end, 501)

n_steps = length(time_span)

sol1 = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))
cros!(sol1, initial, mass, func, jacoby, step(time_span))
animate(sys, sol1, time_span, "Flexia/out/test.mp4"; framerate = 30, limits = (-0.1, 0.6, -0.1, 0.5))

